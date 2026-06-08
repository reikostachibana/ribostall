#!/usr/bin/env python3
import argparse
import datetime
import gzip
import json
import logging
import pickle
import sys
from pathlib import Path

import pandas as pd

from functions_stall_sites import (
    call_stalls,
    codonize_counts_cds,
    consensus_stalls_across_reps,
    consensus_to_long_df,
    filter_tx,
    stalls_to_long_df,
)


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(processName)s  %(message)s",
)


def parse_groups(groups_arg):
    groups = {}
    for block in groups_arg.split(";"):
        block = block.strip()
        if not block:
            continue
        name, reps = block.split(":", 1)
        groups[name] = [rep.strip() for rep in reps.split(",") if rep.strip()]
    return groups


def write_jsonl(df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for rec in df.to_dict(orient="records"):
            json.dump(rec, f)
            f.write("\n")


def empty_consensus_df():
    return pd.DataFrame(columns=["group", "transcript", "tx_id", "gene", "pos_codon"])


def empty_replicate_df():
    return pd.DataFrame(
        columns=["group", "replicate", "gene", "tx_id", "transcript", "pos_codon", "obs", "z"]
    )


def add_option(parser, *option_strings, legacy=None, **kwargs):
    parser.add_argument(*option_strings, **kwargs)
    if legacy:
        legacy_kwargs = dict(kwargs)
        legacy_kwargs["help"] = argparse.SUPPRESS
        legacy_kwargs["default"] = argparse.SUPPRESS
        legacy_kwargs.pop("required", None)
        dest = kwargs.get("dest", option_strings[0].lstrip("-").replace("-", "_"))
        parser.add_argument(legacy, dest=dest, **legacy_kwargs)


def parse_args():
    parser = argparse.ArgumentParser(description="Detect ribosome stall sites from adjusted coverage.")
    parser.add_argument("--pickle", required=True, metavar="PATH", help="Path to adjusted coverage pickle.gz file")
    parser.add_argument("--ribo", required=True, metavar="PATH", help="Path to .ribo file, recorded for provenance")
    parser.add_argument(
        "--alias",
        action="store_true",
        help="Record that APPRIS-style aliasing was used upstream; this script does not reopen the .ribo file",
    )
    parser.add_argument(
        "--groups",
        required=True,
        metavar="GROUPS",
        help="Semicolon-separated group:rep1,rep2 definitions",
    )
    add_option(
        parser,
        "--tx-threshold",
        legacy="--tx_threshold",
        type=float,
        default=1.0,
        metavar="THRESHOLD",
        help="Minimum mean reads/nt in CDS for transcript filtering",
    )
    add_option(
        parser,
        "--tx-min-reps",
        legacy="--tx_min_reps",
        type=int,
        default=2,
        metavar="N",
        help="Minimum number of replicates that must pass the transcript threshold",
    )
    add_option(
        parser,
        "--min-z",
        legacy="--min_z",
        type=float,
        default=1.0,
        metavar="Z",
        help="Minimum z-score for a stall site",
    )
    add_option(
        parser,
        "--min-reads",
        legacy="--min_reads",
        type=int,
        default=2,
        metavar="READS",
        help="Minimum reads for a stall site",
    )
    add_option(
        parser,
        "--trim-edges",
        legacy="--trim_edges",
        type=int,
        default=10,
        metavar="CODONS",
        help="Trim this many codons after start and before stop",
    )
    parser.add_argument(
        "--pseudocount",
        type=float,
        default=0.5,
        metavar="VALUE",
        help="Pseudocount used for log/z-score stall calling",
    )
    add_option(
        parser,
        "--stall-min-reps",
        legacy="--stall_min_reps",
        type=int,
        default=2,
        metavar="N",
        help="Minimum number of replicates that must support a consensus stall site",
    )
    parser.add_argument(
        "--tol",
        type=int,
        default=0,
        metavar="CODONS",
        help="Tolerance window for matching stall sites across replicates",
    )
    add_option(
        parser,
        "--min-sep",
        legacy="--min_sep",
        type=int,
        default=7,
        metavar="CODONS",
        help="Minimum separation between consensus sites",
    )
    parser.add_argument(
        "--out-json",
        default="../ribostall_results/stall_sites.jsonl",
        metavar="PATH",
        help="Consensus stall-site JSONL output",
    )
    parser.add_argument(
        "--out-replicate-json",
        default=None,
        metavar="PATH",
        help="Per-replicate stall-site JSONL output; defaults next to --out-json",
    )
    parser.add_argument(
        "--summary-json",
        default=None,
        metavar="PATH",
        help="Run summary JSON output; defaults next to --out-json",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    groups = parse_groups(args.groups)

    with gzip.open(args.pickle, "rb") as f:
        cov = pickle.load(f)

    missing = sorted({rep for reps in groups.values() for rep in reps if rep not in cov})
    if missing:
        available = ", ".join(sorted(cov.keys()))
        raise SystemExit(
            "The following replicates are missing from the coverage pickle: "
            f"{', '.join(missing)}\nAvailable experiments: {available}"
        )

    filt_tx_dict = {
        group: filter_tx(cov, reps, min_reps=args.tx_min_reps, threshold=args.tx_threshold)
        for group, reps in groups.items()
    }
    filt_tx_set = set.intersection(*(set(v) for v in filt_tx_dict.values())) if filt_tx_dict else set()

    cov_filt = {
        exp: {tx: arr for tx, arr in tx_dict.items() if tx in filt_tx_set}
        for exp, tx_dict in cov.items()
    }

    codon_cov = {
        exp: {tx: codonize_counts_cds(arr) for tx, arr in tx_dict.items()}
        for exp, tx_dict in cov_filt.items()
    }

    stalls = {
        exp: {
            tx: call_stalls(
                arr,
                min_z=args.min_z,
                min_obs=args.min_reads,
                trim_edges=args.trim_edges,
                pseudocount=args.pseudocount,
            )
            for tx, arr in tx_dict.items()
        }
        for exp, tx_dict in codon_cov.items()
    }

    consensus = {
        group: consensus_stalls_across_reps(
            stalls,
            reps,
            min_support=args.stall_min_reps,
            tol=args.tol,
            min_sep=args.min_sep,
        )
        for group, reps in groups.items()
    }

    total_counts = {
        group: sum(len(idxs) for idxs in tx_map.values())
        for group, tx_map in consensus.items()
    }

    print(f"Number of filtered transcripts: {len(filt_tx_set)}")
    print(f"Number of total stall sites per group: {total_counts}")

    total_consensus_sites = sum(total_counts.values())
    consensus_df = consensus_to_long_df(consensus) if total_consensus_sites else empty_consensus_df()
    write_jsonl(consensus_df, args.out_json)
    logging.info(f"Saved consensus stall sites to {args.out_json}")

    rep_to_group = {rep: group for group, reps in groups.items() for rep in reps}
    per_rep_total = sum(len(stall_list) for tx_map in stalls.values() for stall_list in tx_map.values())
    replicate_df = stalls_to_long_df(stalls, rep_to_group=rep_to_group) if per_rep_total else empty_replicate_df()

    replicate_json = args.out_replicate_json
    if replicate_json is None:
        replicate_json = str(Path(args.out_json).with_name("stall_sites_by_replicate.jsonl"))
    write_jsonl(replicate_df, replicate_json)
    logging.info(f"Saved per-replicate stall sites to {replicate_json}")

    summary_json = args.summary_json
    if summary_json is None:
        summary_json = str(Path(args.out_json).with_name("stall_sites_summary.json"))

    summary = {
        "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
        "command": " ".join(sys.argv),
        "parameters": vars(args),
        "groups": groups,
        "filtered_transcripts": len(filt_tx_set),
        "filtered_transcripts_by_group": {g: len(v) for g, v in filt_tx_dict.items()},
        "total_stall_sites": total_counts,
        "output_files": {
            "consensus_stall_sites_jsonl": args.out_json,
            "per_replicate_stall_sites_jsonl": replicate_json,
            "summary_json": summary_json,
        },
    }

    Path(summary_json).parent.mkdir(parents=True, exist_ok=True)
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)
    logging.info(f"Saved run summary to {summary_json}")


if __name__ == "__main__":
    main()
