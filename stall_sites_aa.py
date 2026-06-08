#!/usr/bin/env python3
import argparse
import datetime
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ribopy
from ribopy import Ribo

from functions import get_cds_range_lookup, get_sequence
from functions_AA import (
    AA_CLASS,
    AA_ORDER,
    CLASS_COLORS,
    background_aa_freq,
    count_matrix,
    plot_logo,
    pwm_position_weighted_log2,
    windows_aa,
)


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(processName)s  %(message)s",
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


def load_stall_sites_jsonl(path):
    consensus = defaultdict(lambda: defaultdict(list))
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            rec = json.loads(line)
            consensus[rec["group"]][rec["transcript"]].append(int(rec["pos_codon"]))
    return {group: dict(tx_map) for group, tx_map in consensus.items()}


def finite_weight_matrix(weight_matrix):
    return weight_matrix.replace([np.inf, -np.inf], 0).fillna(0)


def tx_alias(tx):
    tx = str(tx)
    parts = tx.split("|")
    return parts[4] if len(parts) > 4 else tx


def motif_usable_key(tx, cds_range, sequence):
    candidates = [str(tx)]
    alias = tx_alias(tx)
    if alias != tx:
        candidates.append(alias)

    for candidate in candidates:
        if candidate not in sequence:
            continue
        if candidate in cds_range:
            return candidate
        candidate_alias = tx_alias(candidate)
        if candidate_alias != candidate and candidate_alias in cds_range:
            return candidate
    return None


def filter_tx_map_for_motif(tx_map, cds_range, sequence):
    filtered = defaultdict(list)
    skipped = []
    for tx, positions in tx_map.items():
        key = motif_usable_key(tx, cds_range, sequence)
        if key is None:
            skipped.append(tx)
            continue
        filtered[key].extend(int(p) for p in positions)
    return dict(filtered), skipped


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate amino-acid motif logos from stall-site JSONL output."
    )
    parser.add_argument(
        "--stall-sites",
        required=True,
        metavar="PATH",
        help="Consensus stall-site JSONL from stall_sites.py",
    )
    parser.add_argument("--ribo", required=True, metavar="PATH", help="Path to .ribo file")
    parser.add_argument("--reference", required=True, metavar="PATH", help="Reference FASTA file")
    parser.add_argument(
        "--alias",
        action="store_true",
        help="Use APPRIS-style human/mouse aliasing when opening the .ribo file",
    )
    add_option(
        parser,
        "--flank-left",
        legacy="--flank_left",
        type=int,
        default=10,
        metavar="CODONS",
        help="Number of amino acids to include left/upstream of the stall codon",
    )
    add_option(
        parser,
        "--flank-right",
        legacy="--flank_right",
        type=int,
        default=6,
        metavar="CODONS",
        help="Number of amino acids to include right/downstream of the stall codon",
    )
    add_option(
        parser,
        "--psite-offset",
        legacy="--psite_offset",
        type=int,
        default=0,
        metavar="CODONS",
        help="Codon offset to add before extracting motif windows",
    )
    parser.add_argument(
        "--pseudocount",
        type=float,
        default=0.5,
        metavar="VALUE",
        help="Pseudocount for amino-acid enrichment calculation",
    )
    parser.add_argument(
        "--out-dir",
        default="../ribostall_results/aa_motif",
        metavar="PATH",
        help="Output directory for motif PNG, matrices, and summary JSON",
    )
    parser.add_argument(
        "--out-png",
        default=None,
        metavar="PATH",
        help="Motif logo PNG path; defaults to OUT_DIR/motif_logo.png",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.alias:
        ribo_object = Ribo(args.ribo, alias=ribopy.api.alias.apris_human_alias)
    else:
        ribo_object = Ribo(args.ribo)

    consensus = load_stall_sites_jsonl(args.stall_sites)
    if not consensus:
        raise SystemExit(f"No stall sites found in {args.stall_sites}")

    cds_range = get_cds_range_lookup(ribo_object)
    sequence = get_sequence(ribo_object, args.reference, alias=args.alias)

    motif_pwm_dir = out_dir / "motif_pwm"
    motif_counts_dir = out_dir / "motif_counts"
    motif_pwm_dir.mkdir(exist_ok=True)
    motif_counts_dir.mkdir(exist_ok=True)

    counts_by_group = {}
    weights_by_group = {}
    windows_by_group = {}
    skipped_groups = []
    skipped_transcripts_by_group = {}

    for group, tx_map in consensus.items():
        n_sites = sum(len(v) for v in tx_map.values())
        if n_sites == 0:
            skipped_groups.append({"group": group, "reason": "no stall sites"})
            continue

        motif_tx_map, skipped_tx = filter_tx_map_for_motif(tx_map, cds_range, sequence)
        if skipped_tx:
            skipped_transcripts_by_group[group] = skipped_tx
            logging.warning(
                f"[{group}] skipped {len(skipped_tx)} transcripts missing from CDS/sequence lookup"
            )

        win = windows_aa(
            motif_tx_map,
            cds_range,
            sequence,
            flank_left=args.flank_left,
            flank_right=args.flank_right,
            psite_offset_codons=args.psite_offset,
        )
        if not win:
            skipped_groups.append({"group": group, "reason": "no usable amino-acid windows"})
            continue

        counts = count_matrix(
            win,
            AA_ORDER,
            flank_left=args.flank_left,
            flank_right=args.flank_right,
        )
        bg = background_aa_freq(motif_tx_map.keys(), cds_range, sequence, AA_ORDER)
        weights = pwm_position_weighted_log2(counts, bg, pseudocount=args.pseudocount)
        weights = finite_weight_matrix(weights)

        counts_by_group[group] = counts
        weights_by_group[group] = weights
        windows_by_group[group] = len(win)

        counts.to_csv(motif_counts_dir / f"{group}_raw_counts.csv")
        weights.to_csv(motif_pwm_dir / f"{group}_pwm_log2_enrichment.csv")

    if not weights_by_group:
        summary_path = out_dir / "motif_summary.json"
        with open(summary_path, "w") as f:
            json.dump(
                {
                    "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
                    "command": " ".join(sys.argv),
                    "parameters": vars(args),
                    "groups": list(consensus.keys()),
                    "skipped_groups": skipped_groups,
                    "skipped_transcripts_by_group": skipped_transcripts_by_group,
                    "output_files": {},
                },
                f,
                indent=2,
            )
        raise SystemExit(f"No motif could be made. See {summary_path}")

    ymax = max(
        weights.loc[:, pos][weights.loc[:, pos] > 0].sum()
        for weights in weights_by_group.values()
        for pos in weights.columns
    )
    ymin = max(
        abs(weights.loc[:, pos][weights.loc[:, pos] < 0].sum())
        for weights in weights_by_group.values()
        for pos in weights.columns
    )
    ylim = max(float(ymax), float(ymin), 1e-6)

    groups_to_plot = list(weights_by_group.keys())
    fig, axes = plt.subplots(
        1,
        len(groups_to_plot),
        figsize=(5 * len(groups_to_plot), 5),
        sharey=True,
    )
    if len(groups_to_plot) == 1:
        axes = [axes]

    for ax, group in zip(axes, groups_to_plot):
        plot_logo(
            weights_by_group[group],
            title=group,
            aa_class=AA_CLASS,
            ax=ax,
        )
        ax.set_ylim(-ylim, ylim)

    for ax in axes[1:]:
        ax.set_ylabel("")
        ax.tick_params(axis="y", left=False, labelleft=False)
        ax.spines["left"].set_visible(False)

    patches = [mpatches.Patch(color=color, label=cls) for cls, color in CLASS_COLORS.items()]
    fig.legend(handles=patches, loc="lower center", ncol=len(patches), bbox_to_anchor=(0.5, 0.04))
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.17)

    out_png = Path(args.out_png) if args.out_png else out_dir / "motif_logo.png"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=600)
    plt.close(fig)
    logging.info(f"Saved motif image to {out_png}")

    summary_path = out_dir / "motif_summary.json"
    summary = {
        "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
        "command": " ".join(sys.argv),
        "parameters": vars(args),
        "groups": list(consensus.keys()),
        "plotted_groups": groups_to_plot,
        "skipped_groups": skipped_groups,
        "skipped_transcripts_by_group": skipped_transcripts_by_group,
        "stall_sites_by_group": {
            group: sum(len(v) for v in tx_map.values())
            for group, tx_map in consensus.items()
        },
        "usable_windows_by_group": windows_by_group,
        "output_files": {
            "motif_logo_png": str(out_png),
            "motif_pwm_dir": str(motif_pwm_dir),
            "motif_counts_dir": str(motif_counts_dir),
            "summary_json": str(summary_path),
        },
    }
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    logging.info(f"Saved motif summary to {summary_path}")


if __name__ == "__main__":
    main()
