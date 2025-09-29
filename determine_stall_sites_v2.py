import argparse
import gzip
import pickle
import matplotlib.pyplot as plt
import ribopy
from ribopy import Ribo
from functions_stall_sites import filter_tx, codonize_counts_cds, call_stalls, consensus_stalls_across_reps
from functions_AA import translate_cds_nt_to_aa, windows_aa, count_matrix, background_aa_freq, pwm_position_weighted_log2, plot_logo, CODON2AA, AA_ORDER, AA_CLASS
from functions import get_sequence, get_cds_range_lookup

def main():
    parser = argparse.ArgumentParser(
        description="Detect ribosome stall sites"
    )
    parser.add_argument("--pickle", required=True, help="Path to coverage pickle.gz file")
    parser.add_argument("--ribo", required=True, help="Path to ribo file")
    parser.add_argument("--tx_threshold", type=float, default=1.0,
                        help="Minimum reads/nt (in CDS) for filtering transcripts")
    parser.add_argument("--tx_min_reps", type=int, default=2,
                        help="Minimum number of replicates passing threshold for filtering transcripts")
    parser.add_argument("--min_z", type=float, default=1.0,
                        help="Minimum z-score to pass as stall site")
    parser.add_argument("--min_reads", type=int, default=2,
                        help="Minimum number of reads to pass as stall site")
    parser.add_argument("--trim_edges", type=int, default=10,
                        help="Trim x number of codons after start and before stop")
    parser.add_argument("--pseudocount", type=float, default=0.5,
                        help="Pseudocount for stall calling")
    parser.add_argument("--stall_min_reps", type=int, default=2,
                        help="Minimum number of replicates that must support a site")
    parser.add_argument("--tol", type=int, default=0,
                        help="Tolerance window for matching sites across reps (same units as indices)")
    parser.add_argument("--min_sep", type=int, default=7,
                        help="Minimum separation between consensus sites; prefer downstream when closer than this")
    parser.add_argument("--motif", action="store_true", help="Plot motif")
    parser.add_argument("--reference", help="Reference file path")
    parser.add_argument("--flank-left", type=int, default=10, help="Motif")
    parser.add_argument("--flank-right", type=int, default=6, help="Motif")
    parser.add_argument("--psite-offset", type=int, default=0, help="Motif")
    parser.add_argument("--out-motif", default="../motif.png", help="Motif")

    args = parser.parse_args()

    # Load coverage
    with gzip.open(args.pickle, "rb") as f:
        cov = pickle.load(f)

    # Load ribo object (adjust alias to your organism as needed)
    ribo_object = Ribo(args.ribo, alias=ribopy.api.alias.apris_human_alias)

    # Rep groups
    groups = {
        "kidney": ["kidney_rep1", "kidney_rep2", "kidney_rep3"],
        "lung":   ["lung_rep1", "lung_rep2", "lung_rep3"],
        "liver":  ["liver_rep1", "liver_rep2", "liver_rep3"],
    }

    # (Optional) quick sanity check that all reps exist in coverage
    missing = [r for rs in groups.values() for r in rs if r not in cov]
    if missing:
        print("Warning: the following replicates are missing from coverage:", ", ".join(missing))

    # Filter transcripts
    filt_tx_dict = {
        group: filter_tx(cov, reps, min_reps=args.tx_min_reps, threshold=args.tx_threshold)
        for group, reps in groups.items()
    }
    # Intersection across tissues
    filt_tx_set = set.intersection(*(set(v) for v in filt_tx_dict.values())) if filt_tx_dict else set()

    # Keep only filtered transcripts in coverage
    cov_filt = {
        exp: {tx: arr for tx, arr in tx_dict.items() if tx in filt_tx_set}
        for exp, tx_dict in cov.items()
    }

    # Codonize counts
    codon_cov = {
        exp: {tx: codonize_counts_cds(arr) for tx, arr in tx_dict.items()}
        for exp, tx_dict in cov_filt.items()
    }

    # Identify stall sites per experiment
    stalls = {
        exp: {
            tx: call_stalls(
                arr,
                min_z=args.min_z,
                min_obs=args.min_reads,
                trim_edges=args.trim_edges,
                pseudocount=args.pseudocount
            )
            for tx, arr in tx_dict.items()
        }
        for exp, tx_dict in codon_cov.items()
    }

    # Consensus stalls per tissue group
    consensus = {
        group: consensus_stalls_across_reps(
            stalls,
            reps,
            min_support=args.stall_min_reps,
            tol=args.tol,
            min_sep=args.min_sep
        )
        for group, reps in groups.items()
    }

    # Print
    print(f"Number of filtered transcripts: {len(filt_tx_set)}")
    total_counts = {
        group: sum(len(idxs) for idxs in stalls.values())
        for group, stalls in consensus.items()
    }
    print(f"Number of total stall sites per group: {total_counts}")

    # Motif
    if args.motif:
        reference_file_path = args.reference
        cds_range = get_cds_range_lookup(ribo_object)
        sequence = get_sequence(ribo_object, reference_file_path, alias = ribopy.api.alias.apris_human_alias)   
        def compute_W_for_group(g):
            stalls = consensus[g]
            win = windows_aa(consensus[g], cds_range, sequence,
                        flank_left=args.flank_left, flank_right=args.flank_right, psite_offset_codons=args.psite_offset)
            counts = count_matrix(win, AA_ORDER, flank_left=args.flank_left, flank_right=args.flank_right)
            bg = background_aa_freq(consensus[g].keys(), cds_range, sequence, AA_ORDER)
            return pwm_position_weighted_log2(counts, bg, pseudocount=args.pseudocount)

        # compute all, then unify y-limits for fair visual comparison
        W_by_group = {g: compute_W_for_group(g) for g in groups.keys()}

        # per-position heights (sum over amino acids)
        ymax = max(
            W.loc[:, pos][W.loc[:, pos] > 0].sum()
            for g, W in W_by_group.items()
            for pos in W.columns
        )
        ymin = max(
            abs(W.loc[:, pos][W.loc[:, pos] < 0].sum())
            for g, W in W_by_group.items()
            for pos in W.columns
        )

        # plot side-by-side
        fig, axes = plt.subplots(1, len(groups.keys()), figsize=(5*len(groups.keys()), 5), sharey=True)
        if len(groups.keys()) == 1:
            axes = [axes]

        for ax, g in zip(axes, groups.keys()):
            plt.sca(ax)
            plot_logo(W_by_group[g],
                    title=f"{g.capitalize()}: EPA-centered AA motif",
                    aa_class=AA_CLASS)
            ax.set_ylim(-ymin, ymax)   # same scale across panels

        plt.tight_layout()
        fig.savefig(args.out_motif, dpi=600)
    

if __name__ == "__main__":
    main()