import pandas as pd
import numpy as np
from collections import Counter, defaultdict

# --- mappings ---
CODON2AA = {
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R",
    "AAT":"N","AAC":"N","GAT":"D","GAC":"D","TGT":"C","TGC":"C",
    "GAA":"E","GAG":"E","CAA":"Q","CAG":"Q","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
    "CAT":"H","CAC":"H","ATT":"I","ATC":"I","ATA":"I",
    "TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "AAA":"K","AAG":"K","ATG":"M","TTT":"F","TTC":"F",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","TGG":"W",
    "TAT":"Y","TAC":"Y","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TAA":"*","TAG":"*","TGA":"*"
}

AA_CLASS = {
    # Acidic
    "D":"acidic","E":"acidic",
    # Basic
    "K":"basic","R":"basic","H":"basic",
    # Polar (uncharged)
    "S":"polar","T":"polar","N":"polar","Q":"polar","Y":"polar","C":"polar",
    # Hydrophobic (nonpolar)
    "A":"hydrophobic","V":"hydrophobic","L":"hydrophobic","I":"hydrophobic",
    "M":"hydrophobic","F":"hydrophobic","W":"hydrophobic","P":"hydrophobic","G":"hydrophobic",
    # Neutral / stop (drop later)
    "*":"neutral"
}

AA_ORDER = [a for a in "ACDEFGHIKLMNPQRSTVWY"]  # no stop

def translate_cds_nt_to_aa(cds_nt: str) -> str:
    cds_nt = cds_nt.upper().replace("U","T")
    aas = []
    for i in range(0, len(cds_nt) - (len(cds_nt) % 3), 3):
        cod = cds_nt[i:i+3]
        aas.append(CODON2AA.get(cod, "X"))  # X for unknown/ambiguous
    return "".join(aas)

def windows_aa(consensus_group: dict, cds_range: dict, sequence: dict,
               flank_left=10, flank_right=10, psite_offset_codons=0):
    """
    Build AA windows of size (flank_left + 1 + flank_right) around each stall
    aligned at the P-site index (i + psite_offset_codons).
    Returns: list of AA lists, each length = window_len.
    """
    W = flank_left + 1 + flank_right
    win_list = []
    for tx, idx_list in consensus_group.items():
        # harmonize key for cds_range
        key = tx if tx in cds_range else tx.split("|")[4]
        start, stop = cds_range[key]
        cds_nt = sequence[tx][start:stop]
        aa_seq = translate_cds_nt_to_aa(cds_nt)
        Lcod = len(aa_seq)
        for i in idx_list:
            center = i + psite_offset_codons  # codon index for P-site
            lo = center - flank_left
            hi = center + flank_right + 1
            if lo < 0 or hi > Lcod:
                continue
            win = list(aa_seq[lo:hi])
            if "*" in win:    # optional: drop windows containing stop
                continue
            if "X" in win:    # drop ambiguous
                continue
            win_list.append(win)
    return win_list  # list of lists

def count_matrix(win_list, aa_order=AA_ORDER, flank_left=10, flank_right=10):
    """rows=AA, cols=relative positions (-flank_left..+flank_right)."""
    if not win_list:
        return pd.DataFrame(0, index=aa_order, columns=[])
    W = len(win_list[0])
    assert W == flank_left + 1 + flank_right, "Window length != flank sizes"
    cols = list(range(-flank_left, flank_right + 1))
    counts = pd.DataFrame(0, index=aa_order, columns=cols, dtype=int)
    for win in win_list:
        for j, aa in enumerate(win):
            if aa in counts.index:
                counts.loc[aa, cols[j]] += 1
    return counts


def background_aa_freq(transcripts: dict, cds_range: dict, sequence: dict,
                       aa_order=AA_ORDER):
    """
    Background AA frequency across all callable CDS codons of the same transcripts
    (uniform over all positions). Returns Series summing to 1.
    """
    bg_counts = Counter()
    for tx in transcripts:
        key = tx if tx in cds_range else tx.split("|")[4]
        start, stop = cds_range[key]
        aa_seq = translate_cds_nt_to_aa(sequence[tx][start:stop])
        for aa in aa_seq:
            if aa in aa_order:
                bg_counts[aa] += 1
    bg = pd.Series({aa: bg_counts.get(aa, 0) for aa in aa_order}, dtype=float)
    bg = (bg + 1e-6) / (bg.sum() + 1e-6 * len(bg))  # tiny pseudocount
    return bg

def pwm_position_weighted_log2(counts_pos, bg_freq, pseudocount=0.5):
    """Ignore columns with zero total counts to avoid pseudocount artifacts."""
    nonempty = counts_pos.sum(axis=0) > 0
    counts = counts_pos.loc[:, nonempty].copy()
    probs = (counts + pseudocount).div((counts + pseudocount).sum(axis=0), axis=1)
    lo = np.log2(probs.div(bg_freq, axis=0))
    W = probs * lo
    return W  # only non-empty positions



def plot_logo(weight_mat, title="", aa_class=None, ax=None, ylim=None):
    import logomaker, matplotlib.pyplot as plt
    df = weight_mat.T.copy()  # index=positions, columns=AAs

    # Color map by AA class
    CLASS_COLORS = {
        "acidic":"#D62728",      # red
        "basic":"#1F77B4",       # blue
        "hydrophobic":"#2CA02C", # green
        "polar":"#9467BD",       # purple
        "neutral":"#7F7F7F",     # gray
    }
    if aa_class is None:
        aa_class = {aa:"neutral" for aa in df.columns}
    color_scheme = {aa: CLASS_COLORS.get(aa_class.get(aa,"neutral"), "#000000")
                    for aa in df.columns}

    if ax is None:
        ax = plt.gca()

    # Force all letters (positive & negative) to use class-based colors
    logomaker.Logo(
        df,
        color_scheme=color_scheme,
        shade_below=False,   # don't gray out negatives
        fade_below=False,
        flip_below=True,     # flip letters for negative values
        vpad=0.02,
        baseline_width=0,
        ax=ax
    )

    ax.set_title(title, pad=6)
    ax.set_xlabel("Relative position (codons; 0 = P-site)")
    ax.set_ylabel("p · log2(p / bg)")

    # set y-limit if specified
    if ylim is not None:
        ax.set_ylim(-ylim, ylim)   # symmetric so depletion shows clearly

    return ax