import ribopy
from ribopy import Ribo
from Fasta import FastaFile
from ribopy.core.get_gadgets import get_region_boundaries, get_reference_names
import pandas as pd

def get_sequence(ribo_object, reference_file_path, alias=False):
    """
    Retrieves transcript sequences from a reference FASTA file.
    """
    fasta = FastaFile(reference_file_path)

    sequence_dict = {}

    for e in fasta:
        header = str(e.header)
        seq = e.sequence

        sequence_dict[header] = seq

        parts = header.split("|")
        if len(parts) > 4:
            sequence_dict[parts[4]] = seq

    return sequence_dict

def get_cds_range_lookup(ribo_object):
    """
    Create a dict of gene to CDS ranges.
    
    Parameters:
        ribo_object (Ribo): The Ribo object containing ribosome profiling data.

    Returns:
        dict: A dictionary mapping transcript identifiers to the start and stop positions of CDS.
    """
    names = get_reference_names(ribo_object._handle)
    if ribo_object.alias is not None:
        names = map(ribo_object.alias.get_alias, names)
    
    boundaries = get_region_boundaries(ribo_object._handle)
    cds_ranges = [boundary[1] for boundary in boundaries]
    boundary_lookup = dict(zip(list(names), cds_ranges))

    return boundary_lookup

def apris_human_alias(x):
    return x.split("|")[4]

from typing import Dict, Tuple, Literal
import pandas as pd

CODON_NT = 3  # one codon in nt

def get_offset(
    ribo_object,
    exp: str,
    mmin: int,
    mmax: int,
    landmark: Literal["start","stop"],
    search_window: Tuple[int,int] | None = None,
    return_site: Literal["P","A"] = "P",
) -> Dict[int,int]:
    """
    Calibrate offsets (from 5' end) per read length using metagene profiles.

    landmark: which reference to align metagene to ("start" or "stop").
              Using "start" typically yields P-site offsets,
              using "stop" typically yields A-site offsets.
    return_site: which site you want the final offsets for ("P" or "A").
                 We'll convert between A and P by ±3 nt as needed.
    """
    # sensible default peak windows
    if search_window is None:
        search_window = (-25, -10) if landmark == "start" else (-60, -30)

    mg = ribo_object.get_metagene(
        site_type=landmark,
        experiments=exp,
        range_lower=mmin,
        range_upper=mmax,
        sum_lengths=False,
        sum_references=True,
    )

    def _all_intlike(vals):
        try:
            pd.Series(vals).astype(int)
            return True
        except Exception:
            return False

    # Ensure columns are integer positions
    if isinstance(mg.columns, pd.MultiIndex):
        pos_level = None
        for lvl in range(mg.columns.nlevels):
            vals = mg.columns.get_level_values(lvl)
            if _all_intlike(vals):
                pos_level = lvl
                break
        if pos_level is None:
            raise ValueError("Couldn't identify a position level in metagene columns.")
        pos_vals = mg.columns.get_level_values(pos_level).astype(int)
        mg = mg.groupby(pos_vals, axis=1).sum()
    else:
        if not _all_intlike(mg.columns):
            mg = mg.T
        mg.columns = mg.columns.astype(int)

    mg = mg.loc[:, sorted(mg.columns)]

    # Group by read length
    if isinstance(mg.index, pd.MultiIndex):
        length_level = mg.index.nlevels - 1
        if not _all_intlike(mg.index.get_level_values(length_level)):
            for lvl in range(mg.index.nlevels):
                if _all_intlike(mg.index.get_level_values(lvl)):
                    length_level = lvl
                    break
        grouped = mg.groupby(level=length_level)
    else:
        if not _all_intlike(mg.index):
            raise ValueError("Metagene rows are not indexed by read length.")
        mg.index = mg.index.astype(int)
        grouped = ((int(L), mg.loc[[L]]) for L in mg.index.unique())

    lo, hi = search_window
    offsets = {}
    for L, block in grouped:
        prof = block.sum(axis=0)  # 1D profile across positions
        windowed = prof.loc[(prof.index >= lo) & (prof.index <= hi)]
        if windowed.empty or (windowed.max() == 0):
            continue
        peak_pos = int(windowed.idxmax())          # 5' end peak position (relative to landmark)
        base_offset = -peak_pos + 1                # offset from 5' end to the *landmark-defined* site
        offsets[int(L)] = base_offset

    # Convert to the requested site
    inferred_site = "P" if landmark == "start" else "A"
    if return_site != inferred_site:
        # A -> P: subtract 3;  P -> A: add 3
        delta = -CODON_NT if (inferred_site == "A" and return_site == "P") else CODON_NT
        for L in list(offsets.keys()):
            offsets[L] = int(offsets[L] + delta)

    return offsets