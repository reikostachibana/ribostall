import numpy as np
import pandas as pd
import bisect
from typing import Dict, List, Iterable, Any

def filter_tx(cov_by_exp: dict, reps: list[str], min_reps: int = 2, threshold: float = 1.0):
    """
    Keep transcripts where at least `min_reps` of the given replicates
    have mean coverage per-nt > threshold.
    """
    common = set.intersection(*(set(cov_by_exp[r].keys()) for r in reps))
    keep = []
    for tx in common:
        n_pass = sum(np.asarray(cov_by_exp[r][tx], float).mean() > threshold for r in reps)
        if n_pass >= min_reps:
            keep.append(tx)
    return keep

def codonize_counts_cds(x_nt: np.ndarray, frame: int = 0):
    """
    Codon-sum a CDS-only nt-coverage vector.
    x_nt: 1D array of per-nt counts covering ONLY the CDS (5'UTR and 3'UTR removed).
    frame: which nt within the codon the indexing should start from (0, 1, or 2).
           Use 0 if x_nt[0] corresponds to the first nt of the start codon, etc.

    Returns:
      x_codon: float array of length = number of full codons
      map_codon_to_nt: list of (lo, hi) nt index slices in the original x_nt
    """
    assert x_nt.ndim == 1, "x_nt must be 1D"
    frame = int(frame) % 3
    # Align start by frame and trim tail to multiple of 3
    start = frame
    usable_len = ((len(x_nt) - start) // 3) * 3
    if usable_len <= 0:
        return np.zeros(0, dtype=float), []
    stop = start + usable_len
    cds_slice = x_nt[start:stop]                # length is multiple of 3
    x3 = cds_slice.reshape(-1, 3)               # (codons, 3)
    x_codon = x3.sum(axis=1).astype(float)
    return x_codon

def global_z_log(x: np.ndarray, pseudocount: float = 0.5) -> np.ndarray:
    """Transcript-wide z on log2(x+pc)."""
    x = np.asarray(x, dtype=float)
    v = np.log2(x + pseudocount)
    sd = v.std()
    if sd == 0:
        return np.zeros_like(v)
    return (v - v.mean()) / (sd + 1e-12)

def call_stalls(
    x_codon: np.ndarray,
    min_z: float = 1.0,
    min_obs: int = 2,
    trim_edges: int = 10,
    pseudocount: float = 0.5,
):
    """
    Keep any codon with global z >= min_z and obs >= min_obs (no local filters).
    Returns list of dicts: {index, obs, z}.
    """
    x = np.asarray(x_codon, dtype=float)
    n = x.size
    if n == 0:
        return []
    z = global_z_log(x, pseudocount=pseudocount)
    lo, hi = trim_edges, n - trim_edges
    if hi <= lo:
        return []
    cand = np.flatnonzero((z >= min_z) & (x >= min_obs) & (np.arange(n) >= lo) & (np.arange(n) < hi))
    return [dict(index=int(i), obs=float(x[i]), z=float(z[i])) for i in cand]

def _indices_from_stalls(stalls):
    return sorted({int(d["index"]) for d in stalls})

import numpy as np
import bisect
from typing import Dict, List

def consensus_stalls_across_reps(
    stalls_by_exp: Dict[str, dict],
    reps: List[str],
    *,
    min_support: int = 2,
    tol: int = 0,
    min_sep: int = 7
):
    """
    Compute consensus stall indices per transcript across replicate experiments,
    allowing a tolerance window and preferring the downstream site when two
    candidates are closer than `min_sep`.
    """

    def _indices_from_stalls(stalls_for_tx: Any) -> List[int]:
        """
        Accepts multiple shapes:
          - list[int]
          - list[tuple]            -> take first element as index
          - list[dict]             -> try common keys for index
          - dict with "indices"    -> use that list
        Returns sorted unique indices.
        """
        if not stalls_for_tx:
            return []

        # Case: dict with "indices"
        if isinstance(stalls_for_tx, dict) and "indices" in stalls_for_tx:
            idxs = stalls_for_tx["indices"]
            return sorted(set(int(i) for i in idxs))

        # Case: iterable (list/tuple) of elements
        if isinstance(stalls_for_tx, Iterable) and not isinstance(stalls_for_tx, (str, bytes, dict)):
            first = next(iter(stalls_for_tx), None)
            if first is None:
                return []

            # list of dicts
            if isinstance(first, dict):
                # try common field names for the codon index
                candidate_keys = ["idx", "index", "pos", "position", "codon", "codon_idx", "codon_index"]
                key = next((k for k in candidate_keys if k in first), None)
                if key is None:
                    raise TypeError(
                        f"Don't know which key holds the stall index in dicts. "
                        f"Expected one of {candidate_keys}, got keys: {list(first.keys())}"
                    )
                idxs = [int(d[key]) for d in stalls_for_tx if d is not None]
                return sorted(set(idxs))

            # list of tuples -> take first element
            if isinstance(first, (list, tuple)):
                idxs = [int(x[0]) for x in stalls_for_tx]
                return sorted(set(idxs))

            # list of ints
            if isinstance(first, (int, np.integer)):
                return sorted(set(int(x) for x in stalls_for_tx))

            # list of strings (just in case)
            if isinstance(first, str):
                return sorted(set(int(x) for x in stalls_for_tx))

        # Fallback: single int
        if isinstance(stalls_for_tx, (int, np.integer)):
            return [int(stalls_for_tx)]

        raise TypeError(f"Unsupported stalls_for_tx shape/type: {type(stalls_for_tx)}")

    # transcripts present in all requested replicates
    common_txs = set.intersection(*(set(stalls_by_exp[r].keys()) for r in reps))
    out: Dict[str, List[int]] = {}

    for tx in common_txs:
        per_rep_idx = {r: _indices_from_stalls(stalls_by_exp[r][tx]) for r in reps}

        # union of all candidate positions
        candidates = sorted(set().union(*per_rep_idx.values()))
        consensus: List[int] = []

        for c in candidates:
            # count support within ±tol, and collect the actual hits
            support = 0
            supporting_hits: List[int] = []
            for r in reps:
                arr = per_rep_idx[r]
                j = bisect.bisect_left(arr, c - tol)
                hit = None
                while j < len(arr) and arr[j] <= c + tol:
                    hit = arr[j]
                    break
                if hit is not None:
                    support += 1
                    supporting_hits.append(hit)

            if support >= min_support:
                # representative index across reps (merge with median if tol>0)
                rep_idx = int(np.median(supporting_hits)) if (tol > 0 and supporting_hits) else int(c)

                if consensus:
                    # downstream preference: if too close to last kept site, replace it
                    if rep_idx - consensus[-1] < min_sep:
                        consensus[-1] = rep_idx
                    else:
                        consensus.append(rep_idx)
                else:
                    consensus.append(rep_idx)

        out[tx] = sorted(consensus)

    return out

def parse_key(k: str):
    s = str(k)
    parts = s.split("|")
    tx_id  = parts[0] if len(parts) > 0 else None
    gene   = parts[5] if len(parts) > 5 else None  # based on your key pattern
    return s, tx_id, gene

def consensus_to_long_df(consensus: dict) -> pd.DataFrame:
    rows = []
    for grp, tx_map in consensus.items():
        for tx_key, positions in tx_map.items():
            tx_str, tx_id, gene = parse_key(tx_key)
            for p in positions:
                rows.append({
                    "group": grp,
                    "transcript": tx_str,
                    "tx_id": tx_id,
                    "gene": gene,
                    "pos_codon": int(p),
                })
    return pd.DataFrame(rows).sort_values(["group","gene","tx_id","pos_codon"])