#!/usr/bin/env python3
import argparse
import gzip
import logging
import multiprocessing as mp
import pickle
import sys
import time
from typing import Dict, Tuple, Iterable, List

import numpy as np
import ribopy
from ribopy import Ribo

# --- import your helpers ---
from functions import get_cds_range_lookup, get_psite_offset

# =========================
# Logging
# =========================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(processName)s  %(message)s",
)

# =========================
# Worker globals
# =========================
_RIBO = None
_CDS_RANGE: Dict[str, Tuple[int, int]] = {}
_TRANSCRIPTS: List[str] = []

def _safe_window(arr: np.ndarray, lo: int, hi: int) -> np.ndarray:
    """Return arr[lo:hi] with zero-padding if the window runs off the array."""
    # lo/hi are signed Python ints
    n = hi - lo
    if n <= 0:
        return np.zeros(0, dtype=arr.dtype)
    left = max(0, -lo)
    right = max(0, hi - int(len(arr)))
    lo_c, hi_c = max(0, lo), min(int(len(arr)), hi)
    core = arr[lo_c:hi_c]
    if left or right:
        return np.concatenate(
            (np.zeros(left, dtype=arr.dtype), core, np.zeros(right, dtype=arr.dtype))
        )
    return core

def _init_worker(ribo_path: str, use_alias: bool,
                 cds_range: Dict[str, Tuple[int, int]],
                 transcripts: List[str]) -> None:
    """Open the .ribo once per worker and stash shared metadata."""
    global _RIBO, _CDS_RANGE, _TRANSCRIPTS
    if use_alias:
        _RIBO = Ribo(ribo_path, alias=ribopy.api.alias.apris_human_alias)
    else:
        _RIBO = Ribo(ribo_path)
    # cast CDS ranges to signed ints defensively
    _CDS_RANGE = {t: (int(s), int(e)) for t, (s, e) in cds_range.items()}
    _TRANSCRIPTS = list(transcripts)

def _preallocate_output(transcripts: Iterable[str]) -> Dict[str, np.ndarray]:
    """Pre-allocate a zero array per transcript sized to its CDS window."""
    out: Dict[str, np.ndarray] = {}
    for t in transcripts:
        start, stop = _CDS_RANGE[t]
        win = max(0, int(stop) - int(start))
        out[t] = np.zeros(win, dtype=np.int64)
    return out

def _add_length_into_out(exp: str, L: int, ps: int,
                         out: Dict[str, np.ndarray],
                         batch: Iterable[str] = None) -> None:
    """Read coverage for a single length L once, then accumulate into 'out'."""
    cov_all = _RIBO.get_coverage(experiment=exp, range_lower=int(L), range_upper=int(L))
    to_iter = _TRANSCRIPTS if batch is None else batch
    ps_i = int(ps)
    for t in to_iter:
        raw = cov_all.get(t)
        if raw is None:
            continue
        start, stop = _CDS_RANGE[t]
        start_i, stop_i = int(start), int(stop)
        if stop_i <= start_i:
            continue
        lo = start_i - ps_i
        hi = stop_i - ps_i
        seg = _safe_window(raw, lo, hi)
        # Defensive: enforce exact length match
        if seg.shape[0] != out[t].shape[0]:
            fix = np.zeros_like(out[t])
            m = min(fix.shape[0], seg.shape[0])
            fix[:m] = seg[:m]
            seg = fix
        # In-place add
        out[t] += seg.astype(np.int64, copy=False)

def _process_experiment(exp: str, min_len: int, max_len: int,
                        offsets: Dict[int, int],
                        batch_size: int = 0) -> Tuple[str, Dict[str, np.ndarray]]:
    """
    Compute CDS-aligned, length-summed coverage for an experiment by:
      - loading each read length once (bulk per-length I/O),
      - applying that length's P-site offset,
      - accumulating into per-transcript arrays.
    Optional batching of transcripts to bound memory.
    """
    t0 = time.time()
    logging.info(f"[{exp}] start: lengths {min_len}..{max_len}")

    out = _preallocate_output(_TRANSCRIPTS)
    lengths = list(range(int(min_len), int(max_len) + 1))

    # Process lengths; one HDF5 call per L
    if batch_size and batch_size > 0:
        for L in lengths:
            ps = int(offsets.get(L, 0))
            for i in range(0, len(_TRANSCRIPTS), batch_size):
                batch = _TRANSCRIPTS[i:i + batch_size]
                _add_length_into_out(exp, L, ps, out, batch=batch)
    else:
        for L in lengths:
            ps = int(offsets.get(L, 0))
            _add_length_into_out(exp, L, ps, out, batch=None)

    dt = time.time() - t0
    logging.info(f"[{exp}] done in {dt:.2f}s")
    return exp, out

def parse_args():
    p = argparse.ArgumentParser(
        description="CDS-aligned coverage via bulk per-length I/O, parallelized by experiment."
    )
    p.add_argument("--ribo", required=True, help="Path to .ribo file")
    p.add_argument("--min-len", type=int, required=True, help="Minimum read length (inclusive)")
    p.add_argument("--max-len", type=int, required=True, help="Maximum read length (inclusive)")
    p.add_argument("--site-type", help="Site type for p-site offset")
    p.add_argument("--alias", action="store_true",
                   help="Use apris_human_alias (set if your .ribo uses mouse/human aliasing)")
    p.add_argument("--procs", type=int, default=1, help="Number of parallel worker processes (experiments run in parallel)")
    p.add_argument("--batch-size", type=int, default=0,
                   help="Optional transcript batch size to bound memory (0 = process all transcripts at once per length)")
    p.add_argument("--out", default="coverage_bulk_perlen_perexp.pkl.gz",
                   help="Output pickle.gz path")
    return p.parse_args()

def main():
    args = parse_args()

    # Single handle to gather metadata and compute offsets up front
    ribo0 = Ribo(args.ribo, alias=ribopy.api.alias.apris_human_alias) if args.alias else Ribo(args.ribo)
    experiments = list(ribo0.experiments)
    transcripts = list(ribo0.transcript_names)

    logging.info(f"Experiments: {experiments}")
    logging.info(f"Transcripts: {len(transcripts)} total")
    logging.info(f"Lengths: {args.min_len}..{args.max_len}")
    logging.info(f"Processes: {args.procs}, batch_size: {args.batch_size}")

    # CDS ranges (dict: transcript -> (start, stop)), cast to signed ints
    cds_range = get_cds_range_lookup(ribo0)
    cds_range = {t: (int(s), int(e)) for t, (s, e) in cds_range.items()}

    # Precompute per-experiment offsets (dict of dict: exp -> {L -> offset})
    exp_offsets: Dict[str, Dict[int, int]] = {}
    for exp in experiments:
        od = get_psite_offset(ribo0, exp, args.min_len, args.max_len, args.site_type)  # user-provided
        # cast keys/values to plain ints to avoid uint overflows later
        exp_offsets[exp] = {int(L): int(o) for L, o in od.items()}

    all_coverage_dict: Dict[str, Dict[str, np.ndarray]] = {}

    if len(experiments) == 0:
        logging.warning("No experiments found in the .ribo file.")
        with gzip.open(args.out, "wb") as f:
            pickle.dump(all_coverage_dict, f)
        return 0

    # Parallelize by experiment
    with mp.Pool(
        processes=int(args.procs),
        initializer=_init_worker,
        initargs=(args.ribo, bool(args.alias), cds_range, transcripts),
    ) as pool:
        tasks = [
            (exp, int(args.min_len), int(args.max_len), exp_offsets[exp], int(args.batch_size))
            for exp in experiments
        ]
        # starmap keeps order; imap_unordered also fine if you want streaming
        for exp, out_dict in pool.starmap(_process_experiment, tasks):
            all_coverage_dict[exp] = out_dict

    # Save
    with gzip.open(args.out, "wb") as f:
        pickle.dump(all_coverage_dict, f)
    logging.info(f"Saved coverage to {args.out}")

    return 0

if __name__ == "__main__":
    sys.exit(main())