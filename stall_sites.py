{
  "timestamp": "2026-06-08T15:24:07",
  "command": "stall_sites.py --pickle ../celegans/adj_coverage.pkl.gz --ribo ../celegans/all.ribo --groups CHX:CHX_rep1,CHX_rep2,CHX_rep3;noCHX:noCHX_rep1,noCHX_rep2,noCHX_rep3",
  "parameters": {
    "pickle": "../celegans/adj_coverage.pkl.gz",
    "ribo": "../celegans/all.ribo",
    "alias": false,
    "tx_threshold": 1.0,
    "groups": "CHX:CHX_rep1,CHX_rep2,CHX_rep3;noCHX:noCHX_rep1,noCHX_rep2,noCHX_rep3",
    "tx_min_reps": 2,
    "min_z": 1.0,
    "min_reads": 2,
    "trim_edges": 10,
    "pseudocount": 0.5,
    "stall_min_reps": 2,
    "tol": 0,
    "min_sep": 7,
    "out_json": "../ribostall_results/stall_sites.jsonl",
    "motif": false,
    "reference": null,
    "flank_left": 10,
    "flank_right": 6,
    "psite_offset": 0,
    "out_png": "../ribostall_results/motif.png",
    "out_csv": "../ribostall_results/motif_csv"
  },
  "groups": {
    "CHX": [
      "CHX_rep1",
      "CHX_rep2",
      "CHX_rep3"
    ],
    "noCHX": [
      "noCHX_rep1",
      "noCHX_rep2",
      "noCHX_rep3"
    ]
  },
  "filtered_transcripts": 287,
  "total_stall_sites": {
    "CHX": 8116,
    "noCHX": 7945
  },
  "output_files": {
    "stall_sites_jsonl": "../ribostall_results/stall_sites.jsonl"
  }
}
