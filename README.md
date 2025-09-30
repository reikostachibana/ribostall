Clone GitHub repository:
```
git clone https://github.com/reikostachibana/ribostall
```

## Create gzipped pickle file of coverage data with applied P-site offset.

`adj_coverage.py` creates a gzipped pickle file of coverage data with applied P-site offset.

Input:
* Ribo file path - REQUIRED
* Minimum read length: Inclusive
* Maximum read length: Inclusive
* Site type: "start" or "stop" for P-site offset
* Alias: DO NOT USE DUE TO UNSOLVED ERROR
* Processes: Number of parallel worker processes (experiments run in parallel)
* Batch-size: USE 0
* Output file

Example input:
```
python adj_coverage.py --ribo "../bxc/bxc_disome.ribo" --site-type "stop" \
--min-len 57 --max-len 65 --procs 9 --out "../bxc/cov_di.pkl.gz"
```
Example output:
```
2025-09-29 14:10:03,451  INFO  MainProcess  Experiments: ['kidney_rep1', 'kidney_rep2', 'kidney_rep3', 'liver_rep1', 'liver_rep2', 'liver_rep3', 'lung_rep1', 'lung_rep2', 'lung_rep3']
2025-09-29 14:10:03,451  INFO  MainProcess  Transcripts: 21568 total
2025-09-29 14:10:03,451  INFO  MainProcess  Lengths: 57..65
2025-09-29 14:10:03,451  INFO  MainProcess  Processes: 9, batch_size: 0
2025-09-29 14:10:20,026  INFO  SpawnPoolWorker-1  [kidney_rep1] start: lengths 57..65
2025-09-29 14:10:20,026  INFO  SpawnPoolWorker-2  [kidney_rep2] start: lengths 57..65
2025-09-29 14:10:20,026  INFO  SpawnPoolWorker-3  [kidney_rep3] start: lengths 57..65
2025-09-29 14:10:20,027  INFO  SpawnPoolWorker-6  [liver_rep3] start: lengths 57..65
2025-09-29 14:10:20,026  INFO  SpawnPoolWorker-5  [liver_rep2] start: lengths 57..65
2025-09-29 14:10:20,027  INFO  SpawnPoolWorker-7  [lung_rep1] start: lengths 57..65
2025-09-29 14:10:20,026  INFO  SpawnPoolWorker-4  [liver_rep1] start: lengths 57..65
2025-09-29 14:10:20,027  INFO  SpawnPoolWorker-8  [lung_rep2] start: lengths 57..65
2025-09-29 14:10:20,114  INFO  SpawnPoolWorker-9  [lung_rep3] start: lengths 57..65
2025-09-29 14:10:29,574  INFO  SpawnPoolWorker-4  [liver_rep1] done in 9.55s
2025-09-29 14:10:29,647  INFO  SpawnPoolWorker-3  [kidney_rep3] done in 9.62s
2025-09-29 14:10:29,652  INFO  SpawnPoolWorker-1  [kidney_rep1] done in 9.63s
2025-09-29 14:10:29,664  INFO  SpawnPoolWorker-2  [kidney_rep2] done in 9.64s
2025-09-29 14:10:30,774  INFO  SpawnPoolWorker-8  [lung_rep2] done in 10.75s
2025-09-29 14:10:30,780  INFO  SpawnPoolWorker-6  [liver_rep3] done in 10.75s
2025-09-29 14:10:30,786  INFO  SpawnPoolWorker-7  [lung_rep1] done in 10.76s
2025-09-29 14:10:30,802  INFO  SpawnPoolWorker-5  [liver_rep2] done in 10.77s
2025-09-29 14:10:30,850  INFO  SpawnPoolWorker-9  [lung_rep3] done in 10.74s
2025-09-29 14:11:57,981  INFO  MainProcess  Saved coverage to ../bxc/cov_di.pkl.gz
```

## Get stall sites

Stall sites must pass
* z-score threshold
* If multiple stall sites within ``--min_sep``, then the most downstream is taken
* Minimum replicates ``--stall_min_reps``

Check filtering:
```
python stall_sites.py --pickle "../bxc/cov_di.pkl.gz" --ribo "../bxc/bxc_disome.ribo" \
--groups "kidney:kidney_rep1,kidney_rep2,kidney_rep3;liver:liver_rep1,liver_rep2,liver_rep3;lung:lung_rep1,lung_rep2,lung_rep3" \
--tx_threshold 0.3 --min_z 1.0
```

Motif analysis:
```
python stall_sites.py --pickle "../bxc/cov_di.pkl.gz" --ribo "../bxc/bxc_disome.ribo" \
--groups "kidney:kidney_rep1,kidney_rep2,kidney_rep3;liver:liver_rep1,liver_rep2,liver_rep3;lung:lung_rep1,lung_rep2,lung_rep3" \
--tx_threshold 0.0 --min_z 0.0
--motif --reference "../reference_files/appris_mouse_v2_selected.fa.gz"
--flank-left 20 --flank-right 10
```
```
Number of filtered transcripts: 11031
Number of total stall sites per group: {'kidney': 50786, 'liver': 22308, 'lung': 20020}
```

Output:
* JSON file of stall sites. The indices are codons (not nucleotides), where 0 is the start codon "ATG."
* PNG
