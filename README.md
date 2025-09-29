1. Clone GitHub repository:
```
git clone https://github.com/reikostachibana/ribostall
```

2. Create gzipped pickle file of {transcript: coverage array}.
```
python adj_coverage_v2.py --ribo "../bxc/bxc_disome.ribo" --site-type "stop" --min-len 57 --max-len 65 --procs 9 --out "../bxc/cov_di.pkl.gz"
```

3. 
