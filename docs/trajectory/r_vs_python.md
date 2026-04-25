# R vs Python Parity

Validation of pytscan against the original R TSCAN package.

## Results

| Metric | Value |
|--------|-------|
| R vs Python Pearson \|r\| | **0.9896** ✅ |
| Python vs ground truth \|r\| | **0.9921** ✅ |
| R vs ground truth \|r\| | **0.9977** ✅ |
| Best K agreement | **3 = 3** ✅ |

## Algorithm Equivalence

| R TSCAN | pytscan | Notes |
|---------|---------|-------|
| mclust VVV | GaussianMixture(covariance_type='full') | BIC-optimal K |
| igraph MST | scipy.sparse.csgraph.minimum_spanning_tree | Euclidean distances |
| TSCANorder | TSCANorder | Cell projection onto MST path |

## Reproduce

```bash
# Step 1: Run R TSCAN
conda activate CMAP
Rscript tests/r_reference/run_tscan_r.R

# Step 2: Run pytscan comparison
conda activate omicverse
python tests/compare_r_python.py
```

## Notebook

See [comparison_R_vs_Python.ipynb](../../notebooks/comparison_R_vs_Python.ipynb)
for full visual comparison with 6-panel figure.
