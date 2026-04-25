# Differential Expression Along Pseudotime

Test for genes that change significantly along the trajectory.

## Usage

```python
from pytscan import preprocess, exprmclust, TSCANorder, difftest

# Run trajectory first
X_pca = preprocess(adata, n_pcs=10)
mobj = exprmclust(X_pca, clusternum=range(3, 7))
pt_df = TSCANorder(mobj)

# Test all genes
deg_df = difftest(adata, pt_df)
print(f"Significant genes (q<0.05): {(deg_df.qvalue < 0.05).sum()}")
print(deg_df.head(10))

# Test specific genes
deg_df = difftest(adata, pt_df, genes=['CD3D', 'CD14', 'FCGR3A'])
```

## Output

| Column | Description |
|--------|-------------|
| gene | Gene name |
| pvalue | Linear regression p-value |
| qvalue | BH-corrected FDR |
| slope | Direction of change along pseudotime |

## Method

Fits linear model: `log2(expr+1) ~ pseudotime` per gene.
P-values corrected with Benjamini-Hochberg FDR.
