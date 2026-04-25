# preprocess

```python
pytscan.preprocess(adata, n_highly_variable=500, log_transform=True,
                   scale=True, n_pcs=50, use_counts_layer=True, random_state=0)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| adata | AnnData | required | Input data with raw counts |
| n_highly_variable | int | 500 | Top variable genes |
| log_transform | bool | True | log2(CPM/10+1) normalization |
| scale | bool | True | Z-score scale genes |
| n_pcs | int | 50 | PCA components |
| use_counts_layer | bool | True | Use layers['counts'] |
| random_state | int | 0 | Random seed |

## Returns

`np.ndarray` of shape `(n_cells, n_pcs)`

## Example

```python
from pytscan import preprocess
X_pca = preprocess(adata, n_pcs=10, n_highly_variable=200)
print(X_pca.shape)  # (n_cells, 10)
```
