# exprmclust

```python
pytscan.exprmclust(X, clusternum=range(2,10), reduce=True, n_pcs=2, random_state=0)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| X | np.ndarray | required | PCA matrix (n_cells, n_dims) |
| clusternum | int or list | range(2,10) | K values to test |
| reduce | bool | True | Use first n_pcs components |
| n_pcs | int | 2 | PCs for clustering |
| random_state | int | 0 | Random seed |

## Returns

dict with keys:

| Key | Type | Description |
|-----|------|-------------|
| clusterid | pd.Series | Cluster assignment (1-indexed) |
| clucenter | np.ndarray | Cluster centroids |
| MSTtree | nx.Graph | Minimum spanning tree |
| pcareduceres | np.ndarray | PCA coords used |
| bic_scores | dict | BIC per K tested |
| best_k | int | Optimal K selected |

## Example

```python
from pytscan import preprocess, exprmclust

X_pca = preprocess(adata, n_pcs=10)
mobj = exprmclust(X_pca, clusternum=range(3, 7), n_pcs=10)
print(f"Best K: {mobj['best_k']}")
print(f"BIC scores: {mobj['bic_scores']}")
```
