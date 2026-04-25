# Basic TSCAN Tutorial

Linear trajectory inference on synthetic data.

## Quick Start

```python
import numpy as np
import anndata as ad
from pytscan import preprocess, exprmclust, TSCANorder, plotmclust

# Load data (must have layers['counts'])
adata = ad.read_h5ad("your_data.h5ad")

# Step 1: Preprocess
X_pca = preprocess(adata, n_pcs=50, n_highly_variable=500)

# Step 2: Cluster + MST
mobj = exprmclust(X_pca, clusternum=range(2, 10))
print(f"Best K: {mobj['best_k']}")

# Step 3: Pseudotime
pt_df = TSCANorder(mobj)
print(pt_df.head())

# Visualize
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7, 6))
plotmclust(mobj, ax=ax)
plt.show()
```

## Key Parameters

### preprocess()
| Parameter | Default | Description |
|-----------|---------|-------------|
| n_highly_variable | 500 | Top variable genes to select |
| n_pcs | 50 | Number of PCA components |
| log_transform | True | Apply log2(CPM/10+1) |
| use_counts_layer | True | Use adata.layers['counts'] |

### exprmclust()
| Parameter | Default | Description |
|-----------|---------|-------------|
| clusternum | range(2,10) | K values to test (BIC selects best) |
| n_pcs | 2 | PCs used for clustering |
| reduce | True | Use reduced PCA space |

### TSCANorder()
| Parameter | Default | Description |
|-----------|---------|-------------|
| MSTorder | None | Custom path (auto if None) |
| startcluster | None | Start cluster for path |
| orderonly | False | Return indices only |

## Notes

- Input AnnData must contain `layers['counts']` (raw counts)
- TSCAN uses log2(CPM/10+1), different from scanpy's log1p
- BIC criterion selects optimal K automatically
