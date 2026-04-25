# plotmclust

```python
pytscan.plotmclust(mclustobj, cell_labels=None, MSTorder=None, show_mst=True, ax=None)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| mclustobj | dict | required | Output of exprmclust() |
| cell_labels | array | None | Labels for coloring (uses cluster if None) |
| MSTorder | list | None | Highlight trajectory path |
| show_mst | bool | True | Draw MST edges |
| ax | matplotlib Axes | None | Axes to plot on |

## Returns

`matplotlib.Axes`

## Example

```python
from pytscan import plotmclust
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Color by cluster
plotmclust(mobj, ax=axes[0])

# Color by cell type + show trajectory
plotmclust(mobj, cell_labels=adata.obs['celltype'],
           MSTorder=[2,1,3], ax=axes[1])

plt.tight_layout()
plt.savefig('trajectory.png', dpi=150)
plt.show()
```
