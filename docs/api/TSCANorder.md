# TSCANorder

```python
pytscan.TSCANorder(mclustobj, MSTorder=None, startcluster=None, orderonly=False)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| mclustobj | dict | required | Output of exprmclust() |
| MSTorder | list | None | Custom cluster path (auto if None) |
| startcluster | int | None | Start cluster for path selection |
| orderonly | bool | False | Return cell indices only |

## Returns

`pd.DataFrame` with columns:

| Column | Description |
|--------|-------------|
| cell_index | Original cell index |
| State | Cluster assignment |
| Pseudotime | Normalized pseudotime [0, 1] |

Or `list` of cell indices if `orderonly=True`.

## Example

```python
from pytscan import TSCANorder

pt_df = TSCANorder(mobj)
print(pt_df.head())

# Custom path
pt_df = TSCANorder(mobj, MSTorder=[2, 1, 3], startcluster=2)

# Indices only
order = TSCANorder(mobj, orderonly=True)
```
