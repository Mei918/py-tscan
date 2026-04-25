# difftest

```python
pytscan.difftest(adata, pseudotime_df, genes=None, use_counts_layer=True)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| adata | AnnData | required | Expression data |
| pseudotime_df | pd.DataFrame | required | Output of TSCANorder() |
| genes | list | None | Genes to test (all if None) |
| use_counts_layer | bool | True | Use layers['counts'] |

## Returns

`pd.DataFrame` sorted by pvalue:

| Column | Description |
|--------|-------------|
| gene | Gene name |
| pvalue | Linear regression p-value |
| qvalue | BH-corrected FDR |
| slope | Effect direction along pseudotime |

## Example

```python
from pytscan import difftest

# Test all genes
deg_df = difftest(adata, pt_df)
sig = deg_df[deg_df.qvalue < 0.05]
print(f"Significant: {len(sig)} genes")

# Test specific genes
deg_df = difftest(adata, pt_df, genes=['CD3D', 'CD14'])
```
