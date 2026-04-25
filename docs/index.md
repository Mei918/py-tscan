# pytscan Documentation

Python reimplementation of [TSCAN](https://github.com/zji90/TSCAN) for trajectory inference in single-cell RNA-seq data.

> **Reference**: Ji, Z. & Ji, H. (2016). TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis. *Nucleic Acids Research*, 44(13), e117. [doi:10.1093/nar/gkw430](https://doi.org/10.1093/nar/gkw430)

## Installation

```bash
pip install pytscan
```

## Overview

TSCAN orders cells along a pseudotime trajectory in three steps:

1. **Preprocess** — log2(CPM/10+1) normalization → top variable genes → PCA
2. **Cluster** — Gaussian Mixture Model (BIC-optimal K) + MST on cluster centers
3. **Pseudotime** — project cells onto MST path → distance-along-path

## Tutorials

| Notebook | Description |
|----------|-------------|
| [Basic TSCAN](trajectory/tscan_basic.md) | Linear trajectory on synthetic data |
| [R vs Python Parity](trajectory/r_vs_python.md) | Validation against R TSCAN |
| [Differential Expression](trajectory/difftest.md) | DE genes along pseudotime |

## API Reference

| Function | Description |
|----------|-------------|
| [preprocess](api/preprocess.md) | Normalize + PCA |
| [exprmclust](api/exprmclust.md) | GMM clustering + MST |
| [TSCANorder](api/TSCANorder.md) | Pseudotime ordering |
| [difftest](api/difftest.md) | DE along pseudotime |
| [plotmclust](api/plotmclust.md) | Visualization |

## Validation

| Metric | Result |
|--------|--------|
| R vs Python \|r\| | **0.9896** |
| Python vs truth \|r\| | **0.9921** |
| Significant DE genes | **197/200** |
| Tests passing | **15/15** |
