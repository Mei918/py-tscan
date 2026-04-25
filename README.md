# pytscan

[![PyPI version](https://badge.fury.io/py/pytscan.svg)](https://pypi.org/project/pytscan/)
[![Tests](https://img.shields.io/badge/tests-15%20passed-brightgreen)](https://github.com/Mei918/py-tscan)
[![R parity](https://img.shields.io/badge/R%20parity-%7Cr%7C%3D0.99-blue)](https://github.com/Mei918/py-tscan)

Python reimplementation of [TSCAN](https://github.com/zji90/TSCAN) — trajectory inference via Gaussian mixture clustering and minimum spanning tree.

> **Reference**: Ji, Z. & Ji, H. (2016). TSCAN. *Nucleic Acids Research*, 44(13), e117. [doi:10.1093/nar/gkw430](https://doi.org/10.1093/nar/gkw430)

## Install

\`\`\`bash
pip install pytscan
\`\`\`

## Quick Start

\`\`\`python
from pytscan import preprocess, exprmclust, TSCANorder, plotmclust

X_pca = preprocess(adata, n_pcs=10)
mobj = exprmclust(X_pca, clusternum=range(2,10))
pt_df = TSCANorder(mobj)
plotmclust(mobj)
\`\`\`

## Validation

| Metric | Result |
|--------|--------|
| R vs Python \|r\| | **0.9896** |
| Python vs truth \|r\| | **0.9921** |
| Significant DE genes | **197/200** |
| Tests passing | **15/15** |

## Algorithm

1. **Preprocess**: log2(CPM/10+1) → top variable genes → PCA
2. **Cluster**: GaussianMixture BIC-optimal K (equivalent to R mclust VVV)
3. **MST**: minimum spanning tree on cluster centers
4. **Pseudotime**: project cells onto MST path

## License

MIT
