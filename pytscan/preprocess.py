import numpy as np
from scipy.sparse import issparse
from sklearn.decomposition import PCA
import anndata as ad


def preprocess(adata, n_highly_variable=500, log_transform=True, scale=True, n_pcs=50, use_counts_layer=True, random_state=0):
    if use_counts_layer and "counts" in adata.layers:
        X = adata.layers["counts"]
    else:
        X = adata.X
    if issparse(X):
        X = X.toarray()
    X = X.astype(np.float64)
    if log_transform:
        lib_size = X.sum(axis=1, keepdims=True)
        lib_size = np.where(lib_size == 0, 1, lib_size)
        X = np.log2(X / lib_size * 10 + 1)
    gene_var = X.var(axis=0)
    top_genes = np.argsort(gene_var)[::-1][:n_highly_variable]
    X = X[:, top_genes]
    if scale:
        gene_mean = X.mean(axis=0)
        gene_std = X.std(axis=0)
        gene_std = np.where(gene_std == 0, 1, gene_std)
        X = (X - gene_mean) / gene_std
    n_pcs_actual = min(n_pcs, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n_pcs_actual, random_state=random_state)
    return pca.fit_transform(X)
