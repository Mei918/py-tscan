import numpy as np
from scipy.sparse import issparse
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
import anndata as ad


def preprocess(data, takelog=True, logbase=2, pseudocount=1,
               minexpr_value=1, minexpr_percent=0.5, cvcutoff=1, clusternum=None):
    if isinstance(data, ad.AnnData):
        X = data.layers["counts"] if "counts" in data.layers else data.X
        if issparse(X):
            X = X.toarray()
        mat = X.T.astype(np.float64)
    else:
        mat = np.array(data, dtype=np.float64)
    if takelog:
        mat = np.log(mat + pseudocount) / np.log(logbase)
    expr_filter = (mat > minexpr_value).mean(axis=1) > minexpr_percent
    row_mean = mat.mean(axis=1)
    row_std = mat.std(axis=1, ddof=1)
    cv = np.where(row_mean != 0, row_std / row_mean, 0)
    cv_filter = cv > cvcutoff
    mat = mat[expr_filter & cv_filter, :]
    if clusternum is not None:
        Z = linkage(mat, method="complete", metric="euclidean")
        cluster = fcluster(Z, clusternum, criterion="maxclust")
        mat = np.array([mat[cluster == k].mean(axis=0) for k in range(1, clusternum+1)])
    return mat