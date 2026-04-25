import numpy as np
import pandas as pd
import pytest
import anndata as ad
import scipy.sparse as sp
from pytscan import exprmclust, TSCANorder, preprocess
from pytscan.core import _project_point_to_segment, _get_longest_path
import networkx as nx


@pytest.fixture(scope="module")
def linear_adata():
    np.random.seed(42)
    n_cells, n_genes = 300, 200
    true_pt = np.linspace(0, 1, n_cells)
    X = np.outer(true_pt, np.random.randn(n_genes)*10) + np.random.randn(n_cells, n_genes)*0.5
    X = (X - X.min() + 1) * 50
    adata = ad.AnnData(X=X)
    adata.layers["counts"] = X.astype(int)
    adata.obs["true_pseudotime"] = true_pt
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    return adata


def test_preprocess_shape(linear_adata):
    X = preprocess(linear_adata, n_pcs=10)
    assert X.shape == (300, 10)

def test_preprocess_no_nan(linear_adata):
    X = preprocess(linear_adata)
    assert not np.any(np.isnan(X))

def test_preprocess_reproducible(linear_adata):
    X1 = preprocess(linear_adata, random_state=0)
    X2 = preprocess(linear_adata, random_state=0)
    np.testing.assert_array_almost_equal(X1, X2)

def test_preprocess_sparse(linear_adata):
    adata2 = linear_adata.copy()
    adata2.X = sp.csr_matrix(adata2.X)
    X = preprocess(adata2, n_pcs=10)
    assert X.shape[0] == 300

def test_exprmclust_keys(linear_adata):
    X = preprocess(linear_adata, n_pcs=10, random_state=0)
    mobj = exprmclust(X, clusternum=[3], random_state=0)
    for k in ["clusterid","clucenter","MSTtree","pcareduceres","model","bic_scores","best_k"]:
        assert k in mobj

def test_exprmclust_1indexed(linear_adata):
    X = preprocess(linear_adata, n_pcs=10, random_state=0)
    mobj = exprmclust(X, clusternum=[4], random_state=0)
    assert mobj["clusterid"].min() >= 1

def test_exprmclust_is_tree(linear_adata):
    X = preprocess(linear_adata, n_pcs=10, random_state=0)
    mobj = exprmclust(X, clusternum=[4], random_state=0)
    assert nx.is_tree(mobj["MSTtree"])

def test_exprmclust_bic(linear_adata):
    X = preprocess(linear_adata, n_pcs=10, random_state=0)
    mobj = exprmclust(X, clusternum=list(range(2,7)), random_state=0)
    best_bic = mobj["bic_scores"][mobj["best_k"]]
    for k, bic in mobj["bic_scores"].items():
        assert best_bic <= bic

def test_tscanorder_columns(linear_adata):
    X = preprocess(linear_adata, n_pcs=10, random_state=0)
    mobj = exprmclust(X, clusternum=[4], random_state=0)
    df = TSCANorder(mobj)
    for col in ["cell_index","State","Pseudotime"]:
        assert col in df.columns

def test_tscanorder_all_cells(linear_adata):
    X = preprocess(linear_adata, n_pcs=10, random_state=0)
    mobj = exprmclust(X, clusternum=[4], random_state=0)
    df = TSCANorder(mobj)
    assert len(df) == 300

def test_tscanorder_range(linear_adata):
    X = preprocess(linear_adata, n_pcs=10, random_state=0)
    mobj = exprmclust(X, clusternum=[4], random_state=0)
    df = TSCANorder(mobj)
    assert df["Pseudotime"].min() >= 0.0
    assert df["Pseudotime"].max() <= 1.0

def test_pseudotime_correlation(linear_adata):
    from scipy.stats import pearsonr
    X = preprocess(linear_adata, n_pcs=10, random_state=0)
    mobj = exprmclust(X, clusternum=list(range(3,7)), n_pcs=10, random_state=0)
    df = TSCANorder(mobj)
    true_pt = linear_adata.obs["true_pseudotime"].values
    pred_pt = np.zeros(300)
    for _, row in df.iterrows():
        pred_pt[int(row["cell_index"])] = row["Pseudotime"]
    r, _ = pearsonr(true_pt, pred_pt)
    assert abs(r) >= 0.90, f"Correlation {abs(r):.4f} < 0.90"

def test_project_segment_midpoint():
    a, b = np.array([0.,0.]), np.array([2.,0.])
    proj, t = _project_point_to_segment(np.array([1.,1.]), a, b)
    np.testing.assert_array_almost_equal(proj, [1.,0.])
    assert abs(t - 0.5) < 1e-9

def test_project_segment_clamp():
    a, b = np.array([0.,0.]), np.array([1.,0.])
    _, t = _project_point_to_segment(np.array([-1.,0.]), a, b)
    assert t == 0.0
    _, t = _project_point_to_segment(np.array([2.,0.]), a, b)
    assert t == 1.0

def test_difftest(linear_adata):
    from pytscan import difftest
    X = preprocess(linear_adata, n_pcs=10, random_state=0)
    mobj = exprmclust(X, clusternum=[4], random_state=0)
    df = TSCANorder(mobj)
    result = difftest(linear_adata, df, genes=["gene_0","gene_1","gene_2"])
    assert "pvalue" in result.columns
    assert "qvalue" in result.columns
    assert result["pvalue"].is_monotonic_increasing
