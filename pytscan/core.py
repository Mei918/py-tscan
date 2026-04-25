import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial.distance import cdist, euclidean
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse import csr_matrix
from sklearn.mixture import GaussianMixture
from typing import Optional, Union, List, Dict, Tuple
import warnings


def exprmclust(X, clusternum=None, modelNames="VVV", reduce=True, n_pcs=2, random_state=0):
    if clusternum is None:
        clusternum = list(range(2, 10))
    if reduce and X.shape[1] > n_pcs:
        X_use = X[:, :n_pcs]
    else:
        X_use = X
    if isinstance(clusternum, int):
        clusternum = [clusternum]
    best_bic = np.inf
    best_model = None
    best_k = None
    bic_scores = {}
    for k in clusternum:
        if k >= X_use.shape[0]:
            continue
        try:
            gm = GaussianMixture(n_components=k, covariance_type="full", n_init=5, random_state=random_state, max_iter=300)
            gm.fit(X_use)
            bic = gm.bic(X_use)
            bic_scores[k] = bic
            if bic < best_bic:
                best_bic = bic
                best_model = gm
                best_k = k
        except Exception as e:
            warnings.warn(f"GMM k={k}: {e}")
            continue
    if best_model is None:
        raise ValueError("GMM fitting failed for all tested cluster numbers.")
    labels = best_model.predict(X_use) + 1
    clusterid = pd.Series(labels, name="cluster")
    clucenter = np.array([X_use[labels == (k+1)].mean(axis=0) for k in range(best_k)])
    dist_matrix = cdist(clucenter, clucenter, metric="euclidean")
    mst_sparse = minimum_spanning_tree(csr_matrix(dist_matrix))
    mst_dense = mst_sparse.toarray()
    G = nx.Graph()
    for i in range(best_k):
        G.add_node(i + 1)
    rows, cols = np.where(mst_dense > 0)
    for r, c in zip(rows, cols):
        G.add_edge(r + 1, c + 1, weight=mst_dense[r, c])
    return {"clusterid": clusterid, "clucenter": clucenter, "MSTtree": G, "pcareduceres": X_use, "model": best_model, "bic_scores": bic_scores, "best_k": best_k}


def _get_longest_path(G, startcluster=None):
    leaves = [n for n, d in G.degree() if d == 1]
    if len(leaves) == 0:
        return list(G.nodes())
    if len(leaves) == 1:
        return list(nx.dfs_preorder_nodes(G, leaves[0]))
    if startcluster is not None and startcluster in G.nodes():
        best_path = []
        for end in leaves:
            if end == startcluster:
                continue
            try:
                path = nx.shortest_path(G, startcluster, end)
                if len(path) > len(best_path):
                    best_path = path
            except nx.NetworkXNoPath:
                continue
    else:
        best_path = []
        for i, s in enumerate(leaves):
            for end in leaves[i+1:]:
                try:
                    path = nx.shortest_path(G, s, end)
                    if len(path) > len(best_path):
                        best_path = path
                except nx.NetworkXNoPath:
                    continue
    if not best_path:
        start = startcluster if (startcluster and startcluster in G.nodes()) else leaves[0]
        best_path = list(nx.dfs_preorder_nodes(G, start))
    return best_path


def _project_point_to_segment(p, a, b):
    ab = b - a
    ab_norm_sq = np.dot(ab, ab)
    if ab_norm_sq < 1e-12:
        return a.copy(), 0.0
    t = np.clip(np.dot(p - a, ab) / ab_norm_sq, 0.0, 1.0)
    return a + t * ab, t


def _project_cells_to_path(X_use, clucenter, clusterid, path_clusters):
    n_cells = X_use.shape[0]
    path_centers = np.array([clucenter[c - 1] for c in path_clusters])
    seg_lengths = np.array([euclidean(path_centers[i], path_centers[i+1]) for i in range(len(path_centers)-1)])
    cum_lengths = np.concatenate([[0], np.cumsum(seg_lengths)])
    pseudotimes = np.zeros(n_cells)
    for cell_idx in range(n_cells):
        cell_pt = X_use[cell_idx]
        best_pt = 0.0
        best_dist = np.inf
        for seg_idx in range(len(path_centers) - 1):
            proj_pt, t = _project_point_to_segment(cell_pt, path_centers[seg_idx], path_centers[seg_idx+1])
            dist = euclidean(cell_pt, proj_pt)
            pt_along = cum_lengths[seg_idx] + t * seg_lengths[seg_idx]
            if dist < best_dist:
                best_dist = dist
                best_pt = pt_along
        pseudotimes[cell_idx] = best_pt
    return pseudotimes


def TSCANorder(mclustobj, MSTorder=None, startcluster=None, orderonly=False):
    G = mclustobj["MSTtree"]
    clucenter = mclustobj["clucenter"]
    clusterid = mclustobj["clusterid"]
    X_use = mclustobj["pcareduceres"]
    if MSTorder is None:
        MSTorder = _get_longest_path(G, startcluster)
    cell_pseudotimes = _project_cells_to_path(X_use, clucenter, clusterid, MSTorder)
    order = np.argsort(cell_pseudotimes)
    ordered_pt = cell_pseudotimes[order]
    pt_min, pt_max = ordered_pt.min(), ordered_pt.max()
    if pt_max > pt_min:
        ordered_pt_norm = (ordered_pt - pt_min) / (pt_max - pt_min)
    else:
        ordered_pt_norm = ordered_pt * 0.0
    if orderonly:
        return order.tolist()
    return pd.DataFrame({"cell_index": order, "State": clusterid.values[order], "Pseudotime": ordered_pt_norm})


def difftest(adata, pseudotime_df, genes=None, use_counts_layer=True):
    from scipy import stats
    from scipy.sparse import issparse
    from statsmodels.stats.multitest import multipletests
    if use_counts_layer and "counts" in adata.layers:
        X = adata.layers["counts"]
    else:
        X = adata.X
    if issparse(X):
        X = X.toarray()
    cell_order = pseudotime_df["cell_index"].values
    X_ordered = X[cell_order, :]
    pt = pseudotime_df["Pseudotime"].values
    gene_names = list(adata.var_names)
    if genes is not None:
        gene_idx = [gene_names.index(g) for g in genes if g in gene_names]
        gene_names = [gene_names[i] for i in gene_idx]
        X_ordered = X_ordered[:, gene_idx]
    results = []
    for i, gene in enumerate(gene_names):
        expr_log = np.log2(X_ordered[:, i].astype(float) + 1)
        slope, intercept, r, p, se = stats.linregress(pt, expr_log)
        results.append({"gene": gene, "pvalue": p, "slope": slope})
    df = pd.DataFrame(results)
    if len(df) > 0:
        _, qvalues, _, _ = multipletests(df["pvalue"].values, method="fdr_bh")
        df["qvalue"] = qvalues
        df = df.sort_values("pvalue").reset_index(drop=True)
    return df


def plotmclust(mclustobj, cell_labels=None, MSTorder=None, show_mst=True, ax=None):
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    X_use = mclustobj["pcareduceres"]
    clusterid = mclustobj["clusterid"].values
    clucenter = mclustobj["clucenter"]
    G = mclustobj["MSTtree"]
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
    color_source = cell_labels if cell_labels is not None else clusterid
    unique_labels = np.unique(color_source)
    cmap = cm.get_cmap("tab20", len(unique_labels))
    colors = [cmap(i) for i, lbl in enumerate(unique_labels) for _ in np.where(color_source == lbl)[0]]
    color_map = {lbl: cmap(i) for i, lbl in enumerate(unique_labels)}
    colors = [color_map[lbl] for lbl in color_source]
    ax.scatter(X_use[:, 0], X_use[:, 1], c=colors, s=20, alpha=0.7, linewidths=0)
    if show_mst:
        for u, v in G.edges():
            cu, cv = clucenter[u-1], clucenter[v-1]
            ax.plot([cu[0], cv[0]], [cu[1], cv[1]], "k-", lw=2, alpha=0.8, zorder=5)
    if MSTorder is not None:
        path_centers = np.array([clucenter[c-1] for c in MSTorder])
        ax.plot(path_centers[:, 0], path_centers[:, 1], "r-", lw=3, alpha=0.9, zorder=6)
    ax.scatter(clucenter[:, 0], clucenter[:, 1], c="white", s=120, edgecolors="black", lw=2, zorder=7)
    for i, (cx, cy) in enumerate(clucenter):
        ax.text(cx, cy, str(i+1), ha="center", va="center", fontsize=9, fontweight="bold", zorder=8)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"TSCAN clustering (k={mclustobj['best_k']})")
    return ax
