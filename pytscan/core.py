import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse import csr_matrix
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
import warnings


def _elbow_pcadim(sdev):
    x = np.arange(1, len(sdev)+1, dtype=float)
    best_rss = np.inf
    best_i = 2
    for i in range(2, 11):
        x2 = np.maximum(0, x - i)
        A = np.column_stack([np.ones_like(x), x, x2])
        coef, _, _, _ = np.linalg.lstsq(A, sdev, rcond=None)
        pred = A @ coef
        rss = np.sum((sdev - pred) ** 2)
        if rss < best_rss:
            best_rss = rss
            best_i = i
    return best_i + 1


def exprmclust(data, clustermethod="mclust", clusternum=None, modelNames="VVV", reduce=True, cluster=None, random_state=12345, n_pcs=None):
    np.random.seed(random_state)
    if clusternum is None:
        clusternum = list(range(2, 10))
    if isinstance(clusternum, int):
        clusternum = [clusternum]
    if reduce:
        data = np.array(data, dtype=np.float64)
        tmpdata = scale(data.T).T
        tmpdata = np.where(np.isnan(tmpdata), 0, tmpdata)
        cells_x_genes = tmpdata.T
        n_pc = min(20, cells_x_genes.shape[0]-1, cells_x_genes.shape[1])
        pca = PCA(n_components=n_pc, random_state=random_state)
        pca.fit(cells_x_genes)
        sdev = pca.singular_values_[:20] / np.sqrt(max(cells_x_genes.shape[0]-1, 1))
        pcadim = _elbow_pcadim(sdev[:min(20, len(sdev))])
        pcadim = min(pcadim, n_pc)
        rotation = pca.components_[:pcadim].T
        pcareduceres = cells_x_genes @ rotation
    else:
        pcareduceres = np.array(data, dtype=np.float64)
    if cluster is not None:
        clusterid_arr = np.array(cluster, dtype=int)
        clunum = len(np.unique(clusterid_arr))
        bic_scores = {}
        best_k = clunum
    elif clustermethod == "mclust":
        clusternum = [k for k in clusternum if k > 1]
        best_bic = np.inf
        best_model = None
        best_k = None
        bic_scores = {}
        for k in clusternum:
            if k >= pcareduceres.shape[0]:
                continue
            try:
                gm = GaussianMixture(n_components=k, covariance_type="full", n_init=5, random_state=random_state, max_iter=300)
                gm.fit(pcareduceres)
                bic = gm.bic(pcareduceres)
                bic_scores[k] = bic
                if bic < best_bic:
                    best_bic = bic
                    best_model = gm
                    best_k = k
            except Exception as e:
                warnings.warn(f"GMM k={k}: {e}")
        if best_model is None:
            raise ValueError("GMM fitting failed.")
        clusterid_arr = best_model.predict(pcareduceres) + 1
        clunum = best_k
    else:
        from sklearn.cluster import KMeans
        k = clusternum[0] if isinstance(clusternum, list) else clusternum
        km = KMeans(n_clusters=k, random_state=random_state, n_init=10)
        clusterid_arr = km.fit_predict(pcareduceres) + 1
        clunum = k
        bic_scores = {}
        best_k = k
    clucenter = np.array([pcareduceres[clusterid_arr==cid].mean(axis=0) for cid in range(1, clunum+1)])
    dist_matrix = cdist(clucenter, clucenter, metric="euclidean")
    mst_sparse = minimum_spanning_tree(csr_matrix(dist_matrix))
    mst_dense = mst_sparse.toarray()
    G = nx.Graph()
    for i in range(clunum):
        G.add_node(i+1)
    rows, cols = np.where(mst_dense > 0)
    for r, c in zip(rows, cols):
        G.add_edge(r+1, c+1, weight=mst_dense[r, c])
    adjmat = np.zeros((clunum, clunum))
    for r, c in zip(rows, cols):
        adjmat[r, c] = 1
        adjmat[c, r] = 1
    return {"pcareduceres": pcareduceres, "MSTtree": G, "clusterid": pd.Series(clusterid_arr, name="cluster"), "clucenter": clucenter, "adjmat": adjmat, "best_k": clunum, "bic_scores": bic_scores}


def _get_mst_order(G, clutable, startcluster=None):
    alldeg = dict(G.degree())
    leaves = [n for n, d in alldeg.items() if d == 1]
    if len(leaves) < 2:
        return list(G.nodes())
    if startcluster is not None:
        allcomb = [(startcluster, e) for e in leaves if e != startcluster]
    else:
        allcomb = [(s, e) for i, s in enumerate(leaves) for e in leaves if s < e]
    best_path = None
    best_score = (-1, -1)
    for s, e in allcomb:
        try:
            path = nx.shortest_path(G, s, e)
            n_nodes = len(path)
            n_cells = sum(clutable.get(n, 0) for n in path)
            score = (n_nodes, n_cells)
            if score > best_score:
                best_score = score
                best_path = path
        except nx.NetworkXNoPath:
            continue
    return best_path


def TSCANorder(mclustobj, MSTorder=None, startcluster=None, orderonly=False, flip=False, listbranch=False, divide=True):
    clucenter = mclustobj["clucenter"]
    clusterid = mclustobj["clusterid"].values
    pcareduceres = mclustobj["pcareduceres"]
    G = mclustobj["MSTtree"]
    adjmat = mclustobj.get("adjmat", None)
    clunum = clucenter.shape[0]
    n_cells = len(clusterid)
    cell_names = np.array([str(i) for i in range(n_cells)])
    if adjmat is None:
        adjmat = np.zeros((clunum, clunum))
        for u, v in G.edges():
            adjmat[u-1, v-1] = 1
            adjmat[v-1, u-1] = 1
    clutable = {k: int((clusterid==k).sum()) for k in range(1, clunum+1)}
    orderinMST = 1
    if MSTorder is None:
        MSTorder = _get_mst_order(G, clutable, startcluster)
        if flip:
            MSTorder = list(reversed(MSTorder))
    else:
        if divide:
            edge_in_mst = all(adjmat[MSTorder[i]-1, MSTorder[i+1]-1]==1 for i in range(len(MSTorder)-1))
            orderinMST = 1 if edge_in_mst else 0
        else:
            orderinMST = 0

    def internal_order(internalorder, MSTinout):
        tscan_order = []
        for i in range(len(internalorder)-1):
            currentcluid = internalorder[i]
            nextcluid = internalorder[i+1]
            currentclucenter = clucenter[currentcluid-1]
            nextclucenter = clucenter[nextcluid-1]
            current_mask = clusterid == currentcluid
            current_idx = np.where(current_mask)[0]
            currentreduceres = pcareduceres[current_mask]
            if MSTinout:
                connectcluid_curr = [j+1 for j in range(clunum) if adjmat[currentcluid-1, j]==1]
            else:
                connectcluid_curr = [nextcluid] if i==0 else [nextcluid, internalorder[i-1]]
            cludist_curr = np.column_stack([np.sum((currentreduceres - clucenter[x-1])**2, axis=1) for x in connectcluid_curr])
            if cludist_curr.ndim == 1:
                cludist_curr = cludist_curr.reshape(-1, 1)
            mindistid_curr = np.argmin(cludist_curr, axis=1)
            next_pos = connectcluid_curr.index(nextcluid)
            edgecell_curr_idx = current_idx[mindistid_curr == next_pos]
            difvec = nextclucenter - currentclucenter
            tmppos = pcareduceres[edgecell_curr_idx] @ difvec
            tscan_order.extend(edgecell_curr_idx[np.argsort(tmppos)].tolist())
            next_mask = clusterid == nextcluid
            next_idx = np.where(next_mask)[0]
            nextreduceres = pcareduceres[next_mask]
            if MSTinout:
                connectcluid_next = [j+1 for j in range(clunum) if adjmat[nextcluid-1, j]==1]
            else:
                connectcluid_next = [currentcluid] if i==len(internalorder)-2 else [currentcluid, internalorder[i+2]]
            cludist_next = np.column_stack([np.sum((nextreduceres - clucenter[x-1])**2, axis=1) for x in connectcluid_next])
            if cludist_next.ndim == 1:
                cludist_next = cludist_next.reshape(-1, 1)
            mindistid_next = np.zeros(len(next_idx), dtype=int) if cludist_next.shape[1]==1 else np.argmin(cludist_next, axis=1)
            curr_pos = connectcluid_next.index(currentcluid)
            edgecell_next_idx = next_idx[mindistid_next == curr_pos]
            tmppos2 = pcareduceres[edgecell_next_idx] @ difvec
            tscan_order.extend(edgecell_next_idx[np.argsort(tmppos2)].tolist())
        if orderonly:
            return tscan_order
        return pd.DataFrame({"sample_name": cell_names[tscan_order], "State": clusterid[tscan_order], "Pseudotime": np.arange(1, len(tscan_order)+1)})

    if not orderinMST:
        return internal_order(MSTorder, 0)
    else:
        return internal_order(MSTorder, 1)


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
    cell_order = pseudotime_df["cell_index"].values if "cell_index" in pseudotime_df.columns else np.arange(len(pseudotime_df))
    X_ordered = X[cell_order, :]
    pt = pseudotime_df["Pseudotime"].values.astype(float)
    if pt.max() > 1:
        pt = (pt - pt.min()) / (pt.max() - pt.min())
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


def orderscore(subpopulation, orders):
    subinfo = dict(zip(subpopulation.iloc[:, 0], subpopulation.iloc[:, 1]))
    def score_one(order):
        scoreorder = np.array([subinfo[c] for c in order if c in subinfo], dtype=float)
        if len(scoreorder) < 2:
            return 0.0
        optscoreorder = np.sort(scoreorder)
        n = len(scoreorder)
        def pairsum(arr):
            total = 0.0
            for i in range(n-1):
                total += np.sum(arr[i+1:] - arr[i])
            return total
        optscore = pairsum(optscoreorder)
        return 0.0 if optscore == 0 else pairsum(scoreorder) / optscore
    return np.array([score_one(o) for o in orders])


def singlegeneplot(geneexpr, order, cell_size=2, ax=None, k=3):
    import matplotlib.pyplot as plt
    from scipy.interpolate import UnivariateSpline
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))
    pt = order["Pseudotime"].values
    if "cell_index" in order.columns:
        expr = np.array(geneexpr)[order["cell_index"].values]
    else:
        expr = np.array(geneexpr)
    states = order["State"].values
    unique_states = np.unique(states)
    cmap = plt.cm.get_cmap("tab10", len(unique_states))
    color_map = {s: cmap(i) for i, s in enumerate(unique_states)}
    for s in unique_states:
        mask = states == s
        ax.scatter(pt[mask], expr[mask], color=color_map[s], s=cell_size*10, alpha=0.7, label=f"State {s}")
    try:
        sort_idx = np.argsort(pt)
        spl = UnivariateSpline(pt[sort_idx], expr[sort_idx], k=min(k,3), s=len(pt))
        pt_smooth = np.linspace(pt.min(), pt.max(), 200)
        ax.plot(pt_smooth, spl(pt_smooth), color="black", lw=1.5)
    except Exception:
        pass
    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("Expression")
    ax.legend(title="State")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return ax


def genedynamics(geneexpr, order, k=3, ax=None):
    import matplotlib.pyplot as plt
    from scipy.interpolate import UnivariateSpline
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))
    pt = order["Pseudotime"].values
    if "cell_index" in order.columns:
        expr = np.array(geneexpr)[order["cell_index"].values]
    else:
        expr = np.array(geneexpr)
    ax.scatter(pt, expr, color="royalblue", alpha=0.5, s=20)
    try:
        sort_idx = np.argsort(pt)
        spl = UnivariateSpline(pt[sort_idx], expr[sort_idx], k=min(k,3), s=len(pt))
        pt_smooth = np.linspace(pt.min(), pt.max(), 200)
        ax.plot(pt_smooth, spl(pt_smooth), color="orange", lw=2)
    except Exception:
        pass
    ax.set_xlabel("pseudotime")
    ax.set_ylabel("gene expression")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return ax


def guided_MST(dr, clu, ax=None):
    import matplotlib.pyplot as plt
    clu = np.array(clu).astype(str)
    unique_clu = np.unique(clu)
    clucenter = np.array([dr[clu==c].mean(axis=0) for c in unique_clu])
    dist_matrix = cdist(clucenter, clucenter, metric="euclidean")
    mst_sparse = minimum_spanning_tree(csr_matrix(dist_matrix))
    mst_dense = mst_sparse.toarray()
    G = nx.Graph()
    for c in unique_clu:
        G.add_node(c)
    rows, cols = np.where(mst_dense > 0)
    for r, c in zip(rows, cols):
        G.add_edge(unique_clu[r], unique_clu[c], weight=mst_dense[r, c])
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    pos = {c: clucenter[i, :2] for i, c in enumerate(unique_clu)}
    nx.draw(G, pos=pos, with_labels=True, node_color="white", edgecolors="black", node_size=800, ax=ax)
    ax.set_title("guided MST")
    return {"MSTtree": G, "clucenter": clucenter, "unique_clu": unique_clu}


def guided_tscan(dr, clu, cluorder):
    clu = np.array(clu).astype(str)
    cluorder = [str(c) for c in cluorder]
    unique_clu = np.unique(clu)
    cm = {c: dr[clu==c].mean(axis=0) for c in unique_clu}
    left = {c: [] for c in unique_clu}
    right = {c: [] for c in unique_clu}
    for i, c in enumerate(cluorder):
        cell_idx = np.where(clu==c)[0]
        cells_dr = dr[clu==c]
        if i == 0:
            right[c] = cell_idx.tolist()
        elif i == len(cluorder)-1:
            left[c] = cell_idx.tolist()
        else:
            leftdist = np.sum((cells_dr - cm[cluorder[i-1]])**2, axis=1)
            rightdist = np.sum((cells_dr - cm[cluorder[i+1]])**2, axis=1)
            left[c] = cell_idx[leftdist <= rightdist].tolist()
            right[c] = cell_idx[leftdist > rightdist].tolist()
    ord_indices = []
    for i in range(len(cluorder)-1):
        c_curr, c_next = cluorder[i], cluorder[i+1]
        difvec = cm[c_next] - cm[c_curr]
        seg_idx = np.array(right[c_curr] + left[c_next])
        if len(seg_idx) == 0:
            continue
        proj = dr[seg_idx] @ difvec
        ord_indices.extend(seg_idx[np.argsort(proj)].tolist())
    return np.array(ord_indices)


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
    for i, center in enumerate(clucenter):
        ax.text(center[0], center[1], str(i+1), ha="center", va="center", fontsize=9, fontweight="bold", zorder=8)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"TSCAN clustering (k={mclustobj[chr(39)+chr(98)+chr(101)+chr(115)+chr(116)+chr(95)+chr(107)+chr(39)]})")
    return ax
