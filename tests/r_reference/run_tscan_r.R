library(TSCAN)
set.seed(42)

# 生成与 Python 完全相同的测试数据
n_cells <- 300
n_genes <- 200
true_pt <- seq(0, 1, length.out=n_cells)

set.seed(42)
basis <- matrix(rnorm(n_genes) * 10, nrow=1)
X <- matrix(rnorm(n_cells * n_genes, sd=0.5), n_cells, n_genes)
X <- X + outer(true_pt, as.vector(basis))
X <- X - min(X) + 1
X <- X * 50

# TSCAN 需要 genes x cells
expr_mat <- t(X)
rownames(expr_mat) <- paste0("gene_", 0:(n_genes-1))
colnames(expr_mat) <- paste0("cell_", 0:(n_cells-1))

cat("Running exprmclust...\n")
mobj <- exprmclust(expr_mat, clusternum=3:6, modelNames="VVV", reduce=TRUE)
cat("Best K:", length(unique(mobj$clusterid)), "\n")

cat("Running TSCANorder...\n")
pt_result <- TSCANorder(mobj, orderonly=FALSE)

# 保存结果
out_dir <- "tests/r_reference/results"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# 伪时间结果
write.csv(pt_result, file.path(out_dir, "r_pseudotime.csv"), row.names=FALSE)

# 聚类结果
clusterid_df <- data.frame(
  cell = names(mobj$clusterid),
  cluster = as.integer(mobj$clusterid)
)
write.csv(clusterid_df, file.path(out_dir, "r_clusterid.csv"), row.names=FALSE)

# PCA坐标
write.csv(as.data.frame(mobj$pcareduceres),
          file.path(out_dir, "r_pca.csv"), row.names=FALSE)

# 聚类中心
write.csv(as.data.frame(mobj$clucenter),
          file.path(out_dir, "r_clucenter.csv"), row.names=FALSE)

cat("Results saved to:", out_dir, "\n")
cat("Pseudotime range:", range(pt_result$Pseudotime), "\n")
cat("Done!\n")
