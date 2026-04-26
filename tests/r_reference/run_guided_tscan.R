library(TSCAN)
library(Seurat)

# 下载 pbmc3k
options(timeout=300)
if (!file.exists("pbmc3k_final.rds")) {
    cat("Downloading pbmc3k...\n")
    download.file("https://www.dropbox.com/s/63gnlw45jf7cje8/pbmc3k_final.rds?dl=1",
                  "pbmc3k_final.rds", mode="wb")
}

s <- readRDS("pbmc3k_final.rds")
cat("Loaded pbmc3k:", ncol(s), "cells\n")

# 提取数据
dr <- Embeddings(s, reduction="pca")[,1:10]
clu <- as.character(Idents(s))
names(clu) <- colnames(s)
mat <- as.matrix(s@assays$RNA@data)

# guided_tscan
ord <- guided_tscan(dr=dr, clu=clu, cluorder=c("Naive CD4 T","Memory CD4 T"))
cat("Pseudotime cells:", length(ord), "\n")
cat("Head:", head(ord), "\n")

# 保存图
library(ggplot2)
p <- genedynamics(mat["IL32",], ord)
ggsave("tests/r_reference/results/r_IL32_genedynamics.png", p, width=8, height=5)
cat("IL32 genedynamics saved\n")

# 保存数据
write.csv(dr, "tests/r_reference/results/r_pbmc_pca.csv")
write.csv(data.frame(cell=names(clu), cluster=clu),
          "tests/r_reference/results/r_pbmc_clu.csv", row.names=FALSE)
write.csv(data.frame(cell=ord, stringsAsFactors=FALSE),
          "tests/r_reference/results/r_pbmc_ord.csv", row.names=FALSE)
write.csv(data.frame(cell=colnames(mat), IL32=as.numeric(mat["IL32",])),
          "tests/r_reference/results/r_pbmc_IL32.csv", row.names=FALSE)
cat("All saved!\n")
