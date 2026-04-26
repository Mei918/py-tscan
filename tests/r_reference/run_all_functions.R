library(TSCAN)
library(ggplot2)
set.seed(42)

data(lpsdata)
procdata <- preprocess(lpsdata)
lpsmclust <- exprmclust(procdata, clusternum=2:6)
lpsorder <- TSCANorder(lpsmclust, orderonly=FALSE)
lpsorder_only <- TSCANorder(lpsmclust, orderonly=TRUE)

# orderscore
subpopulation <- data.frame(
    cell = colnames(procdata),
    sub = ifelse(grepl("Unstimulated", colnames(procdata)), 0, 1),
    stringsAsFactors = FALSE
)
order1 <- TSCANorder(lpsmclust, orderonly=TRUE)
order2 <- TSCANorder(lpsmclust, c(1,2,3), orderonly=TRUE)
scores <- orderscore(subpopulation, list(order1, order2))
cat("orderscore:", scores, "\n")

# plotmclust
png("tests/r_reference/results/r_plotmclust.png", width=800, height=600)
plotmclust(lpsmclust)
dev.off()

# singlegeneplot
STAT2expr <- log2(lpsdata["STAT2",] + 1)
p <- singlegeneplot(STAT2expr, lpsorder)
ggsave("tests/r_reference/results/r_singlegeneplot.png", p, width=8, height=4)

# 保存数据
write.csv(procdata, "tests/r_reference/results/r_procdata.csv")
write.csv(lpsorder, "tests/r_reference/results/r_order.csv")
write.csv(
    data.frame(name=names(STAT2expr), expr=as.numeric(STAT2expr)),
    "tests/r_reference/results/r_STAT2expr.csv", row.names=FALSE
)
write.csv(subpopulation, "tests/r_reference/results/r_subpop.csv", row.names=FALSE)
write.csv(
    data.frame(order1=order1, stringsAsFactors=FALSE),
    "tests/r_reference/results/r_order_only.csv", row.names=FALSE
)
cat("scores:", scores, "\n")
cat("All done!\n")
