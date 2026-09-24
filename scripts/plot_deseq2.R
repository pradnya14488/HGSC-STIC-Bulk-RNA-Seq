library(DESeq2)
library(ggplot2)
library(pheatmap)

dds <- readRDS("results/dds.rds")

res <- read.csv(
  "results/deseq2_results.csv",
  row.names = 1
)

vsd <- vst(dds)

dir.create("plots", showWarnings = FALSE)

################################################
# PCA
################################################

tiff(
  "plots/PCA_plot.tiff",
  width = 2400,
  height = 2000,
  res = 300
)

plotPCA(
  vsd,
  intgroup = "condition"
)

dev.off()

################################################
# MA Plot
################################################

tiff(
  "plots/MA_plot.tiff",
  width = 2400,
  height = 2000,
  res = 300
)

plotMA(
  results(dds),
  ylim = c(-10,10)
)

dev.off()

################################################
# Volcano Plot
################################################

res$Significant <- "No"

res$Significant[
  which(
    res$padj < 0.05 &
    abs(res$log2FoldChange) > 1
  )
] <- "Yes"

tiff(
  "plots/Volcano_plot.tiff",
  width = 2400,
  height = 2000,
  res = 300
)

ggplot(
  res,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = Significant
  )
) +
  geom_point(size = 1.5) +
  theme_bw()

dev.off()

################################################
# Heatmap
################################################

topgenes <- head(
  rownames(
    res[order(res$padj),]
  ),
  50
)

mat <- assay(vsd)[topgenes,]

tiff(
  "plots/Heatmap.tiff",
  width = 2400,
  height = 2400,
  res = 300
)

pheatmap(
  mat,
  scale = "row"
)

dev.off()

################################################
# Boxplot
################################################

norm_counts <- counts(
  dds,
  normalized = TRUE
)

tiff(
  "plots/Boxplot.tiff",
  width = 2400,
  height = 2000,
  res = 300
)

boxplot(
  log2(norm_counts + 1),
  las = 2,
  col = "lightblue"
)

dev.off()
