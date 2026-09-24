library(ggplot2)
library(ggrepel)

res <- read.csv(
  "results/significant_genes_annotated.csv"
)

res$neglog10 <- -log10(res$padj)

top_genes <- res[
  order(res$padj),
][1:15, ]

tiff(
  "plots/Volcano_plot_labeled.tiff",
  width = 2800,
  height = 2400,
  res = 300
)

ggplot(
  res,
  aes(
    x = log2FoldChange,
    y = neglog10
  )
) +
  geom_point(
    color = "grey60",
    size = 2
  ) +
  geom_point(
    data = top_genes,
    color = "red",
    size = 3
  ) +
  geom_text_repel(
    data = top_genes,
    aes(label = SYMBOL),
    size = 4
  ) +
  theme_bw() +
  xlab("Log2 Fold Change") +
  ylab("-log10 adjusted P-value") +
  ggtitle("HGSC STIC+ vs STIC-")

dev.off()
