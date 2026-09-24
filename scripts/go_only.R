library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

sig <- read.csv(
  "results/significant_genes_annotated.csv",
  stringsAsFactors = FALSE
)

sig <- sig[!is.na(sig$SYMBOL), ]

gene.df <- bitr(
  sig$SYMBOL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego <- enrichGO(
  gene = gene.df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

tiff(
  "plots/GO_dotplot.tiff",
  width = 3000,
  height = 2400,
  res = 300
)

print(
  dotplot(
    ego,
    showCategory = 20
  )
)

dev.off()

tiff(
  "plots/GO_barplot.tiff",
  width = 3000,
  height = 2400,
  res = 300
)

print(
  barplot(
    ego,
    showCategory = 20
  )
)

dev.off()
