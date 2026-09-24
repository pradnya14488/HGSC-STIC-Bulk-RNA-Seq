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

################################################
# GO Biological Process
################################################

ego <- enrichGO(
  gene = gene.df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

write.csv(
  as.data.frame(ego),
  "results/GO_enrichment.csv",
  row.names = FALSE
)

################################################
# KEGG
################################################

ekegg <- enrichKEGG(
  gene = gene.df$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH"
)

write.csv(
  as.data.frame(ekegg),
  "results/KEGG_enrichment.csv",
  row.names = FALSE
)

################################################
# GO dotplot
################################################

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

################################################
# KEGG dotplot
################################################

tiff(
  "plots/KEGG_dotplot.tiff",
  width = 3000,
  height = 2400,
  res = 300
)

print(
  dotplot(
    ekegg,
    showCategory = 20
  )
)

dev.off()
