library(DESeq2)
library(pheatmap)

dds <- readRDS("results/dds.rds")

vsd <- vst(dds)

sig <- read.csv(
  "results/significant_genes_annotated.csv",
  stringsAsFactors = FALSE
)

mat <- assay(vsd)

# Remove version numbers from matrix rownames
mat_ids <- sub("\\..*", "", rownames(mat))

# keep top 30 significant genes
top <- sig[order(sig$padj), ]
top <- top[1:min(30, nrow(top)), ]

idx <- match(top$ENSEMBL, mat_ids)

keep <- !is.na(idx)

mat2 <- mat[idx[keep], ]

rownames(mat2) <- top$SYMBOL[keep]

# remove rows with missing values
mat2 <- mat2[complete.cases(mat2), ]

tiff(
  "plots/Heatmap_gene_symbols.tiff",
  width = 3000,
  height = 3000,
  res = 300
)

pheatmap(
  mat2,
  scale = "row",
  fontsize_row = 10,
  fontsize_col = 12
)

dev.off()
