library(org.Hs.eg.db)
library(AnnotationDbi)

res <- read.csv(
  "results/deseq2_results.csv",
  row.names = 1,
  check.names = FALSE
)

ensembl_ids <- sub("\\..*", "", rownames(res))

annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids,
  columns = c(
    "SYMBOL",
    "GENENAME"
  ),
  keytype = "ENSEMBL"
)

res$ENSEMBL <- ensembl_ids

merged <- merge(
  res,
  annot,
  by.x = "ENSEMBL",
  by.y = "ENSEMBL",
  all.x = TRUE
)

write.csv(
  merged,
  "results/deseq2_results_annotated.csv",
  row.names = FALSE
)

sig <- subset(
  merged,
  !is.na(padj) & padj < 0.05
)

write.csv(
  sig,
  "results/significant_genes_annotated.csv",
  row.names = FALSE
)

cat("Annotated files written\n")
