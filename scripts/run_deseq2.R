library(DESeq2)

counts <- read.csv(
  "counts/gene_counts.csv",
  row.names = 1,
  check.names = FALSE
)

colnames(counts) <- sub(".sorted.bam","",colnames(counts))

metadata <- read.csv(
  "metadata/sample_metadata.csv",
  row.names = 1
)

counts <- counts[, rownames(metadata)]

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c(
    "condition",
    "STIC_Yes",
    "STIC_No"
  )
)

res <- res[order(res$padj),]

write.csv(
  as.data.frame(res),
  file = "results/deseq2_results.csv"
)

saveRDS(
  dds,
  file = "results/dds.rds"
)

write.csv(
  counts(dds, normalized = TRUE),
  file = "results/normalized_counts.csv"
)
