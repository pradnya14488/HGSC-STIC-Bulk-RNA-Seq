library(Rsubread)

bam_files <- list.files(
  "bam",
  pattern = "sorted.bam$",
  full.names = TRUE
)

fc <- featureCounts(
  files = bam_files,
  annot.ext = "data/reference/gencode.v48.annotation.gtf",
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  isPairedEnd = TRUE,
  nthreads = 26
)

write.csv(
  fc$counts,
  file = "counts/gene_counts.csv",
  quote = FALSE
)

saveRDS(
  fc,
  file = "counts/featureCounts.rds"
)
