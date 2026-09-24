# Zifo-HGSC
_The Zifo bioinformatics team is supporting the biologist client by providing complete script used for running RNAseq pipeline and completing downstream data analysis.  ._

---

## Overview

Commands to run for Bulk RNASeq Analysis

### In bash

```bash
# Create the conda environment
cd /mnt/active_storage/pradnya
mkdir RNASeqMock
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p /mnt/active_storage/pradnya/miniconda3
export PATH=/mnt/active_storage/pradnya/miniconda3/bin:$PATH
export PATH=/mnt/active_storage/pradnya/miniconda3/bin:$PATH
echo 'export PATH=/mnt/active_storage/pradnya/miniconda3/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
conda --version
conda init bash
source ~/.bashrc

# Create the conda env
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main

conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r  
codna tos view

conda create -n rnaseq \  
python=3.10 \  
hisat2 \  
samtools \  
r-base \  
r-essentials \  
bioconductor-deseq2 \  
bioconductor-rsubread \  
bioconductor-annotationdbi \  
bioconductor-org.hs.eg.db \  
-c conda-forge \  
-c bioconda \  
-y  

# Activate the conda env
# Locate the working dir /mnt/active_storage/pradnya
conda activate rnaseq
# Verify the package versions in command line in conda env
hisat2 --version  
samtools --version  
R --version

# Check the package version inside r, Start R
R
# Check package version
library(DESeq2)  
library(Rsubread)  
library(AnnotationDbi)  
library(org.Hs.eg.db)

sessionInfo()

# quit the R
q()




mkdir HGSC_STIC_4samples
cd HGSC_STIC_4samples

mkdir -p \  
data/raw_fastq \  
data/reference \  
hisat2_index \  
qc/fastqc \  
qc/multiqc \  
alignment \  
bam \  
counts \  
plots \  
results \  
logs \  
scripts \  
metadata

# Download the reference genome
cd /data/reference
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz

# Download the annotation:
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz

# Prepare the reference genome
# Uncompress the files
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v48.annotation.gtf.gz
ls -lh

# Verify the HISAT2
hisat2 --version
# Create the index directory
mkdir -p /mnt/active_storage/pradnya/RNASeqMock_GSE102094/HGSC_STIC_4samples/hisat2_index
ls -ld /mnt/active_storage/pradnya/RNASeqMock_GSE102094/HGSC_STIC_4samples/hisat2_index
# Check availability of disk space
df -h /mnt/active_storage
# Build the HISAT2 Index
cd /mnt/active_storage/pradnya/RNASeqMock
# Run the index build
hisat2-build \  
-p 16 \  
data/reference/GRCh38.primary_assembly.genome.fa \  
hisat2_index/grch38
ls -lh

# Package installation
conda config --add channels conda-forge  
conda config --set channel_priority strict
conda install -c conda-forge -c bioconda fastqc multiqc -y
## Save the HISAT2 Index Build
cd /mnt/active_storage/pradnya/RNASeqMock/hisat2_index  
ls -lh *.ht2
```
# Download the metadata

```
mkdir -p /mnt/active_storage/pradnya/RNASeqMock_GSE102094/metadata  
cd /mnt/active_storage/pradnya/RNASeqMock_GSE102094/metadata
 

wget -O index.html \  
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102094/
grep -i "matrix\|soft\|family" index.html

wget -O GSE102094_family.soft.gz \  
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102094/soft/GSE102094_family.soft.gz

ls -lh
gunzip GSE102094_family.soft.gz
ls -lh
# Extract the key metadata
grep '^!Sample_geo_accession' GSE102094_family.soft | head
grep '^!Sample_title' GSE102094_family.soft | head
grep 'stic:' GSE102094_family.soft | head -50
grep 'stic:' GSE102094_family.soft | sort | uniq -c
```
### 
### Create the metadata extraction script extract_samples.py
```python
import re

samples = []

with open("GSE102094_family.soft") as f:
    gsm = None
    title = None
    stic = None

    for line in f:

        if line.startswith("^SAMPLE"):
            if gsm:
                samples.append([gsm, title, stic])

            gsm = line.strip().split("=")[1].strip()
            title = None
            stic = None

        elif line.startswith("!Sample_title"):
            title = line.split("=",1)[1].strip()

        elif "stic:" in line:
            stic = line.split("stic:")[1].strip()

    if gsm:
        samples.append([gsm, title, stic])

with open("sample_metadata.tsv","w") as out:
    out.write("GSM\tTITLE\tSTIC\n")

    for s in samples:
        out.write(
            f"{s[0]}\t{s[1]}\t{s[2]}\n"
        )
```
### Run the Python script
```bash
python extract_samples.py
head sample_metadata.tsv  
grep "Yes" sample_metadata.tsv | head
grep "No" sample_metadata.tsv | head
# Count samples
cut -f3 sample_metadata.tsv | sort | uniq -c
# 147 No
# 1 STIC
# 139 Yes


# Extract GSM-SRX Relationship for detecting less file size 
grep '^!Sample_relation' GSE102094_family.soft | head -20
grep -B 20 'GSM2719894' GSE102094_family.soft | head -50
grep -B 20 'GSM2719895' GSE102094_family.soft | head -50  
```
### Extract four samples extract_first4.py
```python
import re

keep = {
    "GSM2719894",
    "GSM2719895",
    "GSM2719896",
    "GSM2719897"
}

gsm=None
title=None
stic=None
srx=None

with open("GSE102094_family.soft") as f:

    for line in f:

        if line.startswith("^SAMPLE"):
            if gsm in keep:
                print(gsm, title, stic, srx, sep="\t")

            gsm=line.split("=")[1].strip()
            title=None
            stic=None
            srx=None

        elif line.startswith("!Sample_title"):
            title=line.split("=",1)[1].strip()

        elif "stic:" in line:
            stic=line.split("stic:")[1].strip()

        elif "Sample_relation" in line and "SRX" in line:
            m=re.search(r"SRX\d+", line)
            if m:
                srx=m.group(0)

    if gsm in keep:
        print(gsm, title, stic, srx, sep="\t")
        
```
#### Run the script
```bash
python extract_first4.py

# Extract the SRX

for SRX in SRX3050119 SRX3050120 SRX3050121 SRX3050122
do
    echo "==== $SRX ===="
    wget -qO- "https://www.ncbi.nlm.nih.gov/sra/?term=$SRX" | grep -o 'SRR[0-9]*' | head
done

for SRR in SRR5884094 SRR5884095 SRR5884096 SRR5884097
do
    echo "==== $SRR ===="
    wget -qO- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=run_accession,scientific_name,library_layout,read_count,base_count"
done
```
### Download the samples
```bash
# Create sample metadata
nana sample_metadata.csv
sample,condition  
SRR5884094,STIC_No  
SRR5884096,STIC_No  
SRR5884095,STIC_Yes  
SRR5884097,STIC_Yes

# Verify
cat metadata/sample_metadata.csv  

# Create download script
cd data/raw_fastq

nano download_fastq.sh
#!/bin/bash
wget -c --no-check-certificate https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/004/SRR5884094/SRR5884094_1.fastq.gz
wget -c --no-check-certificate https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/004/SRR5884094/SRR5884094_2.fastq.gz

wget -c --no-check-certificate https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/005/SRR5884095/SRR5884095_1.fastq.gz
wget -c --no-check-certificate https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/005/SRR5884095/SRR5884095_2.fastq.gz

wget -c --no-check-certificate https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/006/SRR5884096/SRR5884096_1.fastq.gz
wget -c --no-check-certificate https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/006/SRR5884096/SRR5884096_2.fastq.gz

wget -c --no-check-certificate https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/007/SRR5884097/SRR5884097_1.fastq.gz
wget -c --no-check-certificate https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/007/SRR5884097/SRR5884097_2.fastq.gz

echo "Download complete"
#
# Make the script executable
chmod +x download_fastq.sh
# Start downloads
nohup ./download_fastq.sh > download.log 2>&1 &

# Check in dir /mnt/active_storage/pradnya/RNASeqMock_GSE102094/HGSC_STIC_4samples/data/raw_fastq

cd /mnt/active_storage/pradnya/RNASeqMock_GSE102094/HGSC_STIC_4samples/data/raw_fastq
ps -ef | grep wget
ls -lh *.fastq.gz

# Monitor 
tail -f download.log
```
### FastQC Analysis
```bash
# Verify all FASTQ files
cd /mnt/active_storage/pradnya/RNASeqMock_GSE102094/HGSC_STIC_4samples/data/raw_fastq

for f in *.fastq.gz
do
    gzip -t "$f" && echo "$f OK"
done

# FastQC using 26 cores
cd /mnt/active_storage/pradnya/RNASeqMock_GSE102094/HGSC_STIC_4samples

fastqc \
-t 26 \
-o qc/fastqc \
data/raw_fastq/*.fastq.gz

# Check FastQC completion
ls qc/fastqc/*.html | wc -l  
# 8
ls qc/fastqc/*.zip | wc -l
# 8
# Generate MultiQC
multiqc qc/fastqc \  
-o qc/multiqc
# Verify
ls -lh qc/multiqc
# 8
# Create alignment script
nano scripts/align_all.sh
#!/bin/bash

THREADS=16

for SAMPLE in \
SRR5884094 \
SRR5884095 \
SRR5884096 \
SRR5884097
do

echo "Processing ${SAMPLE}"

hisat2 \
-p ${THREADS} \
-x hisat2_index/grch38 \
-1 data/raw_fastq/${SAMPLE}_1.fastq.gz \
-2 data/raw_fastq/${SAMPLE}_2.fastq.gz \
2> logs/${SAMPLE}_hisat2.log \
| samtools view -@ ${THREADS} -bS - \
| samtools sort -@ ${THREADS} -o bam/${SAMPLE}.sorted.bam -

samtools index \
-@ ${THREADS} \
bam/${SAMPLE}.sorted.bam

samtools flagstat \
bam/${SAMPLE}.sorted.bam \
> bam/${SAMPLE}.flagstat.txt

echo "${SAMPLE} completed"

done

# Make the script executable
chmod +x scripts/align_all.sh
# Run
nohup scripts/align_all.sh > logs/align_all.log 2>&1 &  
# Monitor
tail -f logs/align_all.log
ls -lh bam/
ls bam/*.sorted.bam 2>/dev/null | wc -l  
# 4
ls bam/*.flagstat.txt 2>/dev/null | wc -l
# 4
ls bam/*.sorted.bam | wc -l
ls bam/*.flagstat.txt | wc -l
```
### Generate Counts Using Rsubread
```bash
# Create the script:
nano scripts/run_featurecounts.R  

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

# Run the script
Rscript scripts/run_featurecounts.R
# Verify
ls -lh counts
```
### Create Metadata File
```bash
nano metadata/sample_metadata.csv
sample,condition  
SRR5884094,STIC_No  
SRR5884096,STIC_No  
SRR5884095,STIC_Yes  
SRR5884097,STIC_Yes

# Verify
cat metadata/sample_metadata.csv
```
### DESeq2 Analysis
```bash
nano scripts/run_deseq2.R
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

# Run the script
Rscript scripts/run_deseq2.R
# Verify the results
ls -lh results
ls -lh counts  
# Check the top genes: 
head counts/gene_counts.csv
```
### Generate the Plots
```bash
conda install -c conda-forge r-pheatmap -y

conda install -c conda-forge \  
r-pheatmap \  
r-ggplot2 \  
r-rcolorbrewer \  
r-ggrepel \  
-y



nano scripts/plot_deseq2.R
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

# Run the script
Rscript scripts/plot_deseq2.R
ls -lh plots

# Create significant gene list
R
res <- read.csv(
  "results/deseq2_results.csv"
)

sig <- subset(
  res,
  !is.na(padj) &
  padj < 0.05
)

nrow(sig)

write.csv(
  sig,
  "results/significant_genes.csv",
  row.names = FALSE
)

head(sig)
q()

ls -lh results
ls -lh plots
```
### Plots with the annotions
```bash
conda install -c bioconda bioconductor-biomart -y

# For annotations of genes from ENSEMBL ID
nano scripts/annotate_genes_orgdb.R
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

# Run the script
Rscript scripts/annotate_genes_orgdb.R

# Verify
ls -lh results/*annotated*
head results/significant_genes_annotated.csv

# Volcano plot with gene labesl
nano scripts/volcano_labeled.R

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

# Run 
Rscript scripts/volcano_labeled.R

# heatmap
nano scripts/heatmap_symbols.R
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

# Run
Rscript scripts/heatmap_symbols.R

# Check the new files
ls -lh plots/*gene*
ls plots/

# To see which genes are on the heatmap
cut -d',' -f8 results/significant_genes_annotated.csv | head -30

# In R
R
sig <- read.csv("results/significant_genes_annotated.csv")  

sig[order(sig$padj), c("SYMBOL","GENENAME","log2FoldChange","padj")][1:20,]

sig <- read.csv(
"results/significant_genes_annotated.csv"
)

top20 <- sig[
order(sig$padj),
][1:20,]

write.csv(
top20,
"results/top20_significant_genes.csv",
row.names = FALSE
)
```
### Pathway Enrichment Analysis
```bash
# Create separate environment for the Pathway enrichment analysis
conda create -n enrichment \
-c conda-forge \
-c bioconda \
r-base=4.5 \
bioconductor-clusterprofiler \
bioconductor-org.hs.eg.db \
bioconductor-enrichplot \
bioconductor-dose \
-y

# Enrichment script
nano scripts/run_enrichment.R

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
# Run Enrichment
cd /mnt/active_storage/pradnya/RNASeqMock_GSE102094/HGSC_STIC_4samples
Rscript scripts/run_enrichment.R  
ls -lh results/*enrichment*

# GO enrichment
nano scripts/go_only.R

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

# Run the script
Rscript scripts/go_only.R
cat results/GO_enrichment.csv
head -20 results/GO_enrichment.csv
```

