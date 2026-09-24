# HGSC-STIC Bulk RNA-Seq Analysis Pipeline

## Project Overview

This repository contains a complete bulk RNA-seq workflow developed for the analysis of High-Grade Serous Ovarian Carcinoma (HGSC) samples with and without associated Serous Tubal Intraepithelial Carcinoma (STIC).

The workflow covers:

- Data retrieval from GEO/SRA
- FASTQ quality assessment
- Reference genome preparation
- Genome alignment
- Read quantification
- Differential expression analysis
- Gene annotation
- Functional enrichment analysis
- Visualization

---

## Dataset

### Study

Ducie J, Dao F, Considine M, Olvera N, et al.

*Molecular analysis of high-grade serous ovarian carcinoma with and without associated serous tubal intra-epithelial carcinoma.*

Nature Communications (2017)

GEO SuperSeries: **GSE102094**

### Samples Used

| GEO Sample | HGSC ID | STIC Status | SRR |
|-----------|----------|-------------|------|
| GSM2719894 | HGSC_40187 | STIC Negative | SRR5884094 |
| GSM2719895 | HGSC_40190 | STIC Positive | SRR5884095 |
| GSM2719896 | HGSC_40193 | STIC Negative | SRR5884096 |
| GSM2719897 | HGSC_40197 | STIC Positive | SRR5884097 |

---

## Software

All software was installed through Conda.

### Core Software

- Miniconda
- HISAT2
- SAMtools
- FastQC
- MultiQC
- Rsubread
- DESeq2
- clusterProfiler
- org.Hs.eg.db

---

## Reference Resources

### Genome

GRCh38 Primary Assembly

### Annotation

GENCODE Release 48

---

## Workflow

```text
FASTQ
 ↓
FastQC
 ↓
MultiQC
 ↓
HISAT2 Alignment
 ↓
SAMtools Sorting & Indexing
 ↓
featureCounts
 ↓
DESeq2
 ↓
Gene Annotation
 ↓
GO Enrichment Analysis
```

---

## Quality Control

Quality assessment was performed using:

- FastQC
- MultiQC

Metrics assessed included:

- Per-base quality
- GC content
- Adapter contamination
- Sequence duplication levels

---

## Alignment

Reads were aligned against the GRCh38 human reference genome using HISAT2.

Representative alignment statistics:

- Overall mapping rate > 90%
- Properly paired reads > 90%

---

## Quantification

Read quantification was performed using featureCounts from the Rsubread package.

Output:

```text
gene_counts.csv
featureCounts.rds
```

---

## Differential Expression Analysis

Differential expression analysis was performed using DESeq2.

Experimental design:

```r
design = ~ condition
```

Conditions:

- STIC_Positive
- STIC_Negative

Generated files:

```text
dds.rds
deseq2_results.csv
deseq2_results_annotated.csv
significant_genes_annotated.csv
```

---

## Visualization

The workflow generates:

- PCA Plot
- MA Plot
- Volcano Plot
- Labeled Volcano Plot
- Heatmap
- Heatmap with Gene Symbols
- Boxplot

All figures are exported as high-resolution TIFF files.

---

## Functional Enrichment Analysis

Gene Ontology enrichment was performed using clusterProfiler.

Enriched biological processes included:

- B cell mediated immunity
- Immunoglobulin mediated immune response
- Lymphocyte mediated immunity
- Adaptive immune response

Generated files:

```text
GO_enrichment.csv
GO_dotplot.tiff
GO_barplot.tiff
```

---

## Repository Contents

### Python Scripts

- extract_samples.py
- extract_first4.py

### Bash Scripts

- download_fastq.sh
- align_all.sh

### R Scripts

- run_featurecounts.R
- run_deseq2.R
- plot_deseq2.R
- annotate_genes_orgdb.R
- volcano_labeled.R
- heatmap_symbols.R
- run_enrichment.R
- go_only.R

---

## Author

Pradnya Kamble

Zifo RnD Solutions