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
