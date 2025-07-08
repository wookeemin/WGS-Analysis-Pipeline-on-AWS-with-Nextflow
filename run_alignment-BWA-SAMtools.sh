#!/bin/bash

# User definition
THREADS=32
SAMPLE="sample"  # prefix of file name
REF="ref.fa"     # Reference genomics
FASTQ1="${SAMPLE}_1.fastq.gz"
FASTQ2="${SAMPLE}_2.fastq.gz"
OUTPUT="aligned.sorted.bam"

# 1. execution of alignment
echo "[INFO] Running BWA alignment..."
bwa mem -t ${THREADS} ${REF} ${FASTQ1} ${FASTQ2} | \
  samtools view -@ 8 -Sb - > aligned.bam

# 2. Alig & Indexing
echo "[INFO] Sorting BAM file..."
samtools sort -@ 8 aligned.bam -o ${OUTPUT}

echo "[INFO] Indexing BAM file..."
samtools index ${OUTPUT}
#3. intermediate file remove (optional)
rm aligned.bam

# 4. Print summary output
echo "[INFO] Alignment stats:"
samtools flagstat ${OUTPUT}
