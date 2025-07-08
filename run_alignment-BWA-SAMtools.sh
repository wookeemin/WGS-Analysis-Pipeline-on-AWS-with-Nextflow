#!/bin/bash

# 사용자 설정
THREADS=32
SAMPLE="sample"  # 파일명 접두어
REF="ref.fa"     # 기준 유전체
FASTQ1="${SAMPLE}_1.fastq.gz"
FASTQ2="${SAMPLE}_2.fastq.gz"
OUTPUT="aligned.sorted.bam"

# 1. 정렬 수행
echo "[INFO] Running BWA alignment..."
bwa mem -t ${THREADS} ${REF} ${FASTQ1} ${FASTQ2} | \
  samtools view -@ 8 -Sb - > aligned.bam

# 2. 정렬 및 인덱싱
echo "[INFO] Sorting BAM file..."
samtools sort -@ 8 aligned.bam -o ${OUTPUT}

echo "[INFO] Indexing BAM file..."
samtools index ${OUTPUT}
#3. 중간 파일 제거 (선택)
rm aligned.bam

# 4. 요약 통계 출력
echo "[INFO] Alignment stats:"
samtools flagstat ${OUTPUT}
