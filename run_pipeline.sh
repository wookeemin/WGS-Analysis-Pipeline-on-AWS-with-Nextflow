#!/bin/bash
# Usage: ./run_pipeline.sh [elprep|gatk] /path/to/input.bam /path/to/output_dir

TOOL=$1
INPUT_BAM=$2
OUTPUT_DIR=$3

if [ -z "$TOOL" ] || [ -z "$INPUT_BAM" ] || [ -z "$OUTPUT_DIR" ]; then
  echo "Usage: $0 [elprep|gatk] /path/to/input.bam /path/to/output_dir"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

if [ "$TOOL" == "elprep" ]; then
  echo "Running elPrep pipeline..."
  elprep sfm "$INPUT_BAM" "$OUTPUT_DIR/output.bam" --mark-duplicates --sorting-order coordinate
elif [ "$TOOL" == "gatk" ]; then
  echo "Running GATK pipeline..."
  gatk MarkDuplicates -I "$INPUT_BAM" -O "$OUTPUT_DIR/dedup.bam" -M "$OUTPUT_DIR/metrics.txt"
  gatk SortSam -I "$OUTPUT_DIR/dedup.bam" -O "$OUTPUT_DIR/sorted.bam" --SORT_ORDER coordinate
else
  echo "Invalid tool specified. Use 'elprep' or 'gatk'."
  exit 1
fi

echo "Pipeline execution complete."
