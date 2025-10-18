#!/bin/bash

# Input/output directories
RAW_DIR="data/raw_fastq"
TRIMMED_DIR="data/trimmed_fastq"
REPORT_DIR="results/qc/fastp_reports"

# Create output directories
mkdir -p "$TRIMMED_DIR"
mkdir -p "$REPORT_DIR"

# Loop through all _1.fastq.gz files
for R1 in "$RAW_DIR"/*_1.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="$RAW_DIR/${SAMPLE}_2.fastq.gz"

    echo "Processing $SAMPLE..."

    fastp \
      --in1 "$R1" --in2 "$R2" \
      --out1 "$TRIMMED_DIR/${SAMPLE}_1.trimmed.fastq.gz" \
      --out2 "$TRIMMED_DIR/${SAMPLE}_2.trimmed.fastq.gz" \
      --detect_adapter_for_pe \
      --trim_front1 12 --trim_front2 12 \
      --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
      --length_required 30 \
      --thread 4 \
      --html "$REPORT_DIR/${SAMPLE}_fastp.html" \
      --json "$REPORT_DIR/${SAMPLE}_fastp.json"
done

