#!/bin/bash

INPUT_DIR="data/raw_fastq"
OUTPUT_DIR="results/qc/fastqc_reports"

mkdir -p "$OUTPUT_DIR"

# Run FastQC in parallel on all .fastq.gz files
for file in "$INPUT_DIR"/*.fastq.gz; do
    echo "Running FastQC on $file"
    fastqc "$file" -o "$OUTPUT_DIR" --threads 4
done

