#!/bin/bash

# Set number of threads per sample
THREADS=8

# Read sample list and run 1 at a time in parallel
cat samples.txt | parallel -j 1 --line-buffer '
  SAMPLE={}
  mkdir -p results/salmon_quant/${SAMPLE}

  salmon quant \
    -i references/gencode/salmon_index \
    -l A \
    -1 data/trimmed_fastq/${SAMPLE}_1.trimmed.fastq.gz \
    -2 data/trimmed_fastq/${SAMPLE}_2.trimmed.fastq.gz \
    -p '"$THREADS"' \
    --validateMappings \
    -o results/salmon_quant/${SAMPLE} \
    2>&1 | tee results/salmon_quant/${SAMPLE}/salmon_stdout.log

'

