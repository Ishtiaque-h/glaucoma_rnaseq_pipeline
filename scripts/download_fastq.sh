#!/bin/bash

INPUT_CSV=$1
OUTDIR="data/raw_fastq"
TMPDIR="tmp"
mkdir -p $OUTDIR
mkdir -p $TMPDIR

# Extract SRR IDs and download FASTQs
cut -d',' -f1 "$INPUT_CSV" | grep SRR | while read SRA_ID; do
    echo "Downloading $SRA_ID ..."
    
    # Use custom temp folder to avoid disk-limit errors
    fasterq-dump "$SRA_ID" -O "$OUTDIR" --split-files --threads 4 --temp "$TMPDIR"

    # Compress the output FASTQ files
    pigz "$OUTDIR/${SRA_ID}_1.fastq"
    pigz "$OUTDIR/${SRA_ID}_2.fastq"

    # Clean up temp files after each download
    rm -rf "$TMPDIR"/*
done
