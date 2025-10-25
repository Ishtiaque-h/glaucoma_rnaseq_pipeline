#!/bin/bash
#SBATCH -J salmon_index
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=512G
#SBATCH -t 72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ishtiaq@udel.edu

set -euo pipefail
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

PROJECT="$HOME/glaucoma_rnaseq_pipeline"
TRIM="$PROJECT/data/trimmed_fast"
IDX="$PROJECT/data/reference/salmon_index"
OUT="$PROJECT/results/salmon"
SIF="$PROJECT/containers/salmon_1.10.3.sif" 

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TRIM/samples.txt")
R1="$TRIM/${SAMPLE}_1.trimmed.fastq.gz"
R2="$TRIM/${SAMPLE}_2.trimmed.fastq.gz"

WORKDIR="/work/ishtiaq/${SLURM_JOB_ID}_salmon"
mkdir -p "$WORKDIR"
cp -r "$IDX" "$WORKDIR/"
cp "$R1" "$WORKDIR/"
cp "$R2" "$WORKDIR/"

singularity exec -e --bind "$WORKDIR:$WORKDIR" "$SIF" salmon quant \
  -i "$WORKDIR/salmon_index" \
  -l A \
  -1 "$WORKDIR/$(basename "$R1")" -2 "$WORKDIR/$(basename "$R2")" \
  -p "${SLURM_CPUS_PER_TASK:-8}" \
  --validateMappings --seqBias --gcBias \
  --numBootstraps 30 \
  -o "$WORKDIR/${SAMPLE}"

rsync -av "$WORKDIR/${SAMPLE}" "$OUT/"
rm -rf "$WORKDIR"
