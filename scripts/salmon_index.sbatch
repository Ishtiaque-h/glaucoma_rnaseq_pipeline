#!/bin/bash
#SBATCH -J salmon_index
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ishtiaq@udel.edu

set -euo pipefail
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

PROJECT="$HOME/glaucoma_rnaseq_pipeline"
REF="$PROJECT/data/reference"
REF_FASTA="$REF/mmus_transcripts.fa.gz"
IDX="$REF/salmon_index"
SIF="$PROJECT/containers/salmon_1.10.3.sif" 

# Build index inside the container, writing to the bound host paths
singularity exec -e --bind "$PROJECT:$PROJECT" "$SIF" \
  salmon index -t "$REF_FASTA" -i "$IDX" -k 31

echo "[DONE] Index built at $IDX"

