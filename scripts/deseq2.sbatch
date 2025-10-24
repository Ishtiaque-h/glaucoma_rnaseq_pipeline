#!/bin/bash
#SBATCH -J deseq2
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ishtiaq@udel.edu

set -euo pipefail

singularity exec --bind /work/,/scratch  /work/ishtiaq/glaucoma_rnaseq_pipeline/containers/bioc_3.21.sif Rscript scripts/tximport_deseq2.R
