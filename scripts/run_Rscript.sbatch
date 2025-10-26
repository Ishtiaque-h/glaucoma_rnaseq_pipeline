#!/bin/bash
#SBATCH -J Rjob
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ishtiaq@udel.edu

set -euo pipefail
module purge

# Load container runtime (adjust to your cluster’s module name)
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

export PROJECT="/work/ishtiaq/glaucoma_rnaseq_pipeline"

# Point to your SIF image (adjust the filename to the actual one you have)
IMAGE="$PROJECT/containers/bioc_3.21.sif"

# Bind PROJECT into the container so paths resolve the same inside & outside
# Add /scratch if you want to write there as well
BIND_OPTS="--bind $PROJECT:$PROJECT"

# Run the R script passed as argument(s)
# Example submission:
#   sbatch scripts/run_Rscript.sbatch scripts/02_export_matrices.R
apptainer exec --cleanenv $BIND_OPTS "$IMAGE" Rscript "$@"

echo "Done!"
