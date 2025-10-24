#!/bin/bash
#SBATCH -J salmon_quant
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=512G
#SBATCH -t 72:00:00
#SBATCH --tmp=100000                 # ~100 GB local scratch requested
#SBATCH --array=1-20                 # submit with %N throttling at the CLI
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ishtiaq@udel.edu

set -euo pipefail
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

PROJECT="$HOME/glaucoma_rnaseq_pipeline"
TRIM="$PROJECT/data/trimmed_fastq"
IDX_SRC="$PROJECT/data/reference/salmon_index"
OUT="$PROJECT/results/salmon"
SIF="$PROJECT/containers/salmon_1.10.3.sif"     # <- your updated SIF

# --------------------------
# Pick fastest scratch space
# --------------------------
if [ -n "${SLURM_TMPDIR:-}" ] && [ -d "$SLURM_TMPDIR" ] && [ -w "$SLURM_TMPDIR" ]; then
  JOB_SCRATCH="$SLURM_TMPDIR"
elif [ -d "/scratch" ] && [ -w "/scratch" ]; then
  JOB_SCRATCH="/scratch/${USER}/${SLURM_JOB_ID}"
  mkdir -p "$JOB_SCRATCH"
else
  JOB_SCRATCH="/work/${USER}/${SLURM_JOB_ID}_salmon"
  mkdir -p "$JOB_SCRATCH"
fi
echo "[INFO] Using scratch: $JOB_SCRATCH"

# Cleanup on success or failure (safe with set -u)
trap '[[ -n "${STAGE:-}" ]] && rm -rf "$STAGE"' EXIT

# Keep container caches off /home
export SINGULARITY_TMPDIR="$JOB_SCRATCH"
export SINGULARITY_CACHEDIR="$JOB_SCRATCH/.sif_cache"
mkdir -p "$SINGULARITY_CACHEDIR"

# --------------------------
# Resolve sample
# --------------------------
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TRIM/samples.txt")
R1="$TRIM/${SAMPLE}_1.trimmed.fastq.gz"
R2="$TRIM/${SAMPLE}_2.trimmed.fastq.gz"

# --------------------------
# Stage-in to scratch
# --------------------------
STAGE="$JOB_SCRATCH/${SAMPLE}"
mkdir -p "$STAGE"

# Copy reads (compressed) to scratch
cp "$R1" "$R2" "$STAGE/"

# Copy index to scratch (fastest). If scratch is tight, comment the next two lines
# and set IDX="$IDX_SRC" to read from /work instead.
cp -r "$IDX_SRC" "$STAGE/salmon_index"
IDX="$STAGE/salmon_index"

# If not copying index, use:
# IDX="$IDX_SRC"

# --------------------------
# Run Salmon on scratch
# --------------------------
singularity exec -e --bind "$JOB_SCRATCH:$JOB_SCRATCH" "$SIF" salmon quant \
  -i "$IDX" \
  -l A \
  -1 "$STAGE/$(basename "$R1")" \
  -2 "$STAGE/$(basename "$R2")" \
  -p "${SLURM_CPUS_PER_TASK:-8}" \
  --validateMappings \
  --seqBias --gcBias \
  --numBootstraps 30 \
  -o "$STAGE/${SAMPLE}"

# --------------------------
# Stage-out 
# --------------------------
rsync -a "$STAGE/${SAMPLE}" "$OUT/"


echo "[DONE] ${SAMPLE} → $OUT/${SAMPLE}"
