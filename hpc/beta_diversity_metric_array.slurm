#!/usr/bin/env bash
#SBATCH --job-name=blackice-beta-metric
#SBATCH --partition=compute
#SBATCH --account=CHANGE_ME
#SBATCH --array=0-4
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/beta-metric-%A_%a.out
#SBATCH --error=logs/beta-metric-%A_%a.err

set -euo pipefail
set -x
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$ROOT"
mkdir -p logs
METRICS=(bray_curtis jaccard_binary quantitative_jaccard_ruzička euclidean_hellinger aitchison)
INDEX="${SLURM_ARRAY_TASK_ID:-0}"
METRIC="${METRICS[$INDEX]}"
SCRATCH="${SLURM_TMPDIR:-/tmp/${USER}/blackice-beta-${SLURM_JOB_ID:-local}-${INDEX}}"
mkdir -p "$SCRATCH"
cleanup() { status=$?; rm -rf "$SCRATCH"; if (( status != 0 )); then echo "FAILED metric=$METRIC" >&2; fi; exit "$status"; }
trap cleanup EXIT INT TERM
module purge || true
# module load R/4.3
if command -v conda >/dev/null 2>&1; then source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_NAME:-glmm-hpc}"; fi
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
Rscript analysis/beta_diversity/03_compute_distances.R --config "${BETA_CONFIG:-config/beta_diversity.yml}" --metric "$METRIC"
