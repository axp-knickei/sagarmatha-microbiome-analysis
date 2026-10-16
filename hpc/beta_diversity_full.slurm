#!/usr/bin/env bash
#SBATCH --job-name=blackice-beta-full
#SBATCH --partition=compute
#SBATCH --account=CHANGE_ME
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/beta-full-%j.out
#SBATCH --error=logs/beta-full-%j.err
##SBATCH --mail-type=FAIL
##SBATCH --mail-user=user@example.org

set -euo pipefail
set -x
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$ROOT"
mkdir -p logs
SCRATCH="${SLURM_TMPDIR:-/tmp/${USER}/blackice-beta-${SLURM_JOB_ID:-local}}"
mkdir -p "$SCRATCH"
cleanup() { status=$?; rm -rf "$SCRATCH"; if (( status != 0 )); then echo "FAILED beta-diversity full workflow" >&2; fi; exit "$status"; }
trap cleanup EXIT INT TERM
module purge || true
# module load R/4.3
if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${CONDA_ENV_NAME:-glmm-hpc}"
fi
export BETA_CONFIG="${BETA_CONFIG:-config/beta_diversity.yml}"
export BETA_SEED="${BETA_SEED:-20260711}"
echo "host=$(hostname) job=${SLURM_JOB_ID:-local} commit=$(git rev-parse HEAD) seed=$BETA_SEED"
bash scripts/validate_beta_diversity_environment.sh
bash scripts/run_beta_diversity_hpc.sh
