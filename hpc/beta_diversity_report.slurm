#!/usr/bin/env bash
#SBATCH --job-name=blackice-beta-report
#SBATCH --partition=compute
#SBATCH --account=CHANGE_ME
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/beta-report-%j.out
#SBATCH --error=logs/beta-report-%j.err

set -euo pipefail
set -x
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$ROOT"
mkdir -p logs
SCRATCH="${SLURM_TMPDIR:-/tmp/${USER}/blackice-beta-report-${SLURM_JOB_ID:-local}}"
mkdir -p "$SCRATCH"
cleanup() { status=$?; rm -rf "$SCRATCH"; if (( status != 0 )); then echo "FAILED beta-diversity report stage" >&2; fi; exit "$status"; }
trap cleanup EXIT INT TERM
module purge || true
# module load R/4.3
if command -v conda >/dev/null 2>&1; then source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_NAME:-glmm-hpc}"; fi
CONFIG="${BETA_CONFIG:-config/beta_diversity.yml}"
Rscript analysis/beta_diversity/04_run_ordination.R --config "$CONFIG"
Rscript analysis/beta_diversity/05_run_permanova.R --config "$CONFIG"
Rscript analysis/beta_diversity/06_run_longitudinal_analysis.R --config "$CONFIG"
Rscript analysis/beta_diversity/07_make_figures.R --config "$CONFIG"
Rscript analysis/beta_diversity/08_build_report.R --config "$CONFIG"
