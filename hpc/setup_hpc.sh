#!/bin/bash
# ============================================
# setup_hpc.sh — Run ONCE when first accessing HPC
# ============================================
# This script:
#   1. Clones your GitHub repo
#   2. Creates Conda environment
#   3. Verifies setup
#
# Usage: bash hpc/setup_hpc.sh
# ============================================

set -euo pipefail

# --- Configuration ---
GITHUB_URL="${GITHUB_URL:-}"  # Supply the approved public repository URL when cloning.
PROJECT_DIR="${PROJECT_DIR:-$HOME/projects}"          # or path to 1TB project dir
REPO_DIR_NAME="${REPO_DIR_NAME:-sagarmatha-microbiome-analysis}"
CONDA_ENV="${CONDA_ENV_NAME:-glmm-hpc}"

echo "=== Setting up HPC environment ==="

# --- Create project directory ---
mkdir -p "$PROJECT_DIR"

# --- Clone repository, or use the current checkout if already inside it. ---
if git rev-parse --show-toplevel >/dev/null 2>&1; then
    cd "$(git rev-parse --show-toplevel)"
    echo "Using existing repository checkout."
elif [ -d "$PROJECT_DIR/$REPO_DIR_NAME" ]; then
    echo "Repo exists, pulling latest..."
    cd "$PROJECT_DIR/$REPO_DIR_NAME"
    git pull --ff-only
else
    echo "Cloning repository..."
    cd "$PROJECT_DIR"
    : "${GITHUB_URL:?Set GITHUB_URL to the approved public repository URL}"
    git clone "$GITHUB_URL" "$REPO_DIR_NAME"
    cd "$REPO_DIR_NAME"
fi

# --- Create Conda environment ---
if conda env list | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
    echo "Environment '$CONDA_ENV' already exists."
else
    echo "Creating Conda environment '$CONDA_ENV'..."
    conda env create -f environment.yml
fi

# --- Verify ---
mkdir -p data/raw data/processed logs results/hpc_per_mag_lmertest tmp

echo ""
echo "=== Setup Complete ==="
echo "Repo:       $(pwd)"
echo "Branch:     $(git branch --show-current)"
echo "Environment: $CONDA_ENV"
echo "R version:  $(conda run -n $CONDA_ENV R --version | head -1)"
echo ""
echo "Next steps:"
echo "  1. Place sample_info.csv and the selected abundance CSV in: data/raw/"
echo "     Or set BLACKICE_RAW_DIR=/path/to/raw_data before running scripts/01_data_prep.R."
echo "  2. Prepare data: conda run -n $CONDA_ENV Rscript scripts/01_data_prep.R"
echo "  3. Submit a small OpenPBS test: qsub -v RUN_MODE=dryrun hpc/job_per_mag_lmertest.pbs"
echo "  4. Submit the full run after the dry run succeeds: qsub hpc/job_per_mag_lmertest.pbs"
