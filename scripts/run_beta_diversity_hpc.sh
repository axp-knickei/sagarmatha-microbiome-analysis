#!/usr/bin/env bash
# Purpose: Apply scheduler-aware thread limits and delegate beta-diversity work
# to the canonical local runner.
# Inputs: BETA_CONFIG/BETA_THREADS and arguments accepted by the local runner.
# Outputs: The configured beta-diversity output tree.
# Dependencies: Bash, the local runner, and an activated workflow environment.
# Usage: bash scripts/run_beta_diversity_hpc.sh
# Workflow: Called by beta-diversity scheduler wrappers.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-${BETA_THREADS:-1}}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
CONFIG="${BETA_CONFIG:-config/beta_diversity.yml}"
exec bash scripts/run_beta_diversity_local.sh --config "$CONFIG" "$@"
