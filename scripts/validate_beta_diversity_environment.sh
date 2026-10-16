#!/usr/bin/env bash
# Purpose: Fail early when required beta-diversity executables or R packages are
# missing from the active environment.
# Inputs: The current PATH and installed R/Python runtime packages.
# Outputs: Validation messages on stdout/stderr; no project files are written.
# Dependencies: Bash, Rscript, and python3.
# Usage: bash scripts/validate_beta_diversity_environment.sh
# Workflow: Environment preflight used locally and by HPC update/run jobs.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
command -v Rscript >/dev/null || { echo "FAILED: Rscript not found" >&2; exit 1; }
command -v python3 >/dev/null || { echo "FAILED: python3 not found" >&2; exit 1; }
Rscript -e 'required <- c("vegan","ape","permute","ggplot2","scales","yaml","digest","sessioninfo","testthat","lme4"); missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]; if(length(missing)) stop("Missing R packages: ", paste(missing, collapse=", ")); cat("R environment valid\n")'
python3 -c 'import csv, hashlib, pathlib; print("Python standard-library utilities valid")'
echo "COMPLETED_ENVIRONMENT_VALIDATION"
