#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

Rscript scripts/hpc_run_per_mag_lmertest.R --max-mags "${MAX_MAGS:-10}" "$@"
