#!/usr/bin/env bash
# Purpose: Orchestrate the staged beta-diversity workflow from validation
# through tests and reporting on a local machine.
# Inputs: A beta-diversity YAML file; defaults to config/beta_diversity.yml.
# Outputs: Workflow-specific files under the configured output directory.
# Dependencies: Bash, Rscript, environment.yml packages, and analysis/beta_diversity/.
# Usage: bash scripts/run_beta_diversity_local.sh --config config/beta_diversity.yml
# Workflow: Canonical local runner; supports --force, --metric, and --from-stage.

set -euo pipefail

CONFIG="config/beta_diversity.yml"
FORCE=0
METRIC=""
FROM_STAGE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --metric) METRIC="$2"; shift 2 ;;
    --from-stage) FROM_STAGE="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
ARGS=(--config "$CONFIG")
[[ "$FORCE" == 1 ]] && ARGS+=(--force)

run_stage() {
  local number="$1" script="$2"
  if (( number >= FROM_STAGE )); then
    echo "[$(date --iso-8601=seconds)] START $script"
    Rscript "$script" "${ARGS[@]}"
  fi
}

run_stage 1 analysis/beta_diversity/01_validate_inputs.R
run_stage 2 analysis/beta_diversity/02_prepare_data.R
if (( 3 >= FROM_STAGE )); then
  DIST_ARGS=("${ARGS[@]}")
  [[ -n "$METRIC" ]] && DIST_ARGS+=(--metric "$METRIC")
  Rscript analysis/beta_diversity/03_compute_distances.R "${DIST_ARGS[@]}"
fi
run_stage 4 analysis/beta_diversity/04_run_ordination.R
run_stage 5 analysis/beta_diversity/05_run_permanova.R
run_stage 6 analysis/beta_diversity/06_run_longitudinal_analysis.R
run_stage 7 analysis/beta_diversity/07_make_figures.R
run_stage 8 analysis/beta_diversity/08_build_report.R

if [[ "$CONFIG" == *test* ]]; then
  Rscript analysis/beta_diversity/tests/testthat.R
fi
echo "[$(date --iso-8601=seconds)] COMPLETED beta-diversity workflow"
