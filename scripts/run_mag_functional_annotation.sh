#!/usr/bin/env bash
# Purpose: Run the staged MAG functional annotation workflow locally or on HPC.
# Inputs: Production or smoke-test YAML configuration.
# Outputs: Numbered QC, tables, figures, intermediates, and RESULTS_INDEX files.
# Dependencies: Bash, Rscript, environment.yml packages, and analysis/mag_functional_annotation/.
# Usage: bash scripts/run_mag_functional_annotation.sh --config config/mag_functional_annotation.yml
# Workflow: Canonical orchestrator; supports restartable --from-stage/--to-stage and intentional --force.

set -euo pipefail

CONFIG="config/mag_functional_annotation.yml"
FROM_STAGE=0
TO_STAGE=10
FORCE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --from-stage) FROM_STAGE="$2"; shift 2 ;;
    --to-stage) TO_STAGE="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
test -f "$CONFIG" || { echo "Configuration not found: $CONFIG" >&2; exit 1; }
[[ "$FROM_STAGE" =~ ^([0-9]|10)$ ]] || { echo "--from-stage must be 0-10" >&2; exit 2; }
[[ "$TO_STAGE" =~ ^([0-9]|10)$ ]] || { echo "--to-stage must be 0-10" >&2; exit 2; }
(( FROM_STAGE <= TO_STAGE )) || { echo "--from-stage cannot exceed --to-stage" >&2; exit 2; }

RSCRIPT_BIN="${MAGFUNC_RSCRIPT:-Rscript}"
ARGS=(--config "$CONFIG")
[[ "$FORCE" == 1 ]] && ARGS+=(--force)

stages=(
  analysis/mag_functional_annotation/00_validate_inputs.R
  analysis/mag_functional_annotation/01_prepare_catalogue.R
  analysis/mag_functional_annotation/02_catalogue_analyses.R
  analysis/mag_functional_annotation/03_longitudinal_analyses.R
  analysis/mag_functional_annotation/04_make_publication_outputs.R
  analysis/mag_functional_annotation/05_build_results_index.R
  analysis/mag_functional_annotation/06_validate_selected_candidates.R
  analysis/mag_functional_annotation/07_explore_pathway_unknowns.R
  analysis/mag_functional_annotation/08_decompose_shared_unique_carriers.R
  analysis/mag_functional_annotation/09_analyze_shared_component_drivers.R
  analysis/mag_functional_annotation/10_audit_candidate_annotations.R
)

for stage in "${!stages[@]}"; do
  if (( stage >= FROM_STAGE && stage <= TO_STAGE )); then
    echo "[$(date --iso-8601=seconds)] START stage $stage: ${stages[$stage]}"
    "$RSCRIPT_BIN" "${stages[$stage]}" "${ARGS[@]}"
    echo "[$(date --iso-8601=seconds)] END stage $stage"
  fi
done
echo "[$(date --iso-8601=seconds)] COMPLETED MAG functional annotation workflow"
