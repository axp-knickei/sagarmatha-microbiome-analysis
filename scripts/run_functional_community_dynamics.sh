#!/usr/bin/env bash
# Purpose: Run the independent functional-community-dynamics workflow.
# Inputs: MAG abundance, metadata, and Stage-03 predicted functional-potential files.
# Outputs: Tables 26-30, Figures 16-19, reproducibility files, and a dated report.
# Usage: bash scripts/run_functional_community_dynamics.sh --config config/functional_community_dynamics.yml

set -euo pipefail
CONFIG="config/functional_community_dynamics.yml"
FROM_STAGE=0
TO_STAGE=4
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
[[ "$FROM_STAGE" =~ ^[0-4]$ ]] || { echo "--from-stage must be 0-4" >&2; exit 2; }
[[ "$TO_STAGE" =~ ^[0-4]$ ]] || { echo "--to-stage must be 0-4" >&2; exit 2; }
(( FROM_STAGE <= TO_STAGE )) || { echo "--from-stage cannot exceed --to-stage" >&2; exit 2; }

RSCRIPT_BIN="${FUNCDYN_RSCRIPT:-Rscript}"
ARGS=(--config "$CONFIG")
[[ "$FORCE" == 1 ]] && ARGS+=(--force)
stages=(
  analysis/functional_community_dynamics/00_validate_inputs.R
  analysis/functional_community_dynamics/01_prepare_profiles.R
  analysis/functional_community_dynamics/02_run_multivariate_analysis.R
  analysis/functional_community_dynamics/03_analyze_longitudinal_stability.R
  analysis/functional_community_dynamics/04_make_figures_and_report.R
)
for stage in "${!stages[@]}"; do
  if (( stage >= FROM_STAGE && stage <= TO_STAGE )); then
    echo "[$(date --iso-8601=seconds)] START stage $stage: ${stages[$stage]}"
    "$RSCRIPT_BIN" "${stages[$stage]}" "${ARGS[@]}"
    echo "[$(date --iso-8601=seconds)] END stage $stage"
  fi
done
echo "[$(date --iso-8601=seconds)] COMPLETED functional community dynamics workflow"
