#!/usr/bin/env bash
# Purpose: Build the manuscript-facing functional synthesis figure and tables.
# Inputs: Completed MAG-functional-annotation and functional-dynamics results.
# Outputs: Main Figure 1, Main Tables 1-2, panel data, and a completion marker.
# Usage: bash scripts/run_manuscript_functional_synthesis.sh \
#   --config config/manuscript_functional_synthesis.yml

set -euo pipefail

CONFIG="config/manuscript_functional_synthesis.yml"
FORCE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
test -f "$CONFIG" || { echo "Configuration not found: $CONFIG" >&2; exit 1; }

RSCRIPT_BIN="${MANSYN_RSCRIPT:-Rscript}"
ARGS=(--config "$CONFIG")
[[ "$FORCE" == 1 ]] && ARGS+=(--force)

echo "[$(date --iso-8601=seconds)] START manuscript functional synthesis"
"$RSCRIPT_BIN" analysis/manuscript_functional_synthesis/00_build_main_figure_and_tables_public.R "${ARGS[@]}"
"$RSCRIPT_BIN" analysis/manuscript_functional_synthesis/tests/validate_outputs_public.R --config "$CONFIG"
echo "[$(date --iso-8601=seconds)] COMPLETED manuscript functional synthesis"
