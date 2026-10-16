#!/usr/bin/env bash

# Purpose: Build manuscript-facing NetMoss main figure and tables.
# Inputs: Preserved completed NetMoss result tables, matrices, and manifest.
# Outputs: Four-panel figure, two main tables, panel data, captions, and index.
# Usage: bash scripts/run_netmoss_manuscript_synthesis.sh [--config PATH] [--force]

set -euo pipefail

CONFIG="config/netmoss_manuscript_synthesis.yml"
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

read -r -a RSCRIPT_CMD <<<"${NETMOSS_MANUSCRIPT_RSCRIPT:-Rscript}"
ARGS=(--config "$CONFIG")
[[ "$FORCE" == 1 ]] && ARGS+=(--force)

echo "[$(date --iso-8601=seconds)] START NetMoss manuscript synthesis"
"${RSCRIPT_CMD[@]}" scripts/validate_netmoss_manuscript_environment.R
"${RSCRIPT_CMD[@]}" analysis/netmoss_manuscript_synthesis/00_build_main_figure_and_tables_public.R "${ARGS[@]}"
"${RSCRIPT_CMD[@]}" analysis/netmoss_manuscript_synthesis/tests/validate_outputs_public.R --config "$CONFIG"
echo "[$(date --iso-8601=seconds)] COMPLETED NetMoss manuscript synthesis"
