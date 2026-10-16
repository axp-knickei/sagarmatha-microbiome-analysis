#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$ROOT"
SMOKE_ROOT="${1:-$(mktemp -d /tmp/blackice_funcdyn_smoke.XXXXXX)}"
R_BIN="${FUNCDYN_RSCRIPT:-Rscript}"
SMOKE_CONFIG="$($R_BIN analysis/functional_community_dynamics/05_create_smoke_test_data.R --output-dir "$SMOKE_ROOT" | tail -1)"
FUNCDYN_RSCRIPT="$R_BIN" FUNCDYN_PERMUTATIONS=99 bash scripts/run_functional_community_dynamics.sh --config "$SMOKE_CONFIG"
$R_BIN analysis/functional_community_dynamics/run_tests.R
required=(
  "$SMOKE_ROOT/results/01_prepared_profiles/Table_26_functional_profile_input_and_filtering_audit.csv"
  "$SMOKE_ROOT/results/02_multivariate_analysis/Table_27_global_PERMANOVA_and_dispersion_results.csv"
  "$SMOKE_ROOT/results/02_multivariate_analysis/Table_28_supported_functional_timepoint_contrasts.csv"
  "$SMOKE_ROOT/results/03_longitudinal_stability/Table_29_participant_distance_reference_and_adjacent_transitions.csv"
  "$SMOKE_ROOT/results/03_longitudinal_stability/Table_30_matched_taxonomic_functional_stability_models.csv"
  "$SMOKE_ROOT/results/03_longitudinal_stability/Table_31_participant_level_transition_sensitivity.csv"
  "$SMOKE_ROOT/results/figures/Figure_16_functional_PCoA_participant_trajectories.pdf"
  "$SMOKE_ROOT/results/figures/Figure_17_distance_to_earliest_observation.pdf"
  "$SMOKE_ROOT/results/figures/Figure_18_taxonomic_vs_functional_transitions.pdf"
  "$SMOKE_ROOT/results/figures/Figure_19_exploratory_T1_T10_recovery.pdf"
  "$SMOKE_ROOT/results/reports/2026-08-12_functional_community_dynamics_interpretation.md"
  "$SMOKE_ROOT/results/00_quality_control/04_make_figures_and_report.complete"
  "$SMOKE_ROOT/results/RESULTS_INDEX.md"
)
for path in "${required[@]}"; do test -s "$path" || { echo "Missing smoke output: $path" >&2; exit 1; }; done
$R_BIN analysis/functional_community_dynamics/tests/validate-smoke-outputs.R "$SMOKE_ROOT/results"
echo "Functional community dynamics smoke test passed: $SMOKE_ROOT"
