#!/usr/bin/env bash
# End-to-end synthetic smoke test; never reads or writes scientific results.

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$ROOT"
SMOKE_ROOT="${1:-$(mktemp -d /tmp/blackice_magfunc_smoke.XXXXXX)}"
R_BIN="${MAGFUNC_RSCRIPT:-Rscript}"

SMOKE_CONFIG="$($R_BIN analysis/mag_functional_annotation/06_create_smoke_test_data.R --output-dir "$SMOKE_ROOT" | tail -1)"
MAGFUNC_RSCRIPT="$R_BIN" bash scripts/run_mag_functional_annotation.sh --config "$SMOKE_CONFIG"
$R_BIN analysis/mag_functional_annotation/run_tests.R

required=(
  "$SMOKE_ROOT/results/00_quality_control/Table_02_full_catalogue_qc_summary.csv"
  "$SMOKE_ROOT/results/01_functional_convergence/Table_03_functional_convergence_all_ontologies.csv.gz"
  "$SMOKE_ROOT/results/02_unknown_by_cluster_size/Table_04_unknown_annotation_by_cluster_size.csv"
  "$SMOKE_ROOT/results/02_unknown_by_cluster_size/unknown_annotation_model_status.csv"
  "$SMOKE_ROOT/results/03_taxonomic_distribution/Table_05_taxonomic_function_distribution_all_ontologies.csv.gz"
  "$SMOKE_ROOT/results/04_longitudinal_functional_potential/Table_06_longitudinal_functional_models.csv.gz"
  "$SMOKE_ROOT/results/05_ascent_descent_hysteresis/Table_07_equal_altitude_paired_contrasts.csv.gz"
  "$SMOKE_ROOT/results/figures/Figure_06_equal_altitude_paired_contrasts.pdf"
  "$SMOKE_ROOT/results/reproducibility/Figure_03_plot_data.csv"
  "$SMOKE_ROOT/results/reproducibility/Figure_04_plot_data.csv"
  "$SMOKE_ROOT/results/reproducibility/Figure_05_selected_models.csv"
  "$SMOKE_ROOT/results/reproducibility/Figure_06_plot_data.csv"
  "$SMOKE_ROOT/results/RESULTS_INDEX.md"
  "$SMOKE_ROOT/results/06_targeted_candidate_validation/Table_08_candidate_carrier_MAG_taxonomy.csv"
  "$SMOKE_ROOT/results/06_targeted_candidate_validation/Table_09_leave_one_participant_out_models.csv"
  "$SMOKE_ROOT/results/06_targeted_candidate_validation/Table_10_candidate_posthoc_timepoint_contrasts.csv"
  "$SMOKE_ROOT/results/06_targeted_candidate_validation/Table_11_captured_MAG_fraction_sensitivity.csv"
  "$SMOKE_ROOT/results/07_pathway_unknown_exploration/Table_12_candidate_function_co_carriage.csv.gz"
  "$SMOKE_ROOT/results/07_pathway_unknown_exploration/Table_13_unknown_singleton_catalogue_audit.csv"
  "$SMOKE_ROOT/results/07_pathway_unknown_exploration/Table_14_novelty_fragmentation_evidence_readiness.csv"
  "$SMOKE_ROOT/results/figures/Figure_07_candidate_participant_trajectories.pdf"
  "$SMOKE_ROOT/results/figures/Figure_08_candidate_taxonomic_contributions.pdf"
  "$SMOKE_ROOT/results/figures/Figure_09_candidate_leave_one_participant_out.pdf"
  "$SMOKE_ROOT/results/figures/Figure_10_unknown_singleton_catalogue_audit.pdf"
  "$SMOKE_ROOT/results/08_shared_unique_carrier_decomposition/Table_15_candidate_carrier_set_membership.csv"
  "$SMOKE_ROOT/results/08_shared_unique_carrier_decomposition/Table_16_sample_component_decomposition.csv.gz"
  "$SMOKE_ROOT/results/08_shared_unique_carrier_decomposition/Table_17_component_temporal_models.csv"
  "$SMOKE_ROOT/results/figures/Figure_11_candidate_carrier_set_overlap.pdf"
  "$SMOKE_ROOT/results/figures/Figure_12_shared_unique_component_trajectories.pdf"
  "$SMOKE_ROOT/results/00_quality_control/08_decompose_shared_unique_carriers.complete"
  "$SMOKE_ROOT/results/09_shared_component_followup/Table_18_shared_component_timepoint_contrasts.csv"
  "$SMOKE_ROOT/results/09_shared_component_followup/Table_19_shared_component_leave_one_participant_out.csv"
  "$SMOKE_ROOT/results/09_shared_component_followup/Table_20_shared_component_MAG_drivers.csv.gz"
  "$SMOKE_ROOT/results/09_shared_component_followup/Table_21_shared_component_genus_drivers.csv.gz"
  "$SMOKE_ROOT/results/09_shared_component_followup/Table_22_shared_component_driver_concentration.csv"
  "$SMOKE_ROOT/results/figures/Figure_13_shared_component_timepoint_contrasts.pdf"
  "$SMOKE_ROOT/results/figures/Figure_14_shared_component_leave_one_participant_out.pdf"
  "$SMOKE_ROOT/results/figures/Figure_15_shared_component_genus_drivers.pdf"
  "$SMOKE_ROOT/results/00_quality_control/09_analyze_shared_component_drivers.complete"
  "$SMOKE_ROOT/results/10_candidate_annotation_audit/Table_23_candidate_representative_protein_clusters.csv"
  "$SMOKE_ROOT/results/10_candidate_annotation_audit/Table_24_candidate_carrier_MAG_concordance.csv"
  "$SMOKE_ROOT/results/10_candidate_annotation_audit/Table_25_candidate_curation_summary.csv"
  "$SMOKE_ROOT/results/00_quality_control/10_audit_candidate_annotations.complete"
  "$SMOKE_ROOT/results/reproducibility/10_audit_candidate_annotations_session_info.txt"
  "$SMOKE_ROOT/results/FOLLOWUP_RESULTS_INDEX.md"
)
for path in "${required[@]}"; do
  test -s "$path" || { echo "Missing smoke-test output: $path" >&2; exit 1; }
done
$R_BIN analysis/mag_functional_annotation/tests/validate-smoke-outputs.R "$SMOKE_ROOT/results"
MAGFUNC_RSCRIPT="$R_BIN" bash scripts/run_mag_functional_annotation.sh \
  --config "$SMOKE_CONFIG" --from-stage 9 --to-stage 9 --force
grep -q 'Table 25' "$SMOKE_ROOT/results/FOLLOWUP_RESULTS_INDEX.md" || {
  echo "Stage 9 rerun replaced the completed Stage 10 follow-up index." >&2
  exit 1
}
echo "MAG functional annotation smoke test passed: $SMOKE_ROOT"
