#!/usr/bin/env Rscript

# Purpose: Validate semantic contracts that file-existence checks cannot catch.
# Inputs: The synthetic MAG functional annotation results directory.
# Outputs: No files; exits nonzero when required model fields or ordering fail.
# Dependencies: R and data.table.
# Usage: Rscript analysis/mag_functional_annotation/tests/validate-smoke-outputs.R /tmp/.../results
# Workflow: Called by run_smoke_test.sh after the complete synthetic workflow.

source("analysis/mag_functional_annotation/R/common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Expected one smoke-test results directory.", call. = FALSE)
result_dir <- args[[1]]

models <- read_table(file.path(
  result_dir,
  "04_longitudinal_functional_potential",
  "Table_06_longitudinal_functional_models.csv.gz"
))
required_model_columns <- c(
  "term", "ontology", "omnibus_chisq", "omnibus_df", "omnibus_p_value",
  "omnibus_p_adj_bh", "singular", "status"
)
missing_model_columns <- setdiff(required_model_columns, names(models))
if (length(missing_model_columns) > 0L) {
  stop("Table 06 is missing columns: ", paste(missing_model_columns, collapse = ", "), call. = FALSE)
}
fitted_models <- models[status == "fitted"]
if (nrow(fitted_models) == 0L) stop("Smoke Table 06 contains no fitted models.", call. = FALSE)
if (any(!is.finite(fitted_models$omnibus_df) | fitted_models$omnibus_df <= 0)) {
  stop("Fitted Table 06 models must retain positive finite omnibus degrees of freedom.", call. = FALSE)
}

hysteresis <- read_table(file.path(
  result_dir,
  "05_ascent_descent_hysteresis",
  "Table_07_equal_altitude_paired_contrasts.csv.gz"
))
required_hysteresis_columns <- c(
  "term", "ontology", "contrast_id", "n_pairs", "primary_contrast",
  "p_adj_bh", "wilcoxon_p_adj_bh", "status"
)
missing_hysteresis_columns <- setdiff(required_hysteresis_columns, names(hysteresis))
if (length(missing_hysteresis_columns) > 0L) {
  stop("Table 07 is missing columns: ", paste(missing_hysteresis_columns, collapse = ", "), call. = FALSE)
}
if (any(hysteresis$primary_contrast) && !isTRUE(hysteresis$primary_contrast[[1]])) {
  stop("Table 07 must place the prespecified primary contrast before supporting contrasts.", call. = FALSE)
}

figure_three_data <- read_table(file.path(result_dir, "reproducibility", "Figure_03_plot_data.csv"))
if (!setequal(unique(figure_three_data$model), c("annotation_row_fraction", "unknown_or_unannotated_fraction"))) {
  stop("Figure 03 must retain observed data for both annotation-status definitions.", call. = FALSE)
}
if (any(!is.finite(figure_three_data$probability)) || any(figure_three_data$probability < 0 | figure_three_data$probability > 1)) {
  stop("Figure 03 observed fractions must be finite and fall between zero and one.", call. = FALSE)
}
annotation_model_status <- read_table(file.path(result_dir, "02_unknown_by_cluster_size", "unknown_annotation_model_status.csv"))
if (!setequal(annotation_model_status$model, c("annotation_probability", "unknown_probability")) ||
    any(annotation_model_status$status != "fitted")) {
  stop("Both cluster-size annotation models must report fitted status.", call. = FALSE)
}

results_index <- read_table(file.path(result_dir, "RESULTS_INDEX.csv"))
figure_three_index <- results_index[result_number == "Figure 03"]
if (nrow(figure_three_index) == 0L ||
    any(!grepl("Observed annotation-record", figure_three_index$description, fixed = TRUE))) {
  stop("The results index must describe Figure 03 as observed size-group fractions.", call. = FALSE)
}

figure_five_models <- read_table(file.path(result_dir, "reproducibility", "Figure_05_selected_models.csv"))
if (anyDuplicated(figure_five_models$profile_signature)) {
  stop("Figure 05 must not repeat identical longitudinal mean profiles.", call. = FALSE)
}
if (any(figure_five_models$singular %in% TRUE)) {
  stop("Figure 05 must not highlight singular longitudinal models.", call. = FALSE)
}

figure_six_data <- read_table(file.path(result_dir, "reproducibility", "Figure_06_plot_data.csv"))
if (uniqueN(figure_six_data$contrast_id) != 3L) {
  stop("Figure 06 must include top-ranked effects from all three prespecified contrasts.", call. = FALSE)
}
if (any(figure_six_data[, .N, by = contrast_id]$N > 6L)) {
  stop("Figure 06 must include no more than six effects per contrast.", call. = FALSE)
}

candidate_carriers <- read_table(file.path(
  result_dir, "06_targeted_candidate_validation", "Table_08_candidate_carrier_MAG_taxonomy.csv"
))
if (!setequal(unique(candidate_carriers$candidate_id), c("KO:K00001", "EC:1.1.1.1"))) {
  stop("Table 08 must contain both configured synthetic candidates.", call. = FALSE)
}
if (anyNA(candidate_carriers$MAG) || anyNA(candidate_carriers$phylum)) {
  stop("Every synthetic candidate carrier must retain exact MAG and taxonomy identifiers.", call. = FALSE)
}

lopo <- read_table(file.path(
  result_dir, "06_targeted_candidate_validation", "Table_09_leave_one_participant_out_models.csv"
))
if (uniqueN(lopo$excluded_subject) != 6L || any(lopo$status != "fitted")) {
  stop("The synthetic leave-one-participant-out analysis must fit all six exclusions.", call. = FALSE)
}
if (any(!grepl("selection-aware", lopo$inference_scope, fixed = TRUE))) {
  stop("LOPO results must preserve their selection-aware interpretation scope.", call. = FALSE)
}

posthoc <- read_table(file.path(
  result_dir, "06_targeted_candidate_validation", "Table_10_candidate_posthoc_timepoint_contrasts.csv"
))
if (any(posthoc[, .N, by = candidate_id]$N != 45L) ||
    any(!is.finite(posthoc$p_adj_bh_within_candidate))) {
  stop("Each synthetic candidate must have 45 finite pairwise T1-T10 contrasts.", call. = FALSE)
}

sensitivity <- read_table(file.path(
  result_dir, "06_targeted_candidate_validation", "Table_11_captured_MAG_fraction_sensitivity.csv"
))
required_scenarios <- c(
  "primary_unnormalized_refit", "captured_fraction_adjusted",
  "row_sum_normalized_sensitivity"
)
if (any(sensitivity[, .(scenarios = uniqueN(scenario)), by = candidate_id]$scenarios != 3L) ||
    !setequal(unique(sensitivity$scenario), required_scenarios)) {
  stop("Each candidate must retain all three primary/sensitivity scenarios.", call. = FALSE)
}

co_carriage <- read_table(file.path(
  result_dir, "07_pathway_unknown_exploration", "Table_12_candidate_function_co_carriage.csv.gz"
))
if (nrow(co_carriage) == 0L || any(!grepl("not pathway completeness", co_carriage$analysis_scope, fixed = TRUE))) {
  stop("Functional co-carriage must be non-empty and explicitly separated from pathway completeness.", call. = FALSE)
}

unknown_audit <- read_table(file.path(
  result_dir, "07_pathway_unknown_exploration", "Table_13_unknown_singleton_catalogue_audit.csv"
))
if (!setequal(unique(unknown_audit$scope), c("global_catalogue", "gut_selected_MAGs")) ||
    any(unknown_audit$unknown_or_unannotated_fraction < 0 | unknown_audit$unknown_or_unannotated_fraction > 1)) {
  stop("Unknown-singleton audit must retain both scopes and valid descriptive fractions.", call. = FALSE)
}

readiness <- read_table(file.path(
  result_dir, "07_pathway_unknown_exploration", "Table_14_novelty_fragmentation_evidence_readiness.csv"
))
if (!all(c("protein_length_amino_acids", "MAG_completeness_and_contamination") %in% readiness$evidence_component)) {
  stop("Novelty/fragmentation readiness must state the missing protein and MAG-quality evidence.", call. = FALSE)
}

followup_index <- read_table(file.path(result_dir, "FOLLOWUP_RESULTS_INDEX.csv"))
if (!all(sprintf("Table %02d", 8:25) %in% followup_index$result_number) ||
    !all(sprintf("Figure %02d", 7:15) %in% followup_index$result_number)) {
  stop("The follow-up index must contain Tables 08-25 and Figures 07-15.", call. = FALSE)
}

chart_map <- read_table(file.path(result_dir, "reproducibility", "followup_publication_chart_map.csv"))
if (!setequal(chart_map$figure, sprintf("Figure %02d", 7:15)) ||
    any(!nzchar(chart_map$denominator_or_unit))) {
  stop("Follow-up chart map must document Figures 07-15 and their denominators.", call. = FALSE)
}

membership <- read_table(file.path(
  result_dir, "08_shared_unique_carrier_decomposition", "Table_15_candidate_carrier_set_membership.csv"
))
membership_counts <- membership[, .N, by = membership_class]
expected_membership <- data.table(
  membership_class = c(
    "shared_by_both_candidates", "KO:K00001_only", "EC:1.1.1.1_only"
  ),
  N = c(2L, 2L, 1L)
)
if (!isTRUE(all.equal(
  membership_counts[order(membership_class)],
  expected_membership[order(membership_class)],
  check.attributes = FALSE
))) {
  stop("Synthetic carrier membership must retain shared and both candidate-specific classes.", call. = FALSE)
}

components <- read_table(file.path(
  result_dir, "08_shared_unique_carrier_decomposition", "Table_16_sample_component_decomposition.csv.gz"
))
if (any(components[, .N, by = .(candidate_id, sampleID)]$N != 2L) ||
    max(abs(components$reconstruction_error), na.rm = TRUE) > 1e-10 ||
    max(abs(components[, sum(component_share), by = .(candidate_id, sampleID)]$V1 - 1), na.rm = TRUE) > 1e-10) {
  stop("Synthetic component rows must exactly reconstruct each candidate-sample.", call. = FALSE)
}

component_models <- read_table(file.path(
  result_dir, "08_shared_unique_carrier_decomposition", "Table_17_component_temporal_models.csv"
))
if (nrow(component_models) != 5L || any(component_models$status != "fitted") ||
    any(!is.finite(component_models$p_adj_bh_within_family))) {
  stop("Synthetic decomposition must retain five fitted finite component models.", call. = FALSE)
}

shared_contrasts <- read_table(file.path(
  result_dir, "09_shared_component_followup", "Table_18_shared_component_timepoint_contrasts.csv"
))
shared_q_column <- "p_adj_bh_across_45_timepoint_contrasts"
if (nrow(shared_contrasts) != 45L ||
    !shared_q_column %in% names(shared_contrasts) ||
    any(!is.finite(shared_contrasts[[shared_q_column]]))) {
  stop("Table 18 must contain 45 finite BH-adjusted shared-component contrasts.", call. = FALSE)
}

shared_lopo <- read_table(file.path(
  result_dir, "09_shared_component_followup", "Table_19_shared_component_leave_one_participant_out.csv"
))
if (nrow(shared_lopo) != 6L || any(shared_lopo$status != "fitted") ||
    any(!is.finite(shared_lopo$omnibus_raw_p_value)) ||
    any(shared_lopo$proportion_original_significant_retained < 0 |
      shared_lopo$proportion_original_significant_retained > 1, na.rm = TRUE)) {
  stop("Table 19 must retain six fitted finite selection-aware LOPO diagnostics.", call. = FALSE)
}

mag_drivers <- read_table(file.path(
  result_dir, "09_shared_component_followup", "Table_20_shared_component_MAG_drivers.csv.gz"
))
mag_share_sums <- mag_drivers[, .(share_sum = sum(absolute_change_share)), by = contrast_id]
if (uniqueN(mag_drivers$contrast_id) != 45L || nrow(mag_drivers) != 90L ||
    any(abs(mag_share_sums$share_sum - 1) > 1e-10)) {
  stop("Table 20 must decompose all 45 contrasts across two shared synthetic MAGs.", call. = FALSE)
}

genus_drivers <- read_table(file.path(
  result_dir, "09_shared_component_followup", "Table_21_shared_component_genus_drivers.csv.gz"
))
if (uniqueN(genus_drivers$contrast_id) != 45L ||
    any(!is.finite(genus_drivers$mean_delta))) {
  stop("Table 21 must retain finite genus drivers for all 45 contrasts.", call. = FALSE)
}

driver_concentration <- read_table(file.path(
  result_dir, "09_shared_component_followup", "Table_22_shared_component_driver_concentration.csv"
))
if (nrow(driver_concentration) != 45L ||
    any(driver_concentration$cancellation_fraction < 0 |
      driver_concentration$cancellation_fraction > 1) ||
    any(driver_concentration$top_1_MAG_absolute_change_share >
      driver_concentration$top_3_MAG_absolute_change_share) ||
    any(driver_concentration$top_3_MAG_absolute_change_share >
      driver_concentration$top_5_MAG_absolute_change_share)) {
  stop("Table 22 must retain coherent concentration summaries for all 45 contrasts.", call. = FALSE)
}

candidate_clusters <- read_table(file.path(
  result_dir, "10_candidate_annotation_audit", "Table_23_candidate_representative_protein_clusters.csv"
))
if (!setequal(unique(candidate_clusters$candidate_id), c("KO:K00001", "EC:1.1.1.1")) ||
    anyNA(candidate_clusters$representative_id) || anyNA(candidate_clusters$representative_mag)) {
  stop("Table 23 must retain both synthetic candidates and complete representative joins.", call. = FALSE)
}
same_representatives <- unique(candidate_clusters[
  carries_both_candidate_terms_same_representative %in% TRUE, representative_id
])
if (length(same_representatives) != 1L ||
    candidate_clusters[representative_id %in% same_representatives, uniqueN(candidate_id)] != 2L) {
  stop("Table 23 must detect the one synthetic representative carrying both candidate terms.", call. = FALSE)
}

carrier_concordance <- read_table(file.path(
  result_dir, "10_candidate_annotation_audit", "Table_24_candidate_carrier_MAG_concordance.csv"
))
if (anyNA(carrier_concordance$MAG) || anyNA(carrier_concordance$phylum) ||
    sum(carrier_concordance$candidate_1_carried) != 4L ||
    sum(carrier_concordance$candidate_2_carried) != 3L ||
    sum(carrier_concordance$shared_carrier) != 2L ||
    nrow(carrier_concordance) != 5L ||
    sum(carrier_concordance$both_terms_on_same_representative_cluster) < 1L) {
  stop("Table 24 must exactly reconcile synthetic candidate carrier and overlap evidence.", call. = FALSE)
}

curation_summary <- read_table(file.path(
  result_dir, "10_candidate_annotation_audit", "Table_25_candidate_curation_summary.csv"
))
if (!setequal(curation_summary$candidate_id, c("KO:K00001", "EC:1.1.1.1")) ||
    any(curation_summary$shared_carrier_MAGs != 2L) ||
    any(curation_summary$carrier_union_MAGs != 5L) ||
    any(curation_summary$representative_clusters_carrying_both_terms != 1L) ||
    any(!grepl("pathway", curation_summary$pathway_claim_readiness, ignore.case = TRUE))) {
  stop("Table 25 must reconcile candidate curation and preserve pathway claim boundaries.", call. = FALSE)
}

message("Smoke-test output contracts passed.")
