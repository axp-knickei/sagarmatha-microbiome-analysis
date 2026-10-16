#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1]]
audit <- fread(file.path(root, "01_prepared_profiles", "Table_26_functional_profile_input_and_filtering_audit.csv"))
stopifnot(nrow(audit) == 5L, uniqueN(audit$samples) == 1L, all(audit$sample_order_identical),
          all(audit$metadata_identical), all(audit$aitchison_zero_replacement_passed),
          all(nchar(audit$source_checksum_md5) == 32L), all(nchar(audit$retained_matrix_checksum_md5) == 32L))
global <- fread(file.path(root, "02_multivariate_analysis", "Table_27_global_PERMANOVA_and_dispersion_results.csv"))
stopifnot(nrow(global) == 15L, all(c("permanova_p_value", "dispersion_p_value", "supported_global") %in% names(global)))
pairwise <- fread(file.path(root, "02_multivariate_analysis", "Table_28_supported_functional_timepoint_contrasts.csv"))
if (nrow(pairwise)) stopifnot(all(pairwise[, .N, by = .(layer, metric)]$N == 45L))
longitudinal <- fread(file.path(root, "03_longitudinal_stability", "Table_29_participant_distance_reference_and_adjacent_transitions.csv"))
stopifnot(all(c("distance_to_earliest_observation", "adjacent_observed_transition", "exploratory_t1_recovery") %in% longitudinal$record_type))
stopifnot(uniqueN(longitudinal[record_type == "exploratory_t1_recovery", subject_id]) == 3L)
models <- fread(file.path(root, "03_longitudinal_stability", "Table_30_matched_taxonomic_functional_stability_models.csv"))
stopifnot(nrow(models) == 8L, all(c("taxonomic_functional_transition_coupling", "functional_smoothing_relative_to_taxonomy") %in% models$model_family))
sensitivity <- fread(file.path(root, "03_longitudinal_stability", "Table_31_participant_level_transition_sensitivity.csv"))
stopifnot(
  nrow(sensitivity) == 8L,
  all(c(
    "participant_spearman_transition_coupling",
    "participant_median_relative_movement"
  ) %in% sensitivity$sensitivity_family),
  all(sensitivity$status == "fitted"),
  all(sensitivity$participants >= 3L),
  all(is.finite(sensitivity$p_value))
)
report <- readLines(file.path(root, "reports", "2026-08-12_functional_community_dynamics_interpretation.md"))
stopifnot(
  any(grepl("do not support functional smoothing", report, fixed = TRUE)),
  any(grepl("Table 31", report, fixed = TRUE))
)
cat("Functional community dynamics smoke outputs validated\n")
