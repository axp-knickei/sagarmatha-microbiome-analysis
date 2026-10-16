#!/usr/bin/env Rscript
# Purpose: Record which outcome formulations are implemented, disabled, or
# derivable from the existing multiclass out-of-fold predictions.
# Inputs: Completed LOSO predictions and validation configuration.
# Outputs: comparisons/outcome_formulation_comparison.csv.
# Dependencies: R/common.R; requires RF_VALIDATION_RUN_DIR from a completed run.
# Usage: Normally invoked by run_random_forest_subject_validation.R.
# Workflow: Stage 08 interpretation boundary; it does not fit disabled outcomes.

source("analysis/random_forest_subject_validation/R/common.R")
cfg <- read_validation_config(); run_dir <- Sys.getenv("RF_VALIDATION_RUN_DIR", ""); if (!nzchar(run_dir)) stop("Set RF_VALIDATION_RUN_DIR to a completed smoke run.")
p <- read.csv(file.path(run_dir, "predictions", "loso_predictions.csv"), stringsAsFactors = FALSE)
levels <- unlist(cfg$data$altitude_levels); m <- metrics_from_prediction_table(p, levels)$metrics
out <- tibble(formulation = c("four_class_multiclass", "ordinal_metrics_from_multiclass", "numeric_altitude_regression", "expedition_phase", "within_subject_baseline_change"),
 status = c("implemented", "implemented", "disabled_not_implemented", "disabled_missing_mapping", "disabled_not_implemented"),
 primary_metric = c(m$balanced_accuracy, m$ordinal_mae, NA, NA, NA),
 metric_name = c("balanced_accuracy", "mean_absolute_rank_error", "MAE", "balanced_accuracy", "balanced_accuracy"),
 caution = c("Altitude is confounded with expedition time.", "Uses the same multiclass out-of-fold predictions.", "Requires separate grouped regression fitting.", "No phase variable or documented mapping exists.", "Requires baseline at prediction time; subjects without baseline trajectories cannot support this formulation."))
write.csv(out, file.path(run_dir, "comparisons", "outcome_formulation_comparison.csv"), row.names = FALSE)
