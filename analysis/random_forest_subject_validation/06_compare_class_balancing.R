#!/usr/bin/env Rscript
# Purpose: Compare unweighted, class-weighted, balanced-sampling, and combined
# LOSO Random Forest strategies.
# Inputs: Validation data/config and the shared LOSO fitter.
# Outputs: Class-balance comparison, fold-metric, and fold-setting CSVs under
# RF_VALIDATION_RUN_DIR/comparisons.
# Dependencies: 04_run_loso_random_forest.R and its dependencies.
# Usage: Normally invoked by run_random_forest_subject_validation.R.
# Workflow: Stage 06 sensitivity comparison; it does not replace the primary fit.

source("analysis/random_forest_subject_validation/04_run_loso_random_forest.R")
cfg <- read_validation_config(); run_dir <- Sys.getenv("RF_VALIDATION_RUN_DIR", ""); if (!nzchar(run_dir)) run_dir <- make_run_directory(cfg)
dir.create(file.path(run_dir, "comparisons"), recursive = TRUE, showWarnings = FALSE)
d <- read_long_data(cfg); levels <- unlist(cfg$data$altitude_levels)
strategies <- c("unweighted", "class_weights", "balanced_sampling", "combined")
fits <- setNames(lapply(strategies, function(s) fit_loso_rf(d, cfg, strategy = s)), strategies)
comparison <- bind_rows(lapply(strategies, function(s) {
  bind_cols(tibble(strategy = s), metrics_from_prediction_table(fits[[s]]$predictions, levels)$metrics)
}))
write.csv(comparison, file.path(run_dir, "comparisons", "class_balance_strategy_comparison.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(fits, `[[`, "fold_metrics")), file.path(run_dir, "comparisons", "class_balance_fold_metrics.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(fits, `[[`, "training_settings")), file.path(run_dir, "comparisons", "class_balance_fold_settings.csv"), row.names = FALSE)
