#!/usr/bin/env Rscript
# Purpose: Compare fold-local variance-ranked feature counts with all eligible
# features while retaining subject-disjoint LOSO evaluation.
# Inputs: Validation data/config and the shared LOSO fitter.
# Outputs: Feature-strategy comparison and fold-metric CSVs under
# RF_VALIDATION_RUN_DIR/comparisons.
# Dependencies: 04_run_loso_random_forest.R and its dependencies.
# Usage: Normally invoked by run_random_forest_subject_validation.R.
# Workflow: Stage 07 feature-selection sensitivity comparison.

source("analysis/random_forest_subject_validation/04_run_loso_random_forest.R")
cfg <- read_validation_config(); run_dir <- Sys.getenv("RF_VALIDATION_RUN_DIR", ""); if (!nzchar(run_dir)) run_dir <- make_run_directory(cfg)
dir.create(file.path(run_dir, "comparisons"), recursive = TRUE, showWarnings = FALSE)
d <- read_long_data(cfg); levels <- unlist(cfg$data$altitude_levels)
counts <- if (cfg$analysis$mode == "smoke") c(20, 50, NA) else c(10, 20, 50, 100, NA)
labels <- ifelse(is.na(counts), "all_eligible", paste0("top_", counts))
fits <- setNames(lapply(counts, function(n) {
  fit_loso_rf(d, cfg, top_n = if (is.na(n)) NULL else n)
}), labels)
comparison <- bind_rows(lapply(seq_along(counts), function(i) {
  n <- counts[[i]]
  z <- fits[[labels[[i]]]]
  bind_cols(tibble(
    feature_strategy = ifelse(is.na(n), "all_eligible", paste0("top_", n)),
    ranking_method = ifelse(is.na(n), "none", "training_fold_variance")
  ), metrics_from_prediction_table(z$predictions, levels)$metrics)
}))
write.csv(comparison, file.path(run_dir, "comparisons", "feature_strategy_comparison.csv"), row.names = FALSE)
feature_fold_metrics <- bind_rows(lapply(names(fits), function(label) {
  fits[[label]]$fold_metrics |> mutate(feature_strategy = label)
}))
write.csv(feature_fold_metrics, file.path(run_dir, "comparisons", "feature_strategy_fold_metrics.csv"), row.names = FALSE)
