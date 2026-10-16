#!/usr/bin/env Rscript
# Purpose: Evaluate prediction, probability, performance, runtime, and
# importance-rank stability across configured Random Forest tree counts.
# Inputs: Validation data/config and the shared LOSO fitter.
# Outputs: Tree-stability and hyperparameter-summary CSVs under
# RF_VALIDATION_RUN_DIR/comparisons.
# Dependencies: 04_run_loso_random_forest.R and its dependencies.
# Usage: Normally invoked by run_random_forest_subject_validation.R.
# Workflow: Stage 11 sensitivity analysis; no automatic tree count is selected.

source("analysis/random_forest_subject_validation/04_run_loso_random_forest.R")

cfg <- read_validation_config()
run_dir <- Sys.getenv("RF_VALIDATION_RUN_DIR", "")
if (!nzchar(run_dir)) stop("Set RF_VALIDATION_RUN_DIR.", call. = FALSE)
d <- read_long_data(cfg)
levels <- unlist(cfg$data$altitude_levels)
counts <- as.integer(unlist(cfg$tree_stability[[cfg$analysis$mode]]))
if (!length(counts) || any(!is.finite(counts)) || any(counts < 1L)) stop("Invalid tree-stability sequence.", call. = FALSE)

fits <- lapply(counts, function(n) fit_loso_rf(d, cfg, ntree = n))
reference <- fits[[length(fits)]]$predictions
probability_columns <- paste0("prob_", levels)
sid <- cfg$data$sample_id

stability <- bind_rows(lapply(seq_along(counts), function(i) {
  fit <- fits[[i]]
  prediction <- fit$predictions
  aligned <- match(reference[[sid]], prediction[[sid]])
  if (anyNA(aligned)) stop("Tree-stability predictions do not contain identical samples.", call. = FALSE)
  current_prob <- as.matrix(prediction[aligned, probability_columns, drop = FALSE])
  reference_prob <- as.matrix(reference[, probability_columns, drop = FALSE])
  metric <- metrics_from_prediction_table(prediction, levels)$metrics
  bind_cols(
    tibble(
      trees = counts[[i]],
      elapsed_seconds = sum(fit$fold_metrics$elapsed_seconds),
      class_agreement_with_largest_forest = mean(prediction$predicted[aligned] == reference$predicted),
      mean_absolute_probability_difference = mean(abs(current_prob - reference_prob)),
      maximum_absolute_probability_difference = max(abs(current_prob - reference_prob))
    ),
    metric
  )
}))

importance_rank_stability <- bind_rows(lapply(seq_along(counts), function(i) {
  imp <- fits[[i]]$importance
  score <- intersect(c("MeanDecreaseAccuracy", "MeanDecreaseGini"), names(imp))[[1]]
  imp |>
    group_by(fold_id) |>
    mutate(rank = rank(-.data[[score]], ties.method = "average")) |>
    ungroup() |>
    transmute(trees = counts[[i]], fold_id, feature, rank)
})) |>
  group_by(fold_id, feature) |>
  mutate(reference_rank = rank[trees == max(trees)]) |>
  ungroup() |>
  group_by(trees, fold_id) |>
  summarise(rank_correlation_with_largest_forest = cor(rank, reference_rank, method = "spearman"), .groups = "drop")

write.csv(stability, file.path(run_dir, "comparisons", "tree_count_stability.csv"), row.names = FALSE)
write.csv(importance_rank_stability, file.path(run_dir, "comparisons", "tree_importance_rank_stability.csv"), row.names = FALSE)
write.csv(
  tibble(
    mode = cfg$analysis$mode, mtry = cfg$model$mtry,
    min_node_size = cfg$model$min_node_size, replace = cfg$model$replace,
    sample_fraction = cfg$model$sample_fraction,
    trees = paste(counts, collapse = ";"),
    tuning = "fixed; grouped nested tuning is not implemented",
    selection_rule = "No automatic tree-count selection. Inspect metric, class, probability, importance-rank, and runtime stability."
  ),
  file.path(run_dir, "comparisons", "hyperparameter_summary.csv"),
  row.names = FALSE
)
