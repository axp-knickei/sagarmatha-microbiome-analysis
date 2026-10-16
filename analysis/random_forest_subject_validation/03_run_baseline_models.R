#!/usr/bin/env Rscript
# Purpose: Evaluate majority-class and stratified-random baselines with the same
# subject-level holdouts used for model validation.
# Inputs: Canonical sample metadata and validation settings.
# Outputs: Baseline predictions and comparison CSVs in RF_VALIDATION_RUN_DIR.
# Dependencies: R/common.R and 02_define_grouped_resampling.R.
# Usage: Normally invoked by run_random_forest_subject_validation.R; it can also
# run directly with RF_VALIDATION_RUN_DIR set.
# Workflow: Stage 03; provides the non-model reference for Random Forest results.

source("analysis/random_forest_subject_validation/R/common.R")
source("analysis/random_forest_subject_validation/02_define_grouped_resampling.R")

run_loso_baselines <- function(metadata, cfg) {
  sid <- cfg$data$sample_id; sub <- cfg$data$subject_id; outcome <- cfg$data$outcome
  levels <- unlist(cfg$data$altitude_levels); splits <- make_loso_splits(metadata, sub)
  repeats <- as.integer(cfg$analysis$baseline_repeats[[cfg$analysis$mode]])
  if (!is.finite(repeats) || repeats < 1L) stop("Configured baseline repeats must be positive.", call. = FALSE)
  bind_rows(lapply(seq_along(splits), function(i) {
    sp <- splits[[i]]; train <- metadata[sp$analysis, ]; test <- metadata[sp$assessment, ]
    counts <- table(factor(as.character(train[[outcome]]), levels = levels)); props <- as.numeric(counts / sum(counts))
    majority <- levels[which.max(counts)]
    majority_rows <- tibble(
      model = "majority", repeat_id = 0L, fold_id = sp$id,
      !!sid := test[[sid]], !!sub := test[[sub]],
      observed = as.character(test[[outcome]]), predicted = majority
    )
    random_rows <- bind_rows(lapply(seq_len(repeats), function(r) {
      set.seed(cfg$analysis$seed + i * 10000L + r)
      tibble(
        model = "stratified_random", repeat_id = r, fold_id = sp$id,
        !!sid := test[[sid]], !!sub := test[[sub]],
        observed = as.character(test[[outcome]]),
        predicted = sample(levels, nrow(test), replace = TRUE, prob = props)
      )
    }))
    output <- bind_rows(majority_rows, random_rows)
    for (j in seq_along(levels)) {
      # Both baselines use the training prevalence distribution for probabilistic
      # scoring. The majority baseline differs only in its hard class decision.
      output[[paste0("prob_", levels[[j]])]] <- props[[j]]
    }
    output
  }))
}

summarize_baselines <- function(predictions, levels) {
  probability_columns <- paste0("prob_", levels)
  by_repeat <- predictions |>
    group_by(model, repeat_id) |>
    group_modify(function(.x, .g) {
      probabilities <- .x[, probability_columns, drop = FALSE]
      names(probabilities) <- levels
      classification_metrics(.x$observed, .x$predicted, levels, probabilities)$metrics
    }) |>
    ungroup()
  metric_columns <- setdiff(names(by_repeat)[vapply(by_repeat, is.numeric, logical(1))], "repeat_id")
  summary <- by_repeat |>
    group_by(model) |>
    summarise(
      repeats = n(),
      across(all_of(metric_columns),
        list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ if (n() > 1L) sd(.x, na.rm = TRUE) else NA_real_),
        .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
  list(by_repeat = by_repeat, summary = summary)
}

if (sys.nframe() == 0L) {
  cfg <- read_validation_config(); run_dir <- Sys.getenv("RF_VALIDATION_RUN_DIR", "")
  if (!nzchar(run_dir)) run_dir <- make_run_directory(cfg)
  meta <- sample_metadata(read_long_data(cfg), cfg); p <- run_loso_baselines(meta, cfg); levels <- unlist(cfg$data$altitude_levels)
  summaries <- summarize_baselines(p, levels)
  write.csv(p, file.path(run_dir, "predictions", "baseline_predictions.csv"), row.names = FALSE)
  write.csv(summaries$by_repeat, file.path(run_dir, "comparisons", "baseline_metrics_by_repeat.csv"), row.names = FALSE)
  write.csv(summaries$summary, file.path(run_dir, "comparisons", "baseline_model_comparison.csv"), row.names = FALSE)
  cat(run_dir, "\n")
}
