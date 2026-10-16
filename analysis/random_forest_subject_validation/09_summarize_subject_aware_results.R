#!/usr/bin/env Rscript
# Purpose: Summarize LOSO performance at class and subject levels and assemble
# stable exploratory feature candidates from fold-level outputs.
# Inputs: Completed baseline, LOSO, feature-importance, and selection tables.
# Outputs: Additional prediction, feature-stability, candidate, and report files
# within RF_VALIDATION_RUN_DIR.
# Dependencies: R/common.R and outputs from earlier validation stages.
# Usage: Normally invoked by run_random_forest_subject_validation.R.
# Workflow: Stage 09; translates validated outputs into cautious summaries.

source("analysis/random_forest_subject_validation/R/common.R")

cfg <- read_validation_config()
run_dir <- Sys.getenv("RF_VALIDATION_RUN_DIR", "")
if (!nzchar(run_dir)) stop("Set RF_VALIDATION_RUN_DIR.", call. = FALSE)

overall <- read.csv(file.path(run_dir, "predictions", "loso_overall_metrics.csv"))
fold <- read.csv(file.path(run_dir, "predictions", "loso_fold_metrics.csv"))
predictions <- read.csv(file.path(run_dir, "predictions", "loso_predictions.csv"), check.names = FALSE)
baseline <- read.csv(file.path(run_dir, "comparisons", "baseline_model_comparison.csv"))
levels <- unlist(cfg$data$altitude_levels)

per_class <- metrics_from_prediction_table(predictions, levels)$per_class
write.csv(per_class, file.path(run_dir, "predictions", "loso_per_class_metrics.csv"), row.names = FALSE)

subject_summary <- tibble(
  metric = c(
    "subject_macro_accuracy", "subject_macro_present_class_balanced_accuracy",
    "subject_macro_f1", "subjects", "single_observation_subjects"
  ),
  estimate = c(
    mean(fold$accuracy),
    mean(fold$balanced_accuracy),
    mean(fold$macro_f1),
    nrow(fold),
    sum(fold$n_test == 1L)
  ),
  rule = c(
    "unweighted mean of subject-fold accuracy",
    "unweighted mean of recall averaged only across classes observed for each subject; not comparable to four-class pooled balanced accuracy",
    "unweighted mean of macro F1 across classes observed for each subject",
    "count", "count"
  )
)
write.csv(subject_summary, file.path(run_dir, "predictions", "loso_subject_macro_metrics.csv"), row.names = FALSE)

majority_ba <- baseline$balanced_accuracy_mean[baseline$model == "majority"]
random_ba <- baseline$balanced_accuracy_mean[baseline$model == "stratified_random"]
mode_label <- cfg$analysis$mode
result_label <- if (mode_label == "production") {
  "Production-mode exploratory estimates; not independent confirmation."
} else {
  paste0(tools::toTitleCase(mode_label), "-mode developmental estimates; not production results.")
}

writeLines(c(
  "# Subject-aware validation report", "",
  result_label, "",
  paste0(
    "Pooled LOSO accuracy: ", sprintf("%.3f", overall$accuracy),
    "; four-class balanced accuracy: ", sprintf("%.3f", overall$balanced_accuracy),
    "; macro F1: ", sprintf("%.3f", overall$macro_f1), "."
  ),
  paste0(
    "Subject-macro accuracy: ", sprintf("%.3f", mean(fold$accuracy)),
    ". Subject-fold present-class balanced accuracy ranges from ",
    sprintf("%.3f", min(fold$balanced_accuracy)), " to ",
    sprintf("%.3f", max(fold$balanced_accuracy)), "."
  ),
  paste0(
    "Majority baseline balanced accuracy: ", sprintf("%.3f", majority_ba),
    "; stratified-random mean balanced accuracy: ", sprintf("%.3f", random_ba), "."
  ),
  "",
  "A subject-fold balanced accuracy averages recall only over classes observed for that subject. It is descriptive and is not a four-class performance estimate.",
  "A fold with one observation can have perfect descriptive accuracy without providing strong evidence of generalization.",
  "Altitude is perfectly confounded with timepoint/date; these results evaluate prediction among expedition-linked categories and cannot establish altitude-specific biology."
), file.path(run_dir, "reports", "subject_aware_validation_report.md"))

writeLines(c(
  "# Scientific interpretation", "",
  result_label, "",
  "The defensible framing is exploratory four-class prediction among altitude/time-confounded expedition categories, with ordinal error reported from the same out-of-subject predictions.",
  "Phase prediction is disabled because no authoritative phase mapping is available.",
  "Numeric regression and within-subject baseline-change prediction are disabled and not implemented. Baseline-change would require a baseline measurement from each new subject and answers a different prediction question.",
  "Feature importance is predictive model behavior, not association, statistical significance, mechanism, or causality."
), file.path(run_dir, "reports", "scientific_interpretation.md"))

writeLines(c(
  "# Limitations", "",
  "- Only 13 independent subjects; sample counts and altitude coverage differ markedly.",
  "- Altitude is perfectly linked to expedition timepoint and sampling date.",
  "- Four-class imbalance is severe.",
  "- Hyperparameters are fixed; full grouped nested tuning is not implemented.",
  "- Model-strategy comparisons reuse the outer LOSO predictions and must not be reported as an unbiased post-selection estimate.",
  "- Native Random Forest importance can favor prevalent or high-variance predictors and distribute importance across correlated MAGs.",
  "- Taxonomy and expedition-phase mappings are absent.",
  "- Numeric-altitude regression and baseline-change prediction are intentionally disabled."
), file.path(run_dir, "reports", "limitations.md"))

capture.output(sessionInfo(), file = file.path(run_dir, "reports", "session_info.txt"))

imp_path <- file.path(run_dir, "features", "feature_importance_by_fold.csv")
frequency_path <- file.path(run_dir, "features", "feature_selection_frequency.csv")
if (file.exists(imp_path)) {
  imp <- read.csv(imp_path, check.names = FALSE)
  score_candidates <- intersect(c("MeanDecreaseAccuracy", "MeanDecreaseGini"), names(imp))
  if (!length(score_candidates)) stop("No supported Random Forest importance column found.", call. = FALSE)
  score <- score_candidates[[1]]
  ranked <- imp |>
    group_by(fold_id) |>
    mutate(importance_rank = rank(-.data[[score]], ties.method = "average")) |>
    ungroup()
  stability <- ranked |>
    group_by(feature) |>
    summarise(
      importance_folds = n_distinct(fold_id),
      median_native_importance = median(.data[[score]], na.rm = TRUE),
      native_importance_iqr = IQR(.data[[score]], na.rm = TRUE),
      median_importance_rank = median(importance_rank, na.rm = TRUE),
      importance_rank_iqr = IQR(importance_rank, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(median_importance_rank, desc(median_native_importance))
  write.csv(stability, file.path(run_dir, "features", "feature_importance_stability.csv"), row.names = FALSE)
  frequency <- if (file.exists(frequency_path)) read.csv(frequency_path) else tibble(feature = character())
  candidates <- stability |>
    left_join(frequency, by = "feature") |>
    slice_head(n = 20) |>
    mutate(
      taxonomy = NA_character_,
      interpretation = "exploratory predictive feature; native importance and ranks are not biological effect estimates"
    )
  write.csv(candidates, file.path(run_dir, "features", "candidate_MAG_summary.csv"), row.names = FALSE)
}
