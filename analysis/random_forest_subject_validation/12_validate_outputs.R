#!/usr/bin/env Rscript
# Purpose: Verify completeness, subject separation, probability consistency,
# metric agreement, and figure coverage before a validation run is accepted.
# Inputs: All required audit, prediction, comparison, feature, report, and figure
# outputs in RF_VALIDATION_RUN_DIR.
# Outputs: reports/output_validation.csv and a failing exit status when any
# integrity check fails.
# Dependencies: R/common.R and completed stages 01-11.
# Usage: Normally invoked as the final orchestrated validation stage.
# Workflow: Stage 12 quality gate; it does not modify model results.

source("analysis/random_forest_subject_validation/R/common.R")

cfg <- read_validation_config()
run_dir <- Sys.getenv("RF_VALIDATION_RUN_DIR", "")
if (!nzchar(run_dir)) stop("Set RF_VALIDATION_RUN_DIR.", call. = FALSE)
levels <- unlist(cfg$data$altitude_levels)
required <- c(
  "audit/dataset_summary.csv",
  "audit/data_integrity_checks.csv",
  "predictions/loso_predictions.csv",
  "predictions/loso_fold_metrics.csv",
  "predictions/loso_overall_metrics.csv",
  "predictions/loso_confusion_matrix.csv",
  "predictions/grouped_cv_predictions.csv",
  "comparisons/baseline_model_comparison.csv",
  "comparisons/class_balance_strategy_comparison.csv",
  "comparisons/feature_strategy_comparison.csv",
  "comparisons/tree_count_stability.csv",
  "features/feature_eligibility_by_fold.csv",
  "features/feature_importance_by_fold.csv",
  "reports/analysis_audit.md",
  "reports/subject_aware_validation_report.md",
  "reports/scientific_interpretation.md",
  "reports/limitations.md",
  "reports/session_info.txt",
  "reports/run_manifest.yml"
)
missing <- required[!file.exists(file.path(run_dir, required))]
if (length(missing)) stop("Run is missing required outputs: ", paste(missing, collapse = ", "), call. = FALSE)

sid <- cfg$data$sample_id
subject <- cfg$data$subject_id
predictions <- read.csv(file.path(run_dir, "predictions", "loso_predictions.csv"), check.names = FALSE)
fold_metrics <- read.csv(file.path(run_dir, "predictions", "loso_fold_metrics.csv"))
saved_metrics <- read.csv(file.path(run_dir, "predictions", "loso_overall_metrics.csv"))
confusion <- read.csv(file.path(run_dir, "predictions", "loso_confusion_matrix.csv"))
grouped <- read.csv(file.path(run_dir, "predictions", "grouped_cv_predictions.csv"), check.names = FALSE)
data <- read_long_data(cfg)
metadata <- sample_metadata(data, cfg)

checks <- list(
  loso_prediction_rows = nrow(predictions) == nrow(metadata),
  loso_unique_samples = !anyDuplicated(predictions[[sid]]),
  loso_subjects_complete = setequal(predictions$heldout_subject, metadata[[subject]]),
  loso_fold_count = n_distinct(predictions$fold_id) == n_distinct(metadata[[subject]]),
  confusion_total = sum(confusion$Freq) == nrow(predictions),
  probability_rows_normalized = all(abs(rowSums(predictions[paste0("prob_", levels)]) - 1) < 1e-8),
  no_missing_predictions = !anyNA(predictions$predicted),
  fold_prediction_total = sum(fold_metrics$n_test) == nrow(predictions),
  grouped_each_sample_once_per_repeat = all(table(grouped[[sid]]) == ifelse(cfg$analysis$mode == "smoke", 1L, cfg$grouped_cv$repeats)),
  figure_formats_complete = all(vapply(c("png", "pdf", "svg"), function(extension) {
    length(list.files(file.path(run_dir, "figures", extension), pattern = paste0("[.]", extension, "$"))) >= 11L
  }, logical(1)))
)

recomputed <- metrics_from_prediction_table(predictions, levels)$metrics
common_metrics <- intersect(names(saved_metrics), names(recomputed))
numeric_agreement <- vapply(common_metrics, function(metric) {
  a <- saved_metrics[[metric]][[1]]
  b <- recomputed[[metric]][[1]]
  (is.na(a) && is.na(b)) || isTRUE(all.equal(a, b, tolerance = 1e-12))
}, logical(1))
checks$metric_csv_agreement <- all(numeric_agreement)

validation <- tibble(
  check = names(checks),
  passed = unlist(checks, use.names = FALSE),
  detail = c(
    paste(nrow(predictions), "predictions for", nrow(metadata), "samples"),
    paste(n_distinct(predictions[[sid]]), "unique sample IDs"),
    paste(n_distinct(predictions$heldout_subject), "held-out subjects"),
    paste(n_distinct(predictions$fold_id), "LOSO folds"),
    paste(sum(confusion$Freq), "confusion observations"),
    "probability tolerance 1e-8",
    paste(sum(is.na(predictions$predicted)), "missing predictions"),
    paste(sum(fold_metrics$n_test), "fold prediction total"),
    paste("mode", cfg$analysis$mode),
    "at least 11 figures in each required format",
    paste(common_metrics, collapse = ";")
  )
)
write.csv(validation, file.path(run_dir, "reports", "output_validation.csv"), row.names = FALSE)
if (any(!validation$passed)) {
  stop("Output validation failed: ", paste(validation$check[!validation$passed], collapse = ", "), call. = FALSE)
}
cat("All output-integrity checks passed.\n")
