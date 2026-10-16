# ============================================
# run_per_mag_lmertest.R - Draft Per-MAG LMM Inference
# ============================================
# Draft inference workflow. Fits one lmerTest LMM per MAG and extracts
# Satterthwaite approximate p-values for the time fixed effect.
#
# Rationale:
# - lme4::lmer() does not provide fixed-effect p-values by default.
# - lmerTest extends lme4 by providing fixed-effect p-values using the
#   Satterthwaite degrees-of-freedom approximation.
# - Benjamini-Hochberg FDR correction is used because many MAG/features are
#   tested. This is similar in spirit to microbiome per-feature association
#   workflows such as MaAsLin2, which perform feature-wise association testing
#   and multiple hypothesis correction while supporting repeated measures or
#   mixed-model designs.

library(lme4)
library(lmerTest)
library(here)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

processed_path <- here("data", "processed", "cleaned_data.rds")
output_csv <- here("results", "tables", "per_mag_lmertest_results.csv")
log_path <- here("results", "logs", "per_mag_lmertest_log.txt")
summary_path <- here("results", "logs", "per_mag_lmertest_summary.txt")

too_sparse_nonzero_samples <- 10
too_sparse_nonzero_subjects <- 3
low_variance_threshold <- 1e-10
p_value_method <- "lmerTest_Satterthwaite"

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

log_connection <- file(log_path, open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")

cat("Draft Per-MAG lmerTest LMM Workflow\n")
cat("===================================\n\n")
cat("Processed file:", processed_path, "\n")
cat("Results CSV:", output_csv, "\n")
cat("Log file:", log_path, "\n")
cat("Summary file:", summary_path, "\n")
cat("Formula: Y ~ days_from_baseline_subject + (1 | Subject)\n")
cat("P-value method:", p_value_method, "\n\n")

if (!file.exists(processed_path)) {
  stop("Processed data not found: ", processed_path, call. = FALSE)
}

df <- readRDS(processed_path)
required_cols <- c("MAG", "Y", "days_from_baseline_subject", "Subject", "sampleID", "Timepoint")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Processed data is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

status_for_result <- function(result) {
  if (!result$fit_succeeded) {
    return("model_failed")
  }
  if (isTRUE(result$too_sparse)) {
    return("too_sparse")
  }
  if (isTRUE(result$low_variance)) {
    return("low_variance")
  }
  if (isTRUE(result$singular)) {
    return("lmm_singular")
  }
  "lmm_ok"
}

fit_mag <- function(mag_id, index, total) {
  mag_df <- df[df$MAG == mag_id, , drop = FALSE]
  mag_df <- mag_df[stats::complete.cases(mag_df[, c("Y", "days_from_baseline_subject", "Subject")]), , drop = FALSE]
  mag_df$Subject <- droplevels(factor(mag_df$Subject))

  nonzero <- mag_df$Y != 0 & !is.na(mag_df$Y)
  subject_nonzero <- tapply(nonzero, mag_df$Subject, any)
  y_variance <- stats::var(mag_df$Y, na.rm = TRUE)

  result <- data.frame(
    MAG = mag_id,
    n_rows = nrow(mag_df),
    n_subjects = length(unique(mag_df$Subject)),
    n_sampleID = length(unique(mag_df$sampleID)),
    n_timepoints = length(unique(mag_df$Timepoint)),
    zero_fraction = mean(mag_df$Y == 0, na.rm = TRUE),
    nonzero_sample_count = sum(nonzero),
    subjects_with_nonzero = sum(subject_nonzero, na.rm = TRUE),
    y_variance = y_variance,
    fit_succeeded = FALSE,
    singular = NA,
    random_intercept_variance = NA_real_,
    time_estimate = NA_real_,
    time_std_error = NA_real_,
    time_df = NA_real_,
    time_statistic = NA_real_,
    time_p_value = NA_real_,
    p_value_method = p_value_method,
    warnings = "",
    error = "",
    too_sparse = sum(nonzero) < too_sparse_nonzero_samples ||
      sum(subject_nonzero, na.rm = TRUE) < too_sparse_nonzero_subjects,
    low_variance = is.na(y_variance) || y_variance < low_variance_threshold,
    status = "",
    stringsAsFactors = FALSE
  )

  warning_messages <- character(0)
  error_messages <- character(0)

  cat(sprintf("[%d/%d] Modeling MAG: %s\n", index, total, mag_id))

  fit <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(
        Y ~ days_from_baseline_subject + (1 | Subject),
        data = mag_df,
        control = lme4::lmerControl(optimizer = "bobyqa")
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_messages <<- c(error_messages, conditionMessage(e))
      NULL
    }
  )

  if (!is.null(fit)) {
    result$fit_succeeded <- TRUE
    result$singular <- lme4::isSingular(fit)

    vc <- as.data.frame(lme4::VarCorr(fit))
    subject_intercept <- subset(vc, grp == "Subject" & var1 == "(Intercept)" & is.na(var2))
    if (nrow(subject_intercept) > 0) {
      result$random_intercept_variance <- subject_intercept$vcov[[1]]
    }

    coef_table <- summary(fit)$coefficients
    if ("days_from_baseline_subject" %in% rownames(coef_table)) {
      result$time_estimate <- coef_table["days_from_baseline_subject", "Estimate"]
      result$time_std_error <- coef_table["days_from_baseline_subject", "Std. Error"]
      if ("df" %in% colnames(coef_table)) {
        result$time_df <- coef_table["days_from_baseline_subject", "df"]
      }
      statistic_column <- intersect(c("t value", "z value"), colnames(coef_table))
      if (length(statistic_column) > 0) {
        result$time_statistic <- coef_table["days_from_baseline_subject", statistic_column[[1]]]
      }
      p_column <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(coef_table))
      if (length(p_column) > 0) {
        result$time_p_value <- coef_table["days_from_baseline_subject", p_column[[1]]]
      }
    }
  }

  result$warnings <- paste(unique(warning_messages), collapse = " | ")
  result$error <- paste(unique(error_messages), collapse = " | ")
  result$status <- status_for_result(result)

  cat("  Status:", result$status, "\n")
  if (nzchar(result$warnings)) cat("  Warnings:", result$warnings, "\n")
  if (nzchar(result$error)) cat("  Errors:", result$error, "\n")

  result
}

mag_ids <- sort(unique(df$MAG))
results <- do.call(
  rbind,
  Map(fit_mag, mag_ids, seq_along(mag_ids), MoreArgs = list(total = length(mag_ids)))
)

results$time_p_adj_bh <- stats::p.adjust(results$time_p_value, method = "BH")
results$primary_analysis <- results$status == "lmm_ok"
results$time_p_adj_bh_primary <- NA_real_
results$time_p_adj_bh_primary[results$primary_analysis] <- stats::p.adjust(
  results$time_p_value[results$primary_analysis],
  method = "BH"
)

utils::write.csv(results, output_csv, row.names = FALSE)

sink(type = "message")
sink()
close(log_connection)

summary_connection <- file(summary_path, open = "wt")
sink(summary_connection, split = TRUE)

cat("Draft Per-MAG lmerTest LMM Summary\n")
cat("==================================\n\n")
cat("Results CSV:", output_csv, "\n")
cat("Log file:", log_path, "\n")
cat("MAGs modeled:", nrow(results), "\n")
cat("Successful models:", sum(results$fit_succeeded %in% TRUE, na.rm = TRUE), "\n")
cat("Failed models:", sum(results$status == "model_failed", na.rm = TRUE), "\n")
cat("lmm_ok:", sum(results$status == "lmm_ok", na.rm = TRUE), "\n")
cat("lmm_singular:", sum(results$status == "lmm_singular", na.rm = TRUE), "\n")
cat("too_sparse:", sum(results$status == "too_sparse", na.rm = TRUE), "\n")
cat("low_variance:", sum(results$status == "low_variance", na.rm = TRUE), "\n")
cat("P-value method:", p_value_method, "\n\n")

cat("MAGs with raw p < 0.05:", sum(results$time_p_value < 0.05, na.rm = TRUE), "\n")
cat("MAGs with BH q < 0.05 across all modeled MAGs:", sum(results$time_p_adj_bh < 0.05, na.rm = TRUE), "\n")
cat(
  "MAGs with BH q < 0.05 in primary-analysis set:",
  sum(results$time_p_adj_bh_primary < 0.05, na.rm = TRUE),
  "\n\n"
)

cat("Top 20 MAGs ranked by primary BH q-value, then raw p-value:\n")
ranked <- results[order(
  is.na(results$time_p_adj_bh_primary),
  results$time_p_adj_bh_primary,
  results$time_p_value
), ]
top_cols <- c(
  "MAG", "status", "time_estimate", "time_std_error", "time_df",
  "time_statistic", "time_p_value", "time_p_adj_bh",
  "time_p_adj_bh_primary", "zero_fraction", "nonzero_sample_count",
  "subjects_with_nonzero", "singular"
)
print(utils::head(ranked[, top_cols], 20), row.names = FALSE)

cat("\nImportant warning:\n")
cat("These results are still draft. Response transformation, singular-fit handling, p-value method assumptions, and biological interpretation must be reviewed before treating results as final or scaling to HPC.\n")

cat("\nSaved summary to results/logs/per_mag_lmertest_summary.txt\n")

sink()
close(summary_connection)
