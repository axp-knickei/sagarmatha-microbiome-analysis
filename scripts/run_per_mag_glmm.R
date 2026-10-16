# ============================================
# run_per_mag_glmm.R - Draft Per-MAG LMM Workflow
# ============================================
# Draft final-modeling script. Fits one LMM per MAG without modifying
# config/config.yml or overwriting results/model_fitted.rds.

library(lme4)
library(here)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

processed_path <- here("data", "processed", "cleaned_data.rds")
output_csv <- here("results", "tables", "per_mag_glmm_results.csv")
log_path <- here("results", "logs", "per_mag_glmm_log.txt")
summary_path <- here("results", "logs", "per_mag_glmm_summary.txt")

too_sparse_nonzero_samples <- 10
too_sparse_nonzero_subjects <- 3
low_variance_threshold <- 1e-10

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

log_connection <- file(log_path, open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_connection)
}, add = TRUE)

cat("Draft Per-MAG LMM Workflow\n")
cat("==========================\n\n")
cat("Processed file:", processed_path, "\n")
cat("Results CSV:", output_csv, "\n")
cat("Log file:", log_path, "\n")
cat("Summary file:", summary_path, "\n")
cat("Formula: Y ~ days_from_baseline_subject + (1 | Subject)\n\n")

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
    time_statistic = NA_real_,
    time_p_value = NA_real_,
    p_value_note = "lme4::lmer does not provide p-values by default; p-value method needs a decision.",
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
      lme4::lmer(
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
      statistic_column <- intersect(c("t value", "z value"), colnames(coef_table))
      if (length(statistic_column) > 0) {
        result$time_statistic <- coef_table["days_from_baseline_subject", statistic_column[[1]]]
      }
      p_column <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(coef_table))
      if (length(p_column) > 0) {
        result$time_p_value <- coef_table["days_from_baseline_subject", p_column[[1]]]
        result$p_value_note <- ""
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

if (any(!is.na(results$time_p_value))) {
  results$time_p_adj_fdr <- stats::p.adjust(results$time_p_value, method = "fdr")
} else {
  results$time_p_adj_fdr <- NA_real_
}

utils::write.csv(results, output_csv, row.names = FALSE)

sink(type = "message")
sink()
close(log_connection)

summary_connection <- file(summary_path, open = "wt")
sink(summary_connection, split = TRUE)
on.exit({
  sink()
  close(summary_connection)
}, add = TRUE)

cat("Draft Per-MAG LMM Summary\n")
cat("=========================\n\n")
cat("Results CSV:", output_csv, "\n")
cat("Log file:", log_path, "\n")
cat("MAGs modeled:", nrow(results), "\n")
cat("Successful models:", sum(results$fit_succeeded %in% TRUE, na.rm = TRUE), "\n")
cat("Singular models:", sum(results$singular %in% TRUE, na.rm = TRUE), "\n")
cat("Too-sparse MAGs:", sum(results$status == "too_sparse", na.rm = TRUE), "\n")
cat("Low-variance MAGs:", sum(results$status == "low_variance", na.rm = TRUE), "\n")
cat("Failed models:", sum(results$status == "model_failed", na.rm = TRUE), "\n\n")

if (any(!is.na(results$time_p_value))) {
  cat("Significant MAGs before FDR (p < 0.05):", sum(results$time_p_value < 0.05, na.rm = TRUE), "\n")
  cat("Significant MAGs after FDR (q < 0.05):", sum(results$time_p_adj_fdr < 0.05, na.rm = TRUE), "\n")
} else {
  cat("Significant MAGs before FDR: not available because lme4::lmer does not provide p-values by default.\n")
  cat("Significant MAGs after FDR: not available because p-values are NA.\n")
}

cat("\nStatus counts:\n")
print(table(results$status))

cat("\nImportant note:\n")
cat("This is a draft final-modeling script. P-value method, response transformation, singular-fit handling, and biological interpretation still need review before treating results as final or scaling to HPC.\n")

cat("\nSaved summary to results/logs/per_mag_glmm_summary.txt\n")
