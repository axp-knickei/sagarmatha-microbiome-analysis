# ============================================
# screen_all_mag_models.R - All-MAG Model Screen
# ============================================
# Diagnostic only. Screens per-MAG LM/LMM model structure without changing
# config/config.yml or overwriting results/model_fitted.rds.

library(lme4)
library(here)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

processed_path <- here("data", "processed", "cleaned_data.rds")
output_csv <- here("results", "tables", "all_mag_model_screen.csv")
log_path <- here("results", "logs", "all_mag_model_screen_log.txt")
summary_path <- here("results", "logs", "all_mag_model_screen_summary.txt")

# Use NA to screen all MAGs. Set an integer here for a faster local subset.
max_mags <- NA_integer_

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

cat("All-MAG Model Screen\n")
cat("====================\n\n")
cat("Processed file:", processed_path, "\n")
cat("Output CSV:", output_csv, "\n")
cat("Log file:", log_path, "\n")
cat("Summary file:", summary_path, "\n")
cat("MAG limit:", ifelse(is.na(max_mags), "all", max_mags), "\n\n")

if (!file.exists(processed_path)) {
  stop("Processed data not found: ", processed_path, call. = FALSE)
}

df <- readRDS(processed_path)
required_cols <- c("MAG", "Y", "days_from_baseline_subject", "Subject", "sampleID", "Timepoint")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Processed data is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

mag_ids <- sort(unique(df$MAG))
if (!is.na(max_mags)) {
  mag_ids <- head(mag_ids, max_mags)
}

status_for_result <- function(result) {
  if (!result$lm_succeeded && !result$lmm_succeeded) {
    return("model_failed")
  }
  if (isTRUE(result$too_sparse)) {
    return("too_sparse")
  }
  if (isTRUE(result$low_variance)) {
    return("low_variance")
  }
  if (isTRUE(result$lmm_succeeded) && isTRUE(result$lmm_singular)) {
    return("lmm_singular")
  }
  if (isTRUE(result$lmm_succeeded) && isFALSE(result$lmm_singular)) {
    return("lmm_ok")
  }
  "model_failed"
}

fit_one_mag <- function(mag_id, index, total) {
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
    zero_count = sum(mag_df$Y == 0, na.rm = TRUE),
    zero_fraction = mean(mag_df$Y == 0, na.rm = TRUE),
    nonzero_sample_count = sum(nonzero),
    subjects_with_nonzero = sum(subject_nonzero, na.rm = TRUE),
    y_variance = y_variance,
    lm_succeeded = FALSE,
    lmm_succeeded = FALSE,
    lmm_singular = NA,
    random_intercept_variance = NA_real_,
    lm_time_estimate = NA_real_,
    lm_time_p_value = NA_real_,
    lmm_time_estimate = NA_real_,
    lmm_time_p_value = NA_real_,
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

  cat(sprintf("[%d/%d] Testing MAG: %s\n", index, total, mag_id))

  lm_fit <- tryCatch(
    withCallingHandlers(
      stats::lm(Y ~ days_from_baseline_subject, data = mag_df),
      warning = function(w) {
        warning_messages <<- c(warning_messages, paste("LM:", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_messages <<- c(error_messages, paste("LM:", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(lm_fit)) {
    result$lm_succeeded <- TRUE
    lm_coef <- summary(lm_fit)$coefficients
    if ("days_from_baseline_subject" %in% rownames(lm_coef)) {
      result$lm_time_estimate <- lm_coef["days_from_baseline_subject", "Estimate"]
      result$lm_time_p_value <- lm_coef["days_from_baseline_subject", "Pr(>|t|)"]
    }
  }

  lmm_fit <- tryCatch(
    withCallingHandlers(
      lme4::lmer(
        Y ~ days_from_baseline_subject + (1 | Subject),
        data = mag_df,
        control = lme4::lmerControl(optimizer = "bobyqa")
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, paste("LMM:", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_messages <<- c(error_messages, paste("LMM:", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(lmm_fit)) {
    result$lmm_succeeded <- TRUE
    result$lmm_singular <- lme4::isSingular(lmm_fit)
    vc <- as.data.frame(lme4::VarCorr(lmm_fit))
    subject_intercept <- subset(vc, grp == "Subject" & var1 == "(Intercept)" & is.na(var2))
    if (nrow(subject_intercept) > 0) {
      result$random_intercept_variance <- subject_intercept$vcov[[1]]
    }

    lmm_coef <- summary(lmm_fit)$coefficients
    if ("days_from_baseline_subject" %in% rownames(lmm_coef)) {
      result$lmm_time_estimate <- lmm_coef["days_from_baseline_subject", "Estimate"]
      if ("Pr(>|t|)" %in% colnames(lmm_coef)) {
        result$lmm_time_p_value <- lmm_coef["days_from_baseline_subject", "Pr(>|t|)"]
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

results <- do.call(
  rbind,
  Map(fit_one_mag, mag_ids, seq_along(mag_ids), MoreArgs = list(total = length(mag_ids)))
)

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

cat("All-MAG Model Screen Summary\n")
cat("============================\n\n")
cat("Output CSV:", output_csv, "\n")
cat("Log file:", log_path, "\n")
cat("MAGs screened:", nrow(results), "\n\n")

status_counts <- table(results$status)
status_names <- c("lmm_ok", "lmm_singular", "too_sparse", "low_variance", "model_failed")
for (status_name in status_names) {
  count <- if (status_name %in% names(status_counts)) unname(status_counts[[status_name]]) else 0
  pct <- if (nrow(results) > 0) 100 * count / nrow(results) else NA_real_
  cat(sprintf("%s: %d (%.1f%%)\n", status_name, count, pct))
}

cat("\nFit counts:\n")
cat("LM fits succeeded:", sum(results$lm_succeeded %in% TRUE, na.rm = TRUE), "\n")
cat("LMM fits succeeded:", sum(results$lmm_succeeded %in% TRUE, na.rm = TRUE), "\n")
cat("LMM singular fits:", sum(results$lmm_singular %in% TRUE, na.rm = TRUE), "\n")

cat("\nInterpretation:\n")
cat("This screen is diagnostic only. Use it to decide how to handle per-MAG singularity, sparse MAGs, and low-variance features before scaling to HPC.\n")

cat("\nSaved summary to results/logs/all_mag_model_screen_summary.txt\n")
