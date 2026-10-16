# ============================================
# test_per_mag_models.R - Per-MAG Model Structure Check
# ============================================
# Exploratory only. Does not modify config/config.yml or overwrite the main
# fitted model. Fits LM and LMM candidates for a small MAG subset.

library(lme4)
library(here)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

processed_path <- here("data", "processed", "cleaned_data.rds")
output_csv <- here("results", "tables", "per_mag_model_test.csv")
log_path <- here("results", "logs", "per_mag_model_test_log.txt")

max_mags <- 10

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)

log_connection <- file(log_path, open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_connection)
}, add = TRUE)

cat("Per-MAG Model Structure Test\n")
cat("============================\n\n")
cat("Processed file:", processed_path, "\n")
cat("Output CSV:", output_csv, "\n")
cat("Log file:", log_path, "\n")
cat("MAG subset size:", max_mags, "\n\n")

if (!file.exists(processed_path)) {
  stop("Processed data not found: ", processed_path, call. = FALSE)
}

df <- readRDS(processed_path)
required_cols <- c("MAG", "Y", "days_from_baseline_subject", "Subject", "sampleID", "Timepoint")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Processed data is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

mag_counts <- sort(table(df$MAG), decreasing = TRUE)
selected_mags <- names(mag_counts)[seq_len(min(max_mags, length(mag_counts)))]

cat("Selected MAGs:\n")
print(selected_mags)
cat("\n")

fit_one_mag <- function(mag_id) {
  mag_df <- df[df$MAG == mag_id, , drop = FALSE]
  mag_df <- mag_df[stats::complete.cases(mag_df[, c("Y", "days_from_baseline_subject", "Subject")]), , drop = FALSE]
  mag_df$Subject <- droplevels(factor(mag_df$Subject))

  result <- data.frame(
    MAG = mag_id,
    n_rows = nrow(mag_df),
    n_subjects = length(unique(mag_df$Subject)),
    n_sampleID = length(unique(mag_df$sampleID)),
    n_timepoints = length(unique(mag_df$Timepoint)),
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
    stringsAsFactors = FALSE
  )

  warning_messages <- character(0)
  error_messages <- character(0)

  cat("Testing MAG:", mag_id, "\n")
  cat("  Rows:", result$n_rows, "\n")
  cat("  Subjects:", result$n_subjects, "\n")
  cat("  sampleID:", result$n_sampleID, "\n")
  cat("  Timepoints:", result$n_timepoints, "\n")

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

  cat("  LM succeeded:", result$lm_succeeded, "\n")
  cat("  LMM succeeded:", result$lmm_succeeded, "\n")
  cat("  LMM singular:", result$lmm_singular, "\n")
  cat("  Random-intercept variance:", result$random_intercept_variance, "\n")
  if (nzchar(result$warnings)) cat("  Warnings:", result$warnings, "\n")
  if (nzchar(result$error)) cat("  Errors:", result$error, "\n")
  cat("\n")

  result
}

results <- do.call(rbind, lapply(selected_mags, fit_one_mag))

utils::write.csv(results, output_csv, row.names = FALSE)

cat("Summary\n")
cat("-------\n")
cat("MAGs tested:", nrow(results), "\n")
cat("LM fits succeeded:", sum(results$lm_succeeded), "\n")
cat("LMM fits succeeded:", sum(results$lmm_succeeded), "\n")
cat("LMM singular fits:", sum(results$lmm_singular %in% TRUE, na.rm = TRUE), "\n")
cat("Results saved to results/tables/per_mag_model_test.csv\n")
cat("Log saved to results/logs/per_mag_model_test_log.txt\n")
