# ============================================
# compare_mag_response_transformations.R
# ============================================
# Sensitivity analysis for the draft per-MAG LMM workflow.
#
# This script compares the current raw-abundance response against a
# log10-pseudocount transformed response. It does not modify the main model
# config, does not touch scripts/02_run_glmm.R, and does not overwrite the
# current per-MAG lmerTest result files.
#
# Rationale:
# - Raw microbiome abundance values can be sparse and skewed.
# - A log-pseudocount response can reduce the influence of large values and
#   offers a simple sensitivity check for the time association.
# - The pseudocount is defined once as half of the smallest non-zero Y value
#   across the processed gut dataset, then reported in all output rows.

library(lme4)
library(lmerTest)
library(here)

processed_path <- here("data", "processed", "cleaned_data.rds")
output_csv <- here("results", "tables", "per_mag_transformation_comparison.csv")
log_path <- here("results", "logs", "per_mag_transformation_comparison_log.txt")
summary_path <- here("results", "logs", "per_mag_transformation_comparison_summary.txt")

too_sparse_nonzero_samples <- 10
too_sparse_nonzero_subjects <- 3
low_variance_threshold <- 1e-10
p_value_method <- "lmerTest_Satterthwaite"
top_n <- 20

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

log_connection <- file(log_path, open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")

cat("Per-MAG Response Transformation Sensitivity Analysis\n")
cat("====================================================\n\n")
cat("Processed file:", processed_path, "\n")
cat("Results CSV:", output_csv, "\n")
cat("Log file:", log_path, "\n")
cat("Summary file:", summary_path, "\n")
cat("Model: response ~ days_from_baseline_subject + (1 | Subject)\n")
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

smallest_nonzero_y <- min(df$Y[df$Y > 0], na.rm = TRUE)
if (!is.finite(smallest_nonzero_y)) {
  stop("Cannot define pseudocount because no positive Y values were found.", call. = FALSE)
}
pseudocount <- smallest_nonzero_y / 2

cat("Smallest non-zero Y:", format(smallest_nonzero_y, scientific = TRUE), "\n")
cat("Pseudocount:", format(pseudocount, scientific = TRUE), "\n\n")

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

fit_mag_response <- function(mag_id, response_type, index, total) {
  mag_df <- df[df$MAG == mag_id, , drop = FALSE]
  mag_df <- mag_df[stats::complete.cases(mag_df[, c("Y", "days_from_baseline_subject", "Subject")]), , drop = FALSE]
  mag_df$Subject <- droplevels(factor(mag_df$Subject))

  if (response_type == "raw") {
    mag_df$response_value <- mag_df$Y
  } else if (response_type == "log10_pseudocount") {
    mag_df$response_value <- log10(mag_df$Y + pseudocount)
  } else {
    stop("Unknown response_type: ", response_type, call. = FALSE)
  }

  nonzero <- mag_df$Y != 0 & !is.na(mag_df$Y)
  subject_nonzero <- tapply(nonzero, mag_df$Subject, any)
  response_variance <- stats::var(mag_df$response_value, na.rm = TRUE)

  result <- data.frame(
    MAG = mag_id,
    response_type = response_type,
    pseudocount = if (response_type == "log10_pseudocount") pseudocount else NA_real_,
    n_rows = nrow(mag_df),
    n_subjects = length(unique(mag_df$Subject)),
    n_sampleID = length(unique(mag_df$sampleID)),
    n_timepoints = length(unique(mag_df$Timepoint)),
    fit_succeeded = FALSE,
    status = "",
    singular = NA,
    random_intercept_variance = NA_real_,
    zero_fraction = mean(mag_df$Y == 0, na.rm = TRUE),
    nonzero_sample_count = sum(nonzero),
    subjects_with_nonzero = sum(subject_nonzero, na.rm = TRUE),
    y_variance = response_variance,
    time_estimate = NA_real_,
    time_std_error = NA_real_,
    time_df = NA_real_,
    time_statistic = NA_real_,
    time_p_value = NA_real_,
    time_p_adj_bh = NA_real_,
    p_value_method = p_value_method,
    warnings = "",
    error = "",
    too_sparse = sum(nonzero) < too_sparse_nonzero_samples ||
      sum(subject_nonzero, na.rm = TRUE) < too_sparse_nonzero_subjects,
    low_variance = is.na(response_variance) || response_variance < low_variance_threshold,
    stringsAsFactors = FALSE
  )

  warning_messages <- character(0)
  error_messages <- character(0)

  cat(sprintf("[%d/%d] Modeling %s with %s response\n", index, total, mag_id, response_type))

  fit <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(
        response_value ~ days_from_baseline_subject + (1 | Subject),
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
response_types <- c("raw", "log10_pseudocount")
fit_grid <- expand.grid(
  MAG = mag_ids,
  response_type = response_types,
  stringsAsFactors = FALSE
)

results <- do.call(
  rbind,
  Map(
    fit_mag_response,
    fit_grid$MAG,
    fit_grid$response_type,
    seq_len(nrow(fit_grid)),
    MoreArgs = list(total = nrow(fit_grid))
  )
)

results$time_p_adj_bh <- NA_real_
for (response_type in response_types) {
  response_rows <- results$response_type == response_type
  results$time_p_adj_bh[response_rows] <- stats::p.adjust(
    results$time_p_value[response_rows],
    method = "BH"
  )
}

utils::write.csv(results, output_csv, row.names = FALSE)

sink(type = "message")
sink()
close(log_connection)

summary_connection <- file(summary_path, open = "wt")
sink(summary_connection, split = TRUE)

cat("Per-MAG Response Transformation Sensitivity Summary\n")
cat("===================================================\n\n")
cat("Results CSV:", output_csv, "\n")
cat("Log file:", log_path, "\n")
cat("P-value method:", p_value_method, "\n")
cat("Pseudocount definition: half of the smallest non-zero Y across the processed gut dataset\n")
cat("Smallest non-zero Y:", format(smallest_nonzero_y, scientific = TRUE), "\n")
cat("Pseudocount:", format(pseudocount, scientific = TRUE), "\n\n")

for (response_type in response_types) {
  response_results <- results[results$response_type == response_type, , drop = FALSE]

  cat("Response type:", response_type, "\n")
  cat(strrep("-", nchar("Response type: ") + nchar(response_type)), "\n", sep = "")
  cat("MAGs modeled:", nrow(response_results), "\n")
  cat("Successful models:", sum(response_results$fit_succeeded %in% TRUE, na.rm = TRUE), "\n")
  cat("Failed models:", sum(response_results$status == "model_failed", na.rm = TRUE), "\n")
  cat("lmm_ok:", sum(response_results$status == "lmm_ok", na.rm = TRUE), "\n")
  cat("Singular models:", sum(response_results$status == "lmm_singular", na.rm = TRUE), "\n")
  cat("Too-sparse MAGs:", sum(response_results$status == "too_sparse", na.rm = TRUE), "\n")
  cat("Low-variance MAGs:", sum(response_results$status == "low_variance", na.rm = TRUE), "\n")
  cat("Raw p < 0.05:", sum(response_results$time_p_value < 0.05, na.rm = TRUE), "\n")
  cat("BH q < 0.05:", sum(response_results$time_p_adj_bh < 0.05, na.rm = TRUE), "\n\n")

  cat("Top 20 MAGs by BH q-value:\n")
  ranked <- response_results[order(
    is.na(response_results$time_p_adj_bh),
    response_results$time_p_adj_bh,
    response_results$time_p_value
  ), ]
  top_cols <- c(
    "MAG", "response_type", "status", "time_estimate", "time_std_error",
    "time_df", "time_statistic", "time_p_value", "time_p_adj_bh",
    "zero_fraction", "y_variance", "singular"
  )
  print(utils::head(ranked[, top_cols], top_n), row.names = FALSE)
  cat("\n")
}

raw_ranked <- results[results$response_type == "raw", , drop = FALSE]
raw_ranked <- raw_ranked[order(is.na(raw_ranked$time_p_adj_bh), raw_ranked$time_p_adj_bh, raw_ranked$time_p_value), ]
log_ranked <- results[results$response_type == "log10_pseudocount", , drop = FALSE]
log_ranked <- log_ranked[order(is.na(log_ranked$time_p_adj_bh), log_ranked$time_p_adj_bh, log_ranked$time_p_value), ]

raw_top <- utils::head(raw_ranked$MAG, top_n)
log_top <- utils::head(log_ranked$MAG, top_n)
top_overlap <- intersect(raw_top, log_top)

cat("Top-MAG overlap\n")
cat("---------------\n")
cat("Top-N used:", top_n, "\n")
cat("Raw top MAG count:", length(raw_top), "\n")
cat("Log-pseudocount top MAG count:", length(log_top), "\n")
cat("Overlap count:", length(top_overlap), "\n")
cat("Overlap MAGs:\n")
if (length(top_overlap) > 0) {
  cat(paste(top_overlap, collapse = "\n"), "\n", sep = "")
} else {
  cat("None\n")
}

cat("\nImportant warning:\n")
cat("This is a sensitivity analysis only. Do not treat either response model as final until response transformation, singular-fit handling, sparsity, and biological interpretation are reviewed.\n")

cat("\nSaved summary to results/logs/per_mag_transformation_comparison_summary.txt\n")

sink()
close(summary_connection)
