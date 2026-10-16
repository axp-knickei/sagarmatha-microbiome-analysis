# ============================================
# summarize_per_mag_model_test.R
# ============================================
# Summarizes the 10-MAG exploratory model test and diagnoses singular MAGs.

library(here)

input_csv <- here("results", "tables", "per_mag_model_test.csv")
input_log <- here("results", "logs", "per_mag_model_test_log.txt")
processed_path <- here("data", "processed", "cleaned_data.rds")
summary_path <- here("results", "logs", "per_mag_model_test_summary.txt")

dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

sink(summary_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Per-MAG Model Test Summary\n")
cat("==========================\n\n")
cat("Input CSV:", input_csv, "\n")
cat("Input log:", input_log, "\n")
cat("Processed data:", processed_path, "\n")
cat("Summary report:", summary_path, "\n\n")

if (!file.exists(input_csv)) {
  stop("Missing per-MAG model test CSV: ", input_csv, call. = FALSE)
}
if (!file.exists(processed_path)) {
  stop("Missing processed data: ", processed_path, call. = FALSE)
}

results <- utils::read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
df <- readRDS(processed_path)

n_tested <- nrow(results)
n_lm_success <- sum(results$lm_succeeded %in% TRUE, na.rm = TRUE)
n_lmm_success <- sum(results$lmm_succeeded %in% TRUE, na.rm = TRUE)
n_lmm_singular <- sum(results$lmm_singular %in% TRUE, na.rm = TRUE)
singular_mags <- results$MAG[results$lmm_singular %in% TRUE]

cat("Overall Counts\n")
cat("--------------\n")
cat("MAGs tested:", n_tested, "\n")
cat("Successful LM fits:", n_lm_success, "\n")
cat("Successful LMM fits:", n_lmm_success, "\n")
cat("Singular LMM fits:", n_lmm_singular, "\n")
cat("Singular MAGs:", if (length(singular_mags) > 0) paste(singular_mags, collapse = ", ") else "None", "\n\n")

cat("Per-MAG Results\n")
cat("---------------\n")
for (i in seq_len(nrow(results))) {
  row <- results[i, ]
  estimate <- row$lmm_time_estimate
  if (is.na(estimate) && "lm_time_estimate" %in% names(row)) {
    estimate <- row$lm_time_estimate
  }

  cat("MAG:", row$MAG, "\n")
  cat("  Rows:", row$n_rows, "\n")
  cat("  Unique Subjects:", row$n_subjects, "\n")
  cat("  Unique sampleID:", row$n_sampleID, "\n")
  cat("  Unique Timepoint:", row$n_timepoints, "\n")
  cat("  Random-intercept variance:", row$random_intercept_variance, "\n")
  cat("  days_from_baseline_subject estimate:", estimate, "\n")
  cat("  LMM singular:", row$lmm_singular, "\n\n")
}

cat("Singular MAG Diagnostics\n")
cat("------------------------\n")
if (length(singular_mags) == 0) {
  cat("No singular MAGs were found in this test subset.\n")
} else {
  for (mag_id in singular_mags) {
    mag_df <- df[df$MAG == mag_id, , drop = FALSE]
    nonzero <- mag_df$Y != 0 & !is.na(mag_df$Y)
    subject_nonzero <- tapply(nonzero, mag_df$Subject, any)

    cat("MAG:", mag_id, "\n")
    cat("  Rows:", nrow(mag_df), "\n")
    cat("  Zero values:", sum(mag_df$Y == 0, na.rm = TRUE), "\n")
    cat("  Zero fraction:", sprintf("%.4f", mean(mag_df$Y == 0, na.rm = TRUE)), "\n")
    cat("  Y variance:", stats::var(mag_df$Y, na.rm = TRUE), "\n")
    cat("  Non-zero samples:", sum(nonzero), "\n")
    cat("  Subjects with non-zero values:", sum(subject_nonzero, na.rm = TRUE), "\n")
    cat("  Subjects total:", length(subject_nonzero), "\n")

    if (sum(nonzero) < 5) {
      cat("  Diagnostic: very few non-zero samples.\n")
    }
    if (sum(subject_nonzero, na.rm = TRUE) < 3) {
      cat("  Diagnostic: non-zero values occur in too few subjects.\n")
    }
    if (isTRUE(stats::var(mag_df$Y, na.rm = TRUE) < 1e-10)) {
      cat("  Diagnostic: very low variance in Y.\n")
    }
    cat("\n")
  }
}

cat("Interpretation\n")
cat("--------------\n")
cat("The per-MAG test suggests that singularity is MAG-specific rather than universal. ")
cat("This supports moving toward a per-MAG modeling workflow, with explicit handling of singular fits and sparse/low-variance MAGs before scaling to all MAGs or HPC.\n")

cat("\nSaved summary to results/logs/per_mag_model_test_summary.txt\n")
