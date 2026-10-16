# ============================================
# diagnose_response_scale_models.R
# ============================================
# Response-scale and residual diagnostics for the draft per-MAG LMM workflow.
#
# This script compares raw abundance and log10-pseudocount response scales with
# global distribution plots, per-MAG trajectory plots, and residual diagnostics
# for top-ranked MAGs from each response scale. It does not modify the main
# config, scripts/02_run_glmm.R, previous result files, or raw data.

library(ggplot2)
library(lme4)
library(lmerTest)
library(here)

processed_path <- here("data", "processed", "cleaned_data.rds")
comparison_path <- here("results", "tables", "per_mag_transformation_comparison.csv")
summary_path <- here("results", "logs", "response_scale_diagnostics_summary.txt")
diagnostics_csv <- here("results", "tables", "response_scale_top_mag_diagnostics.csv")
figure_dir <- here("results", "figures", "response_scale")

top_n <- 10
p_value_method <- "lmerTest_Satterthwaite"

dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(diagnostics_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(processed_path)) {
  stop("Processed data not found: ", processed_path, call. = FALSE)
}
if (!file.exists(comparison_path)) {
  stop("Transformation comparison file not found: ", comparison_path, call. = FALSE)
}

df <- readRDS(processed_path)
comparison <- read.csv(comparison_path, stringsAsFactors = FALSE)

required_processed_cols <- c("MAG", "Y", "days_from_baseline_subject", "Subject", "sampleID")
missing_processed_cols <- setdiff(required_processed_cols, names(df))
if (length(missing_processed_cols) > 0) {
  stop("Processed data is missing required columns: ", paste(missing_processed_cols, collapse = ", "), call. = FALSE)
}

required_comparison_cols <- c(
  "MAG", "response_type", "status", "singular", "zero_fraction",
  "nonzero_sample_count", "subjects_with_nonzero", "y_variance",
  "time_estimate", "time_p_value", "time_p_adj_bh"
)
missing_comparison_cols <- setdiff(required_comparison_cols, names(comparison))
if (length(missing_comparison_cols) > 0) {
  stop("Comparison file is missing required columns: ", paste(missing_comparison_cols, collapse = ", "), call. = FALSE)
}

smallest_nonzero_y <- min(df$Y[df$Y > 0], na.rm = TRUE)
if (!is.finite(smallest_nonzero_y)) {
  stop("Cannot define pseudocount because no positive Y values were found.", call. = FALSE)
}
pseudocount <- smallest_nonzero_y / 2
df$log10_pseudocount_response <- log10(df$Y + pseudocount)

safe_filename <- function(x) {
  safe <- gsub("[^A-Za-z0-9._-]+", "_", x)
  safe <- gsub("_+", "_", safe)
  safe <- gsub("^_|_$", "", safe)
  ifelse(nchar(safe) == 0, "unnamed_MAG", safe)
}

save_plot <- function(plot, filename, width = 8, height = 5) {
  ggplot2::ggsave(
    filename = file.path(figure_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 150
  )
  ggplot2::ggsave(
    filename = file.path(figure_dir, sub("\\.[^.]+$", ".svg", filename)),
    plot = plot,
    width = width,
    height = height,
    device = svglite::svglite
  )
}

global_summary <- data.frame(
  metric = c(
    "min_raw_Y",
    "max_raw_Y",
    "mean_raw_Y",
    "median_raw_Y",
    "zero_fraction",
    "smallest_nonzero_Y",
    "pseudocount",
    "min_log10_response",
    "max_log10_response",
    "mean_log10_response",
    "median_log10_response"
  ),
  value = c(
    min(df$Y, na.rm = TRUE),
    max(df$Y, na.rm = TRUE),
    mean(df$Y, na.rm = TRUE),
    stats::median(df$Y, na.rm = TRUE),
    mean(df$Y == 0, na.rm = TRUE),
    smallest_nonzero_y,
    pseudocount,
    min(df$log10_pseudocount_response, na.rm = TRUE),
    max(df$log10_pseudocount_response, na.rm = TRUE),
    mean(df$log10_pseudocount_response, na.rm = TRUE),
    stats::median(df$log10_pseudocount_response, na.rm = TRUE)
  )
)

raw_hist <- ggplot(df, aes(x = Y)) +
  geom_histogram(bins = 80, fill = "#4C78A8", color = "white", linewidth = 0.1) +
  labs(
    title = "Global Raw MAG Abundance Distribution",
    x = "Raw Y",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)

log_hist <- ggplot(df, aes(x = log10_pseudocount_response)) +
  geom_histogram(bins = 80, fill = "#F58518", color = "white", linewidth = 0.1) +
  labs(
    title = "Global log10(Y + pseudocount) Distribution",
    subtitle = paste("pseudocount =", format(pseudocount, scientific = TRUE)),
    x = "log10(Y + pseudocount)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)

raw_density <- ggplot(df, aes(x = Y)) +
  geom_density(fill = "#4C78A8", alpha = 0.35, linewidth = 0.7, na.rm = TRUE) +
  labs(
    title = "Global Raw MAG Abundance Density",
    x = "Raw Y",
    y = "Density"
  ) +
  theme_minimal(base_size = 12)

log_density <- ggplot(df, aes(x = log10_pseudocount_response)) +
  geom_density(fill = "#F58518", alpha = 0.35, linewidth = 0.7, na.rm = TRUE) +
  labs(
    title = "Global log10(Y + pseudocount) Density",
    subtitle = paste("pseudocount =", format(pseudocount, scientific = TRUE)),
    x = "log10(Y + pseudocount)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12)

save_plot(raw_hist, "global_raw_Y_histogram.png")
save_plot(log_hist, "global_log10_pseudocount_histogram.png")
save_plot(raw_density, "global_raw_Y_density.png")
save_plot(log_density, "global_log10_pseudocount_density.png")

rank_response <- function(response_type) {
  response_rows <- comparison[comparison$response_type == response_type, , drop = FALSE]
  response_rows[order(
    is.na(response_rows$time_p_adj_bh),
    response_rows$time_p_adj_bh,
    response_rows$time_p_value
  ), , drop = FALSE]
}

raw_ranked <- rank_response("raw")
log_ranked <- rank_response("log10_pseudocount")
top_raw <- utils::head(raw_ranked$MAG, top_n)
top_log <- utils::head(log_ranked$MAG, top_n)
selected_mags <- sort(unique(c(top_raw, top_log)))

comparison_subset <- comparison[comparison$MAG %in% selected_mags, , drop = FALSE]
comparison_wide <- reshape(
  comparison_subset,
  idvar = "MAG",
  timevar = "response_type",
  direction = "wide"
)

get_col <- function(data, name) {
  if (name %in% names(data)) data[[name]] else rep(NA, nrow(data))
}

diagnostics <- data.frame(
  MAG = comparison_wide$MAG,
  top_ranked_raw = comparison_wide$MAG %in% top_raw,
  top_ranked_log10_pseudocount = comparison_wide$MAG %in% top_log,
  raw_time_p_value = get_col(comparison_wide, "time_p_value.raw"),
  raw_time_p_adj_bh = get_col(comparison_wide, "time_p_adj_bh.raw"),
  log_time_p_value = get_col(comparison_wide, "time_p_value.log10_pseudocount"),
  log_time_p_adj_bh = get_col(comparison_wide, "time_p_adj_bh.log10_pseudocount"),
  raw_time_estimate = get_col(comparison_wide, "time_estimate.raw"),
  log_time_estimate = get_col(comparison_wide, "time_estimate.log10_pseudocount"),
  raw_status = get_col(comparison_wide, "status.raw"),
  log_status = get_col(comparison_wide, "status.log10_pseudocount"),
  zero_fraction = get_col(comparison_wide, "zero_fraction.raw"),
  nonzero_sample_count = get_col(comparison_wide, "nonzero_sample_count.raw"),
  subjects_with_nonzero = get_col(comparison_wide, "subjects_with_nonzero.raw"),
  raw_y_variance = get_col(comparison_wide, "y_variance.raw"),
  log_y_variance = get_col(comparison_wide, "y_variance.log10_pseudocount"),
  notes = "",
  raw_fit_succeeded = NA,
  raw_singular_refit = NA,
  raw_residual_variance = NA_real_,
  raw_residual_mean = NA_real_,
  raw_residual_sd = NA_real_,
  raw_residual_extreme = NA,
  raw_refit_warnings = "",
  raw_refit_error = "",
  log_fit_succeeded = NA,
  log_singular_refit = NA,
  log_residual_variance = NA_real_,
  log_residual_mean = NA_real_,
  log_residual_sd = NA_real_,
  log_residual_extreme = NA,
  log_refit_warnings = "",
  log_refit_error = "",
  stringsAsFactors = FALSE
)

note_for_row <- function(raw_status, log_status) {
  notes <- character(0)
  if (raw_status == "lmm_singular" || log_status == "lmm_singular") {
    notes <- c(notes, "singular fit in at least one response scale")
  }
  if (raw_status == "too_sparse" || log_status == "too_sparse") {
    notes <- c(notes, "too sparse in at least one response scale")
  }
  if (raw_status == "low_variance" || log_status == "low_variance") {
    notes <- c(notes, "low variance in at least one response scale")
  }
  if (length(notes) == 0) "none" else paste(notes, collapse = "; ")
}

plot_trajectory <- function(mag_df, mag_id, response_type) {
  if (response_type == "raw") {
    response_col <- "Y"
    y_label <- "Raw Y"
    title <- paste(mag_id, "- raw abundance trajectory")
    suffix <- "trajectory_raw.png"
  } else {
    response_col <- "response_value"
    y_label <- "log10(Y + pseudocount)"
    title <- paste(mag_id, "- log10-pseudocount trajectory")
    suffix <- "trajectory_log10_pseudocount.png"
  }

  ggplot(mag_df, aes(x = days_from_baseline_subject, y = .data[[response_col]], group = Subject)) +
    geom_line(alpha = 0.35, linewidth = 0.45, color = "#6B7280") +
    geom_point(aes(color = Subject), size = 1.8, alpha = 0.85, show.legend = FALSE) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "#1F2937", linewidth = 0.8) +
    labs(
      title = title,
      subtitle = paste("Response type:", response_type),
      x = "Days from subject baseline",
      y = y_label
    ) +
    theme_minimal(base_size = 11)
}

fit_for_residuals <- function(mag_df) {
  warning_messages <- character(0)
  error_messages <- character(0)
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
  list(
    fit = fit,
    warnings = paste(unique(warning_messages), collapse = " | "),
    error = paste(unique(error_messages), collapse = " | ")
  )
}

plot_residuals <- function(residual_df, mag_id, response_type) {
  ggplot(residual_df, aes(x = fitted, y = residual)) +
    geom_hline(yintercept = 0, color = "#6B7280", linewidth = 0.4) +
    geom_point(alpha = 0.75, size = 1.8, color = "#4C78A8") +
    labs(
      title = paste(mag_id, "-", response_type, "residuals vs fitted"),
      x = "Fitted value",
      y = "Residual"
    ) +
    theme_minimal(base_size = 11)
}

plot_qq <- function(residual_df, mag_id, response_type) {
  ggplot(residual_df, aes(sample = residual)) +
    stat_qq(color = "#4C78A8", alpha = 0.75, size = 1.8) +
    stat_qq_line(color = "#1F2937", linewidth = 0.7) +
    labs(
      title = paste(mag_id, "-", response_type, "residual QQ plot"),
      x = "Theoretical quantiles",
      y = "Sample quantiles"
    ) +
    theme_minimal(base_size = 11)
}

for (mag_id in selected_mags) {
  safe_mag <- safe_filename(mag_id)
  mag_df_base <- df[df$MAG == mag_id, , drop = FALSE]
  mag_df_base <- mag_df_base[stats::complete.cases(mag_df_base[, c("Y", "days_from_baseline_subject", "Subject")]), , drop = FALSE]
  mag_df_base$Subject <- droplevels(factor(mag_df_base$Subject))

  for (response_type in c("raw", "log10_pseudocount")) {
    mag_df <- mag_df_base
    if (response_type == "raw") {
      mag_df$response_value <- mag_df$Y
      trajectory_plot <- plot_trajectory(mag_df, mag_id, "raw")
      trajectory_file <- paste0(safe_mag, "_trajectory_raw.png")
      residual_file <- paste0(safe_mag, "_raw_residuals_vs_fitted.png")
      qq_file <- paste0(safe_mag, "_raw_qqplot.png")
      prefix <- "raw"
    } else {
      mag_df$response_value <- log10(mag_df$Y + pseudocount)
      trajectory_plot <- plot_trajectory(mag_df, mag_id, "log10_pseudocount")
      trajectory_file <- paste0(safe_mag, "_trajectory_log10_pseudocount.png")
      residual_file <- paste0(safe_mag, "_log10_pseudocount_residuals_vs_fitted.png")
      qq_file <- paste0(safe_mag, "_log10_pseudocount_qqplot.png")
      prefix <- "log"
    }

    save_plot(trajectory_plot, trajectory_file)

    fit_result <- fit_for_residuals(mag_df)
    row_index <- match(mag_id, diagnostics$MAG)
    diagnostics[row_index, paste0(prefix, "_fit_succeeded")] <- !is.null(fit_result$fit)
    diagnostics[row_index, paste0(prefix, "_refit_warnings")] <- fit_result$warnings
    diagnostics[row_index, paste0(prefix, "_refit_error")] <- fit_result$error

    if (!is.null(fit_result$fit)) {
      residuals <- stats::residuals(fit_result$fit)
      fitted <- stats::fitted(fit_result$fit)
      residual_df <- data.frame(fitted = fitted, residual = residuals)
      residual_sd <- stats::sd(residuals, na.rm = TRUE)
      standardized_residual <- if (is.finite(residual_sd) && residual_sd > 0) residuals / residual_sd else rep(NA_real_, length(residuals))
      extreme_residual <- any(abs(standardized_residual) > 3, na.rm = TRUE)

      diagnostics[row_index, paste0(prefix, "_singular_refit")] <- lme4::isSingular(fit_result$fit)
      diagnostics[row_index, paste0(prefix, "_residual_variance")] <- stats::var(residuals, na.rm = TRUE)
      diagnostics[row_index, paste0(prefix, "_residual_mean")] <- mean(residuals, na.rm = TRUE)
      diagnostics[row_index, paste0(prefix, "_residual_sd")] <- residual_sd
      diagnostics[row_index, paste0(prefix, "_residual_extreme")] <- extreme_residual

      save_plot(plot_residuals(residual_df, mag_id, response_type), residual_file)
      save_plot(plot_qq(residual_df, mag_id, response_type), qq_file)
    }
  }
}

diagnostics$notes <- mapply(note_for_row, diagnostics$raw_status, diagnostics$log_status)
diagnostics <- diagnostics[order(!diagnostics$top_ranked_raw, !diagnostics$top_ranked_log10_pseudocount, diagnostics$MAG), ]
utils::write.csv(diagnostics, diagnostics_csv, row.names = FALSE)

summary_connection <- file(summary_path, open = "wt")
sink(summary_connection, split = TRUE)

cat("Response-Scale Diagnostics Summary\n")
cat("==================================\n\n")
cat("Processed file:", processed_path, "\n")
cat("Transformation comparison file:", comparison_path, "\n")
cat("Diagnostics CSV:", diagnostics_csv, "\n")
cat("Figure directory:", figure_dir, "\n")
cat("P-value method for residual refits:", p_value_method, "\n\n")

cat("Global raw Y distribution:\n")
print(global_summary[global_summary$metric %in% c(
  "min_raw_Y", "max_raw_Y", "mean_raw_Y", "median_raw_Y",
  "zero_fraction", "smallest_nonzero_Y", "pseudocount"
), ], row.names = FALSE)
cat("\n")

cat("Global log10-pseudocount distribution:\n")
print(global_summary[global_summary$metric %in% c(
  "min_log10_response", "max_log10_response",
  "mean_log10_response", "median_log10_response"
), ], row.names = FALSE)
cat("\n")

cat("Global plots created:\n")
cat(file.path(figure_dir, "global_raw_Y_histogram.png"), "\n")
cat(file.path(figure_dir, "global_log10_pseudocount_histogram.png"), "\n")
cat(file.path(figure_dir, "global_raw_Y_density.png"), "\n")
cat(file.path(figure_dir, "global_log10_pseudocount_density.png"), "\n\n")

cat("Selected MAGs for detailed diagnostics:", length(selected_mags), "\n\n")

cat("Top raw-response MAGs:\n")
cat(paste(top_raw, collapse = "\n"), "\n\n", sep = "")

cat("Top log10-pseudocount-response MAGs:\n")
cat(paste(top_log, collapse = "\n"), "\n\n", sep = "")

top_overlap <- intersect(top_raw, top_log)
cat("Overlap between top raw and top log-transformed MAGs:", length(top_overlap), "\n")
if (length(top_overlap) > 0) {
  cat(paste(top_overlap, collapse = "\n"), "\n\n", sep = "")
} else {
  cat("None\n\n")
}

cat("Residual refit status:\n")
cat("Raw refit failures:", sum(!diagnostics$raw_fit_succeeded, na.rm = TRUE), "\n")
cat("Log refit failures:", sum(!diagnostics$log_fit_succeeded, na.rm = TRUE), "\n")
cat("Raw singular refits:", sum(diagnostics$raw_singular_refit %in% TRUE, na.rm = TRUE), "\n")
cat("Log singular refits:", sum(diagnostics$log_singular_refit %in% TRUE, na.rm = TRUE), "\n")
cat("Raw models with extreme residual checks:", sum(diagnostics$raw_residual_extreme %in% TRUE, na.rm = TRUE), "\n")
cat("Log models with extreme residual checks:", sum(diagnostics$log_residual_extreme %in% TRUE, na.rm = TRUE), "\n\n")

cat("Interpretation:\n")
if (length(top_overlap) == 0) {
  cat("The top-ranked MAG lists are strongly transformation-sensitive: the top raw-response and log10-pseudocount sets have no overlap.\n")
} else if (length(top_overlap) < min(length(top_raw), length(top_log)) / 2) {
  cat("The top-ranked MAG lists are partially transformation-sensitive: fewer than half of the top MAGs overlap between response scales.\n")
} else {
  cat("The top-ranked MAG lists show moderate stability across response scales, but residual diagnostics and biological interpretation are still needed.\n")
}
cat("Log transformation reduced singular behavior in the all-MAG comparison, but it also changed MAG ranking substantially.\n")
cat("No biological conclusion should be made yet while rankings remain transformation-sensitive. Review the trajectory and residual plots before deciding which response scale is more appropriate for downstream modeling or HPC scaling.\n\n")

cat("Saved diagnostics CSV to results/tables/response_scale_top_mag_diagnostics.csv\n")
cat("Saved summary to results/logs/response_scale_diagnostics_summary.txt\n")

sink()
close(summary_connection)
