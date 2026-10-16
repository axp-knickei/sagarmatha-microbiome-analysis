# ============================================
# prepare_top_mag_manual_review.R
# ============================================
# Build a manual-review package for the selected top MAG trajectory and
# residual plots. This script only reads existing diagnostic outputs and writes
# review helper files; it does not rerun modeling or modify existing figures.

diagnostics_path <- "results/tables/response_scale_top_mag_diagnostics.csv"
figure_dir <- "results/figures/response_scale"
template_path <- "results/tables/top_mag_manual_review_template.csv"
guide_path <- "results/logs/top_mag_manual_review_guide.md"
index_path <- "results/logs/top_mag_plot_index.md"
summary_path <- "results/logs/top_mag_manual_review_summary.txt"

dir.create(dirname(template_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(guide_path), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(diagnostics_path)) {
  stop("Missing diagnostics table: ", diagnostics_path, call. = FALSE)
}
if (!dir.exists(figure_dir)) {
  stop("Missing figure directory: ", figure_dir, call. = FALSE)
}

diagnostics <- read.csv(diagnostics_path, stringsAsFactors = FALSE)

required_cols <- c(
  "MAG",
  "top_ranked_raw",
  "top_ranked_log10_pseudocount",
  "raw_time_estimate",
  "raw_time_p_value",
  "raw_time_p_adj_bh",
  "raw_status",
  "log_time_estimate",
  "log_time_p_value",
  "log_time_p_adj_bh",
  "log_status",
  "zero_fraction",
  "nonzero_sample_count",
  "subjects_with_nonzero"
)
missing_cols <- setdiff(required_cols, names(diagnostics))
if (length(missing_cols) > 0) {
  stop("Diagnostics table is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

safe_filename <- function(x) {
  safe <- gsub("[^A-Za-z0-9._-]+", "_", x)
  safe <- gsub("_+", "_", safe)
  safe <- gsub("^_|_$", "", safe)
  ifelse(nchar(safe) == 0, "unnamed_MAG", safe)
}

plot_path <- function(mag, suffix) {
  file.path(figure_dir, paste0(safe_filename(mag), suffix))
}

template <- data.frame(
  MAG = diagnostics$MAG,
  top_raw = diagnostics$top_ranked_raw,
  top_log10_pseudocount = diagnostics$top_ranked_log10_pseudocount,
  raw_time_estimate = diagnostics$raw_time_estimate,
  raw_time_p_value = diagnostics$raw_time_p_value,
  raw_time_p_adj_bh = diagnostics$raw_time_p_adj_bh,
  raw_status = diagnostics$raw_status,
  log_time_estimate = diagnostics$log_time_estimate,
  log_time_p_value = diagnostics$log_time_p_value,
  log_time_p_adj_bh = diagnostics$log_time_p_adj_bh,
  log_status = diagnostics$log_status,
  zero_fraction = diagnostics$zero_fraction,
  nonzero_sample_count = diagnostics$nonzero_sample_count,
  subjects_with_nonzero = diagnostics$subjects_with_nonzero,
  raw_plot_trajectory = plot_path(diagnostics$MAG, "_trajectory_raw.png"),
  log_plot_trajectory = plot_path(diagnostics$MAG, "_trajectory_log10_pseudocount.png"),
  raw_plot_residuals_vs_fitted = plot_path(diagnostics$MAG, "_raw_residuals_vs_fitted.png"),
  log_plot_residuals_vs_fitted = plot_path(diagnostics$MAG, "_log10_pseudocount_residuals_vs_fitted.png"),
  raw_plot_qq = plot_path(diagnostics$MAG, "_raw_qqplot.png"),
  log_plot_qq = plot_path(diagnostics$MAG, "_log10_pseudocount_qqplot.png"),
  manual_visual_trend = "",
  manual_outlier_driven = "",
  manual_many_subjects_support = "",
  manual_zero_dominated = "",
  manual_residual_quality_raw = "",
  manual_residual_quality_log = "",
  manual_recommended_action = "",
  manual_notes = "",
  stringsAsFactors = FALSE
)

plot_cols <- c(
  "raw_plot_trajectory",
  "log_plot_trajectory",
  "raw_plot_residuals_vs_fitted",
  "log_plot_residuals_vs_fitted",
  "raw_plot_qq",
  "log_plot_qq"
)

missing_plot_rows <- do.call(
  rbind,
  lapply(plot_cols, function(col) {
    missing <- !file.exists(template[[col]])
    if (!any(missing)) {
      return(NULL)
    }
    data.frame(
      MAG = template$MAG[missing],
      plot_type = col,
      path = template[[col]][missing],
      stringsAsFactors = FALSE
    )
  })
)

write.csv(template, template_path, row.names = FALSE)

guide_lines <- c(
  "# Top MAG Manual Review Guide",
  "",
  "Use this guide with `results/tables/top_mag_manual_review_template.csv` and `results/logs/top_mag_plot_index.md`.",
  "",
  "## Manual Review Fields",
  "",
  "`manual_visual_trend` allowed values:",
  "",
  "- `clear_increase`",
  "- `clear_decrease`",
  "- `weak_trend`",
  "- `no_clear_trend`",
  "- `unclear`",
  "",
  "`manual_outlier_driven`, `manual_many_subjects_support`, and `manual_zero_dominated` allowed values:",
  "",
  "- `yes`",
  "- `no`",
  "- `unclear`",
  "",
  "`manual_residual_quality_raw` and `manual_residual_quality_log` allowed values:",
  "",
  "- `acceptable`",
  "- `questionable`",
  "- `poor`",
  "- `not_reviewed`",
  "",
  "`manual_recommended_action` allowed values:",
  "",
  "- `keep_primary`",
  "- `keep_sensitivity_only`",
  "- `flag_outlier_driven`",
  "- `flag_sparse`",
  "- `exclude_from_primary`",
  "- `needs_more_review`",
  "",
  "## Trajectory Plot Review",
  "",
  "- Check whether the trend is visible across many subjects.",
  "- Check whether one subject or one timepoint drives the apparent trend.",
  "- Check whether zeros dominate the pattern.",
  "- Check whether raw and log plots tell the same or different stories.",
  "- Check whether the direction of effect matches the model estimate.",
  "",
  "## Residuals-vs-Fitted Review",
  "",
  "- Prefer residuals scattered roughly around zero.",
  "- Flag funnel shapes, extreme residuals, or strong fitted-value dependence.",
  "- Compare whether log response improves residual behavior compared with raw response.",
  "",
  "## QQ Plot Review",
  "",
  "- Prefer points roughly following the diagonal.",
  "- Flag strong tail deviation.",
  "- Do not expect perfection, but strong deviation should be noted.",
  "",
  "## Final Action Guidance",
  "",
  "- Use `keep_primary` only if the log response plot looks reasonable, is not outlier-driven, is supported by multiple subjects, and has acceptable residuals.",
  "- Use `keep_sensitivity_only` if the raw or log result is interesting but unstable.",
  "- Use `flag_outlier_driven` if one subject/timepoint drives the signal.",
  "- Use `flag_sparse` if too many zeros dominate the signal.",
  "- Use `exclude_from_primary` if the MAG is too unstable for interpretation.",
  "- Use `needs_more_review` if uncertain.",
  "",
  "Manual review is required before creating the HPC-ready final modeling script."
)
writeLines(guide_lines, guide_path)

index_lines <- c(
  "# Top MAG Plot Index",
  "",
  "Paths are relative to the project root.",
  ""
)

for (i in seq_len(nrow(template))) {
  source <- if (template$top_raw[[i]] && template$top_log10_pseudocount[[i]]) {
    "top raw and top log10-pseudocount"
  } else if (template$top_raw[[i]]) {
    "top raw"
  } else if (template$top_log10_pseudocount[[i]]) {
    "top log10-pseudocount"
  } else {
    "selected diagnostic MAG"
  }

  index_lines <- c(
    index_lines,
    paste0("## ", template$MAG[[i]]),
    "",
    paste0("- Source: ", source),
    paste0("- Raw trajectory: `", template$raw_plot_trajectory[[i]], "`"),
    paste0("- Log10-pseudocount trajectory: `", template$log_plot_trajectory[[i]], "`"),
    paste0("- Raw residuals vs fitted: `", template$raw_plot_residuals_vs_fitted[[i]], "`"),
    paste0("- Log10-pseudocount residuals vs fitted: `", template$log_plot_residuals_vs_fitted[[i]], "`"),
    paste0("- Raw QQ plot: `", template$raw_plot_qq[[i]], "`"),
    paste0("- Log10-pseudocount QQ plot: `", template$log_plot_qq[[i]], "`"),
    ""
  )
}
writeLines(index_lines, index_path)

top_raw_count <- sum(template$top_raw %in% TRUE, na.rm = TRUE)
top_log_count <- sum(template$top_log10_pseudocount %in% TRUE, na.rm = TRUE)
overlap_count <- sum(template$top_raw %in% TRUE & template$top_log10_pseudocount %in% TRUE, na.rm = TRUE)
missing_plot_count <- if (is.null(missing_plot_rows)) 0 else nrow(missing_plot_rows)

summary_lines <- c(
  "Top MAG Manual Review Package Summary",
  "=====================================",
  "",
  paste("Manual review template:", template_path),
  paste("Manual review guide:", guide_path),
  paste("Plot index:", index_path),
  "",
  paste("MAGs prepared for manual review:", nrow(template)),
  paste("MAGs from top raw list:", top_raw_count),
  paste("MAGs from top log10-pseudocount list:", top_log_count),
  paste("Overlap count:", overlap_count),
  paste("Missing plot files:", missing_plot_count),
  "",
  "Manual review columns are intentionally blank until the plots are inspected by eye."
)

if (missing_plot_count > 0) {
  summary_lines <- c(
    summary_lines,
    "",
    "Missing plot files:",
    paste(missing_plot_rows$MAG, missing_plot_rows$plot_type, missing_plot_rows$path, sep = " | ")
  )
}

writeLines(summary_lines, summary_path)

cat("Created manual review template:", template_path, "\n")
cat("Created manual review guide:", guide_path, "\n")
cat("Created plot index:", index_path, "\n")
cat("Created summary:", summary_path, "\n")
cat("MAGs prepared:", nrow(template), "\n")
cat("Missing plot files:", missing_plot_count, "\n")
