# Public renderer: values come from external inputs; see docs/public-rendering.md.
# ============================================
# plot_figure_7_hpc_per_mag_lmertest_results.R
# ============================================
# Figure 7: model status and data adequacy dashboard for the full T1-T10
# per-MAG lmerTest LMM workflow.
#
# This figure documents modelability checks across all MAGs, including status,
# sparsity, response variance, non-zero support, singularity, and failures.

required_packages <- c("ggplot2", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
results_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "tables",
  "per_mag_lmertest_hpc_results.csv"
)
output_dir <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "figures")
output_stem <- "figure_7_hpc_per_mag_lmertest_results"

# External inputs and isolated output destination (run from repository root).
results_path <- Sys.getenv("PUBLIC_LMM_RESULTS", unset = results_path)
output_dir <- Sys.getenv("PUBLIC_RENDER_OUTPUT_DIR", unset = output_dir)

if (!file.exists(results_path)) {
  stop("Missing per-MAG LMM result table: ", results_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c(
  "MAG", "status", "fit_succeeded", "singular", "zero_fraction",
  "response_variance", "n_rows", "n_subjects", "n_timepoints",
  "nonzero_sample_count", "subjects_with_nonzero", "warnings", "errors",
  "primary_analysis"
)
missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns) > 0) {
  stop("Missing required result columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

status_levels <- c("lmm_ok", "lmm_singular", "too_sparse", "low_variance", "model_failed")
status_labels <- c(
  lmm_ok = "lmm_ok",
  lmm_singular = "lmm_singular",
  too_sparse = "too_sparse",
  low_variance = "low_variance",
  model_failed = "model_failed"
)

plot_df <- results |>
  dplyr::mutate(
    status = factor(status, levels = status_levels),
    status_label = factor(status_labels[as.character(status)], levels = status_labels),
    fit_succeeded = fit_succeeded %in% TRUE,
    singular = singular %in% TRUE,
    warning_present = !is.na(warnings) & warnings != "",
    error_present = !is.na(errors) & errors != ""
  )

status_summary <- plot_df |>
  dplyr::count(status, status_label, name = "n_mags", .drop = FALSE) |>
  dplyr::mutate(
    fraction = n_mags / sum(n_mags),
    label = paste0(n_mags, "\n", sprintf("%.1f%%", 100 * fraction))
  )

adequacy_summary <- plot_df |>
  dplyr::group_by(status) |>
  dplyr::summarise(
    n_mags = dplyr::n(),
    fit_succeeded = sum(fit_succeeded, na.rm = TRUE),
    singular = sum(singular, na.rm = TRUE),
    warning_present = sum(warning_present, na.rm = TRUE),
    error_present = sum(error_present, na.rm = TRUE),
    median_zero_fraction = stats::median(zero_fraction, na.rm = TRUE),
    median_response_variance = stats::median(response_variance, na.rm = TRUE),
    median_nonzero_samples = stats::median(nonzero_sample_count, na.rm = TRUE),
    median_subjects_with_nonzero = stats::median(subjects_with_nonzero, na.rm = TRUE),
    .groups = "drop"
  )

support_long <- rbind(
  data.frame(
    MAG = plot_df$MAG,
    status = plot_df$status,
    metric = "Non-zero samples",
    value = plot_df$nonzero_sample_count,
    threshold = 10
  ),
  data.frame(
    MAG = plot_df$MAG,
    status = plot_df$status,
    metric = "Subjects with non-zero",
    value = plot_df$subjects_with_nonzero,
    threshold = 3
  )
)
support_long$metric <- factor(support_long$metric, levels = c("Non-zero samples", "Subjects with non-zero"))

issue_table <- plot_df |>
  dplyr::filter(status != "lmm_ok" | warning_present | error_present) |>
  dplyr::mutate(
    issue = dplyr::case_when(
      status == "too_sparse" ~ "Sparse non-zero support",
      status == "lmm_singular" ~ "Singular random effect",
      status == "low_variance" ~ "Low response variance",
      status == "model_failed" ~ "Model failure",
      warning_present ~ "Warning recorded",
      error_present ~ "Error recorded",
      TRUE ~ "Flagged"
    ),
    mag_short = ifelse(nchar(MAG) > 34, paste0(substr(MAG, 1, 31), "..."), MAG)
  ) |>
  dplyr::select(MAG, mag_short, status, issue, zero_fraction, response_variance, nonzero_sample_count, subjects_with_nonzero)

total_mags <- nrow(plot_df)
primary_mags <- sum(plot_df$primary_analysis %in% TRUE, na.rm = TRUE)
fit_success <- sum(plot_df$fit_succeeded, na.rm = TRUE)
singular_mags <- sum(plot_df$singular, na.rm = TRUE)
warning_count <- sum(plot_df$warning_present, na.rm = TRUE)
error_count <- sum(plot_df$error_present, na.rm = TRUE)

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)

pal <- list(
  ink = "#111827",
  text = "#374151",
  muted = "#6B7280",
  hair = "#D1D5DB",
  faint = "#F8FAFC",
  blue = "#0072B2",
  sky = "#DCEEF8",
  orange = "#D55E00",
  peach = "#FFF2E8",
  red = "#CC3311",
  rose = "#FFF0EE",
  green = "#009E73",
  purple = "#785EF0",
  grey = "#9CA3AF"
)

status_palette <- c(
  lmm_ok = pal$blue,
  lmm_singular = pal$orange,
  too_sparse = pal$red,
  low_variance = pal$purple,
  model_failed = pal$grey
)

base_theme <- ggplot2::theme_minimal(base_size = 9.5, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11.0, color = pal$ink),
    plot.subtitle = ggplot2::element_text(size = 8.0, color = pal$muted, margin = ggplot2::margin(b = 5)),
    axis.title = ggplot2::element_text(size = 8.4, color = pal$text),
    axis.text = ggplot2::element_text(size = 7.5, color = pal$text),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.22),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 8),
    legend.text = ggplot2::element_text(size = 7.5),
    plot.margin = ggplot2::margin(5, 7, 5, 7)
  )

panel_a <- ggplot2::ggplot(status_summary, ggplot2::aes(x = status_label, y = n_mags, fill = status)) +
  ggplot2::geom_col(width = 0.68, color = "white", linewidth = 0.35) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(n_mags > 0, label, "")),
    vjust = -0.15,
    size = 2.7,
    color = pal$text,
    lineheight = 0.9
  ) +
  ggplot2::scale_fill_manual(values = status_palette, drop = FALSE, guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10))) +
  ggplot2::labs(
    title = "A  Model status across all MAGs",
    subtitle = paste0(primary_mags, " of ", total_mags, " MAGs retained as primary lmm_ok analyses"),
    x = NULL,
    y = "MAG count"
  ) +
  base_theme +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    panel.grid.major.x = ggplot2::element_blank()
  )

panel_b <- ggplot2::ggplot(plot_df, ggplot2::aes(x = status, y = zero_fraction, fill = status)) +
  ggplot2::geom_boxplot(width = 0.56, outlier.shape = NA, alpha = 0.82, color = pal$text, linewidth = 0.35) +
  ggplot2::geom_jitter(width = 0.10, size = 1.35, alpha = 0.50, color = pal$text) +
  ggplot2::scale_fill_manual(values = status_palette, drop = FALSE, guide = "none") +
  ggplot2::scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    expand = ggplot2::expansion(mult = c(0.03, 0.08))
  ) +
  ggplot2::labs(
    title = "B  Sparsity by model status",
    subtitle = "Zero fraction summarizes how often a MAG is absent across samples",
    x = NULL,
    y = "Zero fraction"
  ) +
  base_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

panel_c <- ggplot2::ggplot(plot_df, ggplot2::aes(x = status, y = response_variance, fill = status)) +
  ggplot2::geom_boxplot(width = 0.56, outlier.shape = NA, alpha = 0.82, color = pal$text, linewidth = 0.35) +
  ggplot2::geom_jitter(width = 0.10, size = 1.35, alpha = 0.50, color = pal$text) +
  ggplot2::scale_fill_manual(values = status_palette, drop = FALSE, guide = "none") +
  ggplot2::labs(
    title = "C  Response variance by model status",
    subtitle = "Variance is computed on log10(Y + pseudocount), the modeled response scale",
    x = NULL,
    y = "Response variance"
  ) +
  base_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

panel_d <- ggplot2::ggplot(support_long, ggplot2::aes(x = value, fill = status)) +
  ggplot2::geom_histogram(binwidth = 3, boundary = 0, color = "white", linewidth = 0.22, alpha = 0.88) +
  ggplot2::geom_vline(
    data = unique(support_long[, c("metric", "threshold")]),
    ggplot2::aes(xintercept = threshold),
    linetype = "longdash",
    linewidth = 0.52,
    color = pal$red,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_wrap(~metric, scales = "free_x", nrow = 1) +
  ggplot2::scale_fill_manual(values = status_palette, drop = FALSE, name = "Model status") +
  ggplot2::labs(
    title = "D  Non-zero support screens used before primary inference",
    subtitle = "Dashed lines show configured minimum support: >=10 non-zero samples and >=3 non-zero subjects",
    x = "Support count per MAG",
    y = "MAG count"
  ) +
  base_theme +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(face = "bold", color = pal$ink),
    strip.background = ggplot2::element_rect(fill = pal$faint, color = pal$hair, linewidth = 0.35)
  )

draw_issue_table <- function() {
  grid::grid.roundrect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.50, "npc"),
    width = grid::unit(0.98, "npc"),
    height = grid::unit(0.88, "npc"),
    r = grid::unit(0.035, "snpc"),
    gp = grid::gpar(fill = "white", col = pal$hair, lwd = 0.8)
  )
  grid::grid.text(
    "E  Flagged MAG audit trail",
    x = grid::unit(0.04, "npc"),
    y = grid::unit(0.88, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 10.1, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    paste0("Warnings recorded: ", warning_count, "; errors recorded: ", error_count, ". Flagged rows are retained in the result table."),
    x = grid::unit(0.04, "npc"),
    y = grid::unit(0.78, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 7.6, col = pal$muted)
  )

  if (nrow(issue_table) == 0) {
    table_text <- "No flagged MAGs, warnings, or errors recorded."
  } else {
    table_text <- paste(
      sprintf(
        "%-34s  %-20s  zero=%4.1f%%  var=%5.2f  nonzero=%2d  subjects=%2d",
        issue_table$mag_short,
        issue_table$issue,
        100 * issue_table$zero_fraction,
        issue_table$response_variance,
        issue_table$nonzero_sample_count,
        issue_table$subjects_with_nonzero
      ),
      collapse = "\n"
    )
  }

  grid::grid.text(
    table_text,
    x = grid::unit(0.04, "npc"),
    y = grid::unit(0.64, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 6.7, col = pal$text, fontfamily = "mono", lineheight = 0.95)
  )
}

summary_text <- paste0(
  "QC summary: ", fmt_int(fit_success), " of ", fmt_int(total_mags),
  " fits succeeded; ", fmt_int(primary_mags), " primary lmm_ok analyses; ",
  fmt_int(singular_mags), " singular LMM; ", fmt_int(warning_count),
  " warnings; ", fmt_int(error_count),
  paste0(" errors. Subject counts range ", min(results$n_subjects, na.rm = TRUE), "-", max(results$n_subjects, na.rm = TRUE), "; timepoint counts range ", min(results$n_timepoints, na.rm = TRUE), "-", max(results$n_timepoints, na.rm = TRUE), ". Non-zero support and model status are additional adequacy checks.")
)

draw_figure <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  layout <- grid::grid.layout(
    nrow = 4,
    ncol = 2,
    heights = grid::unit(c(0.50, 1.05, 1.0, 0.48), c("in", "null", "null", "in")),
    widths = grid::unit(c(1.0, 1.02), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))

  grid::grid.text(
    "Figure 7 | Model status and data adequacy dashboard",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.88, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 15.5, fontface = "bold", col = pal$ink),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )
  grid::grid.text(
    "Per-MAG models were audited for success, sparsity, response variation, singularity, and failure status before interpretation.",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.20, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 9.0, col = pal$text),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )

  print(panel_a, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panel_b, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(panel_c, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
  grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
  split_layout <- grid::grid.layout(nrow = 2, ncol = 1, heights = grid::unit(c(0.62, 0.38), "null"))
  grid::pushViewport(grid::viewport(layout = split_layout))
  print(panel_d, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  draw_issue_table()
  grid::popViewport(3)

  grid::grid.text(
    paste(strwrap(summary_text, width = 170), collapse = "\n"),
    x = grid::unit(0, "npc"),
    y = grid::unit(0.58, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 6.9, col = pal$muted, lineheight = 0.95),
    vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1:2)
  )

  grid::popViewport()
}

png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

png(filename = png_path, width = 13.0, height = 9.8, units = "in", res = 420, bg = "white")
draw_figure()
dev.off()

pdf(file = pdf_path, width = 13.0, height = 9.8, bg = "white", useDingbats = FALSE)
draw_figure()
dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = svg_path, width = 13.0, height = 9.8, bg = "white")
  draw_figure()
  dev.off()
}

status_summary_path <- file.path(output_dir, paste0(output_stem, "_model_status_summary.csv"))
flagged_path <- file.path(output_dir, paste0(output_stem, "_flagged_MAG_audit_table.csv"))

utils::write.csv(adequacy_summary, status_summary_path, row.names = FALSE)
utils::write.csv(issue_table, flagged_path, row.names = FALSE)

message("Figure 7 model-QC outputs created:")
message(" - ", png_path)
message(" - ", pdf_path)
if (file.exists(svg_path)) message(" - ", svg_path)
message(" - ", status_summary_path)
message(" - ", flagged_path)
