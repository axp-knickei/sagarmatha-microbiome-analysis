# ============================================
# plot_figure_2_class1_analysis.R
# ============================================
# Main Figure 2: Class 1 analysis workflow diagram.
# This script reads the Class 1 LMM outputs and creates a publication-style
# schematic explaining the reproducible analysis process.

required_packages <- c("readr", "dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
results_path <- file.path(
  repo_root,
  "results",
  "class1_window_lmertest",
  "tables",
  "class1_window_lmertest_results.csv"
)
windows_path <- file.path(
  repo_root,
  "results",
  "metadata_expedition_windows",
  "proposed_biological_windows.csv"
)
output_dir <- file.path(repo_root, "results", "class1_window_lmertest", "figures")

if (!file.exists(results_path)) {
  stop("Missing Class 1 results table: ", results_path, call. = FALSE)
}
if (!file.exists(windows_path)) {
  stop("Missing Class 1 window definition table: ", windows_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

class1_results <- readr::read_csv(results_path, show_col_types = FALSE)
class1_windows <- readr::read_csv(windows_path, show_col_types = FALSE)

fmt_int <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
fmt_q <- function(x) {
  if (is.na(x)) {
    return("NA")
  }
  sprintf("%.4f", x)
}

n_mags <- dplyr::n_distinct(class1_results$MAG)
n_windows <- nrow(class1_windows)
n_tests <- nrow(class1_results)
n_primary <- sum(class1_results$primary_analysis %in% TRUE, na.rm = TRUE)
n_singular <- sum(class1_results$status == "lmm_singular", na.rm = TRUE)
n_sparse <- sum(class1_results$status == "too_sparse", na.rm = TRUE)
n_failed <- sum(class1_results$status == "model_failed", na.rm = TRUE)
n_low_variance <- sum(class1_results$status == "low_variance", na.rm = TRUE)
n_nominal <- sum(class1_results$primary_analysis %in% TRUE & class1_results$raw_p_value < 0.05, na.rm = TRUE)
n_fdr <- sum(class1_results$primary_analysis %in% TRUE & class1_results$q_global_class1_bh < 0.05, na.rm = TRUE)
best_q <- suppressWarnings(min(class1_results$q_global_class1_bh[class1_results$primary_analysis %in% TRUE], na.rm = TRUE))
if (!is.finite(best_q)) {
  best_q <- NA_real_
}

window_paths <- paste(class1_windows$forward_timepoint_path, collapse = " | ")

# Nature-style visual choices: restrained colors, high contrast, colorblind-safe
# accents, light rules, and no decorative gradients.
palette <- list(
  ink = "#1f2937",
  muted = "#4b5563",
  faint = "#e5e7eb",
  input = "#f8fafc",
  window = "#eff6ff",
  model = "#f0fdf4",
  status = "#fff7ed",
  fdr = "#fefce8",
  outcome = "#fef2f2",
  blue = "#2563eb",
  green = "#16a34a",
  orange = "#ea580c",
  amber = "#ca8a04",
  red = "#dc2626",
  grey = "#64748b"
)

workflow_nodes <- data.frame(
  step = 1:6,
  y = c(8.85, 7.25, 5.65, 4.05, 2.45, 0.85),
  fill = c(palette$input, palette$window, palette$model, palette$status, palette$fdr, palette$outcome),
  border = c(palette$grey, palette$blue, palette$green, palette$orange, palette$amber, palette$red),
  title = c(
    "Processed longitudinal MAG table",
    "Predefined biological windows",
    "Per MAG x window LMM",
    "Model-status classification",
    "Global BH-FDR correction",
    "Class 1 result summary"
  ),
  body = c(
    paste0("cleaned_data.rds\n", fmt_int(n_mags), " MAGs with repeated samples"),
    paste0(fmt_int(n_windows), " Class 1 windows, predefined before MAG p-values\n", window_paths),
    paste0(
      "One model per eligible MAG-window pair\n",
      "log10(Y + pseudocount) ~ window_order + (1 | Subject)\n",
      "window_order = 0, 1, 2 within each three-timepoint window"
    ),
    paste0(
      fmt_int(n_tests), " attempted tests\n",
      fmt_int(n_primary), " lmm_ok (", fmt_pct(n_primary / n_tests), "); ",
      fmt_int(n_singular), " singular; ", fmt_int(n_sparse), " too sparse"
    ),
    paste0(
      "FDR family = all predefined MAG-window tests\n",
      "Primary claim uses q_global_class1_bh, not per-window q-values"
    ),
    paste0(
      fmt_int(n_nominal), " nominal raw p < 0.05 candidates\n",
      fmt_int(n_fdr), " global FDR-significant tests; best q = ", fmt_q(best_q)
    )
  ),
  stringsAsFactors = FALSE
)

workflow_nodes$x <- 5
workflow_nodes$xmin <- 1.25
workflow_nodes$xmax <- 8.75
workflow_nodes$ymin <- workflow_nodes$y - 0.55
workflow_nodes$ymax <- workflow_nodes$y + 0.55

arrows <- data.frame(
  x = 5,
  xend = 5,
  y = workflow_nodes$ymin[-nrow(workflow_nodes)] - 0.06,
  yend = workflow_nodes$ymax[-1] + 0.06
)

side_notes <- data.frame(
  x = c(9.25, 9.25, 9.25),
  y = workflow_nodes$y[c(2, 4, 6)],
  label = c(
    "Predefined only:\nnot all possible\nthree-timepoint scans",
    paste0("Retained but excluded\nfrom primary claims:\n", fmt_int(n_singular + n_sparse + n_low_variance + n_failed), " flagged tests"),
    "Conservative conclusion:\nno FDR-significant\nClass 1 discoveries"
  ),
  color = c(palette$blue, palette$orange, palette$red),
  stringsAsFactors = FALSE
)

status_df <- data.frame(
  status = c("lmm_ok", "lmm_singular", "too_sparse", "low_variance", "model_failed"),
  count = c(n_primary, n_singular, n_sparse, n_low_variance, n_failed),
  stringsAsFactors = FALSE
) |>
  dplyr::filter(count > 0) |>
  dplyr::mutate(
    fraction = count / sum(count),
    xmin = 1.25 + cumsum(dplyr::lag(fraction, default = 0)) * 7.5,
    xmax = 1.25 + cumsum(fraction) * 7.5,
    ymin = workflow_nodes$ymin[workflow_nodes$step == 4] + 0.08,
    ymax = workflow_nodes$ymin[workflow_nodes$step == 4] + 0.24,
    label_x = (xmin + xmax) / 2
  )

status_colors <- c(
  lmm_ok = palette$green,
  lmm_singular = palette$amber,
  too_sparse = palette$orange,
  low_variance = "#7c3aed",
  model_failed = palette$red
)

figure <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = workflow_nodes,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = workflow_nodes$fill,
    color = workflow_nodes$border,
    linewidth = 0.75
  ) +
  ggplot2::geom_segment(
    data = arrows,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.6,
    color = palette$ink,
    arrow = grid::arrow(length = grid::unit(0.18, "inches"), type = "closed")
  ) +
  ggplot2::geom_text(
    data = workflow_nodes,
    ggplot2::aes(x = xmin + 0.25, y = y + 0.23, label = title),
    hjust = 0,
    size = 4.0,
    fontface = "bold",
    color = palette$ink
  ) +
  ggplot2::geom_text(
    data = workflow_nodes,
    ggplot2::aes(x = xmin + 0.25, y = y - 0.13, label = body),
    hjust = 0,
    size = 2.9,
    lineheight = 0.92,
    color = palette$muted
  ) +
  ggplot2::geom_rect(
    data = status_df,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = status),
    color = "white",
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    data = status_df,
    ggplot2::aes(x = label_x, y = workflow_nodes$ymin[workflow_nodes$step == 4] - 0.10, label = paste0(status, "\n", fmt_int(count))),
    size = 2.35,
    lineheight = 0.85,
    color = palette$ink
  ) +
  ggplot2::geom_curve(
    data = side_notes,
    ggplot2::aes(x = x - 0.15, xend = 8.82, y = y, yend = y, color = color),
    curvature = 0,
    linewidth = 0.5,
    arrow = grid::arrow(length = grid::unit(0.11, "inches"), type = "closed")
  ) +
  ggplot2::geom_label(
    data = side_notes,
    ggplot2::aes(x = x, y = y, label = label, color = color),
    fill = "white",
    label.size = 0.4,
    label.r = grid::unit(0.08, "inches"),
    hjust = 0,
    size = 2.8,
    lineheight = 0.9,
    fontface = "bold"
  ) +
  ggplot2::annotate(
    "text",
    x = 5,
    y = 10.15,
    label = "Main Figure 2. Reproducible Class 1 predefined-window LMM workflow",
    size = 5.0,
    fontface = "bold",
    color = palette$ink
  ) +
  ggplot2::annotate(
    "text",
    x = 5,
    y = 9.75,
    label = "The analysis uses only biologically predefined windows and applies one global BH-FDR correction across all MAG-window tests.",
    size = 3.2,
    color = palette$muted
  ) +
  ggplot2::annotate(
    "text",
    x = 5,
    y = 0.12,
    label = paste0(
      "Interpretation: nominal candidates are hypothesis-generating unless q_global_class1_bh < 0.05. ",
      "Overlapping windows are related secondary stage-local tests, not independent discoveries."
    ),
    size = 2.75,
    color = palette$muted
  ) +
  ggplot2::scale_fill_manual(values = status_colors, guide = "none") +
  ggplot2::scale_color_identity() +
  ggplot2::coord_cartesian(xlim = c(0.35, 11.1), ylim = c(0, 10.35), clip = "off") +
  ggplot2::theme_void(base_size = 11) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.margin = ggplot2::margin(12, 18, 12, 18)
  )

png_path <- file.path(output_dir, "figure_2_class1_analysis_workflow.png")
pdf_path <- file.path(output_dir, "figure_2_class1_analysis_workflow.pdf")
svg_path <- file.path(output_dir, "figure_2_class1_analysis_workflow.svg")

ggplot2::ggsave(png_path, figure, width = 10.8, height = 8.0, dpi = 320, bg = "white")
ggplot2::ggsave(pdf_path, figure, width = 10.8, height = 8.0, bg = "white")
ggplot2::ggsave(svg_path, figure, width = 10.8, height = 8.0, bg = "white", device = svglite::svglite)

message("Main Figure 2 Class 1 workflow diagram created:")
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Inputs:")
message(" - ", results_path)
message(" - ", windows_path)
message("Summary:")
message(" - MAGs: ", fmt_int(n_mags))
message(" - Windows: ", fmt_int(n_windows))
message(" - MAG-window tests: ", fmt_int(n_tests))
message(" - lmm_ok: ", fmt_int(n_primary))
message(" - nominal raw p < 0.05: ", fmt_int(n_nominal))
message(" - FDR-significant tests: ", fmt_int(n_fdr))
message(" - best q_global_class1_bh: ", fmt_q(best_q))
