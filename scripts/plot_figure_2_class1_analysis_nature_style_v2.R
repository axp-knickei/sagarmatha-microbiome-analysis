# ============================================
# plot_figure_2_class1_analysis_nature_style_v2.R
# ============================================
# Main Figure 2: Class 1 predefined-window LMM workflow.
#
# Purpose
# -------
# Create a cleaner, publication-style workflow schematic for the Class 1
# predefined-window LMM analysis. This version fixes the overlap in the
# Model-status classification panel by giving the QC layer more vertical room
# and moving the status bar into its own dedicated inset band.
#
# Design principles
# -----------------
# - restrained Nature-style visual hierarchy
# - larger whitespace and safer text placement
# - no decorative gradients or heavy boxes
# - compact, interpretable status bar that never covers sentence text
# - high-resolution PNG plus vector PDF/SVG output

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

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_pct <- function(x) {
  if (is.na(x) || !is.finite(x)) return("NA")
  sprintf("%.1f%%", 100 * x)
}
fmt_q <- function(x) {
  if (is.na(x) || !is.finite(x)) return("NA")
  sprintf("%.4f", x)
}

safe_divide <- function(num, den) {
  if (is.na(den) || den == 0) return(NA_real_)
  num / den
}

wrap_one <- function(x, width = 80) {
  paste(strwrap(x, width = width), collapse = "\n")
}

# -----------------------------
# Summaries from real output
# -----------------------------
n_mags <- dplyr::n_distinct(class1_results$MAG)
n_windows <- nrow(class1_windows)
n_tests <- nrow(class1_results)
n_primary <- sum(class1_results$primary_analysis %in% TRUE, na.rm = TRUE)
n_singular <- sum(class1_results$status == "lmm_singular", na.rm = TRUE)
n_sparse <- sum(class1_results$status == "too_sparse", na.rm = TRUE)
n_failed <- sum(class1_results$status == "model_failed", na.rm = TRUE)
n_low_variance <- sum(class1_results$status == "low_variance", na.rm = TRUE)
n_flagged <- n_singular + n_sparse + n_low_variance + n_failed
n_nominal <- sum(class1_results$primary_analysis %in% TRUE & class1_results$raw_p_value < 0.05, na.rm = TRUE)
n_fdr <- sum(class1_results$primary_analysis %in% TRUE & class1_results$q_global_class1_bh < 0.05, na.rm = TRUE)
best_q <- suppressWarnings(min(class1_results$q_global_class1_bh[class1_results$primary_analysis %in% TRUE], na.rm = TRUE))
if (!is.finite(best_q)) best_q <- NA_real_

window_paths <- paste(class1_windows$forward_timepoint_path, collapse = "  |  ")
window_paths_wrapped <- wrap_one(window_paths, width = 74)

# -----------------------------
# Visual palette
# -----------------------------
# Restrained, high-contrast, colorblind-aware palette.
pal <- list(
  ink = "#111827",
  text = "#374151",
  muted = "#6B7280",
  hair = "#D7DCE2",
  hair2 = "#ECEFF3",
  paper = "#FFFFFF",
  shadow = "#F3F4F6",
  panel = "#FAFBFC",
  blue = "#0072B2",
  sky = "#EEF6FC",
  green = "#009E73",
  mint = "#EEF9F5",
  orange = "#D55E00",
  peach = "#FFF4EC",
  amber = "#B8860B",
  cream = "#FFFBEA",
  red = "#CC3311",
  rose = "#FFF2F0",
  purple = "#785EF0"
)

# -----------------------------
# Main workflow cards
# -----------------------------
# The QC card is intentionally taller than the other cards. This prevents the
# status bar from touching the explanatory text.
workflow_nodes <- data.frame(
  step = 1:6,
  y = c(8.62, 7.18, 5.76, 4.16, 2.55, 1.14),
  height = c(0.84, 1.08, 1.05, 1.42, 1.00, 1.00),
  fill = c(pal$panel, pal$sky, pal$mint, pal$peach, pal$cream, pal$rose),
  accent = c(pal$muted, pal$blue, pal$green, pal$orange, pal$amber, pal$red),
  kicker = c("INPUT", "DESIGN LOCK", "MODEL", "QC FILTER", "MULTIPLE TESTING", "PRIMARY CLAIM"),
  title = c(
    "Processed longitudinal MAG table",
    "Biologically predefined windows",
    "Per-MAG × window LMM",
    "Model-status classification",
    "Global BH–FDR correction",
    "Class 1 result summary"
  ),
  body = c(
    paste0("cleaned_data.rds · ", fmt_int(n_mags), " MAGs with repeated samples"),
    paste0(fmt_int(n_windows), " Class 1 windows specified before testing\n", window_paths_wrapped),
    "One model per eligible MAG–window pair\nlog10(Y + pseudocount) ~ window_order + (1 | Subject)\nwindow_order = 0, 1, 2 within each three-timepoint window",
    paste0(
      fmt_int(n_tests), " attempted tests · ", fmt_int(n_primary), " primary lmm_ok tests (",
      fmt_pct(safe_divide(n_primary, n_tests)), ")\n",
      "Flagged but retained: ", fmt_int(n_singular), " singular · ", fmt_int(n_sparse),
      " too sparse · ", fmt_int(n_low_variance + n_failed), " other"
    ),
    "One correction family across all predefined MAG–window tests\nPrimary inference uses q_global_class1_bh, not per-window q-values",
    paste0(
      fmt_int(n_nominal), " nominal raw p < 0.05 candidates\n",
      fmt_int(n_fdr), " global FDR-significant tests · best q = ", fmt_q(best_q)
    )
  ),
  stringsAsFactors = FALSE
)

workflow_nodes$xmin <- 0.74
workflow_nodes$xmax <- 7.15
workflow_nodes$ymin <- workflow_nodes$y - workflow_nodes$height / 2
workflow_nodes$ymax <- workflow_nodes$y + workflow_nodes$height / 2
workflow_nodes$kicker_y <- workflow_nodes$ymax - 0.17
workflow_nodes$title_y <- workflow_nodes$ymax - 0.39
workflow_nodes$body_y <- workflow_nodes$ymax - 0.68

arrows <- data.frame(
  x = 3.95,
  xend = 3.95,
  y = workflow_nodes$ymin[-nrow(workflow_nodes)] - 0.07,
  yend = workflow_nodes$ymax[-1] + 0.07
)

# -----------------------------
# Dedicated QC status inset
# -----------------------------
status_df <- data.frame(
  status = c("lmm_ok", "lmm_singular", "too_sparse", "low_variance", "model_failed"),
  label = c("OK", "Singular", "Sparse", "Low variance", "Failed"),
  count = c(n_primary, n_singular, n_sparse, n_low_variance, n_failed),
  fill = c(pal$green, "#E69F00", pal$orange, pal$purple, pal$red),
  stringsAsFactors = FALSE
) |>
  dplyr::filter(count > 0) |>
  dplyr::mutate(
    fraction = count / sum(count),
    xmin = 1.04 + cumsum(dplyr::lag(fraction, default = 0)) * 5.70,
    xmax = 1.04 + cumsum(fraction) * 5.70,
    ymin = workflow_nodes$ymin[workflow_nodes$step == 4] + 0.18,
    ymax = workflow_nodes$ymin[workflow_nodes$step == 4] + 0.32,
    label_x = (xmin + xmax) / 2,
    label_show = ifelse(fraction >= 0.055, paste0(label, "\n", fmt_int(count)), "")
  )

qc_bar_frame <- data.frame(
  xmin = 0.98,
  xmax = 6.82,
  ymin = workflow_nodes$ymin[workflow_nodes$step == 4] + 0.12,
  ymax = workflow_nodes$ymin[workflow_nodes$step == 4] + 0.49
)

status_legend <- status_df |>
  dplyr::mutate(
    legend_label = paste0(label, " ", fmt_int(count)),
    legend_x = c(1.05, 2.40, 3.85, 5.15, 6.15)[seq_len(dplyr::n())],
    legend_y = workflow_nodes$ymin[workflow_nodes$step == 4] + 0.06
  )

# -----------------------------
# Right-side interpretation cards
# -----------------------------
right_panel <- data.frame(
  x = 8.05,
  y = c(7.18, 4.16, 1.14),
  title = c("Why this is reproducible", "What is excluded", "Conservative conclusion"),
  body = c(
    "Only predefined windows are tested. This avoids scanning every possible three-timepoint combination.",
    paste0(fmt_int(n_flagged), " flagged tests are shown in the QC layer but excluded from primary biological claims."),
    "Nominal candidates remain hypothesis-generating unless they pass the global BH–FDR threshold."
  ),
  accent = c(pal$blue, pal$orange, pal$red),
  stringsAsFactors = FALSE
)
right_panel$body <- vapply(right_panel$body, wrap_one, character(1), width = 31)
right_panel$ymin <- right_panel$y - 0.55
right_panel$ymax <- right_panel$y + 0.55

right_arrows <- data.frame(
  x = 7.74,
  xend = 7.17,
  y = right_panel$y,
  yend = right_panel$y,
  accent = right_panel$accent
)

# -----------------------------
# Build figure
# -----------------------------
figure <- ggplot2::ggplot() +
  # page title rule
  ggplot2::annotate("segment", x = 0.72, xend = 9.95, y = 9.34, yend = 9.34, linewidth = 0.35, colour = pal$hair) +

  # card shadows: subtle offset, gives the figure a polished editorial feel
  ggplot2::geom_rect(
    data = workflow_nodes,
    ggplot2::aes(xmin = xmin + 0.03, xmax = xmax + 0.03, ymin = ymin - 0.03, ymax = ymax - 0.03),
    fill = pal$shadow,
    colour = NA
  ) +

  # main cards
  ggplot2::geom_rect(
    data = workflow_nodes,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    colour = pal$hair,
    linewidth = 0.35
  ) +
  ggplot2::geom_segment(
    data = workflow_nodes,
    ggplot2::aes(x = xmin, xend = xmin, y = ymin, yend = ymax, colour = accent),
    linewidth = 2.5,
    lineend = "butt"
  ) +

  # step number circles
  ggplot2::geom_point(
    data = workflow_nodes,
    ggplot2::aes(x = xmin - 0.18, y = y, colour = accent),
    fill = "white",
    shape = 21,
    size = 7.4,
    stroke = 0.85
  ) +
  ggplot2::geom_text(
    data = workflow_nodes,
    ggplot2::aes(x = xmin - 0.18, y = y, label = step),
    size = 2.65,
    fontface = "bold",
    colour = pal$ink
  ) +

  # flow arrows
  ggplot2::geom_segment(
    data = arrows,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.44,
    colour = pal$muted,
    arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed")
  ) +

  # text hierarchy in cards
  ggplot2::geom_text(
    data = workflow_nodes,
    ggplot2::aes(x = xmin + 0.22, y = kicker_y, label = kicker, colour = accent),
    hjust = 0,
    size = 2.00,
    fontface = "bold"
  ) +
  ggplot2::geom_text(
    data = workflow_nodes,
    ggplot2::aes(x = xmin + 0.22, y = title_y, label = title),
    hjust = 0,
    size = 3.26,
    fontface = "bold",
    colour = pal$ink
  ) +
  ggplot2::geom_text(
    data = workflow_nodes,
    ggplot2::aes(x = xmin + 0.22, y = body_y, label = body),
    hjust = 0,
    size = 2.26,
    lineheight = 0.90,
    colour = pal$text
  ) +

  # QC status inset frame and bar: isolated from sentence text
  ggplot2::geom_rect(
    data = qc_bar_frame,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    colour = pal$hair2,
    linewidth = 0.25
  ) +
  ggplot2::geom_text(
    data = qc_bar_frame,
    ggplot2::aes(x = xmin + 0.07, y = ymax - 0.10),
    label = "QC composition of attempted tests",
    hjust = 0,
    size = 1.72,
    fontface = "bold",
    colour = pal$muted
  ) +
  ggplot2::geom_rect(
    data = status_df,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    colour = "white",
    linewidth = 0.22
  ) +
  ggplot2::geom_text(
    data = status_df,
    ggplot2::aes(x = label_x, y = ymin - 0.10, label = label_show),
    size = 1.74,
    lineheight = 0.88,
    colour = pal$text
  ) +
  ggplot2::geom_point(
    data = status_legend,
    ggplot2::aes(x = legend_x, y = legend_y, fill = fill),
    shape = 22,
    size = 2.1,
    colour = "white",
    stroke = 0.2
  ) +
  ggplot2::geom_text(
    data = status_legend,
    ggplot2::aes(x = legend_x + 0.08, y = legend_y, label = legend_label),
    hjust = 0,
    size = 1.56,
    colour = pal$muted
  ) +

  # right-side notes
  ggplot2::geom_segment(
    data = right_arrows,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend, colour = accent),
    linewidth = 0.46,
    arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed")
  ) +
  ggplot2::geom_rect(
    data = right_panel,
    ggplot2::aes(xmin = x + 0.03, xmax = x + 1.93, ymin = ymin - 0.03, ymax = ymax - 0.03),
    fill = pal$shadow,
    colour = NA
  ) +
  ggplot2::geom_rect(
    data = right_panel,
    ggplot2::aes(xmin = x, xmax = x + 1.90, ymin = ymin, ymax = ymax),
    fill = "white",
    colour = pal$hair,
    linewidth = 0.34
  ) +
  ggplot2::geom_segment(
    data = right_panel,
    ggplot2::aes(x = x, xend = x, y = ymin, yend = ymax, colour = accent),
    linewidth = 1.75,
    lineend = "butt"
  ) +
  ggplot2::geom_text(
    data = right_panel,
    ggplot2::aes(x = x + 0.15, y = ymax - 0.23, label = title),
    hjust = 0,
    size = 2.34,
    fontface = "bold",
    colour = pal$ink
  ) +
  ggplot2::geom_text(
    data = right_panel,
    ggplot2::aes(x = x + 0.15, y = ymax - 0.54, label = body),
    hjust = 0,
    size = 1.96,
    lineheight = 0.88,
    colour = pal$text
  ) +

  # title and figure note
  ggplot2::annotate(
    "text",
    x = 0.72,
    y = 9.77,
    hjust = 0,
    label = "Main Figure 2 | Reproducible Class 1 predefined-window LMM workflow",
    size = 4.05,
    fontface = "bold",
    colour = pal$ink
  ) +
  ggplot2::annotate(
    "text",
    x = 0.72,
    y = 9.52,
    hjust = 0,
    label = "Biologically predefined windows, transparent model QC, and one global BH–FDR family across all MAG–window tests.",
    size = 2.45,
    colour = pal$muted
  ) +
  ggplot2::annotate(
    "text",
    x = 0.72,
    y = 0.32,
    hjust = 0,
    label = wrap_one(
      "Interpretation: nominal raw-p candidates are hypothesis-generating unless q_global_class1_bh < 0.05. Overlapping windows are stage-local secondary tests, not independent discoveries.",
      width = 140
    ),
    size = 1.98,
    colour = pal$muted
  ) +

  ggplot2::scale_fill_identity() +
  ggplot2::scale_colour_identity() +
  ggplot2::coord_cartesian(xlim = c(0.30, 10.05), ylim = c(0.10, 10.00), clip = "off") +
  ggplot2::theme_void(base_size = 10, base_family = "sans") +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.margin = ggplot2::margin(12, 18, 12, 18)
  )

# -----------------------------
# Export
# -----------------------------
# Save both the canonical filenames and explicit v2 filenames. This makes it
# easy to replace the earlier figure while keeping a versioned copy.
png_path <- file.path(output_dir, "figure_2_class1_analysis_workflow_nature_style.png")
pdf_path <- file.path(output_dir, "figure_2_class1_analysis_workflow_nature_style.pdf")
svg_path <- file.path(output_dir, "figure_2_class1_analysis_workflow_nature_style.svg")

png_path_v2 <- file.path(output_dir, "figure_2_class1_analysis_workflow_nature_style_v2.png")
pdf_path_v2 <- file.path(output_dir, "figure_2_class1_analysis_workflow_nature_style_v2.pdf")
svg_path_v2 <- file.path(output_dir, "figure_2_class1_analysis_workflow_nature_style_v2.svg")

for (path in c(png_path, png_path_v2)) {
  ggplot2::ggsave(path, figure, width = 7.4, height = 7.15, dpi = 600, bg = "white")
}
for (path in c(pdf_path, pdf_path_v2)) {
  ggplot2::ggsave(path, figure, width = 7.4, height = 7.15, bg = "white", device = grDevices::cairo_pdf)
}

if (requireNamespace("svglite", quietly = TRUE)) {
  for (path in c(svg_path, svg_path_v2)) {
    ggplot2::ggsave(path, figure, width = 7.4, height = 7.15, bg = "white", device = svglite::svglite)
  }
} else {
  message("Optional package 'svglite' not installed; SVG export skipped.")
}

message("Main Figure 2 Nature-style workflow diagram created:")
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", png_path_v2)
message(" - ", pdf_path_v2)
if (requireNamespace("svglite", quietly = TRUE)) {
  message(" - ", svg_path)
  message(" - ", svg_path_v2)
}
message("Inputs:")
message(" - ", results_path)
message(" - ", windows_path)
message("Summary:")
message(" - MAGs: ", fmt_int(n_mags))
message(" - Windows: ", fmt_int(n_windows))
message(" - MAG-window tests: ", fmt_int(n_tests))
message(" - lmm_ok: ", fmt_int(n_primary))
message(" - singular: ", fmt_int(n_singular))
message(" - too sparse: ", fmt_int(n_sparse))
message(" - other flagged: ", fmt_int(n_low_variance + n_failed))
message(" - nominal raw p < 0.05: ", fmt_int(n_nominal))
message(" - FDR-significant tests: ", fmt_int(n_fdr))
message(" - best q_global_class1_bh: ", fmt_q(best_q))
