# ============================================
# plot_figure_2_class1_analysis_nature_style.R
# ============================================
# Main Figure 2: Class 1 predefined-window LMM workflow.
# Design goal: journal-style, compact, readable schematic suitable for a
# Nature-family manuscript draft. The script keeps the scientific content
# identical to the original workflow figure but improves layout hierarchy,
# whitespace, typography, export quality, and color accessibility.

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
  if (is.na(x) || !is.finite(x)) return("NA")
  sprintf("%.4f", x)
}

wrap_text <- function(x, width = 82) {
  vapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n", FUN.VALUE = character(1))
}

safe_divide <- function(num, den) {
  if (is.na(den) || den == 0) return(NA_real_)
  num / den
}

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
window_paths_wrapped <- wrap_text(window_paths, width = 88)

# Restrained, colorblind-aware palette. The main figure uses neutral grey with
# limited accent colors so that the hierarchy remains readable when printed.
pal <- list(
  ink = "#111827",
  text = "#374151",
  muted = "#6B7280",
  hair = "#D1D5DB",
  panel = "#F9FAFB",
  panel2 = "#F3F4F6",
  blue = "#0072B2",
  sky = "#EAF3FA",
  green = "#009E73",
  mint = "#EAF7F2",
  orange = "#D55E00",
  peach = "#FFF2E8",
  yellow = "#F0E442",
  cream = "#FFFCE6",
  red = "#CC3311",
  rose = "#FFF0EE",
  purple = "#785EF0"
)

workflow_nodes <- data.frame(
  step = 1:6,
  x = 3.9,
  y = c(8.40, 7.05, 5.70, 4.35, 3.00, 1.65),
  fill = c(pal$panel, pal$sky, pal$mint, pal$peach, pal$cream, pal$rose),
  accent = c(pal$muted, pal$blue, pal$green, pal$orange, "#B8860B", pal$red),
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
    paste0(fmt_int(n_tests), " attempted tests · ", fmt_int(n_primary), " primary lmm_ok tests (", fmt_pct(safe_divide(n_primary, n_tests)), ")\n", fmt_int(n_flagged), " flagged tests retained for transparency but excluded from primary claims"),
    "One correction family across all predefined MAG–window tests\nPrimary inference uses q_global_class1_bh, not per-window q-values",
    paste0(fmt_int(n_nominal), " nominal raw p < 0.05 candidates\n", fmt_int(n_fdr), " global FDR-significant tests · best q = ", fmt_q(best_q))
  ),
  stringsAsFactors = FALSE
)
workflow_nodes$xmin <- 0.70
workflow_nodes$xmax <- 7.10
workflow_nodes$ymin <- workflow_nodes$y - 0.47
workflow_nodes$ymax <- workflow_nodes$y + 0.47

arrows <- data.frame(
  x = 3.9,
  xend = 3.9,
  y = workflow_nodes$ymin[-nrow(workflow_nodes)] - 0.05,
  yend = workflow_nodes$ymax[-1] + 0.05
)

status_df <- data.frame(
  status = c("lmm_ok", "lmm_singular", "too_sparse", "low_variance", "model_failed"),
  label = c("OK", "Singular", "Sparse", "Low variance", "Failed"),
  count = c(n_primary, n_singular, n_sparse, n_low_variance, n_failed),
  color = c(pal$green, "#E69F00", pal$orange, pal$purple, pal$red),
  stringsAsFactors = FALSE
) |>
  dplyr::filter(count > 0) |>
  dplyr::mutate(
    fraction = count / sum(count),
    xmin = 0.95 + cumsum(dplyr::lag(fraction, default = 0)) * 5.90,
    xmax = 0.95 + cumsum(fraction) * 5.90,
    ymin = workflow_nodes$ymin[workflow_nodes$step == 4] + 0.07,
    ymax = workflow_nodes$ymin[workflow_nodes$step == 4] + 0.19,
    label_x = (xmin + xmax) / 2,
    label_show = ifelse(fraction >= 0.025, paste0(label, "\n", fmt_int(count)), "")
  )

right_panel <- data.frame(
  x = 8.05,
  y = c(7.05, 4.35, 1.65),
  title = c("Why this is reproducible", "What is excluded", "Conservative conclusion"),
  body = c(
    "Only predefined windows are tested; the analysis avoids scanning every possible three-timepoint combination.",
    paste0(fmt_int(n_flagged), " flagged tests are reported in the QC layer but do not support primary biological claims."),
    "Nominal candidates remain hypothesis-generating unless they pass the global BH–FDR threshold."
  ),
  accent = c(pal$blue, pal$orange, pal$red),
  stringsAsFactors = FALSE
)
right_panel$body <- wrap_text(right_panel$body, width = 34)

right_arrows <- data.frame(
  x = 7.70,
  xend = 7.12,
  y = right_panel$y,
  yend = right_panel$y,
  accent = right_panel$accent
)

figure <- ggplot2::ggplot() +
  # subtle page title rule
  ggplot2::annotate("segment", x = 0.70, xend = 9.85, y = 9.12, yend = 9.12, linewidth = 0.35, colour = pal$hair) +
  # node cards
  ggplot2::geom_rect(
    data = workflow_nodes,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = workflow_nodes$fill,
    colour = pal$hair,
    linewidth = 0.35
  ) +
  ggplot2::geom_segment(
    data = workflow_nodes,
    ggplot2::aes(x = xmin, xend = xmin, y = ymin, yend = ymax),
    colour = workflow_nodes$accent,
    linewidth = 2.4,
    lineend = "butt"
  ) +
  # step number circles
  ggplot2::geom_point(
    data = workflow_nodes,
    ggplot2::aes(x = xmin - 0.18, y = y),
    size = 7.5,
    colour = workflow_nodes$accent,
    fill = "white",
    shape = 21,
    stroke = 0.8
  ) +
  ggplot2::geom_text(
    data = workflow_nodes,
    ggplot2::aes(x = xmin - 0.18, y = y, label = step),
    size = 2.7,
    fontface = "bold",
    colour = pal$ink
  ) +
  # connecting arrows
  ggplot2::geom_segment(
    data = arrows,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.45,
    colour = pal$muted,
    arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed")
  ) +
  # text hierarchy inside cards
  ggplot2::geom_text(
    data = workflow_nodes,
    ggplot2::aes(x = xmin + 0.22, y = ymax - 0.15, label = kicker),
    hjust = 0,
    size = 2.05,
    fontface = "bold",
    colour = workflow_nodes$accent
  ) +
  ggplot2::geom_text(
    data = workflow_nodes,
    ggplot2::aes(x = xmin + 0.22, y = y + 0.12, label = title),
    hjust = 0,
    size = 3.35,
    fontface = "bold",
    colour = pal$ink
  ) +
  ggplot2::geom_text(
    data = workflow_nodes,
    ggplot2::aes(x = xmin + 0.22, y = y - 0.18, label = body),
    hjust = 0,
    size = 2.35,
    lineheight = 0.90,
    colour = pal$text
  ) +
  # embedded QC bar
  ggplot2::geom_rect(
    data = status_df,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = status_df$color,
    colour = "white",
    linewidth = 0.25
  ) +
  ggplot2::geom_text(
    data = status_df,
    ggplot2::aes(x = label_x, y = workflow_nodes$ymin[workflow_nodes$step == 4] - 0.14, label = label_show),
    size = 1.95,
    lineheight = 0.88,
    colour = pal$text
  ) +
  # right-side interpretation panel
  ggplot2::geom_segment(
    data = right_arrows,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    colour = right_arrows$accent,
    linewidth = 0.45,
    arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed")
  ) +
  ggplot2::geom_rect(
    data = right_panel,
    ggplot2::aes(xmin = x, xmax = x + 1.90, ymin = y - 0.47, ymax = y + 0.47),
    fill = "white",
    colour = pal$hair,
    linewidth = 0.35
  ) +
  ggplot2::geom_segment(
    data = right_panel,
    ggplot2::aes(x = x, xend = x, y = y - 0.47, yend = y + 0.47),
    colour = right_panel$accent,
    linewidth = 1.7
  ) +
  ggplot2::geom_text(
    data = right_panel,
    ggplot2::aes(x = x + 0.15, y = y + 0.23, label = title),
    hjust = 0,
    size = 2.45,
    fontface = "bold",
    colour = pal$ink
  ) +
  ggplot2::geom_text(
    data = right_panel,
    ggplot2::aes(x = x + 0.15, y = y - 0.13, label = body),
    hjust = 0,
    size = 2.10,
    lineheight = 0.88,
    colour = pal$text
  ) +
  # title and caption
  ggplot2::annotate(
    "text",
    x = 0.70,
    y = 9.55,
    hjust = 0,
    label = "Main Figure 2 | Reproducible Class 1 predefined-window LMM workflow",
    size = 4.15,
    fontface = "bold",
    colour = pal$ink
  ) +
  ggplot2::annotate(
    "text",
    x = 0.70,
    y = 9.28,
    hjust = 0,
    label = "Biologically predefined windows, model-status transparency, and one global BH–FDR family across all MAG–window tests.",
    size = 2.55,
    colour = pal$muted
  ) +
  ggplot2::annotate(
    "text",
    x = 0.70,
    y = 0.58,
    hjust = 0,
    label = "Interpretation: nominal raw-p candidates are hypothesis-generating unless q_global_class1_bh < 0.05. Overlapping windows are stage-local secondary tests, not independent discoveries.",
    size = 2.10,
    colour = pal$muted
  ) +
  ggplot2::coord_cartesian(xlim = c(0.30, 10.05), ylim = c(0.35, 9.75), clip = "off") +
  ggplot2::theme_void(base_size = 10, base_family = "sans") +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.margin = ggplot2::margin(12, 16, 12, 16)
  )

png_path <- file.path(output_dir, "figure_2_class1_analysis_workflow_nature_style.png")
pdf_path <- file.path(output_dir, "figure_2_class1_analysis_workflow_nature_style.pdf")
svg_path <- file.path(output_dir, "figure_2_class1_analysis_workflow_nature_style.svg")

ggplot2::ggsave(png_path, figure, width = 7.2, height = 6.9, dpi = 600, bg = "white")
ggplot2::ggsave(pdf_path, figure, width = 7.2, height = 6.9, bg = "white", device = grDevices::cairo_pdf)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggplot2::ggsave(svg_path, figure, width = 7.2, height = 6.9, bg = "white", device = svglite::svglite)
} else {
  message("Optional package 'svglite' not installed; SVG export skipped.")
}

message("Main Figure 2 Nature-style workflow diagram created:")
message(" - ", png_path)
message(" - ", pdf_path)
if (requireNamespace("svglite", quietly = TRUE)) message(" - ", svg_path)
message("Inputs:")
message(" - ", results_path)
message(" - ", windows_path)
message("Summary:")
message(" - MAGs: ", fmt_int(n_mags))
message(" - Windows: ", fmt_int(n_windows))
message(" - MAG-window tests: ", fmt_int(n_tests))
message(" - lmm_ok: ", fmt_int(n_primary))
message(" - flagged tests: ", fmt_int(n_flagged))
message(" - nominal raw p < 0.05: ", fmt_int(n_nominal))
message(" - FDR-significant tests: ", fmt_int(n_fdr))
message(" - best q_global_class1_bh: ", fmt_q(best_q))
