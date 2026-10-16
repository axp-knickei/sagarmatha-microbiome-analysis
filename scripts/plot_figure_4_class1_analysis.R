# ============================================
# plot_figure_4_class1_analysis.R
# ============================================
# Main Figure 4: per-window model usability for the Class 1 analysis.
# The figure summarizes model status counts per predefined biological window.

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
fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)

status_levels <- c("lmm_ok", "lmm_singular", "too_sparse", "low_variance", "model_failed")
status_palette <- c(
  lmm_ok = "#0072B2",       # Okabe-Ito blue
  lmm_singular = "#E69F00", # Okabe-Ito orange
  too_sparse = "#D55E00",   # Okabe-Ito vermillion
  low_variance = "#CC79A7", # Okabe-Ito reddish purple
  model_failed = "#4B5563"  # neutral dark grey
)

compact_window_names <- c(
  initial_ascent_to_high_altitude = "Initial ascent",
  ascent_and_intermediate_adjustment = "Ascent adjustment",
  high_climb_6500m_exposure = "6500 m climb",
  pre_summit_staging_to_post_summit = "Pre/post summit",
  summit_push_and_descent_recovery_transition = "Summit + descent",
  post_expedition_low_altitude_recovery = "Recovery"
)

window_order <- class1_windows$window_label
window_labels <- setNames(
  paste0(
    "W", seq_along(window_order), "\n",
    gsub(" -> ", "-", class1_windows$forward_timepoint_path), "\n",
    unname(compact_window_names[window_order])
  ),
  window_order
)

summit_windows <- c(
  "pre_summit_staging_to_post_summit",
  "summit_push_and_descent_recovery_transition"
)

# Count statuses per predefined window. Explicitly complete zero-count statuses
# so low_variance/model_failed are represented even when absent.
window_status <- class1_results |>
  dplyr::mutate(
    window_label = factor(window_label, levels = window_order),
    status = factor(status, levels = status_levels)
  ) |>
  dplyr::count(window_label, status, name = "n")

all_combinations <- expand.grid(
  window_label = factor(window_order, levels = window_order),
  status = factor(status_levels, levels = status_levels),
  KEEP.OUT.ATTRS = FALSE
)

window_status <- all_combinations |>
  dplyr::left_join(window_status, by = c("window_label", "status")) |>
  dplyr::mutate(
    n = dplyr::coalesce(n, 0L),
    status = factor(status, levels = status_levels)
  ) |>
  dplyr::group_by(window_label) |>
  dplyr::mutate(
    total = sum(n),
    fraction = n / total,
    ymax = cumsum(n),
    ymin = ymax - n,
    label_y = ymin + n / 2
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    window_index = as.integer(window_label),
    is_summit_related = as.character(window_label) %in% summit_windows
  )

window_totals <- window_status |>
  dplyr::group_by(window_label, window_index, is_summit_related) |>
  dplyr::summarise(
    total = unique(total),
    lmm_ok = sum(n[status == "lmm_ok"]),
    lmm_singular = sum(n[status == "lmm_singular"]),
    too_sparse = sum(n[status == "too_sparse"]),
    low_variance = sum(n[status == "low_variance"]),
    model_failed = sum(n[status == "model_failed"]),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    percent_lmm_ok = lmm_ok / total,
    top_label = paste0("OK ", lmm_ok, "\n", fmt_pct(percent_lmm_ok)),
    sparse_label = paste0("Sparse ", too_sparse),
    summit_note = ifelse(is_summit_related, "summit-related\nwindow", "")
  )

total_tests <- sum(window_totals$total)
total_ok <- sum(window_totals$lmm_ok)
total_singular <- sum(window_totals$lmm_singular)
total_sparse <- sum(window_totals$too_sparse)
total_low_variance <- sum(window_totals$low_variance)
total_failed <- sum(window_totals$model_failed)

summit_ok_min <- min(window_totals$percent_lmm_ok[window_totals$is_summit_related])
summit_support_text <- paste0(
  "Summit-related windows retained ",
  fmt_pct(summit_ok_min),
  " or more lmm_ok tests,\n",
  "supporting that these windows were modelable and not lost to sparsity."
)

audit_text <- paste0(
  "Total tests: ", fmt_int(total_tests),
  "  |  lmm_ok: ", fmt_int(total_ok),
  "  |  lmm_singular: ", fmt_int(total_singular),
  "  |  too_sparse: ", fmt_int(total_sparse),
  "  |  low_variance: ", fmt_int(total_low_variance),
  "  |  model_failed: ", fmt_int(total_failed)
)

bar_label_data <- window_status |>
  dplyr::filter(n > 0, status %in% c("lmm_ok", "too_sparse")) |>
  dplyr::mutate(
    label = ifelse(status == "lmm_ok", paste0(n, " OK"), paste0(n, " sparse")),
    label_color = ifelse(status == "lmm_ok", "white", "white")
  )

summit_band <- data.frame(
  xmin = min(window_totals$window_index[window_totals$is_summit_related]) - 0.52,
  xmax = max(window_totals$window_index[window_totals$is_summit_related]) + 0.52,
  ymin = -28,
  ymax = max(window_totals$total) + 76
)

figure <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = summit_band,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#fff7ed",
    color = "#fed7aa",
    linewidth = 0.45
  ) +
  ggplot2::geom_rect(
    data = window_status,
    ggplot2::aes(
      xmin = window_index - 0.36,
      xmax = window_index + 0.36,
      ymin = ymin,
      ymax = ymax,
      fill = status
    ),
    color = "white",
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    data = bar_label_data,
    ggplot2::aes(x = window_index, y = label_y, label = label),
    size = 3.0,
    fontface = "bold",
    color = "white",
    lineheight = 0.9
  ) +
  ggplot2::geom_text(
    data = window_totals,
    ggplot2::aes(x = window_index, y = total + 18, label = top_label),
    size = 3.2,
    fontface = "bold",
    lineheight = 0.88,
    color = "#111827"
  ) +
  ggplot2::annotate(
    "segment",
    x = summit_band$xmin,
    xend = summit_band$xmax,
    y = max(window_totals$total) + 45,
    yend = max(window_totals$total) + 45,
    color = "#ea580c",
    linewidth = 0.7
  ) +
  ggplot2::annotate(
    "text",
    x = mean(c(summit_band$xmin, summit_band$xmax)),
    y = max(window_totals$total) + 61,
    label = "Summit-related windows",
    size = 3.35,
    fontface = "bold",
    color = "#9a3412"
  ) +
  ggplot2::annotate(
    "text",
    x = 0.5,
    y = max(window_totals$total) + 98,
    label = "Main Figure 4 | Per-window model usability in Class 1 analysis",
    hjust = 0,
    size = 5.1,
    fontface = "bold",
    color = "#111827"
  ) +
  ggplot2::annotate(
    "text",
    x = 0.5,
    y = max(window_totals$total) + 82,
    label = "Stacked bars show model-status composition for each biologically predefined three-timepoint window.",
    hjust = 0,
    size = 3.35,
    color = "#4b5563"
  ) +
  ggplot2::annotate(
    "text",
    x = 0.5,
    y = max(window_totals$total) + 67,
    label = summit_support_text,
    hjust = 0,
    size = 3.0,
    lineheight = 0.92,
    color = "#9a3412"
  ) +
  ggplot2::annotate(
    "text",
    x = 0.5,
    y = -52,
    label = audit_text,
    hjust = 0,
    size = 3.15,
    color = "#374151"
  ) +
  ggplot2::scale_fill_manual(
    values = status_palette,
    breaks = status_levels,
    labels = c("lmm_ok", "lmm_singular", "too_sparse", "low_variance", "model_failed"),
    drop = FALSE,
    name = "Model status"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq_along(window_order),
    labels = unname(window_labels[window_order]),
    limits = c(0.35, length(window_order) + 0.65),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 441, by = 100),
    limits = c(-72, max(window_totals$total) + 112),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    x = NULL,
    y = "MAG-window tests per predefined window",
    caption = paste0(
      "Primary-analysis claims use lmm_ok tests only. Flagged tests are retained for auditability. ",
      "\nThe summit-push interval is figure-informed; this plot evaluates statistical usability, not biological significance."
    )
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(color = "#e5e7eb", linewidth = 0.35),
    axis.text.x = ggplot2::element_text(size = 8.8, color = "#111827", lineheight = 0.86, margin = ggplot2::margin(t = 8)),
    axis.text.y = ggplot2::element_text(size = 9, color = "#4b5563"),
    axis.title.y = ggplot2::element_text(size = 10.4, color = "#374151", margin = ggplot2::margin(r = 10)),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 9.5),
    legend.text = ggplot2::element_text(size = 9),
    legend.key.width = grid::unit(0.42, "inches"),
    plot.caption = ggplot2::element_text(size = 8.5, color = "#4b5563", hjust = 0, margin = ggplot2::margin(t = 12)),
    plot.margin = ggplot2::margin(18, 18, 16, 18)
  )

png_path <- file.path(output_dir, "figure_4_class1_analysis_per_window_usability.png")
pdf_path <- file.path(output_dir, "figure_4_class1_analysis_per_window_usability.pdf")
svg_path <- file.path(output_dir, "figure_4_class1_analysis_per_window_usability.svg")

ggplot2::ggsave(png_path, figure, width = 10.8, height = 7.4, dpi = 360, bg = "white")
ggplot2::ggsave(pdf_path, figure, width = 10.8, height = 7.4, bg = "white")
ggplot2::ggsave(svg_path, figure, width = 10.8, height = 7.4, bg = "white", device = svglite::svglite)

message("Main Figure 4 Class 1 per-window usability chart created:")
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Inputs:")
message(" - ", results_path)
message(" - ", windows_path)
message("Summary:")
message(" - Total tests: ", fmt_int(total_tests))
message(" - lmm_ok: ", fmt_int(total_ok))
message(" - lmm_singular: ", fmt_int(total_singular))
message(" - too_sparse: ", fmt_int(total_sparse))
message(" - low_variance: ", fmt_int(total_low_variance))
message(" - model_failed: ", fmt_int(total_failed))
message(" - Summit-related minimum lmm_ok rate: ", fmt_pct(summit_ok_min))
