# ============================================
# plot_figure_9_class1_analysis.R
# ============================================
# Supplementary/Main Figure 9: ranked global q-values by predefined window.
# This figure shows whether Class 1 evidence is broadly distributed across
# windows or concentrated in only a few top-ranked tests.

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
fmt_q <- function(x) ifelse(is.na(x), "NA", sprintf("%.4f", x))

window_order <- class1_windows$window_label
compact_window_names <- c(
  initial_ascent_to_high_altitude = "Initial ascent",
  ascent_and_intermediate_adjustment = "Ascent adjustment",
  high_climb_6500m_exposure = "6500 m climb",
  pre_summit_staging_to_post_summit = "Pre/post summit",
  summit_push_and_descent_recovery_transition = "Summit + descent",
  post_expedition_low_altitude_recovery = "Recovery"
)
window_labels <- setNames(
  paste0(
    "W", seq_along(window_order), ": ",
    gsub(" -> ", "-", class1_windows$forward_timepoint_path),
    "\n",
    unname(compact_window_names[window_order])
  ),
  window_order
)

window_palette <- c(
  initial_ascent_to_high_altitude = "#0072B2",
  ascent_and_intermediate_adjustment = "#56B4E9",
  high_climb_6500m_exposure = "#009E73",
  pre_summit_staging_to_post_summit = "#D55E00",
  summit_push_and_descent_recovery_transition = "#E69F00",
  post_expedition_low_altitude_recovery = "#CC79A7"
)

plot_data <- class1_results |>
  dplyr::filter(primary_analysis %in% TRUE) |>
  dplyr::filter(!is.na(q_global_class1_bh), !is.na(raw_p_value)) |>
  dplyr::mutate(window_label = factor(window_label, levels = window_order)) |>
  dplyr::group_by(window_label) |>
  dplyr::arrange(q_global_class1_bh, raw_p_value, .by_group = TRUE) |>
  dplyr::mutate(rank_within_window = dplyr::row_number()) |>
  dplyr::ungroup()

if (nrow(plot_data) == 0) {
  stop("No primary Class 1 q-values found for ranked q-value plotting.", call. = FALSE)
}

window_summary <- plot_data |>
  dplyr::group_by(window_label) |>
  dplyr::summarise(
    n_primary = dplyr::n(),
    best_q = min(q_global_class1_bh, na.rm = TRUE),
    nominal = sum(raw_p_value < 0.05, na.rm = TRUE),
    fdr_significant = sum(q_global_class1_bh < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    window_index = as.integer(window_label),
    label = paste0(
      "best q=", fmt_q(best_q),
      "\nnominal=", nominal,
      " | FDR=", fdr_significant
    )
  )

n_primary <- nrow(plot_data)
n_nominal <- sum(plot_data$raw_p_value < 0.05, na.rm = TRUE)
n_fdr <- sum(plot_data$q_global_class1_bh < 0.05, na.rm = TRUE)
best_q <- min(plot_data$q_global_class1_bh, na.rm = TRUE)

# Show the full global q-value landscape. Earlier zoomed versions hid windows
# whose best global q-values were far above 0.05, which weakened the audit goal.
q_display_max <- 1
plot_data_display <- plot_data

window_summary_display <- window_summary |>
  dplyr::mutate(label_y = pmin(best_q + 0.055, 0.92))

figure <- ggplot2::ggplot(
  plot_data_display,
  ggplot2::aes(x = rank_within_window, y = q_global_class1_bh, color = window_label)
) +
  ggplot2::annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 0,
    ymax = 0.05,
    fill = "#fee2e2",
    alpha = 0.48
  ) +
  ggplot2::geom_hline(yintercept = 0.05, color = "#dc2626", linetype = "22", linewidth = 0.58) +
  ggplot2::geom_line(linewidth = 0.55, alpha = 0.55) +
  ggplot2::geom_point(size = 1.25, alpha = 0.72) +
  ggplot2::geom_point(
    data = plot_data |>
      dplyr::group_by(window_label) |>
      dplyr::slice_min(q_global_class1_bh, n = 1, with_ties = FALSE) |>
      dplyr::ungroup(),
    ggplot2::aes(x = rank_within_window, y = q_global_class1_bh),
    shape = 21,
    fill = "white",
    size = 3.1,
    stroke = 1.0,
    inherit.aes = TRUE
  ) +
  ggplot2::geom_text(
    data = window_summary_display,
    ggplot2::aes(x = Inf, y = label_y, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 0,
    size = 2.75,
    lineheight = 0.9,
    color = "#111827"
  ) +
  ggplot2::facet_wrap(
    ~ window_label,
    ncol = 3,
    labeller = ggplot2::labeller(window_label = ggplot2::as_labeller(window_labels))
  ) +
  ggplot2::scale_color_manual(values = window_palette, drop = FALSE, guide = "none") +
  ggplot2::scale_y_continuous(
    limits = c(0, q_display_max),
    breaks = seq(0, q_display_max, by = 0.2),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = "Main Figure 9 | Ranked global Class 1 q-values by predefined window",
    subtitle = paste0(
      "Facets compare whether evidence is broad within a window or concentrated among only the top-ranked tests; ",
      "all q-values are global Class 1 BH-FDR values."
    ),
    x = "Rank within window by q_global_class1_bh",
    y = "q_global_class1_bh",
    caption = paste0(
      "Primary lmm_ok tests: ", fmt_int(n_primary),
      "; nominal raw p < 0.05: ", fmt_int(n_nominal),
      "; global q < 0.05: ", fmt_int(n_fdr),
      "; best q: ", fmt_q(best_q), "."
    )
  ) +
  ggplot2::theme_minimal(base_size = 10.5) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#e5e7eb", linewidth = 0.32),
    strip.text = ggplot2::element_text(face = "bold", size = 8.7, lineheight = 0.88, color = "#111827"),
    strip.background = ggplot2::element_rect(fill = "#f8fafc", color = "#d1d5db", linewidth = 0.35),
    plot.title = ggplot2::element_text(face = "bold", size = 15.5, color = "#111827"),
    plot.subtitle = ggplot2::element_text(size = 10.1, color = "#4b5563", margin = ggplot2::margin(b = 8)),
    axis.title = ggplot2::element_text(size = 10.4, color = "#374151"),
    axis.text = ggplot2::element_text(size = 8.8, color = "#4b5563"),
    panel.spacing = grid::unit(0.85, "lines"),
    plot.caption = ggplot2::element_text(size = 8.3, color = "#4b5563", hjust = 0, margin = ggplot2::margin(t = 10)),
    plot.margin = ggplot2::margin(14, 16, 14, 16)
  )

png_path <- file.path(output_dir, "figure_9_class1_analysis_ranked_q_by_window.png")
pdf_path <- file.path(output_dir, "figure_9_class1_analysis_ranked_q_by_window.pdf")
svg_path <- file.path(output_dir, "figure_9_class1_analysis_ranked_q_by_window.svg")

ggplot2::ggsave(png_path, figure, width = 11.4, height = 7.8, dpi = 360, bg = "white")
ggplot2::ggsave(pdf_path, figure, width = 11.4, height = 7.8, bg = "white")
ggplot2::ggsave(svg_path, figure, width = 11.4, height = 7.8, bg = "white", device = svglite::svglite)

message("Main Figure 9 Class 1 ranked q-value by window figure created:")
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Inputs:")
message(" - ", results_path)
message(" - ", windows_path)
message("Summary:")
message(" - Primary lmm_ok tests: ", fmt_int(n_primary))
message(" - Nominal raw p < 0.05 candidates: ", fmt_int(n_nominal))
message(" - Global Class 1 BH-FDR significant tests: ", fmt_int(n_fdr))
message(" - Best q_global_class1_bh: ", fmt_q(best_q))
message(" - Displayed q-value range: 0 to ", q_display_max)
