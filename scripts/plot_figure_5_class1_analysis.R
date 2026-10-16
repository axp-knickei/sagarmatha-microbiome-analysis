# ============================================
# plot_figure_5_class1_analysis.R
# ============================================
# Main Figure 5: exploratory MAG-window association landscape.
# A volcano-style plot shows effect size versus raw-p evidence for Class 1
# primary-analysis MAG-window tests, without implying FDR-supported discovery.

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
  dplyr::filter(!is.na(estimate_window_order), !is.na(raw_p_value), raw_p_value > 0) |>
  dplyr::mutate(
    window_label = factor(window_label, levels = window_order),
    neg_log10_raw_p = -log10(raw_p_value),
    nominal_raw_p_lt_0_05 = raw_p_value < 0.05,
    fdr_significant = q_global_class1_bh < 0.05,
    point_alpha = ifelse(nominal_raw_p_lt_0_05, 0.92, 0.42),
    point_size = ifelse(nominal_raw_p_lt_0_05, 2.1, 1.25)
  )

n_primary <- nrow(plot_data)
n_nominal <- sum(plot_data$nominal_raw_p_lt_0_05, na.rm = TRUE)
n_fdr <- sum(plot_data$fdr_significant, na.rm = TRUE)
best_q <- suppressWarnings(min(plot_data$q_global_class1_bh, na.rm = TRUE))
if (!is.finite(best_q)) {
  best_q <- NA_real_
}

x_limit <- max(abs(plot_data$estimate_window_order), na.rm = TRUE)
x_limit <- ceiling(x_limit * 10) / 10
if (!is.finite(x_limit) || x_limit == 0) {
  x_limit <- 1
}
y_limit <- max(plot_data$neg_log10_raw_p, na.rm = TRUE)
y_limit <- ceiling(y_limit * 10) / 10

nominal_line <- -log10(0.05)

top_by_window <- plot_data |>
  dplyr::group_by(window_label) |>
  dplyr::slice_min(raw_p_value, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    label = paste0(
      as.character(window_labels[as.character(window_label)]),
      "\nraw p=", signif(raw_p_value, 2),
      "; q=", fmt_q(q_global_class1_bh)
    ),
    hjust = ifelse(estimate_window_order >= 0, 0, 1),
    nudge_x = ifelse(estimate_window_order >= 0, 0.08, -0.08)
  )

window_counts <- plot_data |>
  dplyr::group_by(window_label) |>
  dplyr::summarise(
    n_primary = dplyr::n(),
    nominal = sum(nominal_raw_p_lt_0_05, na.rm = TRUE),
    min_q = min(q_global_class1_bh, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    window_index = as.integer(window_label),
    label = paste0(
      unname(window_labels[as.character(window_label)]),
      "\nprimary ", n_primary,
      " | nominal ", nominal
    )
  )

annotation_text <- paste0(
  fmt_int(n_nominal),
  " nominal candidates; ",
  fmt_int(n_fdr),
  " global Class 1 BH-FDR significant\nbest q = ",
  fmt_q(best_q)
)

fdr_note <- paste0(
  "BH-FDR controls the expected false discovery proportion across the predefined Class 1 test family;\n",
  "it does not convert individual nominal raw-p candidates into confirmed discoveries."
)

figure <- ggplot2::ggplot(plot_data, ggplot2::aes(x = estimate_window_order, y = neg_log10_raw_p)) +
  ggplot2::annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = nominal_line,
    ymax = Inf,
    fill = "#fef2f2",
    alpha = 0.42
  ) +
  ggplot2::geom_hline(
    yintercept = nominal_line,
    linewidth = 0.55,
    linetype = "22",
    color = "#dc2626"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.45,
    color = "#9ca3af"
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(plot_data, !nominal_raw_p_lt_0_05),
    ggplot2::aes(color = window_label),
    size = 1.25,
    alpha = 0.38,
    shape = 16
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(plot_data, nominal_raw_p_lt_0_05),
    ggplot2::aes(color = window_label),
    size = 2.25,
    alpha = 0.92,
    shape = 21,
    fill = "white",
    stroke = 1.05
  ) +
  ggplot2::geom_text(
    data = top_by_window,
    ggplot2::aes(
      x = estimate_window_order + nudge_x,
      y = neg_log10_raw_p,
      label = paste0("W", as.integer(window_label)),
      hjust = hjust
    ),
    size = 2.8,
    fontface = "bold",
    color = "#111827",
    inherit.aes = FALSE
  ) +
  ggplot2::annotate(
    "label",
    x = x_limit * 0.96,
    y = y_limit * 0.94,
    label = annotation_text,
    hjust = 1,
    size = 3.35,
    lineheight = 0.95,
    label.r = grid::unit(0.08, "inches"),
    label.size = 0.35,
    fill = "white",
    color = "#111827"
  ) +
  ggplot2::annotate(
    "text",
    x = x_limit * 0.96,
    y = nominal_line + 0.08,
    label = "raw p = 0.05\nnominal threshold",
    hjust = 1,
    vjust = 0,
    size = 2.85,
    color = "#991b1b",
    lineheight = 0.9
  ) +
  ggplot2::annotate(
    "text",
    x = -x_limit,
    y = -0.46,
    label = fdr_note,
    hjust = 0,
    size = 2.75,
    color = "#4b5563",
    lineheight = 0.9
  ) +
  ggplot2::scale_color_manual(
    values = window_palette,
    breaks = window_order,
    labels = unname(window_labels[window_order]),
    name = "Predefined window"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(-x_limit * 1.08, x_limit * 1.08),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(-0.68, y_limit * 1.08),
    breaks = seq(0, ceiling(y_limit), by = 1),
    expand = ggplot2::expansion(mult = c(0, 0.03))
  ) +
  ggplot2::labs(
    title = "Main Figure 5 | Exploratory MAG-window association landscape",
    subtitle = paste0(
      "Volcano-style display of ", fmt_int(n_primary),
      " primary Class 1 MAG-window LMM tests; outlined points indicate nominal raw p < 0.05."
    ),
    x = "Window-order estimate from LMM",
    y = expression(-log[10]("raw p-value")),
    caption = paste0(
      "Model: log10(Y + pseudocount) ~ window_order + (1 | Subject). ",
      "The plot visualizes raw-p evidence and effect direction; primary significance requires q_global_class1_bh < 0.05."
    )
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#e5e7eb", linewidth = 0.35),
    plot.title = ggplot2::element_text(face = "bold", size = 15.5, color = "#111827"),
    plot.subtitle = ggplot2::element_text(size = 10.3, color = "#4b5563", margin = ggplot2::margin(b = 8)),
    axis.title = ggplot2::element_text(size = 10.6, color = "#374151"),
    axis.text = ggplot2::element_text(size = 9.2, color = "#4b5563"),
    legend.position = "right",
    legend.title = ggplot2::element_text(face = "bold", size = 9.3),
    legend.text = ggplot2::element_text(size = 8.2, lineheight = 0.86),
    legend.key.height = grid::unit(0.48, "inches"),
    plot.caption = ggplot2::element_text(size = 8.3, color = "#4b5563", hjust = 0, margin = ggplot2::margin(t = 10)),
    plot.margin = ggplot2::margin(14, 18, 18, 16)
  )

png_path <- file.path(output_dir, "figure_5_class1_analysis_global_fdr_volcano.png")
pdf_path <- file.path(output_dir, "figure_5_class1_analysis_global_fdr_volcano.pdf")
svg_path <- file.path(output_dir, "figure_5_class1_analysis_global_fdr_volcano.svg")

ggplot2::ggsave(png_path, figure, width = 11.2, height = 7.6, dpi = 360, bg = "white")
ggplot2::ggsave(pdf_path, figure, width = 11.2, height = 7.6, bg = "white")
ggplot2::ggsave(svg_path, figure, width = 11.2, height = 7.6, bg = "white", device = svglite::svglite)

message("Main Figure 5 Class 1 volcano-style association landscape created:")
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Inputs:")
message(" - ", results_path)
message(" - ", windows_path)
message("Summary:")
message(" - Primary lmm_ok tests plotted: ", fmt_int(n_primary))
message(" - Nominal raw p < 0.05 candidates: ", fmt_int(n_nominal))
message(" - Global Class 1 BH-FDR significant tests: ", fmt_int(n_fdr))
message(" - Best q_global_class1_bh: ", fmt_q(best_q))
