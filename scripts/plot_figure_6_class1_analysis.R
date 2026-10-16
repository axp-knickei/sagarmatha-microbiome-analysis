# ============================================
# plot_figure_6_class1_analysis.R
# ============================================
# Main Figure 6: window-wise nominal candidate direction.
# This figure summarizes where nominal raw-p candidate signals concentrate,
# without treating them as FDR-significant discoveries.

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

nominal <- class1_results |>
  dplyr::filter(primary_analysis %in% TRUE, nominal_p_lt_0_05 %in% TRUE) |>
  dplyr::filter(direction %in% c("increasing", "decreasing")) |>
  dplyr::mutate(
    window_label = factor(window_label, levels = window_order),
    direction = factor(direction, levels = c("decreasing", "increasing"))
  )

all_combinations <- expand.grid(
  window_label = factor(window_order, levels = window_order),
  direction = factor(c("decreasing", "increasing"), levels = c("decreasing", "increasing")),
  KEEP.OUT.ATTRS = FALSE
)

direction_counts <- nominal |>
  dplyr::count(window_label, direction, name = "n")

direction_counts <- all_combinations |>
  dplyr::left_join(direction_counts, by = c("window_label", "direction")) |>
  dplyr::mutate(
    n = dplyr::coalesce(n, 0L),
    signed_n = ifelse(direction == "decreasing", -n, n),
    window_index = as.integer(window_label)
  )

window_totals <- direction_counts |>
  dplyr::group_by(window_label, window_index) |>
  dplyr::summarise(
    increasing = sum(n[direction == "increasing"]),
    decreasing = sum(n[direction == "decreasing"]),
    total = increasing + decreasing,
    balance = increasing - decreasing,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    label = paste0("total ", total),
    is_recovery = as.character(window_label) == "post_expedition_low_altitude_recovery"
  )

n_nominal <- sum(window_totals$total)
n_increasing <- sum(window_totals$increasing)
n_decreasing <- sum(window_totals$decreasing)
n_fdr <- sum(class1_results$primary_analysis %in% TRUE & class1_results$q_global_class1_bh < 0.05, na.rm = TRUE)
best_q <- suppressWarnings(min(class1_results$q_global_class1_bh[class1_results$primary_analysis %in% TRUE], na.rm = TRUE))
if (!is.finite(best_q)) {
  best_q <- NA_real_
}

max_count <- max(abs(direction_counts$signed_n), na.rm = TRUE)
x_limit <- max(10, ceiling(max_count / 5) * 5 + 5)

direction_palette <- c(
  decreasing = "#D55E00", # Okabe-Ito vermillion
  increasing = "#0072B2"  # Okabe-Ito blue
)

annotation_text <- paste0(
  fmt_int(n_nominal),
  " nominal candidates: ",
  fmt_int(n_increasing),
  " increasing, ",
  fmt_int(n_decreasing),
  " decreasing\n",
  fmt_int(n_fdr),
  " global Class 1 BH-FDR significant; best q = ",
  fmt_q(best_q)
)

recovery <- window_totals |>
  dplyr::filter(is_recovery)

recovery_text <- paste0(
  "Post-expedition recovery window\nis decreasing-skewed among nominal candidates\n",
  recovery$decreasing,
  " decreasing vs ",
  recovery$increasing,
  " increasing"
)

axis_breaks <- seq(-x_limit, x_limit, by = 5)
axis_labels <- abs(axis_breaks)

figure <- ggplot2::ggplot(direction_counts, ggplot2::aes(y = window_index, x = signed_n)) +
  ggplot2::annotate(
    "rect",
    xmin = -x_limit,
    xmax = 0,
    ymin = 0.4,
    ymax = length(window_order) + 0.6,
    fill = "#fff7ed",
    alpha = 0.65
  ) +
  ggplot2::annotate(
    "rect",
    xmin = 0,
    xmax = x_limit,
    ymin = 0.4,
    ymax = length(window_order) + 0.6,
    fill = "#eff6ff",
    alpha = 0.65
  ) +
  ggplot2::geom_vline(xintercept = 0, color = "#111827", linewidth = 0.45) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = signed_n, yend = window_index, color = direction),
    linewidth = 1.25,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = direction),
    size = 5.2,
    shape = 21,
    fill = "white",
    stroke = 1.35
  ) +
  ggplot2::geom_text(
    data = dplyr::filter(direction_counts, n > 0),
    ggplot2::aes(
      x = signed_n + ifelse(direction == "increasing", 1.05, -1.05),
      label = n,
      hjust = ifelse(direction == "increasing", 0, 1)
    ),
    size = 3.35,
    fontface = "bold",
    color = "#111827"
  ) +
  ggplot2::geom_text(
    data = window_totals,
    ggplot2::aes(x = 0, y = window_index + 0.28, label = label),
    size = 2.85,
    color = "#4b5563",
    inherit.aes = FALSE
  ) +
  ggplot2::annotate(
    "text",
    x = -x_limit * 0.74,
    y = length(window_order) + 0.88,
    label = "decreasing",
    size = 3.25,
    fontface = "bold",
    color = direction_palette[["decreasing"]]
  ) +
  ggplot2::annotate(
    "text",
    x = x_limit * 0.74,
    y = length(window_order) + 0.88,
    label = "increasing",
    size = 3.25,
    fontface = "bold",
    color = direction_palette[["increasing"]]
  ) +
  ggplot2::annotate(
    "label",
    x = -x_limit * 0.96,
    y = 0.72,
    label = annotation_text,
    hjust = 0,
    size = 3.15,
    lineheight = 0.92,
    label.r = grid::unit(0.08, "inches"),
    label.size = 0.35,
    fill = "white",
    color = "#111827"
  ) +
  ggplot2::annotate(
    "label",
    x = x_limit * 0.96,
    y = recovery$window_index,
    label = recovery_text,
    hjust = 1,
    size = 3.05,
    lineheight = 0.9,
    label.r = grid::unit(0.08, "inches"),
    label.size = 0.35,
    fill = "white",
    color = "#374151"
  ) +
  ggplot2::scale_color_manual(values = direction_palette, drop = FALSE, name = "Nominal direction") +
  ggplot2::scale_x_continuous(
    limits = c(-x_limit, x_limit),
    breaks = axis_breaks,
    labels = axis_labels,
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq_along(window_order),
    labels = unname(window_labels[window_order]),
    limits = c(0.35, length(window_order) + 1.05),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = "Main Figure 6 | Window-wise direction of nominal Class 1 candidates",
    subtitle = "Mirrored lollipop chart summarizes nominal raw-p signals without treating them as discoveries.",
    x = "Number of nominal raw p < 0.05 MAG-window candidates",
    y = NULL,
    caption = paste0(
      "Counts include primary lmm_ok tests with raw p < 0.05 only. ",
      "No Class 1 MAG-window test passed global BH-FDR correction, so these are hypothesis-generating patterns."
    )
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(color = "#e5e7eb", linewidth = 0.35),
    plot.title = ggplot2::element_text(face = "bold", size = 15.3, color = "#111827"),
    plot.subtitle = ggplot2::element_text(size = 10.3, color = "#4b5563", margin = ggplot2::margin(b = 8)),
    axis.text.y = ggplot2::element_text(size = 9.2, color = "#111827", lineheight = 0.88),
    axis.text.x = ggplot2::element_text(size = 9.2, color = "#4b5563"),
    axis.title.x = ggplot2::element_text(size = 10.5, color = "#374151", margin = ggplot2::margin(t = 8)),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 9.5),
    legend.text = ggplot2::element_text(size = 9),
    plot.caption = ggplot2::element_text(size = 8.5, color = "#4b5563", hjust = 0, margin = ggplot2::margin(t = 10)),
    plot.margin = ggplot2::margin(14, 22, 14, 18)
  )

png_path <- file.path(output_dir, "figure_6_class1_analysis_candidate_direction_lollipop.png")
pdf_path <- file.path(output_dir, "figure_6_class1_analysis_candidate_direction_lollipop.pdf")
svg_path <- file.path(output_dir, "figure_6_class1_analysis_candidate_direction_lollipop.svg")

ggplot2::ggsave(png_path, figure, width = 10.8, height = 7.2, dpi = 360, bg = "white")
ggplot2::ggsave(pdf_path, figure, width = 10.8, height = 7.2, bg = "white")
ggplot2::ggsave(svg_path, figure, width = 10.8, height = 7.2, bg = "white", device = svglite::svglite)

message("Main Figure 6 Class 1 candidate-direction lollipop chart created:")
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Inputs:")
message(" - ", results_path)
message(" - ", windows_path)
message("Summary:")
message(" - Nominal candidates: ", fmt_int(n_nominal))
message(" - Increasing: ", fmt_int(n_increasing))
message(" - Decreasing: ", fmt_int(n_decreasing))
message(" - Global Class 1 BH-FDR significant tests: ", fmt_int(n_fdr))
message(" - Best q_global_class1_bh: ", fmt_q(best_q))
