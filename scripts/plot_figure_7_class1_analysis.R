# ============================================
# plot_figure_7_class1_analysis.R
# ============================================
# Main Figure 7: candidate trajectory small multiples.
# This figure visualizes subject-level trajectories for top nominal Class 1
# MAG-window candidates. It is diagnostic and hypothesis-generating only.

required_packages <- c("readr", "dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
cleaned_data_path <- file.path(repo_root, "data", "processed", "cleaned_data.rds")
results_path <- file.path(
  repo_root,
  "results",
  "class1_window_lmertest",
  "tables",
  "class1_window_lmertest_results.csv"
)
output_dir <- file.path(repo_root, "results", "class1_window_lmertest", "figures")

if (!file.exists(cleaned_data_path)) {
  stop("Missing processed data: ", cleaned_data_path, call. = FALSE)
}
if (!file.exists(results_path)) {
  stop("Missing Class 1 results table: ", results_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cleaned_data <- readRDS(cleaned_data_path)
class1_results <- readr::read_csv(results_path, show_col_types = FALSE)

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
}
fmt_q <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.3f", x))
}
wrap_label <- function(x, width = 34) {
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

compact_window_names <- c(
  initial_ascent_to_high_altitude = "Initial ascent",
  ascent_and_intermediate_adjustment = "Ascent adjustment",
  high_climb_6500m_exposure = "6500 m climb",
  pre_summit_staging_to_post_summit = "Pre/post summit",
  summit_push_and_descent_recovery_transition = "Summit + descent",
  post_expedition_low_altitude_recovery = "Recovery"
)

top_n <- 16

top_candidates <- class1_results |>
  dplyr::filter(
    primary_analysis %in% TRUE,
    nominal_p_lt_0_05 %in% TRUE,
    !is.na(raw_p_value),
    !is.na(q_global_class1_bh)
  ) |>
  dplyr::arrange(q_global_class1_bh, raw_p_value, dplyr::desc(abs(estimate_window_order))) |>
  dplyr::slice_head(n = top_n) |>
  dplyr::mutate(
    candidate_id = dplyr::row_number(),
    window_short = dplyr::coalesce(unname(compact_window_names[window_label]), window_label),
    panel_label = paste0(
      "Rank ", candidate_id, " | ", window_short, " | ", forward_timepoint_path,
      "\n", wrap_label(MAG, 40),
      "\nraw p=", fmt_p(raw_p_value),
      "  q=", fmt_q(q_global_class1_bh),
      "  ", direction
    )
  )

if (nrow(top_candidates) == 0) {
  stop("No nominal Class 1 candidates found for trajectory plotting.", call. = FALSE)
}

# Build candidate-specific trajectory data because each candidate has its own
# MAG and three-timepoint window. The response scale matches the Class 1 LMM.
trajectory_list <- lapply(seq_len(nrow(top_candidates)), function(i) {
  candidate <- top_candidates[i, , drop = FALSE]
  candidate_timepoints <- strsplit(candidate$timepoints, ",", fixed = TRUE)[[1]]
  candidate_timepoints <- trimws(candidate_timepoints)
  order_lookup <- setNames(seq_along(candidate_timepoints) - 1L, candidate_timepoints)

  cleaned_data |>
    dplyr::filter(MAG == candidate$MAG, Timepoint %in% candidate_timepoints) |>
    dplyr::mutate(
      candidate_id = candidate$candidate_id,
      panel_label = candidate$panel_label,
      window_label = candidate$window_label,
      forward_timepoint_path = candidate$forward_timepoint_path,
      raw_p_value = candidate$raw_p_value,
      q_global_class1_bh = candidate$q_global_class1_bh,
      direction = candidate$direction,
      pseudocount = candidate$pseudocount,
      window_order = as.integer(order_lookup[as.character(Timepoint)]),
      response_log10 = log10(Y + candidate$pseudocount),
      observed_nonzero = Y > 0
    )
})

trajectory_data <- dplyr::bind_rows(trajectory_list) |>
  dplyr::mutate(
    panel_label = factor(panel_label, levels = top_candidates$panel_label),
    Subject = factor(Subject)
  )

mean_data <- trajectory_data |>
  dplyr::group_by(panel_label, candidate_id, window_order) |>
  dplyr::summarise(
    mean_response = mean(response_log10, na.rm = TRUE),
    n_samples = dplyr::n(),
    n_nonzero = sum(observed_nonzero, na.rm = TRUE),
    .groups = "drop"
  )

panel_support <- trajectory_data |>
  dplyr::group_by(panel_label, candidate_id) |>
  dplyr::summarise(
    n_subjects = dplyr::n_distinct(Subject),
    n_samples = dplyr::n(),
    nonzero_fraction = mean(observed_nonzero, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    label = paste0(
      "subjects ", n_subjects,
      " | samples ", n_samples,
      " | nonzero ", sprintf("%.0f%%", 100 * nonzero_fraction)
    ),
    x = 0.03,
    y = Inf
  )

direction_palette <- c(
  decreasing = "#D55E00",
  increasing = "#0072B2",
  flat_or_zero = "#4B5563"
)

figure <- ggplot2::ggplot(
  trajectory_data,
  ggplot2::aes(x = window_order, y = response_log10, group = Subject)
) +
  ggplot2::geom_line(color = "#9ca3af", linewidth = 0.38, alpha = 0.55) +
  ggplot2::geom_point(
    ggplot2::aes(shape = observed_nonzero),
    color = "#6b7280",
    fill = "white",
    size = 1.35,
    stroke = 0.35,
    alpha = 0.72
  ) +
  ggplot2::geom_line(
    data = mean_data,
    ggplot2::aes(x = window_order, y = mean_response, group = panel_label),
    inherit.aes = FALSE,
    color = "#111827",
    linewidth = 1.05
  ) +
  ggplot2::geom_point(
    data = mean_data,
    ggplot2::aes(x = window_order, y = mean_response),
    inherit.aes = FALSE,
    color = "#111827",
    fill = "white",
    shape = 21,
    size = 2.25,
    stroke = 0.55
  ) +
  ggplot2::geom_text(
    data = panel_support,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.25,
    size = 2.35,
    color = "#4b5563"
  ) +
  ggplot2::facet_wrap(~ panel_label, ncol = 4, scales = "free_y") +
  ggplot2::scale_x_continuous(
    breaks = c(0, 1, 2),
    labels = c("0", "1", "2"),
    limits = c(-0.08, 2.08),
    expand = c(0.02, 0.02)
  ) +
  ggplot2::scale_shape_manual(
    values = c(`TRUE` = 21, `FALSE` = 1),
    labels = c(`TRUE` = "nonzero", `FALSE` = "zero"),
    name = "Observation"
  ) +
  ggplot2::labs(
    title = "Main Figure 7 | Subject trajectories for top nominal Class 1 candidates",
    subtitle = paste0(
      "Top ", nrow(top_candidates),
      " nominal MAG-window candidates ranked by global q-value then raw p-value; thick line shows window-order mean."
    ),
    x = "Window order within predefined three-timepoint window",
    y = "log10(Y + pseudocount)",
    caption = paste0(
      "Thin grey lines are individual subjects. These panels are diagnostic: no candidate shown here passed global Class 1 BH-FDR correction. ",
      "Free y-axis scaling is used so within-window trajectory shape can be inspected."
    )
  ) +
  ggplot2::theme_minimal(base_size = 10.5) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#e5e7eb", linewidth = 0.28),
    strip.text = ggplot2::element_text(face = "bold", size = 6.8, lineheight = 0.88, color = "#111827"),
    strip.background = ggplot2::element_rect(fill = "#f8fafc", color = "#d1d5db", linewidth = 0.35),
    plot.title = ggplot2::element_text(face = "bold", size = 15.5, color = "#111827"),
    plot.subtitle = ggplot2::element_text(size = 10.3, color = "#4b5563", margin = ggplot2::margin(b = 8)),
    axis.title = ggplot2::element_text(size = 10.3, color = "#374151"),
    axis.text = ggplot2::element_text(size = 8.4, color = "#4b5563"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 9),
    legend.text = ggplot2::element_text(size = 8.5),
    plot.caption = ggplot2::element_text(size = 8.3, color = "#4b5563", hjust = 0, margin = ggplot2::margin(t = 10)),
    panel.spacing = grid::unit(0.72, "lines"),
    plot.margin = ggplot2::margin(14, 16, 14, 16)
  )

png_path <- file.path(output_dir, "figure_7_class1_analysis_candidate_trajectories.png")
pdf_path <- file.path(output_dir, "figure_7_class1_analysis_candidate_trajectories.pdf")
svg_path <- file.path(output_dir, "figure_7_class1_analysis_candidate_trajectories.svg")

ggplot2::ggsave(png_path, figure, width = 13.2, height = 10.4, dpi = 360, bg = "white")
ggplot2::ggsave(pdf_path, figure, width = 13.2, height = 10.4, bg = "white")
ggplot2::ggsave(svg_path, figure, width = 13.2, height = 10.4, bg = "white", device = svglite::svglite)

message("Main Figure 7 Class 1 candidate trajectory small multiples created:")
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Inputs:")
message(" - ", cleaned_data_path)
message(" - ", results_path)
message("Summary:")
message(" - Candidate panels: ", nrow(top_candidates))
message(" - Unique MAGs shown: ", dplyr::n_distinct(top_candidates$MAG))
message(" - Windows represented: ", dplyr::n_distinct(top_candidates$window_label))
message(" - No plotted candidate is FDR-significant; panels are diagnostic only.")
