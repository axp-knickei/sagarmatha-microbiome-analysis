# ============================================
# plot_main_figure_1_class1_analysis.R
# ============================================
# Main Figure 1: expedition-stage design schematic for the Class 1 analysis.
# This is a metadata/figure-context schematic only; it does not run MAG models.

required_packages <- c("readr", "dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
timepoint_path <- file.path(
  repo_root,
  "results",
  "metadata_expedition_windows",
  "timepoint_altitude_date_summary.csv"
)
window_path <- file.path(
  repo_root,
  "results",
  "metadata_expedition_windows",
  "proposed_biological_windows.csv"
)
output_dir <- file.path(repo_root, "results", "class1_window_lmertest", "figures")

if (!file.exists(timepoint_path)) {
  stop("Missing timepoint summary: ", timepoint_path, call. = FALSE)
}
if (!file.exists(window_path)) {
  stop("Missing Class 1 window definitions: ", window_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

timepoints <- readr::read_csv(timepoint_path, show_col_types = FALSE) |>
  dplyr::mutate(
    x = as.numeric(timepoint_numeric_order),
    sampling_altitude_m = as.numeric(median_altitude),
    sampling_y = 4.2 + sampling_altitude_m / 2000,
    sampling_track_y = 3.15,
    date_label = format(as.Date(median_date), "%b %d")
  )

windows <- readr::read_csv(window_path, show_col_types = FALSE) |>
  dplyr::mutate(
    start_tp = sub(",.*$", "", timepoints),
    end_tp = sub("^.*,", "", timepoints),
    start_x = as.numeric(sub("^T", "", start_tp)),
    end_x = as.numeric(sub("^T", "", end_tp)),
    window_index = dplyr::row_number(),
    window_y = 2.38 - ((window_index - 1) %% 3) * 0.42,
    compact_label = dplyr::case_when(
      window_label == "initial_ascent_to_high_altitude" ~ "Initial ascent",
      window_label == "ascent_and_intermediate_adjustment" ~ "Ascent adjustment",
      window_label == "high_climb_6500m_exposure" ~ "6500 m high climb",
      window_label == "pre_summit_staging_to_post_summit" ~ "Pre/post summit",
      window_label == "summit_push_and_descent_recovery_transition" ~ "Summit + descent",
      window_label == "post_expedition_low_altitude_recovery" ~ "Recovery",
      TRUE ~ window_label
    ),
    window_label_text = paste0(compact_label, "\n", forward_timepoint_path)
  )

# Figure-informed exposure trajectory. Sampling altitudes are metadata values;
# the summit-push peak between T7 and T8 is prior context from the expedition
# figure and is not directly encoded as a sampling altitude.
exposure_trajectory <- dplyr::bind_rows(
  timepoints |>
    dplyr::transmute(
      x,
      altitude_m = sampling_altitude_m,
      point_type = "sampling timepoint"
    ),
  data.frame(
    x = 7.55,
    altitude_m = 8848.86,
    point_type = "figure-informed summit-push exposure",
    stringsAsFactors = FALSE
  )
) |>
  dplyr::arrange(x) |>
  dplyr::mutate(exposure_y = 4.2 + altitude_m / 2000)

summit_peak <- exposure_trajectory |>
  dplyr::filter(point_type == "figure-informed summit-push exposure")

annotation_points <- data.frame(
  x = c(5, 7.55, 8),
  y = c(4.2 + 6500 / 2000, 4.2 + 8848.86 / 2000, 4.2 + 5200 / 2000),
  label = c(
    "T5: 6500 m high-climb\nsampling point, not summit",
    "T7 -> T8: figure-informed\nsummit-push exposure to 8848.86 m",
    "T8: post-summit /\npost-extreme-altitude sample"
  ),
  hjust = c(1.05, 0.5, -0.05),
  vjust = c(-0.35, -0.45, 1.15),
  stringsAsFactors = FALSE
)

evidence_palette <- c(
  metadata_only = "#6b7280",
  figure_informed = "#8b5cf6",
  metadata_plus_figure = "#2563eb",
  figure_informed_critical_interval = "#dc2626",
  figure_informed_plus_sampling_descent = "#ea580c"
)

track_labels <- data.frame(
  x = 0.63,
  y = c(8.72, 3.66, 2.6),
  label = c(
    "Track A: exposure altitude trajectory",
    "Track B: sampling timepoints",
    "Track C: predefined Class 1 windows"
  ),
  stringsAsFactors = FALSE
)

altitude_breaks_m <- c(0, 3650, 5200, 6500, 8848.86)
altitude_breaks_y <- 4.2 + altitude_breaks_m / 2000

figure <- ggplot2::ggplot() +
  ggplot2::annotate(
    "rect",
    xmin = 7,
    xmax = 8,
    ymin = 4.05,
    ymax = 8.85,
    fill = "#fee2e2",
    alpha = 0.55
  ) +
  ggplot2::annotate(
    "rect",
    xmin = 0.6,
    xmax = 10.4,
    ymin = 2.78,
    ymax = 3.72,
    fill = "#f3f4f6",
    alpha = 0.9
  ) +
  ggplot2::annotate(
    "rect",
    xmin = 0.6,
    xmax = 10.4,
    ymin = 0.85,
    ymax = 2.62,
    fill = "#eff6ff",
    alpha = 0.75
  ) +
  ggplot2::geom_line(
    data = exposure_trajectory,
    ggplot2::aes(x = x, y = exposure_y),
    linewidth = 1.05,
    linetype = "22",
    color = "#991b1b"
  ) +
  ggplot2::geom_line(
    data = timepoints,
    ggplot2::aes(x = x, y = sampling_y),
    linewidth = 0.75,
    color = "#111827",
    alpha = 0.72
  ) +
  ggplot2::geom_point(
    data = timepoints,
    ggplot2::aes(x = x, y = sampling_y),
    size = 3.2,
    shape = 21,
    fill = "white",
    color = "#111827",
    stroke = 0.9
  ) +
  ggplot2::geom_point(
    data = summit_peak,
    ggplot2::aes(x = x, y = exposure_y),
    size = 3.6,
    shape = 24,
    fill = "#dc2626",
    color = "#7f1d1d",
    stroke = 0.9
  ) +
  ggplot2::geom_segment(
    data = timepoints,
    ggplot2::aes(x = x, xend = x, y = sampling_track_y - 0.12, yend = sampling_y - 0.12),
    linewidth = 0.25,
    color = "#9ca3af",
    linetype = "dotted"
  ) +
  ggplot2::geom_point(
    data = timepoints,
    ggplot2::aes(x = x, y = sampling_track_y),
    size = 3.2,
    shape = 21,
    fill = "#ffffff",
    color = "#111827",
    stroke = 0.9
  ) +
  ggplot2::geom_text(
    data = timepoints,
    ggplot2::aes(x = x, y = sampling_track_y - 0.23, label = Timepoint),
    size = 3.4,
    fontface = "bold",
    color = "#111827"
  ) +
  ggplot2::geom_text(
    data = timepoints,
    ggplot2::aes(x = x, y = sampling_track_y + 0.22, label = paste0(date_label, "\n", sampling_altitude_m, " m")),
    size = 2.45,
    lineheight = 0.9,
    color = "#374151"
  ) +
  ggplot2::geom_segment(
    data = windows,
    ggplot2::aes(
      x = start_x,
      xend = end_x,
      y = window_y,
      yend = window_y,
      color = evidence_type
    ),
    linewidth = 2.1,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    data = windows,
    ggplot2::aes(x = start_x, y = window_y, color = evidence_type),
    size = 2
  ) +
  ggplot2::geom_point(
    data = windows,
    ggplot2::aes(x = end_x, y = window_y, color = evidence_type),
    size = 2
  ) +
  ggplot2::geom_text(
    data = windows,
    ggplot2::aes(x = (start_x + end_x) / 2, y = window_y - 0.16, label = window_label_text),
    size = 2.35,
    lineheight = 0.86,
    color = "#111827"
  ) +
  ggplot2::geom_text(
    data = annotation_points,
    ggplot2::aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    size = 3.05,
    lineheight = 0.92,
    color = "#7f1d1d",
    fontface = "bold"
  ) +
  ggplot2::geom_text(
    data = track_labels,
    ggplot2::aes(x = x, y = y, label = label),
    hjust = 0,
    size = 3.25,
    fontface = "bold",
    color = "#111827"
  ) +
  ggplot2::annotate(
    "text",
    x = 1,
    y = 8.25,
    label = "Dashed red line = figure-informed exposure trajectory\nBlack points/line = metadata sampling altitude",
    hjust = 0,
    size = 3.05,
    lineheight = 0.95,
    color = "#374151"
  ) +
  ggplot2::scale_color_manual(values = evidence_palette, drop = FALSE, name = "Window evidence type") +
  ggplot2::scale_x_continuous(
    breaks = timepoints$x,
    labels = timepoints$Timepoint,
    limits = c(0.55, 10.45),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_y_continuous(
    breaks = c(altitude_breaks_y, 3.15, 2.38),
    labels = c("0 m", "3650 m", "5200 m", "6500 m", "8848.86 m", "sampling\ntrack", "Class 1\nwindows"),
    limits = c(0.65, 9.05),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = "Main Figure 1. Expedition-stage design schematic for Class 1 analysis",
    subtitle = "Sampling altitude is not the full exposure altitude trajectory; predefined three-timepoint windows target stage-local longitudinal responses.",
    x = NULL,
    y = NULL,
    caption = paste(
      "Schematic uses metadata sampling altitude/date summaries plus figure-informed exposure context.\n",
      "The T7 -> T8 summit-push peak is prior expedition-figure context, not a sampled altitude value.",
      "Class 1 windows were defined before evaluating MAG p-values."
    )
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 15, color = "#111827"),
    plot.subtitle = ggplot2::element_text(size = 10.5, color = "#374151"),
    plot.caption = ggplot2::element_text(size = 8.2, color = "#4b5563", hjust = 0),
    axis.text.y = ggplot2::element_text(size = 8.8, color = "#374151"),
    axis.text.x = ggplot2::element_text(size = 9.2, face = "bold", color = "#111827"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(linewidth = 0.25, color = "#e5e7eb"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold"),
    legend.text = ggplot2::element_text(size = 8.2),
    plot.margin = ggplot2::margin(12, 18, 12, 12)
  )

png_path <- file.path(output_dir, "main_figure_1_class1_analysis.png")
pdf_path <- file.path(output_dir, "main_figure_1_class1_analysis.pdf")
svg_path <- file.path(output_dir, "main_figure_1_class1_analysis.svg")

ggplot2::ggsave(png_path, figure, width = 13.5, height = 8.2, dpi = 320, bg = "white")
ggplot2::ggsave(pdf_path, figure, width = 13.5, height = 8.2, bg = "white")
ggplot2::ggsave(svg_path, figure, width = 13.5, height = 8.2, bg = "white", device = svglite::svglite)

message("Main Figure 1 Class 1 schematic created:")
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Inputs:")
message(" - ", timepoint_path)
message(" - ", window_path)
