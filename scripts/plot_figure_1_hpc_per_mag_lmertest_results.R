# ============================================
# plot_figure_1_hpc_per_mag_lmertest_results.R
# ============================================
# Figure 1: longitudinal sampling design and data structure for the full
# T1-T10 per-MAG lmerTest LMM workflow.
#
# The figure is metadata-only: it visualizes repeated sampling, missingness,
# model time, and expedition altitude context. No MAG abundance transformation
# or statistical testing is performed here.

required_packages <- c("ggplot2", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
input_path <- file.path(repo_root, "data", "processed", "cleaned_data.rds")
output_dir <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "figures")
output_stem <- "figure_1_hpc_per_mag_lmertest_results_20260708_024546"

if (!file.exists(input_path)) {
  stop("Missing processed input data: ", input_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- readRDS(input_path)
required_columns <- c(
  "sampleID", "Subject", "Timepoint", "days_from_baseline_subject",
  "date", "altitude", "gender", "age", "MAG", "Y"
)
missing_columns <- setdiff(required_columns, names(df))
if (length(missing_columns) > 0) {
  stop("Missing required columns in cleaned_data.rds: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

timepoint_levels <- paste0("T", 1:10)

# One row per biological sample. The processed table is long by MAG, so sample
# metadata must be deduplicated before drawing sampling-design summaries.
sample_metadata <- df |>
  dplyr::select(
    sampleID, Subject, Timepoint, date, days_from_baseline_subject,
    altitude, gender, age
  ) |>
  dplyr::distinct() |>
  dplyr::mutate(
    Subject = factor(Subject, levels = sort(unique(Subject))),
    Timepoint = factor(Timepoint, levels = timepoint_levels),
    date_label = ifelse(
      grepl("^\\d{8}$", as.character(date)),
      format(as.Date(as.character(date), format = "%Y%m%d"), "%b %d"),
      as.character(date)
    )
  )

subject_levels <- levels(sample_metadata$Subject)

presence_grid <- expand.grid(
  Subject = subject_levels,
  Timepoint = timepoint_levels,
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(
    sample_metadata |>
      dplyr::count(Subject = as.character(Subject), Timepoint = as.character(Timepoint), name = "sample_count"),
    by = c("Subject", "Timepoint")
  ) |>
  dplyr::mutate(
    sample_count = dplyr::coalesce(sample_count, 0L),
    sample_available = sample_count > 0,
    Subject = factor(Subject, levels = rev(subject_levels)),
    Timepoint = factor(Timepoint, levels = timepoint_levels)
  )

timepoint_summary <- sample_metadata |>
  dplyr::group_by(Timepoint) |>
  dplyr::summarise(
    n_samples = dplyr::n_distinct(sampleID),
    median_days = stats::median(days_from_baseline_subject, na.rm = TRUE),
    min_days = min(days_from_baseline_subject, na.rm = TRUE),
    max_days = max(days_from_baseline_subject, na.rm = TRUE),
    altitude = stats::median(altitude, na.rm = TRUE),
    date_label = paste(unique(date_label), collapse = ", "),
    .groups = "drop"
  ) |>
  dplyr::mutate(Timepoint = factor(Timepoint, levels = timepoint_levels))

subject_summary <- sample_metadata |>
  dplyr::group_by(Subject) |>
  dplyr::summarise(
    n_samples = dplyr::n_distinct(sampleID),
    first_day = min(days_from_baseline_subject, na.rm = TRUE),
    last_day = max(days_from_baseline_subject, na.rm = TRUE),
    .groups = "drop"
  )

total_samples <- dplyr::n_distinct(sample_metadata$sampleID)
total_subjects <- dplyr::n_distinct(sample_metadata$Subject)
total_mags <- dplyr::n_distinct(df$MAG)
observed_slots <- sum(presence_grid$sample_available)
possible_slots <- length(subject_levels) * length(timepoint_levels)
complete_subjects <- sum(subject_summary$n_samples == length(timepoint_levels))
median_samples_per_subject <- stats::median(subject_summary$n_samples)
max_subject_day <- max(sample_metadata$days_from_baseline_subject, na.rm = TRUE)
timeline_breaks <- seq(0, ceiling(max_subject_day / 10) * 10, by = 10)

fmt_int <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
fmt_pct <- function(x) sprintf("%.0f%%", 100 * x)

# Colorblind-friendly palette derived from Okabe-Ito plus restrained neutrals.
pal <- list(
  ink = "#111827",
  text = "#374151",
  muted = "#6B7280",
  hair = "#D1D5DB",
  faint = "#F3F4F6",
  panel = "#FAFAF8",
  blue = "#0072B2",
  sky = "#DCEEF8",
  orange = "#D55E00",
  vermillion = "#CC3311",
  green = "#009E73",
  purple = "#785EF0"
)

base_theme <- ggplot2::theme_minimal(base_size = 9.8, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11.2, color = pal$ink),
    plot.subtitle = ggplot2::element_text(size = 8.4, color = pal$muted, margin = ggplot2::margin(b = 6)),
    plot.caption = ggplot2::element_text(size = 7.2, color = pal$muted, hjust = 0),
    axis.title = ggplot2::element_text(size = 8.5, color = pal$text),
    axis.text = ggplot2::element_text(size = 7.6, color = pal$text),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.25),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 8),
    legend.text = ggplot2::element_text(size = 7.6),
    plot.margin = ggplot2::margin(5, 7, 5, 7)
  )

panel_a <- ggplot2::ggplot(
  presence_grid,
  ggplot2::aes(x = Timepoint, y = Subject, fill = sample_available)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.8, width = 0.92, height = 0.86) +
  ggplot2::geom_text(
    data = timepoint_summary,
    ggplot2::aes(x = Timepoint, y = length(subject_levels) + 1.05, label = paste0(altitude, " m")),
    inherit.aes = FALSE,
    size = 2.25,
    color = pal$muted
  ) +
  ggplot2::scale_fill_manual(
    values = c("TRUE" = pal$blue, "FALSE" = "#ECEFF1"),
    labels = c("TRUE" = "Sample available", "FALSE" = "Missing"),
    name = NULL
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "A  Structured but incomplete T1-T10 sampling",
    subtitle = paste0(
      fmt_int(observed_slots), " of ", fmt_int(possible_slots),
      " possible subject-timepoint slots observed; altitude shown above columns"
    ),
    x = NULL,
    y = "Subject"
  ) +
  base_theme +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(face = "bold", size = 8.3),
    axis.text.y = ggplot2::element_text(size = 7.5),
    legend.position = "none",
    plot.margin = ggplot2::margin(11, 8, 6, 7)
  )

panel_b <- ggplot2::ggplot(
  sample_metadata,
  ggplot2::aes(x = days_from_baseline_subject, y = Subject, group = Subject)
) +
  ggplot2::geom_line(color = pal$hair, linewidth = 0.45) +
  ggplot2::geom_point(
    ggplot2::aes(fill = Timepoint),
    shape = 21,
    size = 2.5,
    color = "white",
    stroke = 0.35
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      T1 = "#4D4D4D", T2 = "#0072B2", T3 = "#56B4E9", T4 = "#009E73",
      T5 = "#E69F00", T6 = "#F0E442", T7 = "#D55E00", T8 = "#CC79A7",
      T9 = "#785EF0", T10 = "#000000"
    ),
    name = "Timepoint"
  ) +
  ggplot2::scale_x_continuous(breaks = timeline_breaks, limits = c(-2, max_subject_day + 4)) +
  ggplot2::labs(
    title = "B  Repeated measures along subject-specific study time",
    subtitle = "Each line follows one participant; points are observed samples used by the LMM design",
    x = "Days from subject baseline",
    y = "Subject"
  ) +
  base_theme +
  ggplot2::theme(
    legend.key.width = grid::unit(0.25, "in"),
    legend.key.height = grid::unit(0.15, "in"),
    panel.grid.major.y = ggplot2::element_blank()
  )

panel_c <- ggplot2::ggplot(timepoint_summary, ggplot2::aes(x = Timepoint)) +
  ggplot2::geom_col(
    ggplot2::aes(y = n_samples * 500),
    fill = pal$sky,
    color = "white",
    width = 0.72
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = altitude, group = 1),
    color = pal$orange,
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = altitude),
    color = pal$orange,
    fill = "white",
    shape = 21,
    size = 2.6,
    stroke = 0.7
  ) +
  ggplot2::geom_text(
    ggplot2::aes(y = n_samples * 500 + 260, label = n_samples),
    size = 2.45,
    fontface = "bold",
    color = pal$blue
  ) +
  ggplot2::scale_y_continuous(
    name = "Altitude (m)",
    limits = c(0, 7000),
    breaks = c(0, 3650, 5200, 6500),
    labels = c("0", "3,650", "5,200", "6,500"),
    sec.axis = ggplot2::sec_axis(~ . / 500, name = "Samples per timepoint")
  ) +
  ggplot2::labs(
    title = "C  Expedition altitude context and sample yield",
    subtitle = "Orange profile shows sampling altitude; blue bars show sample count at each stage",
    x = NULL
  ) +
  base_theme +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(face = "bold", size = 8.3),
    axis.title.y.left = ggplot2::element_text(color = pal$orange),
    axis.title.y.right = ggplot2::element_text(color = pal$blue),
    panel.grid.major.x = ggplot2::element_blank(),
    legend.position = "none"
  )

summary_text <- paste0(
  "Repeated-measures design: ", fmt_int(total_samples), " samples from ",
  fmt_int(total_subjects), " participants across T1-T10, profiled over ",
  fmt_int(total_mags), " MAGs. Median samples per participant = ",
  median_samples_per_subject, "; complete T1-T10 participants = ",
  complete_subjects, ". Mixed models use days from subject baseline and a subject random intercept."
)

draw_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(gp = grid::gpar(fill = "white", col = NA)))
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  layout <- grid::grid.layout(
    nrow = 4,
    ncol = 2,
    heights = grid::unit(c(0.42, 3.2, 3.1, 0.36), c("in", "null", "null", "in")),
    widths = grid::unit(c(1.05, 1.0), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))

  grid::grid.text(
    "Figure 1 | Longitudinal sampling design and data structure",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.88, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 15.5, fontface = "bold", col = pal$ink),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )
  grid::grid.text(
    "Incomplete but structured repeated sampling across expedition time motivates subject-aware mixed modelling.",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.24, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 9.2, col = pal$text),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )

  print(panel_a, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panel_b, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(panel_c, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1:2))

  grid::grid.text(
    summary_text,
    x = grid::unit(0, "npc"),
    y = grid::unit(0.55, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 7.8, col = pal$muted),
    vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1:2)
  )

  grid::popViewport(2)
}

png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

png(filename = png_path, width = 12.6, height = 9.0, units = "in", res = 420, bg = "white")
draw_figure()
dev.off()

pdf(file = pdf_path, width = 12.6, height = 9.0, bg = "white", useDingbats = FALSE)
draw_figure()
dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = svg_path, width = 12.6, height = 9.0, bg = "white")
  draw_figure()
  dev.off()
}

metadata_summary_path <- file.path(output_dir, paste0(output_stem, "_metadata_summary.csv"))
utils::write.csv(
  data.frame(
    metric = c(
      "total_samples",
      "total_subjects",
      "total_mags",
      "observed_subject_timepoint_slots",
      "possible_subject_timepoint_slots",
      "slot_completion_fraction",
      "median_samples_per_subject",
      "complete_T1_T10_subjects"
    ),
    value = c(
      total_samples,
      total_subjects,
      total_mags,
      observed_slots,
      possible_slots,
      observed_slots / possible_slots,
      median_samples_per_subject,
      complete_subjects
    )
  ),
  metadata_summary_path,
  row.names = FALSE
)

message("Figure 1 sampling-design outputs created:")
message(" - ", png_path)
message(" - ", pdf_path)
if (file.exists(svg_path)) message(" - ", svg_path)
message(" - ", metadata_summary_path)
