# ============================================
# plot_figure_2_hpc_per_mag_lmertest_results.R
# ============================================
# Figure 2: MAG abundance distribution, sparsity, and transformation
# justification for the full T1-T10 per-MAG lmerTest LMM workflow.
#
# This figure explains why the primary LMM uses log10(Y + pseudocount):
# the MAG abundance table is zero-heavy, sparse across MAGs, and strongly
# right-skewed on the raw relative-abundance scale.

required_packages <- c("ggplot2", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
input_path <- file.path(repo_root, "data", "processed", "cleaned_data.rds")
results_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "tables",
  "per_mag_lmertest_hpc_results.csv"
)
output_dir <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "figures")
output_stem <- "figure_2_hpc_per_mag_lmertest_results_20260708_024546"

if (!file.exists(input_path)) {
  stop("Missing processed input data: ", input_path, call. = FALSE)
}
if (!file.exists(results_path)) {
  stop("Missing per-MAG LMM result table: ", results_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- readRDS(input_path)
required_columns <- c("sampleID", "Subject", "Timepoint", "MAG", "Y")
missing_columns <- setdiff(required_columns, names(df))
if (length(missing_columns) > 0) {
  stop("Missing required columns in cleaned_data.rds: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("MAG", "status", "primary_analysis") %in% names(results))) {
  stop("The per-MAG result table must contain MAG, status, and primary_analysis columns.", call. = FALSE)
}

smallest_nonzero_y <- min(df$Y[df$Y > 0], na.rm = TRUE)
if (!is.finite(smallest_nonzero_y)) {
  stop("Cannot compute pseudocount because no positive Y values were found.", call. = FALSE)
}
pseudocount <- smallest_nonzero_y / 2

df <- df |>
  dplyr::mutate(
    raw_abundance = Y,
    nonzero = Y > 0,
    log10_abundance = log10(Y + pseudocount)
  )

mag_summary <- df |>
  dplyr::group_by(MAG) |>
  dplyr::summarise(
    n_observations = dplyr::n(),
    zero_fraction = mean(Y == 0, na.rm = TRUE),
    nonzero_count = sum(Y > 0, na.rm = TRUE),
    mean_abundance = mean(Y, na.rm = TRUE),
    median_abundance = stats::median(Y, na.rm = TRUE),
    response_variance = stats::var(log10(Y + pseudocount), na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    results |>
      dplyr::select(MAG, status, primary_analysis, time_p_value, time_p_adj_bh_primary),
    by = "MAG"
  ) |>
  dplyr::mutate(
    model_status = dplyr::case_when(
      status == "lmm_ok" ~ "Modelable",
      status == "too_sparse" ~ "Too sparse",
      status == "low_variance" ~ "Low variance",
      status == "lmm_singular" ~ "Singular LMM",
      TRUE ~ "Other / failed"
    ),
    model_status = factor(
      model_status,
      levels = c("Modelable", "Singular LMM", "Too sparse", "Low variance", "Other / failed")
    )
  )

nonzero_df <- df |>
  dplyr::filter(Y > 0)

zero_fraction_threshold <- 1 - 10 / dplyr::n_distinct(df$sampleID)
subject_threshold_label <- "Primary LMM screen also requires >=3 subjects with non-zero abundance"

total_values <- nrow(df)
zero_values <- sum(df$Y == 0, na.rm = TRUE)
nonzero_values <- sum(df$Y > 0, na.rm = TRUE)
zero_fraction_global <- zero_values / total_values
total_mags <- dplyr::n_distinct(df$MAG)
modelable_mags <- sum(mag_summary$model_status == "Modelable", na.rm = TRUE)
too_sparse_mags <- sum(mag_summary$model_status == "Too sparse", na.rm = TRUE)
singular_mags <- sum(mag_summary$model_status == "Singular LMM", na.rm = TRUE)
median_zero_fraction <- stats::median(mag_summary$zero_fraction, na.rm = TRUE)
median_response_variance <- stats::median(mag_summary$response_variance, na.rm = TRUE)

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
fmt_sci <- function(x) formatC(x, format = "e", digits = 2)

pal <- list(
  ink = "#111827",
  text = "#374151",
  muted = "#6B7280",
  hair = "#D1D5DB",
  faint = "#F3F4F6",
  blue = "#0072B2",
  sky = "#DCEEF8",
  orange = "#D55E00",
  green = "#009E73",
  vermillion = "#CC3311",
  purple = "#785EF0",
  grey = "#9CA3AF"
)

status_palette <- c(
  "Modelable" = pal$blue,
  "Singular LMM" = pal$orange,
  "Too sparse" = pal$vermillion,
  "Low variance" = pal$purple,
  "Other / failed" = pal$grey
)

base_theme <- ggplot2::theme_minimal(base_size = 9.6, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11.1, color = pal$ink),
    plot.subtitle = ggplot2::element_text(size = 8.2, color = pal$muted, margin = ggplot2::margin(b = 6)),
    plot.caption = ggplot2::element_text(size = 7.0, color = pal$muted, hjust = 0),
    axis.title = ggplot2::element_text(size = 8.4, color = pal$text),
    axis.text = ggplot2::element_text(size = 7.6, color = pal$text),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.25),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 8),
    legend.text = ggplot2::element_text(size = 7.5),
    plot.margin = ggplot2::margin(6, 7, 6, 7)
  )

panel_a <- ggplot2::ggplot(nonzero_df, ggplot2::aes(x = raw_abundance)) +
  ggplot2::geom_histogram(bins = 55, fill = pal$blue, color = "white", linewidth = 0.18) +
  ggplot2::geom_vline(xintercept = smallest_nonzero_y, color = pal$orange, linewidth = 0.55, linetype = "longdash") +
  ggplot2::scale_x_log10(labels = function(x) formatC(x, format = "e", digits = 1)) +
  ggplot2::annotate(
    "label",
    x = stats::quantile(nonzero_df$raw_abundance, 0.82, na.rm = TRUE),
    y = Inf,
    label = paste0(
      fmt_pct(zero_fraction_global), " of all MAG values are zero\n",
      fmt_int(nonzero_values), " non-zero values shown"
    ),
    vjust = 1.12,
    hjust = 0.5,
    size = 2.7,
    fill = "white",
    color = pal$text
  ) +
  ggplot2::labs(
    title = "A  Raw non-zero abundance is strongly right-skewed",
    subtitle = "Zeros are excluded from the log-scaled x-axis but counted in the annotation",
    x = "Raw relative abundance Y, non-zero values only (log10 axis)",
    y = "Observation count"
  ) +
  base_theme

panel_b <- ggplot2::ggplot(df, ggplot2::aes(x = log10_abundance)) +
  ggplot2::geom_histogram(bins = 70, fill = pal$green, color = "white", linewidth = 0.15) +
  ggplot2::geom_vline(xintercept = log10(pseudocount), color = pal$orange, linewidth = 0.55, linetype = "longdash") +
  ggplot2::annotate(
    "label",
    x = log10(pseudocount),
    y = Inf,
    label = paste0("Zeros map to\nlog10(pseudocount) = ", sprintf("%.2f", log10(pseudocount))),
    vjust = 1.12,
    hjust = -0.02,
    size = 2.7,
    fill = "white",
    color = pal$text
  ) +
  ggplot2::labs(
    title = "B  Modeling scale after log10(Y + pseudocount)",
    subtitle = paste0("Pseudocount = half of smallest non-zero Y = ", fmt_sci(pseudocount)),
    x = "log10(Y + pseudocount)",
    y = "Observation count"
  ) +
  base_theme

panel_c <- ggplot2::ggplot(mag_summary, ggplot2::aes(x = zero_fraction)) +
  ggplot2::geom_histogram(binwidth = 0.05, boundary = 0, fill = pal$sky, color = "white", linewidth = 0.2) +
  ggplot2::geom_vline(xintercept = median_zero_fraction, color = pal$blue, linewidth = 0.65) +
  ggplot2::geom_vline(xintercept = zero_fraction_threshold, color = pal$vermillion, linewidth = 0.6, linetype = "longdash") +
  ggplot2::annotate(
    "text",
    x = median_zero_fraction,
    y = Inf,
    label = paste0("median = ", fmt_pct(median_zero_fraction)),
    vjust = 1.5,
    hjust = -0.05,
    size = 2.55,
    color = pal$blue
  ) +
  ggplot2::annotate(
    "text",
    x = zero_fraction_threshold,
    y = Inf,
    label = "10 non-zero\nsample screen",
    vjust = 1.5,
    hjust = 1.05,
    size = 2.55,
    color = pal$vermillion
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  ggplot2::labs(
    title = "C  Sparsity varies substantially across MAGs",
    subtitle = subject_threshold_label,
    x = "Zero fraction per MAG",
    y = "MAG count"
  ) +
  base_theme

panel_d <- ggplot2::ggplot(
  mag_summary,
  ggplot2::aes(x = mean_abundance + pseudocount, y = zero_fraction, color = model_status)
) +
  ggplot2::geom_point(alpha = 0.78, size = 1.9) +
  ggplot2::scale_x_log10(labels = function(x) formatC(x, format = "e", digits = 1)) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  ggplot2::scale_color_manual(values = status_palette, drop = FALSE, name = "Model status") +
  ggplot2::labs(
    title = "D  Modelability reflects abundance and sparsity",
    subtitle = paste0(modelable_mags, " of ", total_mags, " MAGs were primary modelable; flagged MAGs remain in audit tables"),
    x = "Mean raw abundance per MAG + pseudocount (log10 axis)",
    y = "Zero fraction per MAG"
  ) +
  base_theme +
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(0.30, "in")
  )

summary_text <- paste0(
  "Data structure: ", fmt_int(total_values), " MAG-by-sample observations from ",
  fmt_int(total_mags), " MAGs and ", fmt_int(dplyr::n_distinct(df$sampleID)),
  " samples. Zeros = ", fmt_int(zero_values), " (", fmt_pct(zero_fraction_global),
  "); non-zero values = ", fmt_int(nonzero_values), ". Smallest non-zero Y = ",
  fmt_sci(smallest_nonzero_y), "; pseudocount = ", fmt_sci(pseudocount),
  ". Median per-MAG response variance on log10 scale = ",
  sprintf("%.2f", median_response_variance), "."
)

draw_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(gp = grid::gpar(fill = "white", col = NA)))
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  layout <- grid::grid.layout(
    nrow = 4,
    ncol = 2,
    heights = grid::unit(c(0.48, 1, 1, 0.34), c("in", "null", "null", "in")),
    widths = grid::unit(c(1, 1), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))

  grid::grid.text(
    "Figure 2 | MAG abundance distribution, sparsity, and transformation justification",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.88, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 15.0, fontface = "bold", col = pal$ink),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )
  grid::grid.text(
    "Sparse, zero-heavy, right-skewed MAG abundance data motivate the practical log10(Y + pseudocount) LMM screening scale.",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.20, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 9.0, col = pal$text),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )

  print(panel_a, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panel_b, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(panel_c, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
  print(panel_d, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 2))

  grid::grid.text(
    summary_text,
    x = grid::unit(0, "npc"),
    y = grid::unit(0.57, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 7.6, col = pal$muted),
    vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1:2)
  )

  grid::popViewport(2)
}

png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

png(filename = png_path, width = 12.4, height = 9.2, units = "in", res = 420, bg = "white")
draw_figure()
dev.off()

pdf(file = pdf_path, width = 12.4, height = 9.2, bg = "white", useDingbats = FALSE)
draw_figure()
dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = svg_path, width = 12.4, height = 9.2, bg = "white")
  draw_figure()
  dev.off()
}

summary_path <- file.path(output_dir, paste0(output_stem, "_abundance_transformation_summary.csv"))
utils::write.csv(
  data.frame(
    metric = c(
      "total_mag_sample_observations",
      "total_mags",
      "total_samples",
      "zero_values",
      "nonzero_values",
      "global_zero_fraction",
      "smallest_nonzero_y",
      "pseudocount",
      "median_zero_fraction_per_mag",
      "median_log10_response_variance_per_mag",
      "modelable_mags",
      "too_sparse_mags",
      "singular_lmm_mags"
    ),
    value = c(
      total_values,
      total_mags,
      dplyr::n_distinct(df$sampleID),
      zero_values,
      nonzero_values,
      zero_fraction_global,
      smallest_nonzero_y,
      pseudocount,
      median_zero_fraction,
      median_response_variance,
      modelable_mags,
      too_sparse_mags,
      singular_mags
    )
  ),
  summary_path,
  row.names = FALSE
)

message("Figure 2 abundance/transformation outputs created:")
message(" - ", png_path)
message(" - ", pdf_path)
if (file.exists(svg_path)) message(" - ", svg_path)
message(" - ", summary_path)
