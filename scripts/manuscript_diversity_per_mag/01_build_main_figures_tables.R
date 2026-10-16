# Public replacement: all result reporting comes from supplied outputs.
#!/usr/bin/env Rscript
# Purpose: Build reader-first manuscript figures and tables for the completed
# alpha-diversity, beta-diversity, and primary per-MAG longitudinal LMM analyses.
# Inputs: Existing result tables, ordination coordinates, model-ready data, and
# exact MAG taxonomy mappings already produced by the canonical workflows.
# Outputs: Three composite main figures, three presentation tables, their
# machine-readable companions, and an inventory of supplementary evidence.
# Dependencies: Packages declared in environment.yml; no scientific pipeline is
# rerun. The trajectory panel refits only the documented per-MAG display models.
# Usage: Rscript scripts/manuscript_diversity_per_mag/01_build_main_figures_tables.R
# Relationship: Manuscript post-processing only; canonical analysis outputs are
# read without modification and remain the scientific source of truth.

required_packages <- c(
  "dplyr", "tidyr", "readr", "ggplot2", "patchwork", "gridExtra",
  "grid", "scales", "svglite", "ragg", "lme4", "yaml", "knitr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    ". Update the glmm-hpc environment from environment.yml.",
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop(flag, " requires a value.", call. = FALSE)
  args[[idx + 1L]]
}
has_flag <- function(flag) flag %in% args

repo_root <- normalizePath(getwd(), mustWork = TRUE)
input_root <- Sys.getenv("PUBLIC_DIVERSITY_INPUT_ROOT", unset = repo_root)
raw_data_root <- Sys.getenv("BLACKICE_RAW_DATA_DIR", file.path(input_root, "data", "raw"))
processed_data_root <- Sys.getenv("BLACKICE_PROCESSED_DATA_DIR", file.path(input_root, "data", "processed"))
output_root <- arg_value(
  "--output-dir",
  Sys.getenv("PUBLIC_DIVERSITY_OUTPUT_DIR", unset = file.path(repo_root, "results", "manuscript_diversity_per_mag"))
)
force <- has_flag("--force")
raster_dpi <- suppressWarnings(as.integer(Sys.getenv("PUBLIC_DIVERSITY_RASTER_DPI", unset = "600")))
if (length(raster_dpi) != 1L || is.na(raster_dpi) || raster_dpi <= 0L) stop("PUBLIC_DIVERSITY_RASTER_DPI must be a positive integer.", call. = FALSE)

paths <- list(
  alpha_input = file.path(input_root, "results", "alpha_diversity_lmm", "tables", "alpha_lmm_input.csv"),
  alpha_effects = file.path(input_root, "results", "alpha_diversity_lmm", "tables", "alpha_lmm_fixed_effects.csv"),
  alpha_predictions = file.path(input_root, "results", "alpha_diversity_lmm", "tables", "alpha_lmm_predictions.csv"),
  alpha_diagnostics = file.path(input_root, "results", "alpha_diversity_lmm", "tables", "alpha_lmm_diagnostics_summary.csv"),
  beta_permanova = file.path(input_root, "results", "beta_diversity", "statistics", "permanova_results.tsv"),
  beta_dispersion = file.path(input_root, "results", "beta_diversity", "statistics", "dispersion_results.tsv"),
  beta_lmm = file.path(input_root, "results", "beta_diversity", "statistics", "longitudinal_lmm_results.tsv"),
  beta_baseline = file.path(input_root, "results", "beta_diversity", "statistics", "distance_to_baseline.tsv"),
  beta_bray_coordinates = file.path(input_root, "results", "beta_diversity", "ordination", "pcoa_bray_curtis_coordinates.tsv"),
  beta_bray_eigenvalues = file.path(input_root, "results", "beta_diversity", "ordination", "pcoa_bray_curtis_eigenvalues.tsv"),
  beta_aitchison_coordinates = file.path(input_root, "results", "beta_diversity", "ordination", "pcoa_aitchison_coordinates.tsv"),
  beta_aitchison_eigenvalues = file.path(input_root, "results", "beta_diversity", "ordination", "pcoa_aitchison_eigenvalues.tsv"),
  mag_results = file.path(input_root, "results", "hpc_per_mag_lmertest", "tables", "per_mag_lmertest_hpc_results.csv"),
  mag_model_data = file.path(processed_data_root, "cleaned_data.rds"),
  mag_taxonomy = file.path(raw_data_root, "GDTB.csv"),
  mag_rules = file.path(input_root, "config", "model_decision_rules.yml")
)
missing_inputs <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required manuscript input(s): ",
    paste(paste0(missing_inputs, "=", unlist(paths[missing_inputs])), collapse = "; "),
    call. = FALSE
  )
}

output_dirs <- file.path(
  output_root,
  c("alpha_diversity", "beta_diversity", "per_mag_lmm", "supplementary")
)
expected_stems <- file.path(
  output_dirs[1:3],
  c("main_figure_alpha_diversity", "main_figure_beta_diversity", "main_figure_per_mag_lmm")
)
existing_outputs <- unlist(lapply(expected_stems, function(x) paste0(x, c(".pdf", ".svg", ".png", ".tiff"))))
if (!force && any(file.exists(existing_outputs))) {
  stop(
    "Manuscript outputs already exist. Review them first, then rerun with --force to replace only this post-processing output directory.",
    call. = FALSE
  )
}


require_columns <- function(data, columns, source_name) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stop(source_name, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

read_csv_checked <- function(path) readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
read_tsv_checked <- function(path) readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)

# Fail on required source structures before writing any analytical output.
validate_required_inputs <- function() {
  require_columns(read_csv_checked(paths$alpha_input), c("sampleID", "subjectID", "Timepoint", "days_from_baseline", "altitude", "Shannon", "Simpson", "Observed_richness"), "Alpha-diversity LMM input")
  require_columns(read_csv_checked(paths$alpha_effects), c("metric", "metric_label", "model_type", "term", "estimate", "standard_error", "p_value", "confidence_interval_low", "confidence_interval_high", "n_samples", "n_subjects"), "Alpha-diversity fixed effects")
  require_columns(read_csv_checked(paths$alpha_predictions), c("metric", "metric_label", "model_type", "time_value", "estimate", "confidence_interval_low", "confidence_interval_high"), "Alpha-diversity predictions")
  require_columns(read_tsv_checked(paths$beta_permanova), c("metric", "hypothesis", "r_squared", "statistic", "p_value", "p_adjusted_bh", "sample_count", "permutations"), "Beta-diversity PERMANOVA")
  require_columns(read_tsv_checked(paths$beta_dispersion), c("metric", "hypothesis", "statistic", "p_value", "p_adjusted_bh", "sample_count"), "Beta-diversity dispersion")
  require_columns(read_tsv_checked(paths$beta_lmm), c("metric", "estimate", "std_error", "conf_low", "conf_high", "sample_count", "subject_count", "singular"), "Beta-diversity longitudinal LMM")
  require_columns(read_tsv_checked(paths$beta_baseline), c("metric", "sample_id", "subject_id", "timepoint", "days_from_subject_baseline", "distance_to_baseline"), "Distance-to-baseline table")
  require_columns(read_tsv_checked(paths$beta_bray_coordinates), c("sample_id", "subject_id", "timepoint", "sample_date", "PCoA1", "PCoA2"), "Bray-Curtis PCoA coordinates")
  require_columns(read_tsv_checked(paths$beta_aitchison_coordinates), c("sample_id", "subject_id", "timepoint", "sample_date", "PCoA1", "PCoA2"), "Aitchison PCoA coordinates")
  require_columns(read_csv_checked(paths$mag_results), c("MAG", "n_rows", "n_subjects", "n_timepoints", "zero_fraction", "nonzero_sample_count", "subjects_with_nonzero", "fit_succeeded", "singular", "time_estimate", "time_std_error", "time_p_value", "status", "primary_analysis", "time_p_adj_bh_primary", "primary_significant", "nominal_p_lt_0_05"), "Per-MAG LMM results")
  require_columns(readRDS(paths$mag_model_data), c("sampleID", "MAG", "Y", "Subject", "Timepoint", "days_from_baseline_subject"), "Per-MAG model-ready data")
  require_columns(read_csv_checked(paths$mag_taxonomy), c("Genome", "label", "phylum", "family", "genus", "species"), "GTDB taxonomy")
  for (key in c("beta_bray_eigenvalues", "beta_aitchison_eigenvalues")) {
    require_columns(read_tsv_checked(paths[[key]]), c("axis", "variance_percent"), key)
  }
}
validate_required_inputs()
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

wrap_panel_text <- function(text, width = 48L) {
  paste(base::strwrap(text, width = width), collapse = "\n")
}

write_table_bundle <- function(
    raw_table,
    display_table,
    stem,
    title,
    subtitle,
    footnote,
    width,
    height,
    column_widths = NULL,
    header_labels = NULL) {
  readr::write_csv(raw_table, paste0(stem, ".csv"), na = "")
  readr::write_tsv(raw_table, paste0(stem, ".tsv"), na = "")
  writeLines(
    c(
      paste0("# ", title),
      "",
      subtitle,
      "",
      knitr::kable(display_table, format = "pipe", align = rep("l", ncol(display_table))),
      "",
      paste0("Note: ", footnote)
    ),
    paste0(stem, ".md")
  )

  grob_table <- display_table
  if (!is.null(header_labels)) {
    if (length(header_labels) != ncol(grob_table)) {
      stop("header_labels must match the number of display-table columns.", call. = FALSE)
    }
    names(grob_table) <- header_labels
  }
  if (!is.null(column_widths) && length(column_widths) != ncol(grob_table)) {
    stop("column_widths must match the number of display-table columns.", call. = FALSE)
  }

  table_theme <- gridExtra::ttheme_minimal(
    base_size = 7.2,
    base_family = "sans",
    core = list(
      fg_params = list(col = "#252525", hjust = 0, x = 0.03),
      bg_params = list(fill = rep(c("#FFFFFF", "#F7F7F7"), length.out = nrow(display_table)), col = NA)
    ),
    colhead = list(
      fg_params = list(col = "#FFFFFF", fontface = "bold", hjust = 0, x = 0.03),
      bg_params = list(fill = "#443F90", col = "#443F90")
    )
  )
  table_grob <- if (nrow(grob_table)) {
    gridExtra::tableGrob(grob_table, rows = NULL, theme = table_theme)
  } else {
    grid::textGrob("No eligible rows in supplied results")
  }
  if (!is.null(column_widths) && nrow(grob_table)) {
    table_grob$widths <- grid::unit(column_widths / sum(column_widths), "npc")
  }
  title_grob <- grid::textGrob(
    wrap_panel_text(title, width = 78L), x = 0, hjust = 0,
    gp = grid::gpar(fontfamily = "sans", fontsize = 12, fontface = "bold", col = "#252525")
  )
  subtitle_grob <- grid::textGrob(
    wrap_panel_text(subtitle, width = 100L), x = 0, hjust = 0,
    gp = grid::gpar(fontfamily = "sans", fontsize = 8, lineheight = 1.05, col = "#555555")
  )
  note_grob <- grid::textGrob(
    wrap_panel_text(footnote, width = 118L), x = 0, hjust = 0,
    gp = grid::gpar(fontfamily = "sans", fontsize = 6.4, lineheight = 1.08, col = "#555555")
  )
  composed <- gridExtra::arrangeGrob(
    title_grob, subtitle_grob, table_grob, note_grob,
    ncol = 1,
    heights = grid::unit.c(
      grid::unit(0.38, "in"), grid::unit(0.42, "in"),
      grid::unit(1, "null"), grid::unit(0.68, "in")
    ),
    padding = grid::unit(0.05, "in")
  )
  ggplot2::ggsave(paste0(stem, ".pdf"), composed, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(stem, ".svg"), composed, width = width, height = height, units = "in", device = svglite::svglite)
  ragg::agg_png(paste0(stem, ".png"), width = width, height = height, units = "in", res = raster_dpi, background = "white")
  grid::grid.draw(composed)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width, height = height, units = "in", res = raster_dpi, compression = "lzw", background = "white")
  grid::grid.draw(composed)
  grDevices::dev.off()
}

save_figure_bundle <- function(plot, stem, width = 7.2, height = 7.6) {
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(stem, ".svg"), plot, width = width, height = height, units = "in", device = svglite::svglite)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, units = "in", dpi = raster_dpi, device = ragg::agg_png, background = "white")
  ggplot2::ggsave(paste0(stem, ".tiff"), plot, width = width, height = height, units = "in", dpi = raster_dpi, device = ragg::agg_tiff, compression = "lzw", background = "white")
}

# Exact roots transcribed from the user-specified Ordnance Survey GeoDataViz
# Toolkit palette sheet (version 0.9.1). Time is also encoded by ordering and
# subject trajectories, so interpretation does not depend on color alone.
gdv <- list(
  purple_dark = "#443F90", purple = "#685BA7", lavender = "#A599CA",
  pink_light = "#F5DDEB", pink = "#F492A5", rose = "#EA6E8A",
  magenta = "#D21C5E", teal = "#089099", orange = "#FD8D3C",
  green = "#008042", gold = "#EDBA50", ink = "#252525",
  text = "#4D4D4D", mid_grey = "#8A8A8A", light_grey = "#D9D9D9",
  faint = "#F4F4F4", white = "#FFFFFF"
)
timepoint_levels <- c("T1", "T2", "Ta", "T3", "Tb", "T4", "T5", "T6", "T7", "T8", "T9", "T10")
timepoint_palette <- stats::setNames(
  grDevices::colorRampPalette(c(gdv$purple_dark, gdv$lavender, gdv$pink, gdv$magenta))(length(timepoint_levels)),
  timepoint_levels
)
direction_palette <- c("Increasing" = gdv$teal, "Decreasing" = gdv$orange, "Other modeled MAGs" = "#BDBDBD")

nature_theme <- function(base_size = 7.3) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1.1, color = gdv$ink, margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.4, color = gdv$text, lineheight = 1.05, margin = ggplot2::margin(b = 5)),
      axis.title = ggplot2::element_text(size = base_size, color = gdv$ink),
      axis.text = ggplot2::element_text(size = base_size - 0.5, color = gdv$text),
      axis.line = ggplot2::element_line(linewidth = 0.35, color = gdv$ink),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, color = gdv$ink),
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 0.1, color = gdv$ink),
      plot.title.position = "plot",
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 0.3),
      legend.text = ggplot2::element_text(size = base_size - 0.6),
      legend.key.height = grid::unit(0.13, "in"),
      legend.box = "horizontal",
      plot.margin = ggplot2::margin(5, 8, 5, 8)
    )
}

panel_annotation <- patchwork::plot_annotation(
  tag_levels = "A",
  theme = ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 11, family = "sans", color = gdv$ink),
    plot.tag.position = c(0.005, 0.995)
  )
)

# -------------------------------------------------------------------------
# Alpha diversity
# -------------------------------------------------------------------------
alpha_input <- read_csv_checked(paths$alpha_input)
alpha_effects <- read_csv_checked(paths$alpha_effects)
alpha_predictions <- read_csv_checked(paths$alpha_predictions)
alpha_diagnostics <- read_csv_checked(paths$alpha_diagnostics)
require_columns(alpha_input, c("sampleID", "subjectID", "Timepoint", "days_from_baseline", "altitude", "Shannon", "Simpson", "Observed_richness"), "Alpha-diversity LMM input")
require_columns(alpha_effects, c("metric", "metric_label", "model_type", "term", "estimate", "standard_error", "p_value", "confidence_interval_low", "confidence_interval_high", "n_samples", "n_subjects"), "Alpha-diversity fixed effects")
require_columns(alpha_predictions, c("metric", "metric_label", "model_type", "time_value", "estimate", "confidence_interval_low", "confidence_interval_high"), "Alpha-diversity predictions")

metric_labels <- c(
  Shannon = "Shannon diversity",
  Simpson = "Simpson diversity",
  Observed_richness = "Observed richness"
)
alpha_input <- alpha_input |>
  dplyr::mutate(Timepoint = factor(Timepoint, levels = timepoint_levels))
alpha_long <- alpha_input |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(names(metric_labels)),
    names_to = "metric", values_to = "metric_value"
  ) |>
  dplyr::mutate(metric_label = unname(metric_labels[metric]))
alpha_main_effects <- alpha_effects |>
  dplyr::filter(model_type == "main", term == "days_from_baseline", metric %in% names(metric_labels)) |>
  dplyr::mutate(
    metric_label = unname(metric_labels[metric]),
    observed_sd = vapply(metric, function(x) stats::sd(alpha_input[[x]], na.rm = TRUE), numeric(1)),
    standardized_change_30d = estimate * 30 / observed_sd,
    standardized_low_30d = confidence_interval_low * 30 / observed_sd,
    standardized_high_30d = confidence_interval_high * 30 / observed_sd
  )

alpha_sample_count <- dplyr::n_distinct(alpha_input$sampleID)
alpha_subject_count <- dplyr::n_distinct(alpha_input$subjectID)
if (nrow(alpha_main_effects) != length(metric_labels) || anyDuplicated(alpha_main_effects$metric)) {
  stop("Required alpha main effects must contain one row for each rendered metric.", call. = FALSE)
}

alpha_coverage <- alpha_input |>
  dplyr::count(Timepoint, name = "n_samples", .drop = FALSE) |>
  dplyr::left_join(
    alpha_input |>
      dplyr::group_by(Timepoint) |>
      dplyr::summarise(altitude = stats::median(altitude, na.rm = TRUE), .groups = "drop"),
    by = "Timepoint"
  )
alpha_panel_coverage <- ggplot2::ggplot(alpha_coverage, ggplot2::aes(Timepoint, n_samples, fill = Timepoint)) +
  ggplot2::geom_col(width = 0.72, color = gdv$ink, linewidth = 0.22, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(label = n_samples), vjust = -0.35, size = 2.25, color = gdv$ink) +
  ggplot2::scale_fill_manual(values = timepoint_palette, drop = FALSE) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.14)), breaks = c(0, 4, 8, 12)) +
  ggplot2::labs(
    title = "Sampling coverage",
    subtitle = wrap_panel_text(
      paste0(alpha_sample_count, " samples from ", alpha_subject_count, " participants; labels show samples per timepoint."),
      width = 47L
    ),
    x = NULL, y = "Samples"
  ) +
  nature_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

make_alpha_metric_panel <- function(metric_name) {
  observed <- dplyr::filter(alpha_long, metric == metric_name)
  prediction <- alpha_predictions |>
    dplyr::filter(metric == metric_name, model_type == "main")
  effect <- dplyr::filter(alpha_main_effects, metric == metric_name)
  ggplot2::ggplot(observed, ggplot2::aes(days_from_baseline, metric_value)) +
    ggplot2::geom_line(ggplot2::aes(group = subjectID), color = gdv$mid_grey, alpha = 0.30, linewidth = 0.32) +
    ggplot2::geom_point(ggplot2::aes(fill = Timepoint), shape = 21, color = gdv$white, stroke = 0.25, size = 1.55, alpha = 0.82) +
    ggplot2::geom_ribbon(
      data = prediction,
      ggplot2::aes(x = time_value, ymin = confidence_interval_low, ymax = confidence_interval_high),
      inherit.aes = FALSE, fill = gdv$lavender, alpha = 0.30
    ) +
    ggplot2::geom_line(
      data = prediction,
      ggplot2::aes(x = time_value, y = estimate),
      inherit.aes = FALSE, color = gdv$purple_dark, linewidth = 0.85
    ) +
    ggplot2::scale_fill_manual(values = timepoint_palette, drop = FALSE, name = "Timepoint") +
    ggplot2::labs(
      title = metric_labels[[metric_name]],
      subtitle = wrap_panel_text(
        paste0("LMM slope = ", scales::number(effect$estimate, accuracy = ifelse(metric_name == "Observed_richness", 0.001, 0.0001)), " per day; p = ", scales::pvalue(effect$p_value, accuracy = 0.001)),
        width = 48L
      ),
      x = "Days from participant baseline", y = metric_labels[[metric_name]]
    ) +
    nature_theme() +
    ggplot2::theme(legend.position = "none")
}

alpha_panel_effect <- ggplot2::ggplot(
  alpha_main_effects,
  ggplot2::aes(standardized_change_30d, stats::reorder(metric_label, standardized_change_30d))
) +
  ggplot2::geom_vline(xintercept = 0, color = gdv$mid_grey, linewidth = 0.45) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = standardized_low_30d, xmax = standardized_high_30d),
    height = 0.16, color = gdv$purple_dark, linewidth = 0.6
  ) +
  ggplot2::geom_point(shape = 21, size = 2.5, fill = gdv$pink, color = gdv$purple_dark, stroke = 0.5) +
  ggplot2::labs(
    title = wrap_panel_text("Standardized temporal effects", width = 38L),
    subtitle = wrap_panel_text(
      "30-day change in observed SD units; intervals are 95% CIs.",
      width = 45L
    ),
    x = "Standardized change per 30 days", y = NULL
  ) +
  nature_theme()

alpha_figure <- (
  (alpha_panel_coverage | alpha_panel_effect) /
    (make_alpha_metric_panel("Shannon") | make_alpha_metric_panel("Simpson")) /
    make_alpha_metric_panel("Observed_richness")
) +
  patchwork::plot_layout(heights = c(0.78, 1, 1.02), guides = "collect") +
  panel_annotation &
  ggplot2::theme(legend.position = "bottom")

alpha_figure_stem <- file.path(output_dirs[[1]], "main_figure_alpha_diversity")
save_figure_bundle(alpha_figure, alpha_figure_stem, width = 7.2, height = 8.5)

alpha_table <- alpha_main_effects |>
  dplyr::transmute(
    metric = metric_label,
    samples = n_samples,
    participants = n_subjects,
    estimate_per_day = estimate,
    ci_low_per_day = confidence_interval_low,
    ci_high_per_day = confidence_interval_high,
    estimate_per_30_days = estimate * 30,
    ci_low_per_30_days = confidence_interval_low * 30,
    ci_high_per_30_days = confidence_interval_high * 30,
    p_value = p_value
  )
alpha_table_display <- alpha_table |>
  dplyr::transmute(
    Metric = metric,
    `n` = samples,
    Participants = participants,
    `Estimate/day` = scales::number(estimate_per_day, accuracy = 0.0001),
    `95% CI/day` = paste0(scales::number(ci_low_per_day, accuracy = 0.0001), " to ", scales::number(ci_high_per_day, accuracy = 0.0001)),
    `Estimate/30 days` = scales::number(estimate_per_30_days, accuracy = 0.001),
    `P value` = scales::pvalue(p_value, accuracy = 0.001)
  )
write_table_bundle(
  alpha_table, alpha_table_display,
  file.path(output_dirs[[1]], "main_table_alpha_diversity"),
  "Table 1 | Longitudinal alpha-diversity model estimates",
  paste0("Random-intercept LMMs; metric-specific counts are tabulated. Input: ", alpha_sample_count, " samples from ", alpha_subject_count, " participants."),
  "Model: diversity metric ~ days from participant baseline + (1 | participant). Confidence intervals are 95%. This panel renders Shannon, Simpson and observed richness only; it does not establish availability or modeling status for other metrics.",
  width = 7.2, height = 3.0,
  column_widths = c(1.55, 0.45, 0.80, 1.00, 1.65, 1.35, 0.80),
  header_labels = c(
    "Metric", "n", "Participants", "Estimate/day",
    "95% CI/day", "Estimate/\n30 days", "P value"
  )
)

# -------------------------------------------------------------------------
# Beta diversity
# -------------------------------------------------------------------------
beta_permanova <- read_tsv_checked(paths$beta_permanova)
beta_dispersion <- read_tsv_checked(paths$beta_dispersion)
beta_lmm <- read_tsv_checked(paths$beta_lmm)
beta_baseline <- read_tsv_checked(paths$beta_baseline)
beta_bray <- read_tsv_checked(paths$beta_bray_coordinates)
beta_bray_eigen <- read_tsv_checked(paths$beta_bray_eigenvalues)
beta_aitchison <- read_tsv_checked(paths$beta_aitchison_coordinates)
beta_aitchison_eigen <- read_tsv_checked(paths$beta_aitchison_eigenvalues)
primary_beta_metrics <- c("bray_curtis", "jaccard_binary", "aitchison")
beta_metric_labels <- c(
  bray_curtis = "Bray-Curtis",
  jaccard_binary = "Binary Jaccard",
  aitchison = "Aitchison"
)
require_columns(beta_permanova, c("metric", "hypothesis", "r_squared", "statistic", "p_value", "p_adjusted_bh", "sample_count", "permutations"), "Beta-diversity PERMANOVA")
require_columns(beta_dispersion, c("metric", "hypothesis", "statistic", "p_value", "p_adjusted_bh", "sample_count"), "Beta-diversity dispersion")
require_columns(beta_lmm, c("metric", "estimate", "std_error", "conf_low", "conf_high", "sample_count", "subject_count", "singular"), "Beta-diversity longitudinal LMM")
require_columns(beta_baseline, c("metric", "sample_id", "subject_id", "timepoint", "days_from_subject_baseline", "distance_to_baseline"), "Distance-to-baseline table")
require_columns(beta_bray, c("sample_id", "subject_id", "timepoint", "sample_date", "PCoA1", "PCoA2"), "Bray-Curtis PCoA coordinates")
require_columns(beta_aitchison, c("sample_id", "subject_id", "timepoint", "sample_date", "PCoA1", "PCoA2"), "Aitchison PCoA coordinates")

prepare_pcoa <- function(data) {
  data |>
    dplyr::mutate(
      timepoint = factor(timepoint, levels = timepoint_levels),
      sample_date = as.Date(sample_date)
    ) |>
    dplyr::arrange(subject_id, sample_date)
}
beta_bray <- prepare_pcoa(beta_bray)
beta_aitchison <- prepare_pcoa(beta_aitchison)

make_pcoa_panel <- function(data, eigen, metric_title, show_timepoint_legend = TRUE) {
  centroids <- data |>
    dplyr::group_by(timepoint) |>
    dplyr::summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2), .groups = "drop") |>
    dplyr::arrange(timepoint)
  axis1 <- eigen$variance_percent[eigen$axis == 1][1]
  axis2 <- eigen$variance_percent[eigen$axis == 2][1]
  ggplot2::ggplot(data, ggplot2::aes(PCoA1, PCoA2)) +
    ggplot2::geom_path(ggplot2::aes(group = subject_id), color = gdv$mid_grey, alpha = 0.28, linewidth = 0.32) +
    ggplot2::geom_point(ggplot2::aes(fill = timepoint), shape = 21, size = 1.65, stroke = 0.25, color = gdv$white, alpha = 0.86) +
    ggplot2::geom_path(data = centroids, color = gdv$ink, linetype = "22", linewidth = 0.55, inherit.aes = FALSE, ggplot2::aes(PCoA1, PCoA2, group = 1)) +
    ggplot2::geom_point(data = centroids, shape = 23, size = 2.2, fill = gdv$gold, color = gdv$ink, stroke = 0.35, inherit.aes = FALSE, ggplot2::aes(PCoA1, PCoA2)) +
    ggplot2::scale_fill_manual(
      values = timepoint_palette,
      drop = FALSE,
      name = "Timepoint",
      guide = if (show_timepoint_legend) ggplot2::guide_legend(nrow = 2, byrow = TRUE) else "none"
    ) +
    ggplot2::labs(
      title = wrap_panel_text(paste0(metric_title, " PCoA"), width = 38L),
      subtitle = wrap_panel_text(
        "Grey lines link repeated samples; gold diamonds show timepoint centroids.",
        width = 43L
      ),
      x = paste0("PCoA1 (", scales::number(axis1, accuracy = 0.1), "%)"),
      y = paste0("PCoA2 (", scales::number(axis2, accuracy = 0.1), "%)")
    ) +
    nature_theme() +
    ggplot2::theme(legend.position = if (show_timepoint_legend) "bottom" else "none")
}

beta_primary <- beta_permanova |>
  dplyr::filter(hypothesis == "timepoint", metric %in% primary_beta_metrics) |>
  dplyr::select(metric, r_squared, statistic, p_value, p_adjusted_bh, sample_count, permutations) |>
  dplyr::left_join(
    beta_dispersion |>
      dplyr::filter(hypothesis == "Groups", metric %in% primary_beta_metrics) |>
      dplyr::select(metric, dispersion_statistic = statistic, dispersion_p_value = p_value, dispersion_p_adjusted_bh = p_adjusted_bh),
    by = "metric"
  ) |>
  dplyr::left_join(
    beta_lmm |>
      dplyr::filter(metric %in% primary_beta_metrics) |>
      dplyr::select(metric, baseline_slope_per_day = estimate, baseline_ci_low = conf_low, baseline_ci_high = conf_high, baseline_samples = sample_count, baseline_subjects = subject_count, baseline_singular = singular),
    by = "metric"
  ) |>
  dplyr::mutate(metric_label = unname(beta_metric_labels[metric]))

beta_panel_effect <- ggplot2::ggplot(
  beta_primary,
  ggplot2::aes(r_squared * 100, stats::reorder(metric_label, r_squared), color = metric_label)
) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = r_squared * 100, yend = stats::reorder(metric_label, r_squared)), linewidth = 0.75, color = gdv$light_grey) +
  ggplot2::geom_point(size = 3.0) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0("BH q=", scales::pvalue(p_adjusted_bh, accuracy = 0.0001))),
    hjust = -0.12, size = 2.15, color = gdv$ink
  ) +
  ggplot2::scale_color_manual(values = c("Bray-Curtis" = gdv$purple_dark, "Binary Jaccard" = gdv$magenta, "Aitchison" = gdv$teal), guide = "none") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.25))) +
  ggplot2::labs(
    title = wrap_panel_text("Within-participant timepoint effects", width = 38L),
    subtitle = wrap_panel_text(
      "Restricted PERMANOVA; effect size is variance attributed to timepoint.",
      width = 43L
    ),
    x = "Timepoint R-squared (%)", y = NULL
  ) +
  nature_theme()

bray_baseline <- beta_baseline |>
  dplyr::filter(metric == "bray_curtis") |>
  dplyr::group_by(subject_id) |>
  dplyr::mutate(days_since_earliest_gut = days_from_subject_baseline - min(days_from_subject_baseline, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::mutate(timepoint = factor(timepoint, levels = timepoint_levels))
bray_summary <- bray_baseline |>
  dplyr::group_by(timepoint) |>
  dplyr::summarise(
    day = stats::median(days_since_earliest_gut, na.rm = TRUE),
    median_distance = stats::median(distance_to_baseline, na.rm = TRUE),
    lower = stats::quantile(distance_to_baseline, 0.25, na.rm = TRUE),
    upper = stats::quantile(distance_to_baseline, 0.75, na.rm = TRUE),
    n = dplyr::n(), .groups = "drop"
  ) |>
  dplyr::arrange(timepoint)
beta_panel_baseline <- ggplot2::ggplot(bray_baseline, ggplot2::aes(days_since_earliest_gut, distance_to_baseline)) +
  ggplot2::geom_line(ggplot2::aes(group = subject_id), color = gdv$mid_grey, alpha = 0.30, linewidth = 0.34) +
  ggplot2::geom_point(ggplot2::aes(fill = timepoint), shape = 21, size = 1.45, color = gdv$white, stroke = 0.2, alpha = 0.78) +
  ggplot2::geom_ribbon(data = bray_summary, ggplot2::aes(x = day, ymin = lower, ymax = upper), inherit.aes = FALSE, fill = gdv$pink_light, alpha = 0.75) +
  ggplot2::geom_line(data = bray_summary, ggplot2::aes(x = day, y = median_distance), inherit.aes = FALSE, color = gdv$magenta, linewidth = 0.85) +
  ggplot2::geom_point(data = bray_summary, ggplot2::aes(x = day, y = median_distance), inherit.aes = FALSE, shape = 21, fill = gdv$white, color = gdv$magenta, size = 1.8, stroke = 0.55) +
  ggplot2::scale_fill_manual(
    values = timepoint_palette,
    drop = FALSE,
    name = "Timepoint",
    guide = "none"
  ) +
  ggplot2::labs(
    title = wrap_panel_text("Bray-Curtis distance from participant baseline", width = 38L),
    subtitle = wrap_panel_text(
      "Grey lines: participants; magenta: median and interquartile range.",
      width = 43L
    ),
    x = "Days since earliest available gut sample", y = "Bray-Curtis distance"
  ) +
  nature_theme() +
  ggplot2::theme(legend.position = "none")

beta_figure <- (
  (make_pcoa_panel(beta_bray, beta_bray_eigen, "Bray-Curtis", show_timepoint_legend = TRUE) |
     make_pcoa_panel(beta_aitchison, beta_aitchison_eigen, "Aitchison", show_timepoint_legend = FALSE)) /
    (beta_panel_effect | beta_panel_baseline)
) +
  patchwork::plot_layout(widths = c(1, 1), heights = c(1, 0.92), guides = "collect") +
  panel_annotation &
  ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")
beta_figure_stem <- file.path(output_dirs[[2]], "main_figure_beta_diversity")
save_figure_bundle(beta_figure, beta_figure_stem, width = 7.2, height = 7.2)

beta_table <- beta_primary |>
  dplyr::transmute(
    metric = metric_label,
    samples = sample_count,
    permutations = permutations,
    timepoint_r_squared = r_squared,
    permanova_p_value = p_value,
    permanova_bh_q = p_adjusted_bh,
    dispersion_statistic = dispersion_statistic,
    dispersion_p_value = dispersion_p_value,
    dispersion_bh_q = dispersion_p_adjusted_bh,
    baseline_slope_per_day = baseline_slope_per_day,
    baseline_ci_low = baseline_ci_low,
    baseline_ci_high = baseline_ci_high
  )
beta_table_display <- beta_table |>
  dplyr::transmute(
    Metric = metric,
    `Timepoint R2` = scales::percent(timepoint_r_squared, accuracy = 0.1),
    `PERMANOVA P` = scales::pvalue(permanova_p_value, accuracy = 0.0001),
    `BH q` = scales::pvalue(permanova_bh_q, accuracy = 0.0001),
    `Dispersion P` = scales::pvalue(dispersion_p_value, accuracy = 0.0001),
    `Dispersion BH q` = scales::pvalue(dispersion_bh_q, accuracy = 0.0001),
    `Distance-to-baseline slope/day (95% CI)` = paste0(
      scales::number(baseline_slope_per_day, accuracy = 0.0001), " (",
      scales::number(baseline_ci_low, accuracy = 0.0001), " to ",
      scales::number(baseline_ci_high, accuracy = 0.0001), ")"
    )
  )
write_table_bundle(
  beta_table, beta_table_display,
  file.path(output_dirs[[2]], "main_table_beta_diversity"),
  "Table 2 | Primary beta-diversity inference",
  "Within-participant timepoint PERMANOVA, dispersion tests, and exploratory distance-to-baseline LMMs.",
  paste0("Reported permutation counts: ", paste(sort(unique(beta_table$permutations)), collapse = ", "), ". P and BH q values are supplied by the upstream analysis; consult its recorded participant restrictions and correction families. This renderer does not refit those tests. Distance slopes have metric-specific units and should not be compared numerically across rows."),
  width = 7.2, height = 3.2,
  column_widths = c(1.15, 0.85, 0.95, 0.60, 0.90, 1.00, 2.25),
  header_labels = c(
    "Metric", "Timepoint\nR2", "PERMANOVA\nP", "BH q",
    "Dispersion\nP", "Dispersion\nBH q",
    "Distance-to-baseline slope/day\n(95% CI)"
  )
)

# -------------------------------------------------------------------------
# Primary per-MAG LMM
# -------------------------------------------------------------------------
mag_results <- read_csv_checked(paths$mag_results)
# readr infers an entirely missing numeric column as logical; preserve missingness.
mag_results$time_p_adj_bh_primary <- as.numeric(mag_results$time_p_adj_bh_primary)
mag_data <- readRDS(paths$mag_model_data)
taxonomy <- read_csv_checked(paths$mag_taxonomy)
rules <- yaml::read_yaml(paths$mag_rules)
require_columns(mag_results, c("MAG", "n_rows", "n_subjects", "n_timepoints", "zero_fraction", "nonzero_sample_count", "subjects_with_nonzero", "fit_succeeded", "singular", "time_estimate", "time_std_error", "time_p_value", "status", "primary_analysis", "time_p_adj_bh_primary", "primary_significant", "nominal_p_lt_0_05"), "Per-MAG LMM results")
require_columns(mag_data, c("sampleID", "MAG", "Y", "Subject", "Timepoint", "days_from_baseline_subject"), "Per-MAG model-ready data")
require_columns(taxonomy, c("Genome", "label", "phylum", "family", "genus", "species"), "GTDB taxonomy")

clean_taxon <- function(x, prefix) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == paste0(prefix, "__")] <- "Unclassified"
  x <- sub(paste0("^", prefix, "__"), "", x)
  gsub("_", " ", x)
}
taxonomy_clean <- taxonomy |>
  dplyr::transmute(
    MAG = Genome,
    mag_label = label,
    phylum = clean_taxon(phylum, "p"),
    family = clean_taxon(family, "f"),
    genus = clean_taxon(genus, "g"),
    species = clean_taxon(species, "s")
  )
mag_annotated <- mag_results |>
  dplyr::left_join(taxonomy_clean, by = "MAG") |>
  dplyr::mutate(
    estimate_30d = time_estimate * 30,
    ci_low_30d = (time_estimate - 1.96 * time_std_error) * 30,
    ci_high_30d = (time_estimate + 1.96 * time_std_error) * 30,
    direction = dplyr::case_when(
      time_estimate > 0 ~ "Increasing",
      time_estimate < 0 ~ "Decreasing",
      TRUE ~ "Other modeled MAGs"
    ),
    result_group = dplyr::case_when(
      primary_analysis %in% TRUE & nominal_p_lt_0_05 %in% TRUE ~ direction,
      TRUE ~ "Other modeled MAGs"
    )
  )
if (any(is.na(mag_annotated$mag_label))) {
  stop("The exact Genome-to-MAG taxonomy join is incomplete.", call. = FALSE)
}

screening <- data.frame(
  stage = factor(
    c("MAGs evaluated", "Primary valid models", "Nominal raw P < 0.05", "BH-FDR q < 0.05"),
    levels = rev(c("MAGs evaluated", "Primary valid models", "Nominal raw P < 0.05", "BH-FDR q < 0.05"))
  ),
  count = c(
    nrow(mag_annotated),
    sum(mag_annotated$primary_analysis %in% TRUE, na.rm = TRUE),
    sum(mag_annotated$primary_analysis %in% TRUE & mag_annotated$nominal_p_lt_0_05 %in% TRUE, na.rm = TRUE),
    sum(mag_annotated$primary_significant %in% TRUE, na.rm = TRUE)
  )
)
primary_count <- sum(mag_annotated$primary_analysis %in% TRUE)
nominal_count <- sum(mag_annotated$primary_analysis %in% TRUE & mag_annotated$nominal_p_lt_0_05 %in% TRUE)
fdr_count <- sum(mag_annotated$primary_significant %in% TRUE)
primary_q <- mag_annotated$time_p_adj_bh_primary[mag_annotated$primary_analysis %in% TRUE]
q_available <- is.finite(primary_q)
minimum_q_text <- if (any(q_available)) scales::pvalue(min(primary_q[q_available]), accuracy = 0.001) else "unavailable"
screening_text <- paste0(nominal_count, " nominal raw-P candidates; ", fdr_count, " supplied primary BH-FDR discovery flags.")
if (any(!q_available) || anyNA(mag_annotated$primary_significant[mag_annotated$primary_analysis %in% TRUE])) {
  screening_text <- paste0(screening_text, " Some primary adjusted values or flags are unavailable; absence of a flag is not evidence of no discovery.")
}

mag_panel_screening <- ggplot2::ggplot(screening, ggplot2::aes(count, stage)) +
  ggplot2::geom_col(ggplot2::aes(fill = stage), width = 0.64, color = gdv$ink, linewidth = 0.25, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(label = count), hjust = -0.15, size = 2.5, color = gdv$ink) +
  ggplot2::scale_fill_manual(values = c(gdv$purple_dark, gdv$purple, gdv$pink, gdv$faint)) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
  ggplot2::labs(
    title = "Feature-wise screening outcome",
    subtitle = wrap_panel_text(
      screening_text,
      width = 43L
    ),
    x = "MAG count", y = NULL
  ) +
  nature_theme()

volcano_data <- mag_annotated |>
  dplyr::filter(primary_analysis %in% TRUE, is.finite(time_p_value), is.finite(estimate_30d)) |>
  dplyr::mutate(neg_log10_p = -log10(pmax(time_p_value, .Machine$double.xmin)))
mag_panel_volcano <- ggplot2::ggplot(volcano_data, ggplot2::aes(estimate_30d, neg_log10_p)) +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = gdv$mid_grey, linetype = "22", linewidth = 0.45) +
  ggplot2::geom_vline(xintercept = 0, color = gdv$light_grey, linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(fill = result_group), shape = 21, size = 1.55, color = gdv$white, stroke = 0.18, alpha = 0.82) +
  ggplot2::scale_fill_manual(
    values = direction_palette,
    name = "Result class",
    guide = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  ) +
  ggplot2::labs(
    title = wrap_panel_text(paste0("Time effects across ", primary_count, " primary models"), width = 38L),
    subtitle = wrap_panel_text(
      "30-day log10 effects; dashed line marks raw P = 0.05.",
      width = 43L
    ),
    x = "LMM estimate per 30 days", y = expression(-log[10](italic(P)))
  ) +
  nature_theme() +
  ggplot2::theme(legend.position = "bottom")

top_mags <- mag_annotated |>
  dplyr::filter(primary_analysis %in% TRUE, is.finite(time_p_value)) |>
  dplyr::arrange(time_p_value) |>
  dplyr::slice_head(n = 10L) |>
  dplyr::mutate(
    forest_label = paste0(mag_label, " | ", genus),
    forest_label = factor(forest_label, levels = rev(forest_label))
  )
mag_panel_forest <- ggplot2::ggplot(top_mags, ggplot2::aes(estimate_30d, forest_label, color = direction)) +
  ggplot2::geom_vline(xintercept = 0, color = gdv$mid_grey, linewidth = 0.45) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_low_30d, xmax = ci_high_30d), height = 0.13, linewidth = 0.55) +
  ggplot2::geom_point(size = 2.3) +
  ggplot2::scale_color_manual(
    values = direction_palette[c("Increasing", "Decreasing")],
    guide = "none"
  ) +
  ggplot2::labs(
    title = "Lowest-P primary time-effect estimates",
    subtitle = wrap_panel_text(
      paste0(nrow(top_mags), " smallest finite raw P values; minimum available primary BH q = ", minimum_q_text, "."),
      width = 43L
    ),
    x = "LMM estimate per 30 days (95% CI)", y = NULL
  ) +
  nature_theme() +
  ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 6.2))

pseudocount <- rules$model_decision_rules$pseudocount$value
if (is.null(pseudocount) || !is.finite(as.numeric(pseudocount))) {
  pseudocount <- min(mag_data$Y[mag_data$Y > 0], na.rm = TRUE) / 2
}
pseudocount <- as.numeric(pseudocount)
trajectory_selection <- mag_annotated |>
  dplyr::filter(primary_analysis %in% TRUE, nominal_p_lt_0_05 %in% TRUE, is.finite(time_p_value)) |>
  dplyr::group_by(direction) |>
  dplyr::arrange(time_p_value, .by_group = TRUE) |>
  dplyr::slice_head(n = 3L) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    panel_label = paste0(
      mag_label, "\n", genus, "\nraw P=",
      scales::pvalue(time_p_value, accuracy = 0.0001),
      "; BH q=", scales::pvalue(time_p_adj_bh_primary, accuracy = 0.001)
    )
  )
trajectory_data <- mag_data |>
  dplyr::filter(MAG %in% trajectory_selection$MAG) |>
  dplyr::inner_join(trajectory_selection |>
                      dplyr::select(MAG, direction, panel_label), by = "MAG") |>
  dplyr::mutate(
    response_value = log10(Y + pseudocount),
    Timepoint = factor(Timepoint, levels = paste0("T", 1:10))
  )

fit_display_model <- function(mag_id) {
  data <- trajectory_data[trajectory_data$MAG == mag_id, , drop = FALSE]
  fit <- lme4::lmer(
    response_value ~ days_from_baseline_subject + (1 | Subject),
    data = data,
    control = lme4::lmerControl(optimizer = "bobyqa")
  )
  prediction <- data.frame(
    days_from_baseline_subject = seq(
      min(data$days_from_baseline_subject, na.rm = TRUE),
      max(data$days_from_baseline_subject, na.rm = TRUE),
      length.out = 80L
    )
  )
  matrix_x <- stats::model.matrix(~ days_from_baseline_subject, prediction)
  fixed <- lme4::fixef(fit)
  covariance <- as.matrix(stats::vcov(fit))
  prediction$estimate <- as.numeric(matrix_x %*% fixed)
  prediction$standard_error <- sqrt(rowSums((matrix_x %*% covariance) * matrix_x))
  metadata <- trajectory_selection[trajectory_selection$MAG == mag_id, , drop = FALSE]
  prediction |>
    dplyr::mutate(
      MAG = mag_id,
      direction = metadata$direction[[1]],
      panel_label = metadata$panel_label[[1]],
      low = estimate - 1.96 * standard_error,
      high = estimate + 1.96 * standard_error
    )
}
trajectory_predictions <- if (nrow(trajectory_selection)) {
  dplyr::bind_rows(lapply(trajectory_selection$MAG, fit_display_model))
} else {
  data.frame(days_from_baseline_subject = numeric(), estimate = numeric(), standard_error = numeric(), MAG = character(), direction = character(), panel_label = character(), low = numeric(), high = numeric())
}
trajectory_levels <- trajectory_selection$panel_label
trajectory_data$panel_label <- factor(trajectory_data$panel_label, levels = trajectory_levels)
trajectory_predictions$panel_label <- factor(trajectory_predictions$panel_label, levels = trajectory_levels)
mag_panel_trajectories <- ggplot2::ggplot(trajectory_data, ggplot2::aes(days_from_baseline_subject, response_value)) +
  ggplot2::geom_line(ggplot2::aes(group = Subject), color = gdv$mid_grey, alpha = 0.28, linewidth = 0.30) +
  ggplot2::geom_point(ggplot2::aes(fill = direction), shape = 21, color = gdv$white, stroke = 0.2, size = 1.2, alpha = 0.78) +
  ggplot2::geom_ribbon(data = trajectory_predictions, ggplot2::aes(x = days_from_baseline_subject, ymin = low, ymax = high, fill = direction), inherit.aes = FALSE, alpha = 0.22) +
  ggplot2::geom_line(data = trajectory_predictions, ggplot2::aes(x = days_from_baseline_subject, y = estimate, color = direction), inherit.aes = FALSE, linewidth = 0.72) +
  ggplot2::facet_wrap(~panel_label, ncol = 2, scales = "free_y") +
  ggplot2::scale_color_manual(values = direction_palette[c("Increasing", "Decreasing")], guide = "none") +
  ggplot2::scale_fill_manual(values = direction_palette[c("Increasing", "Decreasing")], guide = "none") +
  ggplot2::labs(
    title = wrap_panel_text("Representative exploratory trajectories", width = 38L),
    subtitle = wrap_panel_text(
      paste0(nrow(trajectory_selection), " nominal candidates displayed (up to three per direction); grey lines show participants."),
      width = 43L
    ),
    x = "Days from participant baseline", y = expression(log[10](abundance + pseudocount))
  ) +
  nature_theme(base_size = 6.7) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(
      size = 5.8,
      lineheight = 0.92,
      margin = ggplot2::margin(t = 2, r = 2, b = 3, l = 2)
    ),
    strip.clip = "off",
    panel.spacing = grid::unit(0.14, "in")
  )

if (!nrow(trajectory_selection)) {
  mag_panel_trajectories <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = "No eligible nominal candidates in supplied results") +
    ggplot2::theme_void() + ggplot2::labs(title = "Representative trajectories unavailable")
}

mag_figure <- (
  (mag_panel_screening | mag_panel_volcano) /
    (mag_panel_forest | mag_panel_trajectories)
) +
  patchwork::plot_layout(widths = c(0.88, 1.25), heights = c(0.78, 1.42), guides = "collect") +
  panel_annotation &
  ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")
mag_figure_stem <- file.path(output_dirs[[3]], "main_figure_per_mag_lmm")
save_figure_bundle(mag_figure, mag_figure_stem, width = 7.2, height = 8.3)

mag_table <- top_mags |>
  dplyr::transmute(
    rank = dplyr::row_number(),
    MAG = MAG,
    GTDB_label = mag_label,
    genus = genus,
    species = species,
    estimate_per_30_days = estimate_30d,
    ci_low_per_30_days = ci_low_30d,
    ci_high_per_30_days = ci_high_30d,
    raw_p_value = time_p_value,
    primary_bh_q = time_p_adj_bh_primary,
    zero_fraction = zero_fraction,
    subjects_with_nonzero = subjects_with_nonzero
  )
mag_table_display <- mag_table |>
  dplyr::transmute(
    Rank = rank,
    MAG = GTDB_label,
    Genus = genus,
    `Estimate/30 days (95% CI)` = paste0(
      scales::number(estimate_per_30_days, accuracy = 0.001), " (",
      scales::number(ci_low_per_30_days, accuracy = 0.001), " to ",
      scales::number(ci_high_per_30_days, accuracy = 0.001), ")"
    ),
    `Raw P` = scales::pvalue(raw_p_value, accuracy = 0.0001),
    `BH q` = scales::pvalue(primary_bh_q, accuracy = 0.001),
    `Zero samples` = scales::percent(zero_fraction, accuracy = 1),
    `Nonzero participants` = subjects_with_nonzero
  )
write_table_bundle(
  mag_table, mag_table_display,
  file.path(output_dirs[[3]], "main_table_per_mag_lmm"),
  "Table 3 | Highest-ranked exploratory per-MAG temporal signals",
  paste0(nrow(top_mags), " smallest finite raw P values among ", primary_count, " primary valid models; exact MAG identifiers are retained in the machine-readable table."),
  paste0("Model: log10(relative-abundance-like value + pseudocount) ~ days from participant baseline + (1 | participant). Estimates are shown per 30 days. ", screening_text, " Ranking alone does not establish significance or biological causation; interpret supplied adjusted values and diagnostic status."),
  width = 7.2, height = 4.25,
  column_widths = c(0.45, 0.75, 1.20, 2.10, 0.70, 0.55, 0.90, 1.20),
  header_labels = c(
    "Rank", "MAG", "Genus", "Estimate/30 days\n(95% CI)",
    "Raw P", "BH q", "Zero\nsamples", "Nonzero\nparticipants"
  )
)

# -------------------------------------------------------------------------
# Supplementary evidence inventory
# -------------------------------------------------------------------------
supplementary_sources <- data.frame(
  analysis = c("Alpha diversity", "Alpha diversity", "Beta diversity", "Primary per-MAG LMM"),
  source_dir = c(
    file.path(input_root, "results", "alpha_diversity"),
    file.path(input_root, "results", "alpha_diversity_lmm"),
    file.path(input_root, "results", "beta_diversity"),
    file.path(input_root, "results", "hpc_per_mag_lmertest")
  ),
  stringsAsFactors = FALSE
)
supplementary_manifest <- dplyr::bind_rows(lapply(seq_len(nrow(supplementary_sources)), function(i) {
  source_dir <- supplementary_sources$source_dir[[i]]
  files <- list.files(
    source_dir, recursive = TRUE, full.names = TRUE,
    pattern = "\\.(pdf|png|svg|tiff|csv|tsv|md)$", ignore.case = TRUE
  )
  data.frame(
    analysis = rep(supplementary_sources$analysis[[i]], length(files)),
    artifact_type = as.character(ifelse(grepl("\\.(pdf|png|svg|tiff)$", files, ignore.case = TRUE), "figure", ifelse(grepl("\\.(csv|tsv)$", files, ignore.case = TRUE), "table", "report_or_note"))),
    source_path = sub(paste0("^", input_root, "/"), "", files),
    supplementary_role = rep("Retained source evidence; assign final supplementary numbering during manuscript assembly", length(files)),
    stringsAsFactors = FALSE
  )
})) |>
  dplyr::arrange(analysis, artifact_type, source_path)
supplementary_status <- supplementary_sources |>
  dplyr::mutate(
    artifact_count = vapply(source_dir, function(d) length(list.files(d, recursive = TRUE, pattern = "\\.(pdf|png|svg|tiff|csv|tsv|md)$", ignore.case = TRUE)), integer(1)),
    status = ifelse(!dir.exists(source_dir), "missing_optional_directory", ifelse(artifact_count == 0L, "empty_optional_directory", "included")),
    source_dir = sub(paste0("^", input_root, "/"), "", source_dir)
  )
readr::write_csv(supplementary_status, file.path(output_dirs[[4]], "supplementary_inventory_status.csv"))
readr::write_csv(supplementary_manifest, file.path(output_dirs[[4]], "supplementary_evidence_manifest.csv"))
readr::write_tsv(supplementary_manifest, file.path(output_dirs[[4]], "supplementary_evidence_manifest.tsv"))
writeLines(
  c(
    "# Supplementary evidence manifest",
    "",
    "This generated inventory retains every existing figure, table, and report/note from the three source analyses without duplicating or modifying the canonical result artifacts.",
    "",
    knitr::kable(supplementary_manifest, format = "pipe", align = c("l", "l", "l", "l")),
    "", "## Optional supplementary source status", "",
    knitr::kable(supplementary_status, format = "pipe")
  ),
  file.path(output_dirs[[4]], "supplementary_evidence_manifest.md")
)

session_path <- file.path(output_root, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_path)
cat("Manuscript figure/table post-processing completed.\n")
cat("Output directory:", normalizePath(output_root, mustWork = TRUE), "\n")
