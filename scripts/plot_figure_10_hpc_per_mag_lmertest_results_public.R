# Public renderer: values come from external inputs; see docs/public-rendering.md.
# =================================================
# plot_figure_10_hpc_per_mag_lmertest_results.R
# =================================================
# Figure 10: biological interpretation heatmap for top-ranked MAGs from the
# full T1-T10 per-MAG lmerTest LMM workflow.
#
# Selection rules are preserved; significance counts are calculated from supplied results.

required_packages <- c("ggplot2", "dplyr", "grid")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
data_path <- file.path(repo_root, "data", "processed", "cleaned_data.rds")
results_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "tables",
  "per_mag_lmertest_hpc_results.csv"
)
taxonomy_path <- file.path(repo_root, "data", "GDTB.csv")
output_dir <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "figures")
output_stem <- "figure_10_hpc_per_mag_lmertest_results"

# External inputs and isolated output destination (run from repository root).
data_path <- Sys.getenv("PUBLIC_CLEANED_DATA", unset = data_path)
results_path <- Sys.getenv("PUBLIC_LMM_RESULTS", unset = results_path)
taxonomy_path <- Sys.getenv("PUBLIC_TAXONOMY", unset = taxonomy_path)
output_dir <- Sys.getenv("PUBLIC_RENDER_OUTPUT_DIR", unset = output_dir)

if (!file.exists(data_path)) {
  stop("Missing cleaned data file: ", data_path, call. = FALSE)
}
if (!file.exists(results_path)) {
  stop("Missing per-MAG LMM result table: ", results_path, call. = FALSE)
}
if (!file.exists(taxonomy_path)) {
  stop("Missing GTDB-Tk taxonomy table: ", taxonomy_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cleaned <- readRDS(data_path)
results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
taxonomy <- utils::read.csv(taxonomy_path, check.names = FALSE, stringsAsFactors = FALSE)

required_data_columns <- c("sampleID", "MAG", "Y", "Subject", "Timepoint", "days_from_baseline_subject", "altitude")
required_result_columns <- c(
  "MAG", "status", "primary_analysis", "time_estimate", "time_p_value",
  "time_p_adj_bh_primary", "primary_significant", "zero_fraction",
  "n_subjects", "n_timepoints"
)
required_taxonomy_columns <- c("Genome", "phylum", "class", "ordeer", "family", "genus", "species")

missing_data_columns <- setdiff(required_data_columns, names(cleaned))
missing_result_columns <- setdiff(required_result_columns, names(results))
missing_taxonomy_columns <- setdiff(required_taxonomy_columns, names(taxonomy))
if (length(missing_data_columns) > 0) {
  stop("Missing required cleaned-data columns: ", paste(missing_data_columns, collapse = ", "), call. = FALSE)
}
if (length(missing_result_columns) > 0) {
  stop("Missing required result columns: ", paste(missing_result_columns, collapse = ", "), call. = FALSE)
}
if (length(missing_taxonomy_columns) > 0) {
  stop("Missing required taxonomy columns: ", paste(missing_taxonomy_columns, collapse = ", "), call. = FALSE)
}

clean_taxon <- function(x, prefix) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == paste0(prefix, "__")] <- "Unclassified"
  x <- sub(paste0("^", prefix, "__"), "", x)
  x <- gsub("_", " ", x)
  x
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x)))
}

shorten <- function(x, n = 24) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

taxonomy_clean <- taxonomy |>
  dplyr::transmute(
    MAG = Genome,
    phylum = clean_taxon(phylum, "p"),
    class = clean_taxon(class, "c"),
    order = clean_taxon(ordeer, "o"),
    family = clean_taxon(family, "f"),
    genus = clean_taxon(genus, "g"),
    species = clean_taxon(species, "s")
  )

annotated_results <- results |>
  dplyr::left_join(taxonomy_clean, by = "MAG") |>
  dplyr::mutate(
    direction = dplyr::case_when(
      time_estimate > 0 ~ "Increasing",
      time_estimate < 0 ~ "Decreasing",
      TRUE ~ "No directional estimate"
    ),
    taxon_label = dplyr::case_when(
      !is.na(genus) & genus != "Unclassified" ~ genus,
      !is.na(family) & family != "Unclassified" ~ family,
      !is.na(order) & order != "Unclassified" ~ order,
      !is.na(class) & class != "Unclassified" ~ class,
      TRUE ~ "Unclassified"
    )
  )

significant_candidates <- annotated_results |>
  dplyr::filter(primary_significant %in% TRUE, status == "lmm_ok") |>
  dplyr::arrange(time_p_adj_bh_primary, time_p_value)

if (nrow(significant_candidates) > 0) {
  selection_label <- "BH-FDR-significant MAGs"
  selected <- significant_candidates |>
    dplyr::slice_head(n = 20)
} else {
  selection_label <- "Exploratory nominal candidates; no MAG passed BH-FDR q < 0.05"
  selected <- annotated_results |>
    dplyr::filter(primary_analysis %in% TRUE, status == "lmm_ok", time_p_value < 0.05) |>
    dplyr::arrange(time_p_adj_bh_primary, time_p_value, dplyr::desc(abs(time_estimate))) |>
    dplyr::slice_head(n = 20)
}

if (nrow(selected) == 0) {
  stop("No selected MAGs available for Figure 10.", call. = FALSE)
}

selected <- selected |>
  dplyr::mutate(
    mag_rank = dplyr::row_number(),
    mag_display = paste0(sprintf("%02d", mag_rank), "  ", shorten(taxon_label, 22)),
    mag_display = factor(mag_display, levels = rev(mag_display)),
    direction = factor(direction, levels = c("Increasing", "Decreasing", "No directional estimate"))
  )

nonzero_values <- cleaned$Y[is.finite(cleaned$Y) & cleaned$Y > 0]
if (length(nonzero_values) == 0) {
  stop("Cannot derive pseudocount because no positive abundance values were found.", call. = FALSE)
}
pseudocount <- min(nonzero_values, na.rm = TRUE) / 2

sample_meta <- cleaned |>
  dplyr::distinct(sampleID, Subject, Timepoint, days_from_baseline_subject, altitude) |>
  dplyr::arrange(Subject, days_from_baseline_subject, Timepoint, sampleID) |>
  dplyr::mutate(sample_index = dplyr::row_number())

sample_order <- sample_meta$sampleID

heatmap_df <- cleaned |>
  dplyr::filter(MAG %in% selected$MAG) |>
  dplyr::select(sampleID, MAG, Y) |>
  dplyr::left_join(sample_meta, by = "sampleID") |>
  dplyr::left_join(selected |> dplyr::select(MAG, mag_rank, mag_display), by = "MAG") |>
  dplyr::mutate(log10_abundance = log10(Y + pseudocount)) |>
  dplyr::group_by(MAG) |>
  dplyr::mutate(
    mag_mean = mean(log10_abundance, na.rm = TRUE),
    mag_sd = stats::sd(log10_abundance, na.rm = TRUE),
    z_score = ifelse(is.finite(mag_sd) & mag_sd > 0, (log10_abundance - mag_mean) / mag_sd, 0),
    z_score_capped = pmax(pmin(z_score, 2.5), -2.5)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    sampleID = factor(sampleID, levels = sample_order),
    mag_display = factor(mag_display, levels = levels(selected$mag_display))
  )

subject_breaks <- sample_meta |>
  dplyr::group_by(Subject) |>
  dplyr::summarise(
    start = min(sample_index),
    end = max(sample_index),
    center = mean(range(sample_index)),
    .groups = "drop"
  )

timepoint_levels <- paste0("T", seq_len(10))
altitude_levels <- sort(unique(sample_meta$altitude))
altitude_keys <- paste0("Altitude_", altitude_levels)

annotation_df <- dplyr::bind_rows(
  sample_meta |>
    dplyr::transmute(sample_index, annotation = "Subject", fill_key = paste0("Subject_", Subject), label = as.character(Subject)),
  sample_meta |>
    dplyr::transmute(sample_index, annotation = "Timepoint", fill_key = paste0("Timepoint_", Timepoint), label = as.character(Timepoint)),
  sample_meta |>
    dplyr::transmute(sample_index, annotation = "Altitude", fill_key = paste0("Altitude_", altitude), label = paste0(altitude, " m"))
) |>
  dplyr::mutate(annotation = factor(annotation, levels = c("Altitude", "Timepoint", "Subject")))

subject_ids <- unique(sample_meta$Subject)
subject_colors <- rep(c("#F3F4F6", "#D1D5DB"), length.out = length(subject_ids))
names(subject_colors) <- paste0("Subject_", subject_ids)

timepoint_colors <- c(
  "T1" = "#F7FCF0", "T2" = "#E0F3DB", "T3" = "#CCEBC5", "T4" = "#A8DDB5",
  "T5" = "#7BCCC4", "T6" = "#4EB3D3", "T7" = "#2B8CBE", "T8" = "#0868AC",
  "T9" = "#084081", "T10" = "#4B5563"
)
timepoint_colors <- timepoint_colors[intersect(names(timepoint_colors), timepoint_levels)]
names(timepoint_colors) <- paste0("Timepoint_", names(timepoint_colors))

altitude_palette <- grDevices::colorRampPalette(c("#F8FAFC", "#DCEEF8", "#8CC5DF", "#D55E00"))(length(altitude_levels))
names(altitude_palette) <- altitude_keys
annotation_colors <- c(subject_colors, timepoint_colors, altitude_palette)

row_table <- selected |>
  dplyr::transmute(
    mag_display,
    direction,
    beta = sprintf("%+.3f", time_estimate),
    q_value = fmt_p(time_p_adj_bh_primary),
    raw_p = fmt_p(time_p_value),
    zero = paste0(round(100 * zero_fraction), "%"),
    family = shorten(family, 18)
  )

row_table_long <- do.call(
  rbind,
  lapply(seq_len(nrow(row_table)), function(i) {
    data.frame(
      mag_display = row_table$mag_display[[i]],
      metric = c("Direction", "beta/day", "q", "raw p", "zeros", "family"),
      value = c(
        ifelse(row_table$direction[[i]] == "Increasing", "up", "down"),
        row_table$beta[[i]],
        row_table$q_value[[i]],
        row_table$raw_p[[i]],
        row_table$zero[[i]],
        row_table$family[[i]]
      ),
      direction = as.character(row_table$direction[[i]]),
      stringsAsFactors = FALSE
    )
  })
)
row_table_long$mag_display <- factor(row_table_long$mag_display, levels = levels(selected$mag_display))
row_table_long$metric <- factor(row_table_long$metric, levels = c("Direction", "beta/day", "q", "raw p", "zeros", "family"))

selected_output <- selected |>
  dplyr::select(
    mag_rank, MAG, taxon_label, phylum, class, order, family, genus, species,
    direction, time_estimate, time_p_value, time_p_adj_bh_primary, zero_fraction,
    n_subjects, n_timepoints, status
  )

matrix_output <- heatmap_df |>
  dplyr::transmute(
    MAG, sampleID = as.character(sampleID), Subject, Timepoint,
    days_from_baseline_subject, altitude, raw_abundance = Y,
    log10_abundance, z_score, z_score_capped
  )

utils::write.csv(
  selected_output,
  file.path(output_dir, paste0(output_stem, "_selected_mags.csv")),
  row.names = FALSE
)
utils::write.csv(
  matrix_output,
  file.path(output_dir, paste0(output_stem, "_heatmap_values.csv")),
  row.names = FALSE
)

total_mags <- nrow(results)
primary_mags <- sum(results$primary_analysis %in% TRUE, na.rm = TRUE)
nominal_mags <- sum(results$primary_analysis %in% TRUE & results$time_p_value < 0.05, na.rm = TRUE)
fdr_mags <- sum(results$primary_significant %in% TRUE, na.rm = TRUE)
best_q <- min(results$time_p_adj_bh_primary[results$primary_analysis %in% TRUE], na.rm = TRUE)

pal <- list(
  ink = "#111827",
  text = "#374151",
  muted = "#6B7280",
  hair = "#D1D5DB",
  faint = "#F8FAFC",
  blue = "#0072B2",
  orange = "#D55E00",
  red = "#CC3311"
)

base_theme <- ggplot2::theme_minimal(base_size = 9.4, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11.0, color = pal$ink),
    plot.subtitle = ggplot2::element_text(size = 8.0, color = pal$muted, margin = ggplot2::margin(b = 5)),
    axis.title = ggplot2::element_text(size = 8.4, color = pal$text),
    axis.text = ggplot2::element_text(size = 7.3, color = pal$text),
    panel.grid = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 8),
    legend.text = ggplot2::element_text(size = 7.4),
    plot.margin = ggplot2::margin(4, 6, 4, 6)
  )

annotation_plot <- ggplot2::ggplot(annotation_df, ggplot2::aes(x = sample_index, y = annotation, fill = fill_key)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.12) +
  ggplot2::geom_text(
    data = subject_breaks,
    ggplot2::aes(x = center, y = "Subject", label = Subject),
    inherit.aes = FALSE,
    size = 2.0,
    color = pal$text
  ) +
  ggplot2::scale_fill_manual(values = annotation_colors, guide = "none") +
  ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0.5, nrow(sample_meta) + 0.5)) +
  ggplot2::labs(title = "A  Sample ordering and expedition annotations", x = NULL, y = NULL) +
  base_theme +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(face = "bold", size = 7.8),
    plot.title = ggplot2::element_text(size = 9.8)
  )

heatmap_plot <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = sample_index, y = mag_display, fill = z_score_capped)) +
  ggplot2::geom_tile(color = "#F8FAFC", linewidth = 0.08) +
  ggplot2::geom_vline(xintercept = subject_breaks$end + 0.5, color = "white", linewidth = 0.35) +
  ggplot2::scale_fill_gradient2(
    low = pal$blue,
    mid = "#F8FAFC",
    high = pal$orange,
    midpoint = 0,
    limits = c(-2.5, 2.5),
    breaks = c(-2, 0, 2),
    name = "Row z-score"
  ) +
  ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0.5, nrow(sample_meta) + 0.5)) +
  ggplot2::labs(
    title = "B  Row-scaled transformed abundance across subjects and expedition stages",
    subtitle = "Rows are selected MAGs; columns are samples ordered by subject and days from baseline.",
    x = "Samples ordered by Subject and days from baseline",
    y = NULL
  ) +
  base_theme +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 6.9, color = pal$ink),
    legend.key.width = grid::unit(1.15, "cm")
  )

direction_text_colors <- c(
  "Increasing" = pal$blue,
  "Decreasing" = pal$red,
  "No directional estimate" = pal$muted
)

row_table_plot <- ggplot2::ggplot(row_table_long, ggplot2::aes(x = metric, y = mag_display, label = value)) +
  ggplot2::geom_text(
    ggplot2::aes(color = ifelse(metric == "Direction", direction, "Other")),
    size = 2.15,
    hjust = 0.5
  ) +
  ggplot2::scale_color_manual(values = c(direction_text_colors, "Other" = pal$text), guide = "none") +
  ggplot2::labs(title = "C  Model and taxonomy annotations", x = NULL, y = NULL) +
  base_theme +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 6.9, face = "bold", angle = 0),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(size = 9.8),
    plot.margin = ggplot2::margin(4, 3, 4, 0)
  )

title_grob <- grid::textGrob(
  "Figure 10 | Biological interpretation heatmap of top-ranked longitudinal MAGs",
  x = 0,
  y = 0.95,
  just = c("left", "top"),
  gp = grid::gpar(fontsize = 17, fontface = "bold", col = pal$ink)
)

subtitle_grob <- grid::textGrob(
  paste0(
    selection_label,
    ". Displayed abundance is log10(Y + pseudocount), row-scaled for visualization only. ",
    "Primary family: ", primary_mags, " modeled MAGs; nominal raw p < 0.05: ", nominal_mags,
    "; FDR q < 0.05: ", fdr_mags, "; best q = ", sprintf("%.3f", best_q), "."
  ),
  x = 0,
  y = 0.35,
  just = c("left", "top"),
  gp = grid::gpar(fontsize = 9.2, col = pal$muted)
)

footer_grob <- grid::textGrob(
  paste0(
    "Interpretation: this heatmap visualizes structured abundance patterns for selected MAGs; ",
    "it does not create new statistical evidence. Pseudocount = ", formatC(pseudocount, format = "e", digits = 2),
    ". Taxonomy joined from data/GDTB.csv."
  ),
  x = 0,
  y = 0.5,
  just = c("left", "center"),
  gp = grid::gpar(fontsize = 7.8, col = pal$muted)
)

draw_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = 4,
    ncol = 2,
    heights = grid::unit(c(0.65, 0.8, 7.7, 0.42), c("in", "in", "null", "in")),
    widths = grid::unit(c(8.7, 3.6), c("null", "null"))
  )))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  grid::grid.draw(title_grob)
  grid::grid.draw(subtitle_grob)
  grid::popViewport()
  print(annotation_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(heatmap_plot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
  print(row_table_plot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
  grid::pushViewport(grid::viewport(layout.pos.row = 4, layout.pos.col = 1:2))
  grid::grid.draw(footer_grob)
  grid::popViewport(2)
}

png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

grDevices::png(png_path, width = 14, height = 10.2, units = "in", res = 450, type = "cairo")
draw_figure()
grDevices::dev.off()

grDevices::pdf(pdf_path, width = 14, height = 10.2, useDingbats = FALSE)
draw_figure()
grDevices::dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(svg_path, width = 14, height = 10.2)
  draw_figure()
  grDevices::dev.off()
} else {
  warning("Package 'svglite' is not installed; SVG output was skipped.", call. = FALSE)
}

message("Figure 10 saved:")
message("- ", png_path)
message("- ", pdf_path)
if (file.exists(svg_path)) {
  message("- ", svg_path)
}
message("Selected MAG table saved: ", file.path(output_dir, paste0(output_stem, "_selected_mags.csv")))
message("Heatmap value table saved: ", file.path(output_dir, paste0(output_stem, "_heatmap_values.csv")))
