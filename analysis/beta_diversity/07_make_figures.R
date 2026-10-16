#!/usr/bin/env Rscript
# Purpose: Create publication-oriented ordination, distance, longitudinal, and
# dispersion figures from completed beta-diversity stages.
# Inputs: Ordination coordinates, distance matrices, statistical tables, and
# prepared metadata under the configured results directory.
# Outputs: PDF and PNG figures under the configured figures directory.
# Dependencies: Shared plotting helpers, ggplot2, and packages in environment.yml.
# Usage: Rscript analysis/beta_diversity/07_make_figures.R --config config/beta_diversity.yml
# Workflow: Stage 07; consumes stages 03-06 and precedes final reporting.

source("analysis/beta_diversity/utils/io_helpers.R")
source("analysis/beta_diversity/utils/provenance_helpers.R")
source("analysis/beta_diversity/utils/plotting_theme.R")
args <- parse_beta_args(); config <- load_beta_config(args$config); ensure_beta_dirs(config)
require_beta_packages(c("ggplot2", "scales"))

with_stage_log(config, "07_make_figures", {
  coord_files <- list.files(beta_output_dir(config, "ordination"), pattern = "^pcoa_.*_coordinates\\.tsv$", full.names = TRUE)
  coord_files <- coord_files[!grepl("\\.rds_coordinates\\.tsv$", coord_files)]
  summary_coords <- list()
  for (file in coord_files) {
    metric <- sub("_coordinates\\.tsv$", "", sub("^pcoa_", "", basename(file)))
    coord <- utils::read.delim(file, check.names = FALSE); coord$timepoint <- factor(coord$timepoint, levels = unlist(config$columns$ordered_time_levels), ordered = FALSE)
    eig <- utils::read.delim(beta_output_dir(config, "ordination", paste0("pcoa_", metric, "_eigenvalues.tsv")))
    xlab <- sprintf("PCoA1 (%.1f%%)", eig$variance_percent[1]); ylab <- sprintf("PCoA2 (%.1f%%)", eig$variance_percent[2])
    base <- ggplot2::ggplot(coord, ggplot2::aes(PCoA1, PCoA2, color = timepoint)) +
      ggplot2::geom_point(ggplot2::aes(shape = timepoint), size = 2, alpha = 0.85) +
      ggplot2::scale_color_manual(values = beta_palette, drop = FALSE) +
      ggplot2::scale_shape_manual(values = beta_shapes, drop = FALSE) +
      ggplot2::labs(x = xlab, y = ylab, color = "Timepoint", shape = "Timepoint", title = gsub("_", " ", metric)) + theme_beta()
    trajectory <- base + ggplot2::geom_path(ggplot2::aes(group = subject_id), color = "grey45", alpha = 0.35, linewidth = 0.35) +
      ggplot2::geom_point(ggplot2::aes(shape = timepoint), size = 1.8, alpha = 0.9)
    save_beta_plot(base, beta_output_dir(config, "figures", paste0("pcoa_", metric, "_timepoint")))
    save_beta_plot(trajectory, beta_output_dir(config, "figures", paste0("pcoa_", metric, "_trajectories")))
    if ("gender" %in% names(coord) && length(unique(stats::na.omit(coord$gender))) > 1) {
      save_beta_plot(base + ggplot2::facet_wrap(~gender), beta_output_dir(config, "figures", paste0("pcoa_", metric, "_facet_gender")))
    }
    summary_coords[[metric]] <- transform(coord, metric = gsub("_", " ", metric), xlab = xlab, ylab = ylab)
    diag <- utils::read.delim(beta_output_dir(config, "ordination", paste0("pcoa_", metric, "_diagnostics.tsv")))
    scree <- ggplot2::ggplot(head(eig, 10), ggplot2::aes(factor(axis), variance_percent)) + ggplot2::geom_col(fill = "#0072B2", width = 0.7) +
      ggplot2::labs(x = "PCoA axis", y = "Variation represented (%)", subtitle = paste("Correction:", diag$correction[1])) + theme_beta()
    save_beta_plot(scree, beta_output_dir(config, "figures", paste0("pcoa_", metric, "_scree")), 140, 110)
    distance <- as.matrix(readRDS(beta_output_dir(config, "distances", paste0("distance_", metric, ".rds"))))
    order_ids <- coord$sample_id[order(coord$timepoint, coord$subject_id)]
    ordered <- distance[order_ids, order_ids]; long <- as.data.frame(as.table(ordered), stringsAsFactors = FALSE)
    heat <- ggplot2::ggplot(long, ggplot2::aes(Var1, Var2, fill = Freq)) + ggplot2::geom_raster() +
      ggplot2::scale_fill_viridis_c(option = "C", name = "Distance") + ggplot2::labs(x = NULL, y = NULL) +
      theme_beta(7) + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
    save_beta_plot(heat, beta_output_dir(config, "figures", paste0("distance_heatmap_", metric)), 160, 145)
  }
  if (length(summary_coords)) {
    combined <- do.call(rbind, summary_coords[intersect(names(summary_coords), unlist(config$analysis$primary_metrics))])
    composite <- ggplot2::ggplot(combined, ggplot2::aes(PCoA1, PCoA2, color = timepoint, shape = timepoint)) +
      ggplot2::geom_path(ggplot2::aes(group = subject_id), color = "grey60", alpha = 0.25, linewidth = 0.3) +
      ggplot2::geom_point(size = 1.4) + ggplot2::facet_wrap(~metric, scales = "free", nrow = 1) +
      ggplot2::scale_color_manual(values = beta_palette, drop = FALSE) + ggplot2::labs(x = "PCoA axis 1", y = "PCoA axis 2", color = "Timepoint", shape = "Timepoint") + theme_beta(8)
    composite <- composite + ggplot2::scale_shape_manual(values = beta_shapes, drop = FALSE)
    save_beta_plot(composite, beta_output_dir(config, "figures", "primary_multimetric_pcoa"), 270, 105)
  }
  nmds_files <- list.files(beta_output_dir(config, "ordination"), pattern = "^nmds_.*_coordinates\\.tsv$", full.names = TRUE)
  for (file in nmds_files) {
    metric <- sub("_coordinates\\.tsv$", "", sub("^nmds_", "", basename(file)))
    coord <- utils::read.delim(file); coord$timepoint <- factor(coord$timepoint, levels = unlist(config$columns$ordered_time_levels))
    diagnostic <- utils::read.delim(beta_output_dir(config, "ordination", paste0("nmds_", metric, "_diagnostics.tsv")))
    p <- ggplot2::ggplot(coord, ggplot2::aes(NMDS1, NMDS2, color = timepoint, shape = timepoint)) +
      ggplot2::geom_path(ggplot2::aes(group = subject_id), color = "grey55", alpha = 0.35, linewidth = 0.35) +
      ggplot2::geom_point(size = 2) + ggplot2::scale_color_manual(values = beta_palette, drop = FALSE) +
      ggplot2::scale_shape_manual(values = beta_shapes, drop = FALSE) +
      ggplot2::labs(color = "Timepoint", shape = "Timepoint", subtitle = sprintf("2D stress = %.3f", diagnostic$stress[1])) + theme_beta()
    save_beta_plot(p, beta_output_dir(config, "figures", paste0("nmds_", metric)), 180, 140)
  }
  baseline_file <- beta_output_dir(config, "statistics", "distance_to_baseline.tsv")
  if (file.exists(baseline_file)) {
    baseline <- utils::read.delim(baseline_file)
    trend_file <- beta_output_dir(config, "statistics", "longitudinal_lmm_estimated_trends.tsv")
    trend <- if (file.exists(trend_file)) utils::read.delim(trend_file) else NULL
    p <- ggplot2::ggplot(baseline, ggplot2::aes(days_from_subject_baseline, distance_to_baseline, group = subject_id)) +
      ggplot2::geom_line(color = "grey55", alpha = 0.45, linewidth = 0.35) + ggplot2::geom_point(ggplot2::aes(color = timepoint), size = 1.3) +
      ggplot2::facet_wrap(~metric, scales = "free_y") + ggplot2::scale_color_manual(values = beta_palette) +
      ggplot2::labs(x = "Days from subject-specific baseline", y = "Distance to baseline", color = "Timepoint") + theme_beta(8)
    if (!is.null(trend) && "estimate" %in% names(trend)) {
      p <- p + ggplot2::geom_ribbon(data = trend, ggplot2::aes(x = days_from_subject_baseline, ymin = conf_low, ymax = conf_high, group = metric),
                                    inherit.aes = FALSE, fill = "#0072B2", alpha = 0.12) +
        ggplot2::geom_line(data = trend, ggplot2::aes(x = days_from_subject_baseline, y = estimate, group = metric),
                           inherit.aes = FALSE, color = "#0072B2", linewidth = 0.7)
    }
    save_beta_plot(p, beta_output_dir(config, "figures", "distance_to_baseline_trajectories"), 250, 150)
  }
  adjacent_file <- beta_output_dir(config, "statistics", "adjacent_window_distances.tsv")
  if (file.exists(adjacent_file)) {
    adjacent <- utils::read.delim(adjacent_file); adjacent$transition <- paste(adjacent$from_timepoint, adjacent$to_timepoint, sep = " to ")
    p <- ggplot2::ggplot(adjacent, ggplot2::aes(transition, adjacent_distance)) + ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
      ggplot2::geom_jitter(width = 0.12, alpha = 0.55, size = 0.8, color = "#0072B2") + ggplot2::facet_wrap(~metric, scales = "free") +
      ggplot2::labs(x = "Observed adjacent transition", y = "Within-subject distance") + theme_beta(8) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1))
    save_beta_plot(p, beta_output_dir(config, "figures", "adjacent_window_distances"), 260, 170)
  }
  dispersion_file <- beta_output_dir(config, "statistics", "dispersion_results.tsv")
  disp_files <- list.files(beta_output_dir(config, "statistics"), pattern = "^dispersion_distances_.*\\.tsv$", full.names = TRUE)
  if (length(disp_files)) {
    disp <- do.call(rbind, lapply(disp_files, function(f) transform(utils::read.delim(f), metric = sub("\\.tsv$", "", sub("^dispersion_distances_", "", basename(f))))))
    p <- ggplot2::ggplot(disp, ggplot2::aes(timepoint, distance_to_centroid, color = timepoint)) +
      ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.3) + ggplot2::geom_jitter(width = 0.12, size = 0.7, alpha = 0.6) +
      ggplot2::facet_wrap(~metric, scales = "free_y") + ggplot2::scale_color_manual(values = beta_palette, guide = "none") +
      ggplot2::labs(x = "Timepoint", y = "Distance to group centroid") + theme_beta(8)
    save_beta_plot(p, beta_output_dir(config, "figures", "dispersion_by_timepoint"), 250, 150)
  }
  beta_log(config, "07_make_figures", "COMPLETED_FIGURES", paste("ordination_metrics=", length(coord_files)))
})
