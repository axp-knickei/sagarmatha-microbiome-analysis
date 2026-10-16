# Public renderer: values come from external inputs; see docs/public-rendering.md.
# ============================================
# plot_figure_6_hpc_per_mag_lmertest_results.R
# ============================================
# Figure 6: taxonomic identity of time-associated MAGs from the full T1-T10
# per-MAG lmerTest LMM workflow.
#
# Selection rules are preserved; significance counts are calculated from supplied results.

required_packages <- c("ggplot2", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
results_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "tables",
  "per_mag_lmertest_hpc_results.csv"
)
taxonomy_path <- file.path(repo_root, "data", "GDTB.csv")
output_dir <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "figures")
output_stem <- "figure_6_hpc_per_mag_lmertest_results"

# External inputs and isolated output destination (run from repository root).
results_path <- Sys.getenv("PUBLIC_LMM_RESULTS", unset = results_path)
taxonomy_path <- Sys.getenv("PUBLIC_TAXONOMY", unset = taxonomy_path)
output_dir <- Sys.getenv("PUBLIC_RENDER_OUTPUT_DIR", unset = output_dir)

if (!file.exists(results_path)) {
  stop("Missing per-MAG LMM result table: ", results_path, call. = FALSE)
}
if (!file.exists(taxonomy_path)) {
  stop("Missing GTDB-Tk taxonomy table: ", taxonomy_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
taxonomy <- utils::read.csv(taxonomy_path, check.names = FALSE, stringsAsFactors = FALSE)

required_result_columns <- c(
  "MAG", "status", "primary_analysis", "time_estimate", "time_p_value",
  "time_p_adj_bh_primary", "primary_significant", "zero_fraction",
  "nonzero_sample_count", "subjects_with_nonzero"
)
missing_result_columns <- setdiff(required_result_columns, names(results))
if (length(missing_result_columns) > 0) {
  stop("Missing required result columns: ", paste(missing_result_columns, collapse = ", "), call. = FALSE)
}

required_taxonomy_columns <- c("Genome", "label", "domain", "phylum", "class", "ordeer", "family", "genus", "species")
missing_taxonomy_columns <- setdiff(required_taxonomy_columns, names(taxonomy))
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

taxonomy_clean <- taxonomy |>
  dplyr::transmute(
    MAG = Genome,
    gtdb_label = label,
    domain = clean_taxon(domain, "d"),
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
    result_class = dplyr::case_when(
      primary_significant %in% TRUE ~ "BH-FDR q < 0.05",
      primary_analysis %in% TRUE & time_p_value < 0.05 ~ "Nominal raw p < 0.05",
      TRUE ~ "Other modeled MAG"
    )
  )

unmatched_mags <- sum(is.na(annotated_results$phylum))
if (unmatched_mags > 0) {
  warning("Some MAGs did not match data/GDTB.csv: ", unmatched_mags, call. = FALSE)
}

significant_taxonomy <- annotated_results |>
  dplyr::filter(primary_significant %in% TRUE)

nominal_taxonomy <- annotated_results |>
  dplyr::filter(primary_analysis %in% TRUE, time_p_value < 0.05) |>
  dplyr::arrange(time_p_adj_bh_primary, time_p_value)

if (nrow(nominal_taxonomy) == 0) {
  stop("No nominal time-associated MAGs available for Figure 6.", call. = FALSE)
}

top_families <- nominal_taxonomy |>
  dplyr::count(family, sort = TRUE) |>
  dplyr::slice_head(n = 12) |>
  dplyr::pull(family)

family_counts <- nominal_taxonomy |>
  dplyr::filter(family %in% top_families) |>
  dplyr::count(family, direction, name = "n") |>
  dplyr::group_by(family) |>
  dplyr::mutate(total_n = sum(n)) |>
  dplyr::ungroup() |>
  dplyr::mutate(family = factor(family, levels = rev(top_families)))

family_summary <- nominal_taxonomy |>
  dplyr::filter(family %in% top_families) |>
  dplyr::group_by(family) |>
  dplyr::summarise(
    n_mags = dplyr::n(),
    n_increasing = sum(time_estimate > 0, na.rm = TRUE),
    n_decreasing = sum(time_estimate < 0, na.rm = TRUE),
    median_time_estimate = stats::median(time_estimate, na.rm = TRUE),
    best_q = min(time_p_adj_bh_primary, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(family = factor(family, levels = rev(top_families)))

top_heatmap_n <- min(24L, nrow(nominal_taxonomy))
top_mags <- nominal_taxonomy |>
  dplyr::slice_head(n = top_heatmap_n) |>
  dplyr::mutate(
    mag_rank = dplyr::row_number(),
    mag_label = paste0(sprintf("%02d", mag_rank), "  ", ifelse(nchar(MAG) > 30, paste0(substr(MAG, 1, 27), "..."), MAG)),
    direction = factor(direction, levels = c("Increasing", "Decreasing", "No directional estimate"))
  )

taxonomy_levels <- c("phylum", "class", "order", "family", "genus")
taxonomy_heatmap <- do.call(
  rbind,
  lapply(taxonomy_levels, function(level) {
    data.frame(
      mag_label = top_mags$mag_label,
      rank_level = level,
      taxon = top_mags[[level]],
      direction = top_mags$direction,
      time_estimate = top_mags$time_estimate,
      time_p_adj_bh_primary = top_mags$time_p_adj_bh_primary,
      stringsAsFactors = FALSE
    )
  })
)
taxonomy_heatmap$mag_label <- factor(taxonomy_heatmap$mag_label, levels = rev(top_mags$mag_label))
taxonomy_heatmap$rank_level <- factor(taxonomy_heatmap$rank_level, levels = taxonomy_levels)
taxonomy_heatmap$taxon_short <- ifelse(nchar(taxonomy_heatmap$taxon) > 20, paste0(substr(taxonomy_heatmap$taxon, 1, 17), "..."), taxonomy_heatmap$taxon)

total_mags <- nrow(results)
nominal_mags <- nrow(nominal_taxonomy)
fdr_mags <- nrow(significant_taxonomy)
matched_mags <- sum(!is.na(annotated_results$phylum))
top_phylum <- names(sort(table(nominal_taxonomy$phylum), decreasing = TRUE))[1]
top_family <- names(sort(table(nominal_taxonomy$family), decreasing = TRUE))[1]

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)

pal <- list(
  ink = "#111827",
  text = "#374151",
  muted = "#6B7280",
  hair = "#D1D5DB",
  faint = "#F8FAFC",
  blue = "#0072B2",
  sky = "#DCEEF8",
  orange = "#D55E00",
  peach = "#FFF2E8",
  red = "#CC3311",
  green = "#009E73",
  purple = "#785EF0"
)

direction_palette <- c(
  "Increasing" = pal$blue,
  "Decreasing" = pal$red,
  "No directional estimate" = pal$muted
)

base_theme <- ggplot2::theme_minimal(base_size = 9.4, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11.0, color = pal$ink),
    plot.subtitle = ggplot2::element_text(size = 8.0, color = pal$muted, margin = ggplot2::margin(b = 5)),
    axis.title = ggplot2::element_text(size = 8.4, color = pal$text),
    axis.text = ggplot2::element_text(size = 7.4, color = pal$text),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.22),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 8),
    legend.text = ggplot2::element_text(size = 7.4),
    plot.margin = ggplot2::margin(5, 7, 5, 7)
  )

panel_a <- ggplot2::ggplot(family_counts, ggplot2::aes(x = family, y = n, fill = direction)) +
  ggplot2::geom_col(width = 0.78, color = "white", linewidth = 0.25) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = direction_palette, name = "Direction") +
  ggplot2::labs(
    title = "A  Family distribution of nominal time-associated MAGs",
    subtitle = paste0("Top families among ", nominal_mags, " exploratory nominal candidates; ", fdr_mags, " BH-FDR discoveries overall"),
    x = NULL,
    y = "MAG count"
  ) +
  base_theme +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

panel_b <- ggplot2::ggplot(
  family_summary,
  ggplot2::aes(x = n_mags, y = family, size = n_mags, color = median_time_estimate)
) +
  ggplot2::geom_vline(xintercept = 0, color = pal$hair, linewidth = 0.35) +
  ggplot2::geom_point(alpha = 0.90) +
  ggplot2::scale_color_gradient2(
    low = pal$red,
    mid = "#F3F4F6",
    high = pal$blue,
    midpoint = 0,
    name = "Median beta"
  ) +
  ggplot2::scale_size_continuous(range = c(2.4, 8.0), name = "MAG count") +
  ggplot2::labs(
    title = "B  Family-level direction and effect summary",
    subtitle = "Bubble color shows median time-effect estimate on log10(Y + pseudocount) scale",
    x = "Nominal candidate MAGs per family",
    y = NULL
  ) +
  base_theme +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

panel_c <- ggplot2::ggplot(taxonomy_heatmap, ggplot2::aes(x = rank_level, y = mag_label)) +
  ggplot2::geom_tile(fill = pal$faint, color = "white", linewidth = 0.45) +
  ggplot2::geom_text(ggplot2::aes(label = taxon_short), size = 2.2, color = pal$text) +
  ggplot2::geom_point(
    data = top_mags,
    ggplot2::aes(x = 5.75, y = mag_label, color = direction, size = abs(time_estimate)),
    inherit.aes = FALSE,
    alpha = 0.9
  ) +
  ggplot2::scale_color_manual(values = direction_palette, name = "Direction") +
  ggplot2::scale_size_continuous(range = c(1.8, 4.0), name = "|beta|") +
  ggplot2::coord_cartesian(xlim = c(0.5, 6.05), clip = "off") +
  ggplot2::labs(
    title = "C  Taxonomic annotation of top nominal MAGs",
    subtitle = "Rows are top raw-p nominal MAGs; point at right shows direction and effect magnitude",
    x = "GTDB-Tk taxonomy rank",
    y = NULL
  ) +
  base_theme +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 5.8, family = "mono"),
    panel.grid = ggplot2::element_blank()
  )

summary_text <- paste0(
  "Taxonomy source: data/GDTB.csv joined by Genome -> MAG; matched ",
  fmt_int(matched_mags), " of ", fmt_int(total_mags), " result MAGs. ",
  "BH-FDR significant MAGs = ", fmt_int(fdr_mags),
  ", so plotted taxa are exploratory nominal raw-p candidates only. Dominant nominal phylum: ",
  top_phylum, "; top family: ", top_family,
  ". Do not infer function from taxonomy alone."
)

draw_figure <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  layout <- grid::grid.layout(
    nrow = 4,
    ncol = 2,
    heights = grid::unit(c(0.50, 1.0, 1.35, 0.40), c("in", "null", "null", "in")),
    widths = grid::unit(c(1.0, 1.05), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))

  grid::grid.text(
    "Figure 6 | Taxonomic identity of time-associated MAGs",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.88, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 15.5, fontface = "bold", col = pal$ink),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )
  grid::grid.text(
    "Exploratory longitudinal MAG candidates map to GTDB-Tk taxonomic groups, but taxonomy is descriptive and not functional evidence.",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.20, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 9.0, col = pal$text),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )

  print(panel_a, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panel_b, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(panel_c, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1:2))

  grid::grid.text(
    paste(strwrap(summary_text, width = 175), collapse = "\n"),
    x = grid::unit(0, "npc"),
    y = grid::unit(0.58, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 7.2, col = pal$muted, lineheight = 0.95),
    vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1:2)
  )

  grid::popViewport()
}

png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

png(filename = png_path, width = 13.0, height = 10.2, units = "in", res = 420, bg = "white")
draw_figure()
dev.off()

pdf(file = pdf_path, width = 13.0, height = 10.2, bg = "white", useDingbats = FALSE)
draw_figure()
dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = svg_path, width = 13.0, height = 10.2, bg = "white")
  draw_figure()
  dev.off()
}

significant_taxonomy_path <- file.path(output_dir, paste0(output_stem, "_significant_MAG_taxonomy_table.csv"))
nominal_taxonomy_path <- file.path(output_dir, paste0(output_stem, "_exploratory_nominal_MAG_taxonomy_table.csv"))
family_summary_path <- file.path(output_dir, paste0(output_stem, "_family_taxonomy_summary.csv"))

taxonomy_output_cols <- c(
  "MAG", "gtdb_label", "domain", "phylum", "class", "order", "family", "genus", "species",
  "direction", "time_estimate", "time_p_value", "time_p_adj_bh_primary", "primary_significant",
  "zero_fraction", "nonzero_sample_count", "subjects_with_nonzero"
)

utils::write.csv(significant_taxonomy[, taxonomy_output_cols, drop = FALSE], significant_taxonomy_path, row.names = FALSE)
utils::write.csv(nominal_taxonomy[, taxonomy_output_cols, drop = FALSE], nominal_taxonomy_path, row.names = FALSE)
utils::write.csv(family_summary, family_summary_path, row.names = FALSE)

message("Figure 6 taxonomy outputs created:")
message(" - ", png_path)
message(" - ", pdf_path)
if (file.exists(svg_path)) message(" - ", svg_path)
message(" - ", significant_taxonomy_path)
message(" - ", nominal_taxonomy_path)
message(" - ", family_summary_path)
