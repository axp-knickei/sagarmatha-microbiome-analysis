#!/usr/bin/env Rscript

# Purpose: Create numbered Nature-oriented figures, publication data tables,
# and a chart contract map from completed functional analyses.
# Inputs: Stage-01 through stage-03 result tables and intermediates.
# Outputs: Figures 01-06 in PDF/PNG/SVG plus chart-map and plot-data tables.
# Dependencies: R, data.table, ggplot2, ggrepel, gridExtra, scales, svglite.
# Usage: Rscript analysis/mag_functional_annotation/04_make_publication_outputs.R --config config/mag_functional_annotation.yml
# Workflow: Publication rendering only; never refits models or changes raw data.

source("analysis/mag_functional_annotation/R/common.R")
args <- parse_stage_args()
cfg <- read_magfunc_config(args$config)
assert_packages(c("ggrepel", "gridExtra"))
paths <- output_paths(cfg)
if (!file.exists(file.path(paths$qc, "03_longitudinal_analyses.complete"))) stop("Stage 03 completion marker is missing.", call. = FALSE)
guard_stage(file.path(paths$qc, "04_make_publication_outputs.complete"), args, cfg)

figure_registry <- list()
register_figure <- function(number, stem, category, question, description, files) {
  figure_registry[[length(figure_registry) + 1L]] <<- data.table::data.table(
    result_number = paste("Figure", sprintf("%02d", number)),
    category = category,
    question = question,
    artifact_type = "figure",
    path = vapply(files, relative_to_root, character(1)),
    description = description,
    stage = "04_make_publication_outputs"
  )
}

cluster_annotation <- read_object(file.path(paths$intermediate, "cluster_annotation_summary.rds"))
field_coverage <- data.table::fread(file.path(paths$qc, "annotation_field_coverage.csv"))

cluster_sizes <- sort(cluster_annotation$cluster_size)
ccdf <- data.table(
  cluster_size = sort(unique(cluster_sizes)),
  ccdf = vapply(sort(unique(cluster_sizes)), function(x) mean(cluster_sizes >= x), numeric(1))
)
compact_percent <- function(x) {
  percentages <- x * 100
  vapply(percentages, function(value) {
    if (!is.finite(value)) return(NA_character_)
    decimals <- if (value >= 1) 0L else max(1L, ceiling(-log10(value)))
    paste0(formatC(value, format = "f", digits = decimals), "%")
  }, character(1))
}
p1a <- ggplot2::ggplot(ccdf, ggplot2::aes(cluster_size, ccdf)) +
  ggplot2::geom_line(colour = magfunc_palette[["blue"]], linewidth = 0.7) +
  ggplot2::scale_x_log10(labels = scales::label_number()) +
  ggplot2::scale_y_log10(breaks = 10^seq(0, -6, by = -1), labels = compact_percent) +
  ggplot2::labs(title = "a  Cluster-size distribution", subtitle = "Log–log complementary cumulative distribution\n(member sequences are not abundance)", x = "Member protein sequences per cluster", y = "Clusters at or above size") +
  theme_nature()
field_labels <- c(
  seed_ortholog = "Seed ortholog", COG_category = "COG category",
  Description = "Description", PFAMs = "Pfam domain", KEGG_ko = "KEGG ortholog",
  Preferred_name = "Preferred gene name", KEGG_Pathway = "KEGG pathway",
  EC = "EC number", KEGG_Module = "KEGG module", CAZy = "CAZy family"
)
field_coverage[, field_label := unname(field_labels[field])]
field_coverage[is.na(field_label), field_label := gsub("_", " ", field)]
field_coverage[, field_label := factor(field_label, levels = rev(field_label[order(coverage)]))]
p1b <- ggplot2::ggplot(field_coverage, ggplot2::aes(coverage, field_label)) +
  ggplot2::geom_col(fill = magfunc_palette[["gold"]], width = 0.72) +
  ggplot2::geom_text(ggplot2::aes(label = scales::percent(coverage, accuracy = 1)), hjust = 0, nudge_x = 0.012, size = 2.4, colour = magfunc_palette[["ink"]]) +
  ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1.16), breaks = seq(0, 1, 0.25), expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::labs(title = "b  Annotation-field coverage", subtitle = "Field percentages overlap", x = "Representative queries with field", y = NULL) +
  theme_nature()
figure1 <- gridExtra::arrangeGrob(p1a, p1b, ncol = 2, widths = c(1.05, 1))
files1 <- save_publication_plot(figure1, "Figure_01_catalogue_QC_and_annotation_coverage", cfg, width = 7.4, height = 4.1)
register_figure(1, "Figure_01_catalogue_QC_and_annotation_coverage", "quality_control", "What is the catalogue size structure and annotation-field coverage?", "Cluster-size CCDF and representative-level annotation coverage.", files1)

convergence <- read_table(file.path(paths$convergence, "Table_03_functional_convergence_all_ontologies.csv.gz"))
plot_ontologies <- intersect(c("COG", "KO", "KEGG_module", "CAZy"), unique(convergence$ontology))
convergence_plot <- convergence[ontology %in% plot_ontologies & cluster_richness > 0 & member_gene_representation > 0]
convergence_plot[, label := ""]
convergence_plot[order(-cluster_richness), label := ifelse(seq_len(.N) <= 3L, term, ""), by = ontology]
p2 <- ggplot2::ggplot(convergence_plot, ggplot2::aes(cluster_richness, member_gene_representation)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, colour = magfunc_palette[["muted"]], linetype = 2, linewidth = 0.45) +
  ggplot2::geom_point(shape = 21, fill = magfunc_palette[["blue"]], colour = magfunc_palette[["ink"]], stroke = 0.25, alpha = 0.72, size = 1.5) +
  ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 2.4, max.overlaps = 20, min.segment.length = 0, colour = magfunc_palette[["ink"]]) +
  ggplot2::scale_x_log10(labels = scales::label_number()) +
  ggplot2::scale_y_log10(labels = scales::label_number()) +
  ggplot2::facet_wrap(~ontology, scales = "free", ncol = 2) +
  ggplot2::labs(title = "Functional representation across sequence clusters", subtitle = "Dashed line: one member sequence per cluster; distance above it indicates greater sequence representation\nNeither axis is ecological abundance", x = "Distinct representative clusters", y = "Member protein sequences") +
  theme_nature()
files2 <- save_publication_plot(p2, "Figure_02_functional_convergence", cfg, width = 7.2, height = 6.2)
register_figure(2, "Figure_02_functional_convergence", "functional_convergence", "How many sequence families support each predicted function?", "Functional convergence plane for COG, KO, KEGG modules, and CAZy.", files2)

unknown_summary <- data.table::fread(file.path(paths$unknown, "Table_04_unknown_annotation_by_cluster_size.csv"))
unknown_summary <- unknown_summary[!is.na(size_group) & clusters > 0]
unknown_summary[, size_group := factor(as.character(size_group), levels = c("1", "2-5", "6-20", "21-100", ">100"))]
unknown_summary[, size_label := paste0(
  as.character(size_group), "\n(n=", scales::label_number(scale_cut = scales::cut_short_scale())(clusters), ")"
)]
prediction_plot <- data.table::melt(
  unknown_summary,
  id.vars = c("size_group", "size_label", "clusters"),
  measure.vars = c("annotation_row_fraction", "unknown_or_unannotated_fraction"),
  variable.name = "model",
  value.name = "probability"
)
prediction_plot[, size_label := factor(size_label, levels = unique(unknown_summary[order(size_group), size_label]))]
prediction_plot[, label_vjust := ifelse(model == "annotation_row_fraction", -0.9, 1.5)]
prediction_colours <- c(annotation_row_fraction = "#3B6FB6", unknown_or_unannotated_fraction = "#D96C3F")
prediction_labels <- c(annotation_row_fraction = "Annotation record", unknown_or_unannotated_fraction = "Unknown or unannotated")
write_table(prediction_plot, file.path(paths$reproducibility, "Figure_03_plot_data.csv"))
p3 <- ggplot2::ggplot(prediction_plot, ggplot2::aes(size_label, probability, colour = model, group = model))
if (data.table::uniqueN(prediction_plot$size_label) > 1L) {
  p3 <- p3 + ggplot2::geom_line(linewidth = 0.8)
}
p3 <- p3 +
  ggplot2::geom_point(shape = 21, fill = "white", size = 2.2, stroke = 0.75) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::percent(probability, accuracy = 0.1), vjust = label_vjust),
    size = 2.5,
    show.legend = FALSE
  ) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1.06), breaks = seq(0, 1, 0.25), expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_colour_manual(values = prediction_colours, breaks = names(prediction_colours), labels = prediction_labels) +
  ggplot2::labs(title = "Annotation status by protein-cluster size group", subtitle = "Observed catalogue fractions; cluster counts are descriptive and are not biological replicates", x = "Member protein sequences per cluster", y = "Observed fraction of clusters", colour = NULL) +
  theme_nature() + ggplot2::theme(legend.position = "top")
files3 <- save_publication_plot(p3, "Figure_03_unknown_annotation_by_cluster_size", cfg)
register_figure(3, "Figure_03_unknown_annotation_by_cluster_size", "unknown_by_cluster_size", "Are small protein clusters more often unknown?", "Observed annotation-record and unknown/unannotated fractions across prespecified protein-cluster size groups.", files3)

taxonomy <- data.table::fread(resolve_input_paths(cfg)$taxonomy_file)
if ("ordeer" %in% names(taxonomy) && !"order" %in% names(taxonomy)) taxonomy[, order := ordeer]
tax_id <- cfg$columns$taxonomy_mag_id
gut_mags <- read_object(file.path(paths$intermediate, "gut_mag_ids.rds"))
phylum_denominator <- taxonomy[get(tax_id) %in% gut_mags, .(phylum_mags = .N), by = phylum]
cog_tax <- read_table(file.path(paths$taxonomy, "taxonomy_function_distribution_COG.csv.gz"))
cog_phylum <- cog_tax[, .(carrier_mags = sum(carrier_mags)), by = .(term, phylum)]
cog_terms <- sort(unique(convergence[ontology == "COG", term]))
cog_phylum <- merge(data.table::CJ(term = cog_terms, phylum = phylum_denominator$phylum, unique = TRUE), cog_phylum, by = c("term", "phylum"), all.x = TRUE)
cog_phylum[is.na(carrier_mags), carrier_mags := 0L]
cog_phylum <- merge(cog_phylum, phylum_denominator, by = "phylum", all.x = TRUE)
cog_phylum[, carrier_fraction := carrier_mags / phylum_mags]
cog_phylum[, phylum_label := paste0(phylum, " (n=", phylum_mags, ")")]
phylum_levels <- rev(phylum_denominator[order(-phylum_mags), paste0(phylum, " (n=", phylum_mags, ")")])
cog_phylum[, phylum_label := factor(phylum_label, levels = phylum_levels)]
write_table(cog_phylum, file.path(paths$reproducibility, "Figure_04_plot_data.csv"))
p4 <- ggplot2::ggplot(cog_phylum, ggplot2::aes(term, phylum_label, fill = carrier_fraction)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
  ggplot2::scale_fill_gradient(low = "#F4F6F8", high = magfunc_palette[["blue"]], labels = scales::label_percent()) +
  ggplot2::labs(title = "COG functional distribution across gut MAG phyla", subtitle = "Fraction of selected gut MAGs carrying at least one cluster; phylum MAG counts shown in labels", x = "COG category", y = NULL, fill = "Carrier MAGs") +
  theme_nature() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0), legend.position = "top")
files4 <- save_publication_plot(p4, "Figure_04_taxonomic_COG_distribution", cfg, width = 7.2, height = 4.4)
register_figure(4, "Figure_04_taxonomic_COG_distribution", "taxonomic_distribution", "How are broad functional categories distributed across gut MAG phyla?", "Phylum-level COG carrier fractions using exact MAG-to-taxonomy joins.", files4)

models <- read_table(file.path(paths$longitudinal, "Table_06_longitudinal_functional_models.csv.gz"))
summaries <- read_table(file.path(paths$longitudinal, "functional_potential_timepoint_summaries.csv.gz"))
eligible_models <- models[status == "fitted" & is.finite(omnibus_p_value) & singular %in% FALSE][order(omnibus_p_adj_bh, omnibus_p_value)]
profile_rows <- merge(
  summaries,
  eligible_models[, .(ontology, term)],
  by = c("ontology", "term")
)
profile_rows[, timepoint_order := match(as.character(Timepoint), unlist(cfg$analysis$ordered_timepoints))]
profile_signatures <- profile_rows[order(timepoint_order), .(
  profile_signature = paste(sprintf("%.12g", mean_potential), collapse = "|")
), by = .(ontology, term)]
ranked_models <- merge(eligible_models, profile_signatures, by = c("ontology", "term"), all.x = TRUE)
data.table::setorder(ranked_models, omnibus_p_adj_bh, omnibus_p_value, ontology, term, na.last = TRUE)
ranked_models <- ranked_models[!duplicated(profile_signature)]
top_n <- min(6L, nrow(ranked_models))
if (top_n > 0) {
  selected <- ranked_models[seq_len(top_n), .(ontology, term, omnibus_p_adj_bh, significant_omnibus, singular, profile_signature)]
  write_table(selected, file.path(paths$reproducibility, "Figure_05_selected_models.csv"))
  plot_long <- merge(summaries, selected, by = c("ontology", "term"))
  plot_long[, significance_label := ifelse(significant_omnibus, "BH significant", "Not BH significant")]
  plot_long[, feature := paste0(
    ontology, ": ", term, "  [",
    ifelse(significant_omnibus, "BH q=", "not significant; q="),
    formatC(omnibus_p_adj_bh, digits = 2L, format = "g"), "]"
  )]
  plot_long[, Timepoint := factor(Timepoint, levels = unlist(cfg$analysis$ordered_timepoints))]
  longitudinal_colours <- c("BH significant" = "#D96C3F", "Not BH significant" = "#3B6FB6")
  p5 <- ggplot2::ggplot(plot_long, ggplot2::aes(Timepoint, mean_potential, group = 1, colour = significance_label)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(0, mean_potential - 1.96 * standard_error), ymax = mean_potential + 1.96 * standard_error), width = 0.18, alpha = 0.25, linewidth = 0.65) +
    ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::geom_point(shape = 21, fill = "white", size = 1.8, stroke = 0.6) +
    ggplot2::scale_colour_manual(values = longitudinal_colours) +
    ggplot2::facet_wrap(~feature, scales = "free_y", ncol = 2) +
    ggplot2::labs(title = "Top non-redundant longitudinal predicted functional potentials", subtitle = "Mean ± 1.96 SE (descriptive); orange functions pass ontology-specific BH correction\nFree y-scales: compare trajectory shape, not panel height; participant-aware models; Ta/Tb excluded", x = "Ordered expedition timepoint", y = "Predicted functional potential", colour = NULL) +
    theme_nature() + ggplot2::theme(legend.position = "top")
} else {
  p5 <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "No longitudinal models passed the fitted-model gate") + ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) + ggplot2::labs(title = "Longitudinal predicted functional potential") + ggplot2::theme_void(base_family = "sans")
}
files5 <- save_publication_plot(p5, "Figure_05_longitudinal_functional_potential", cfg, width = 7.4, height = 6.4)
register_figure(5, "Figure_05_longitudinal_functional_potential", "longitudinal_functional_potential", "Which predicted functional potentials vary most across T1-T10?", "Top-ranked model-supported trajectories; ranking does not imply significance when adjusted p-values exceed the threshold.", files5)

hysteresis <- read_table(file.path(paths$hysteresis, "Table_07_equal_altitude_paired_contrasts.csv.gz"))
contrast_levels <- vapply(cfg$hysteresis_contrasts, function(x) as.character(x$id), character(1))
hyst_candidates <- hysteresis[status == "fitted" & is.finite(p_value)]
hyst_candidates[, contrast_id := factor(contrast_id, levels = contrast_levels)]
data.table::setorder(hyst_candidates, contrast_id, p_adj_bh, p_value, ontology, term, na.last = TRUE)
hyst_plot <- hyst_candidates[, head(.SD, 6L), by = contrast_id]
hyst_plot[, significant_bh := is.finite(p_adj_bh) & p_adj_bh <= as.numeric(cfg$analysis$fdr_threshold)]
write_table(hyst_plot, file.path(paths$reproducibility, "Figure_06_plot_data.csv"))
if (nrow(hyst_plot) > 0) {
  hyst_plot[, feature := paste0(ontology, ": ", term)]
  hyst_plot[, feature := factor(feature, levels = rev(unique(feature)))]
  contrast_display <- c(
    T3_vs_T6_ascent_post_high_climb = "T3 vs T6: ascent vs post-high climb",
    T4_vs_T6_pre_post_6500m = "T4 vs T6: before vs after 6500 m",
    T7_vs_T8_pre_post_extreme_exposure = "T7 vs T8: pre/post extreme exposure"
  )
  hyst_plot[, contrast_label := paste0(unname(contrast_display[as.character(contrast_id)]), " (n=", n_pairs, " pairs)")]
  contrast_label_levels <- vapply(contrast_levels, function(id) {
    rows <- hyst_plot[as.character(contrast_id) == id]
    paste0(unname(contrast_display[[id]]), " (n=", rows$n_pairs[[1]], " pairs)")
  }, character(1))
  hyst_plot[, contrast_label := factor(contrast_label, levels = contrast_label_levels)]
  significance_note <- if (any(hyst_plot$significant_bh)) {
    "Filled points pass BH correction"
  } else {
    "No paired contrast passes BH correction"
  }
  p6 <- ggplot2::ggplot(hyst_plot, ggplot2::aes(mean_log10_difference, feature, colour = primary_contrast, shape = significant_bh)) +
    ggplot2::geom_vline(xintercept = 0, colour = magfunc_palette[["muted"]], linetype = 2, linewidth = 0.45) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower_95, xmax = upper_95), orientation = "y", width = 0, linewidth = 0.5) +
    ggplot2::geom_point(size = 1.9) +
    ggplot2::facet_wrap(~contrast_label, scales = "free_y", ncol = 1) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = magfunc_palette[["orange"]], `FALSE` = magfunc_palette[["blue"]]), labels = c(`TRUE` = "Primary", `FALSE` = "Supporting")) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "BH significant", `FALSE` = "Not BH significant")) +
    ggplot2::labs(title = "Equal-recorded-altitude paired functional contrasts", subtitle = paste0("Six top-ranked effects per contrast; within-participant 95% confidence intervals\n", significance_note, "; hysteresis-like, not a complete loop"), x = "After minus before (log10 predicted potential)", y = NULL, colour = NULL, shape = NULL) +
    theme_nature() + ggplot2::theme(legend.position = "top")
} else {
  p6 <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "No paired contrasts passed the fitted-model gate") + ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) + ggplot2::labs(title = "Equal-recorded-altitude paired contrasts") + ggplot2::theme_void(base_family = "sans")
}
files6 <- save_publication_plot(p6, "Figure_06_equal_altitude_paired_contrasts", cfg, width = 7.4, height = 8.2)
register_figure(6, "Figure_06_equal_altitude_paired_contrasts", "ascent_descent_hysteresis", "Do predicted functional potentials differ between prespecified equal-altitude phases?", "Paired effect-size forest plot showing six top-ranked functions per primary and supporting contrast, with BH significance encoded explicitly.", files6)

chart_map <- data.table::data.table(
  figure = sprintf("Figure %02d", 1:6),
  analytical_question = c(
    "Catalogue size structure and annotation-field coverage",
    "Functional convergence across sequence clusters",
    "Unknown annotation probability across cluster size",
    "COG carrier distribution across gut MAG phyla",
    "Longitudinal predicted functional potential across T1-T10",
    "Equal-altitude within-participant paired contrasts"
  ),
  family = c("distribution + comparison", "relationship", "relationship + uncertainty", "matrix", "trend", "uncertainty + comparison"),
  chart_type = c("CCDF + horizontal bar", "faceted scatter", "probability curve", "heatmap", "faceted line", "faceted forest plot"),
  palette_policy = c("hard two-root cap", "single-root preferred", "hard two-root cap", "single-root preferred", "single-root preferred", "hard two-root cap"),
  denominator_or_unit = c(
    "unique protein clusters and representative queries",
    "distinct clusters and member protein sequences",
    "unique clusters",
    "selected gut MAGs within phylum",
    "sum of gut MAG values carrying a function; no row renormalization",
    "paired participants; log10 predicted functional potential"
  )
)
write_table(chart_map, file.path(paths$reproducibility, "publication_chart_map.csv"))
register_results(cfg, data.table::rbindlist(figure_registry))
writeLines("Stage 04 complete", file.path(paths$qc, "04_make_publication_outputs.complete"))
write_session_info(cfg, "04_make_publication_outputs")
message("Stage 04 completed.")
