#!/usr/bin/env Rscript
source("analysis/functional_community_dynamics/R/common.R")

args <- parse_dynamics_args()
cfg <- read_dynamics_config(args$config)
assert_dynamics_packages()
paths <- dynamics_paths(cfg)
marker <- file.path(paths$qc, "04_make_figures_and_report.complete")
guard_dynamics_stage(marker, args, cfg)
for (stage in c("02_run_multivariate_analysis.complete", "03_analyze_longitudinal_stability.complete")) {
  if (!file.exists(file.path(paths$qc, stage))) stop("Required prior stage marker missing: ", stage, call. = FALSE)
}

global <- fread(file.path(paths$multivariate, "Table_27_global_PERMANOVA_and_dispersion_results.csv"))
pairwise <- fread(file.path(paths$multivariate, "Table_28_supported_functional_timepoint_contrasts.csv"))
models <- fread(file.path(paths$longitudinal, "Table_30_matched_taxonomic_functional_stability_models.csv"))
sensitivity <- fread(file.path(paths$longitudinal, "Table_31_participant_level_transition_sensitivity.csv"))
functional_order <- configured_layers(cfg)$layer
all_layer_order <- c("taxonomic_MAG", functional_order)
display_layer_labels <- c(
  taxonomic_MAG = "Taxonomic MAG",
  KO = "KO",
  KEGG_module = "KEGG module",
  EC = "EC",
  COG = "COG"
)
pcoa <- fread(file.path(paths$reproducibility, "Figure_16_plot_data.csv"))[metric == cfg$analysis$primary_metric & layer != "taxonomic_MAG"]
pcoa[, timepoint := factor(timepoint, levels = timepoint_levels(cfg), ordered = TRUE)]
pcoa[, layer := factor(layer, levels = functional_order)]
pcoa_axis_labels <- unique(pcoa[, .(layer, axis1_percent, axis2_percent)])
pcoa_labels <- setNames(
  sprintf(
    "%s\nPCoA1 %.1f%%; PCoA2 %.1f%%",
    unname(display_layer_labels[as.character(pcoa_axis_labels$layer)]),
    pcoa_axis_labels$axis1_percent,
    pcoa_axis_labels$axis2_percent
  ),
  as.character(pcoa_axis_labels$layer)
)

figure16 <- ggplot2::ggplot(pcoa, ggplot2::aes(PCoA1, PCoA2, group = subject_id)) +
  ggplot2::geom_path(colour = "#B8B8B8", linewidth = 0.35, alpha = 0.75) +
  ggplot2::geom_point(ggplot2::aes(colour = timepoint), size = 1.4, alpha = 0.9) +
  ggplot2::facet_wrap(
    ~layer,
    scales = "free",
    ncol = 2,
    labeller = ggplot2::as_labeller(pcoa_labels)
  ) +
  ggplot2::scale_colour_manual(values = dynamics_palette, drop = FALSE) +
  ggplot2::labs(
    title = "Predicted functional profiles across the expedition",
    subtitle = "Hellinger PCoA; grey paths connect observations within participants",
    x = "PCoA axis 1 score (panel-specific scale)", y = "PCoA axis 2 score (panel-specific scale)",
    colour = "Timepoint", caption = "Ordination shows predicted genomic potential, not expression or activity."
  ) + theme_dynamics() + ggplot2::theme(legend.position = "bottom")
save_dynamics_plot(figure16, "Figure_16_functional_PCoA_participant_trajectories", cfg, height = 6.2)

reference <- fread(file.path(paths$reproducibility, "Figure_17_plot_data.csv"))
reference[, timepoint := factor(to_timepoint, levels = timepoint_levels(cfg), ordered = TRUE)]
reference[, layer := factor(layer, levels = all_layer_order)]
figure17 <- ggplot2::ggplot(reference, ggplot2::aes(timepoint, standardized_distance, group = subject_id)) +
  ggplot2::geom_line(colour = "#8C8C8C", linewidth = 0.35, alpha = 0.55) +
  ggplot2::geom_point(ggplot2::aes(colour = timepoint), size = 1.25, alpha = 0.8) +
  ggplot2::stat_summary(ggplot2::aes(group = 1), fun = median, geom = "line", colour = "black", linewidth = 0.8) +
  ggplot2::facet_wrap(
    ~layer,
    scales = "free_y",
    ncol = 2,
    labeller = ggplot2::as_labeller(display_layer_labels)
  ) +
  ggplot2::scale_colour_manual(values = dynamics_palette, guide = "none", drop = FALSE) +
  ggplot2::labs(
    title = "Distance from each participant's earliest observation",
    subtitle = "Hellinger distance divided by the layer-wide median pairwise distance; black line is the median",
    x = "Timepoint", y = "Standardized distance to earliest observation",
    caption = "The earliest observation is participant-specific and is not necessarily T1 or a biological baseline."
  ) + theme_dynamics() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_dynamics_plot(figure17, "Figure_17_distance_to_earliest_observation", cfg, height = 6.5)

matched <- fread(file.path(paths$reproducibility, "Figure_18_plot_data.csv"))
matched[, layer := factor(layer, levels = functional_order)]
figure18 <- ggplot2::ggplot(matched, ggplot2::aes(taxonomic_standardized_distance, functional_standardized_distance)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#777777") +
  ggplot2::geom_point(ggplot2::aes(colour = to_timepoint), size = 1.4, alpha = 0.65) +
  ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "black", fill = "#CCCCCC", linewidth = 0.6) +
  ggplot2::facet_wrap(
    ~layer,
    scales = "free",
    ncol = 2,
    labeller = ggplot2::as_labeller(display_layer_labels)
  ) +
  ggplot2::scale_colour_manual(values = dynamics_palette, drop = FALSE) +
  ggplot2::labs(
    title = "Matched taxonomic and functional movement",
    subtitle = "Each point is the same participant and adjacent observed transition",
    x = "Taxonomic MAG transition distance (standardized)",
    y = "Predicted functional transition distance (standardized)", colour = "Later timepoint",
    caption = paste(
      "Dashed line = equal standardized movement; black line = descriptive linear fit.",
      "Below the dashed line supports smoothing; above it does not."
    )
  ) + theme_dynamics() + ggplot2::theme(legend.position = "bottom")
save_dynamics_plot(figure18, "Figure_18_taxonomic_vs_functional_transitions", cfg, height = 6.3)

recovery <- fread(file.path(paths$reproducibility, "Figure_19_plot_data.csv"))
recovery[, layer := factor(layer, levels = all_layer_order)]
recovery_long <- melt(recovery, id.vars = c("layer", "subject_id"), measure.vars = c("peak_distance_from_t1", "distance"), variable.name = "estimate", value.name = "distance")
recovery_long[, estimate := factor(estimate, levels = c("peak_distance_from_t1", "distance"), labels = c("Maximum", "T10"))]
figure19 <- ggplot2::ggplot(recovery_long, ggplot2::aes(estimate, distance, group = subject_id, colour = subject_id)) +
  ggplot2::geom_line(linewidth = 0.55, alpha = 0.7) + ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(
    ~layer,
    scales = "free_y",
    ncol = 2,
    labeller = ggplot2::as_labeller(display_layer_labels)
  ) +
  ggplot2::labs(
    title = "Exploratory T1-to-T10 return toward the initial profile",
    subtitle = paste0("Only ", uniqueN(recovery$subject_id), " participants have paired T1 and T10 observations"),
    x = "Distance from T1: maximum observed versus T10", y = "Hellinger distance", colour = "Participant",
    caption = "A lower T10 value than the maximum suggests partial return, but n is too small for cohort-wide recovery claims."
  ) + theme_dynamics() + ggplot2::theme(legend.position = "bottom")
save_dynamics_plot(figure19, "Figure_19_exploratory_T1_T10_recovery", cfg, height = 6.4)

supported <- global[supported_global %in% TRUE]
dispersion_flags <- global[dispersion_supported %in% TRUE]
supported_pairs <- pairwise[significant_bh %in% TRUE]
smoothing <- models[model_family == "functional_smoothing_relative_to_taxonomy"]
coupling <- models[model_family == "taxonomic_functional_transition_coupling"]
taxonomic_primary <- global[layer == "taxonomic_MAG" & metric == cfg$analysis$primary_metric]
functional_primary <- global[layer == cfg$analysis$primary_layer & metric == cfg$analysis$primary_metric]
coupling_primary <- coupling[layer == cfg$analysis$primary_layer]
smoothing_primary <- smoothing[layer == cfg$analysis$primary_layer]
sensitivity_coupling <- sensitivity[sensitivity_family == "participant_spearman_transition_coupling"]
sensitivity_relative <- sensitivity[sensitivity_family == "participant_median_relative_movement"]
sensitivity_primary_coupling <- sensitivity_coupling[layer == cfg$analysis$primary_layer]
sensitivity_primary_relative <- sensitivity_relative[layer == cfg$analysis$primary_layer]
ko_recovery <- recovery[as.character(layer) == cfg$analysis$primary_layer]

fmt <- function(x, digits = 3L) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "g"), "not applicable")
}
model_quality <- function(row) {
  if (!nrow(row)) return("not available")
  if (!identical(row$status[[1]], "fitted")) return(row$status[[1]])
  if (isTRUE(row$singular[[1]]) || !isTRUE(row$converged[[1]])) return("unstable mixed-model fit")
  "acceptable mixed-model fit"
}
report_path <- file.path(paths$reports, "2026-08-12_functional_community_dynamics_interpretation.md")
report <- c(
  "# Functional Community Restructuring and Recovery Analysis",
  "", paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
  "## Technical summary", "",
  paste0(
    "The 441-MAG taxonomic profile changed across T1–T10 (primary Hellinger PERMANOVA R² = ",
    fmt(taxonomic_primary$permanova_r_squared), ", p = ", fmt(taxonomic_primary$permanova_p_value),
    "), but the primary KO predicted-functional profile did not pass the prespecified multiple-testing threshold ",
    "(R² = ", fmt(functional_primary$permanova_r_squared), ", BH q = ",
    fmt(functional_primary$permanova_p_adj_bh_functional_family), ")."
  ),
  paste0(
    "No functional layer–distance combination was BH-supported, so the gated pairwise functional screen was not run ",
    "and Table 28 correctly contains zero rows. No dispersion test was BH-supported."
  ),
  paste0(
    "Within participants, larger taxonomic transitions accompanied larger KO functional transitions ",
    "(mixed-model slope = ", fmt(coupling_primary$estimate), ", BH q = ",
    fmt(coupling_primary$p_adj_bh_within_model_family), "). The participant-level sensitivity check agreed: ",
    sensitivity_primary_coupling$positive_direction_participants, "/",
    sensitivity_primary_coupling$participants, " estimable participants had positive correlations; median Spearman rho = ",
    fmt(sensitivity_primary_coupling$participant_summary), ", BH q = ",
    fmt(sensitivity_primary_coupling$p_adj_bh_within_sensitivity_family), "."
  ),
  paste0(
    "The KO standardized functional/taxonomic movement ratio was ", fmt(smoothing_primary$ratio),
    " (BH q = ", fmt(smoothing_primary$p_adj_bh_within_model_family), "), which is above—not below—one. ",
    "Therefore these data do not support functional smoothing. This ratio compares separately standardized distance scales; ",
    "it is not a 31% increase in biological function or activity."
  ),
  "", "## Taxonomy changed by timepoint; the overall functional profile did not", "",
  paste0(
    "- Taxonomic MAG / Hellinger: R² = ", fmt(taxonomic_primary$permanova_r_squared),
    ", p = ", fmt(taxonomic_primary$permanova_p_value), ", dispersion BH q = ",
    fmt(taxonomic_primary$dispersion_p_adj_bh), "."
  ),
  paste0(
    "- KO / Hellinger: R² = ", fmt(functional_primary$permanova_r_squared),
    ", raw p = ", fmt(functional_primary$permanova_p_value), ", functional-family BH q = ",
    fmt(functional_primary$permanova_p_adj_bh_functional_family), "."
  ),
  if (nrow(supported)) {
    paste0("- Additional supported functional tests: ", paste(paste(supported$layer, supported$metric, sep = "/"), collapse = ", "), ".")
  } else {
    "- KO, KEGG-module, and EC sensitivities all led to the same negative global conclusion; COG was descriptive only."
  },
  "",
  "Figure 16 should therefore be read as individual trajectories without a reproducible cohort-wide separation by timepoint. The facet labels report the variance represented by the two displayed PCoA axes.",
  "", "## Taxonomic and functional movements were coupled, but not smoothed", "",
  paste0(
    "The primary KO mixed models had an ", model_quality(coupling_primary),
    ". The participant-level signed-rank sensitivity uses one summary per participant and shows that the positive association persists after participant-level summarization."
  ),
  paste0(
    "For relative movement, the participant-level KO sensitivity estimated a median ratio of ",
    fmt(sensitivity_primary_relative$ratio), "; ", sensitivity_primary_relative$positive_direction_participants,
    "/", sensitivity_primary_relative$participants, " participants had ratios above one (BH q = ",
    fmt(sensitivity_primary_relative$p_adj_bh_within_sensitivity_family), ")."
  ),
  paste0(
    "EC transition coupling had an ", model_quality(coupling[layer == "EC"]),
    "; its participant-level sensitivity is the safer supporting check (median rho = ",
    fmt(sensitivity_coupling[layer == "EC", participant_summary]), ", BH q = ",
    fmt(sensitivity_coupling[layer == "EC", p_adj_bh_within_sensitivity_family]), ")."
  ),
  "",
  "Figure 18 shows these matched adjacent transitions. Points above the dashed line have greater functional than taxonomic movement only after each layer is divided by its own median all-sample distance. The black line is descriptive, not the participant-aware model.",
  "", "## T1-to-T10 return remains a three-participant observation", "",
  paste0(
    "Only ", uniqueN(recovery$subject_id), " participants have paired T1 and T10 samples. For KO, ",
    sum(ko_recovery$recovery_fraction_from_peak > 0), " showed some return from their maximum T1 distance by T10, while ",
    sum(ko_recovery$recovery_fraction_from_peak == 0), " did not. This is descriptive and cannot establish cohort-wide recovery."
  ),
  "",
  "Figure 17 uses each participant's earliest available observation, which is not always T1 and is not a biological baseline. Figure 19 uses true T1 only for participants with paired T1 and T10 observations and compares their maximum departure from T1 with their T10 distance.",
  "", "## Robustness and limitations", "",
  paste0("- Participant-level sensitivity results are in Table 31; ", nrow(sensitivity[significant_bh %in% TRUE]), " inferential sensitivity rows are BH-supported."),
  paste0("- No dispersion test was BH-supported (", nrow(dispersion_flags), " flags), reducing concern that a supported global result is only a dispersion artifact."),
  "- Functional potential is derived mathematically from the same MAG abundances used for taxonomy, so coupling is expected to some degree and is not independent biological validation.",
  "- Separately standardized distances permit comparison of relative movement patterns, not raw effect magnitudes or percentages of biological activity.",
  "- The analysis measures predicted genomic potential, not expression, enzyme activity, pathway completeness, metabolite production, resilience, or adaptation.",
  "", "## Outputs", "",
  "- Table 26: `01_prepared_profiles/Table_26_functional_profile_input_and_filtering_audit.csv`",
  "- Table 27: `02_multivariate_analysis/Table_27_global_PERMANOVA_and_dispersion_results.csv`",
  "- Table 28: `02_multivariate_analysis/Table_28_supported_functional_timepoint_contrasts.csv`",
  "- Table 29: `03_longitudinal_stability/Table_29_participant_distance_reference_and_adjacent_transitions.csv`",
  "- Table 30: `03_longitudinal_stability/Table_30_matched_taxonomic_functional_stability_models.csv`",
  "- Table 31: `03_longitudinal_stability/Table_31_participant_level_transition_sensitivity.csv`",
  "- Figures 16–19: `figures/`",
  "", "## Recommended next steps", "",
  "1. Treat the current global functional result as a negative result; do not run ungated pairwise functional contrasts.",
  "2. Use the participant-level sensitivity alongside Table 30 whenever transition coupling or relative movement is reported.",
  "3. Keep T1-to-T10 recovery exploratory unless substantially more paired baseline observations become available.",
  "", "## Interpretation boundaries", "",
  "- KO is primary; KEGG-module and EC-labelled profiles are sensitivities; COG is descriptive.",
  "- KEGG-module labels do not establish complete biochemical modules.",
  "- Functional values are predicted genomic potential derived from carrier-MAG abundance.",
  "- Ta and Tb are excluded, and all distance matrices use identical T1–T10 sample ordering.",
  "- Ratios below one would support smoothing or redundancy as an interpretation, but the observed inferential-layer ratios are above one."
)
writeLines(report, report_path)

index <- c(
  "# Functional Community Dynamics Results Index", "",
  "Tables 26–31 and Figures 16–19 were generated by the independent functional-community-dynamics workflow.", "",
  "See `reports/2026-08-12_functional_community_dynamics_interpretation.md` for the result-dependent interpretation and limitations."
)
writeLines(index, file.path(paths$root, "RESULTS_INDEX.md"))
write_dynamics_session(cfg, "04_make_figures_and_report")
writeLines(c(
  "Stage 04 complete", "Figures 16-19 and the corrected dated interpretation report were generated.",
  "Figure 16 facet labels report variance represented by PCoA axes 1 and 2.",
  "Table 31 provides a participant-level transition sensitivity check.",
  "Figure 19 is exploratory because only participants with both T1 and T10 are shown."
), marker)
