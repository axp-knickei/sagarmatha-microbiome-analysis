#!/usr/bin/env Rscript
# Purpose: Build a reader-first manuscript figure and two compact main tables
# from externally supplied completed MAG-functional-annotation and functional-community-dynamics
# results. This script does not rerun discovery models or modify source results.
# Inputs: Completed production outputs from the two workflows.
# Outputs: Main Figure 1, Main Tables 1-2, panel data, an index, and a marker.
# Usage: Rscript analysis/manuscript_functional_synthesis/00_build_main_figure_and_tables_public.R \
#   --config config/manuscript_functional_synthesis.yml

source("analysis/manuscript_functional_synthesis/R/common.R")

required_packages <- c("data.table", "ggplot2", "gridExtra", "grid", "yaml", "svglite", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

args <- parse_synthesis_args()
cfg <- read_synthesis_config(args$config)
inputs <- resolve_synthesis_inputs(cfg)
paths <- synthesis_paths(cfg)
marker <- file.path(paths$qc, "00_build_main_figure_and_tables.complete")
guard_synthesis(marker, args, cfg)

mag_root <- inputs$mag_functional_results_dir
dyn_root <- inputs$functional_dynamics_results_dir
files <- c(
  mag_marker = file.path(mag_root, "00_quality_control", "10_audit_candidate_annotations.complete"),
  dyn_marker = file.path(dyn_root, "00_quality_control", "04_make_figures_and_report.complete"),
  catalogue = file.path(mag_root, "00_quality_control", "Table_02_full_catalogue_qc_summary.csv"),
  unknown = file.path(mag_root, "02_unknown_by_cluster_size", "Table_04_unknown_annotation_by_cluster_size.csv"),
  singleton = file.path(mag_root, "02_unknown_by_cluster_size", "singleton_unknown_enrichment_model.csv"),
  longitudinal_models = file.path(mag_root, "04_longitudinal_functional_potential", "Table_06_longitudinal_functional_models.csv.gz"),
  trajectories = file.path(mag_root, "06_targeted_candidate_validation", "candidate_participant_trajectories.csv.gz"),
  carrier_overlap = file.path(mag_root, "reproducibility", "Figure_11_plot_data.csv"),
  component_models = file.path(mag_root, "08_shared_unique_carrier_decomposition", "Table_17_component_temporal_models.csv"),
  shared_contrasts = file.path(mag_root, "09_shared_component_followup", "Table_18_shared_component_timepoint_contrasts.csv"),
  genus_drivers = file.path(mag_root, "09_shared_component_followup", "Table_21_shared_component_genus_drivers.csv.gz"),
  curation = file.path(mag_root, "10_candidate_annotation_audit", "Table_25_candidate_curation_summary.csv"),
  global = file.path(dyn_root, "02_multivariate_analysis", "Table_27_global_PERMANOVA_and_dispersion_results.csv"),
  transitions = file.path(dyn_root, "reproducibility", "Figure_18_plot_data.csv"),
  stability_models = file.path(dyn_root, "03_longitudinal_stability", "Table_30_matched_taxonomic_functional_stability_models.csv"),
  participant_sensitivity = file.path(dyn_root, "03_longitudinal_stability", "Table_31_participant_level_transition_sensitivity.csv"),
  recovery = file.path(dyn_root, "reproducibility", "Figure_19_plot_data.csv")
)
require_files(files)

catalogue <- read_result_table(files[["catalogue"]])
unknown <- read_result_table(files[["unknown"]])
singleton <- read_result_table(files[["singleton"]])
longitudinal_models <- read_result_table(files[["longitudinal_models"]])
trajectories <- read_result_table(files[["trajectories"]])
carrier_overlap <- read_result_table(files[["carrier_overlap"]])
component_models <- read_result_table(files[["component_models"]])
shared_contrasts <- read_result_table(files[["shared_contrasts"]])
genus_drivers <- read_result_table(files[["genus_drivers"]])
curation <- read_result_table(files[["curation"]])
global <- read_result_table(files[["global"]])
transitions <- read_result_table(files[["transitions"]])
stability_models <- read_result_table(files[["stability_models"]])
participant_sensitivity <- read_result_table(files[["participant_sensitivity"]])
recovery <- read_result_table(files[["recovery"]])

timepoint_levels <- unlist(cfg$analysis$ordered_timepoints, use.names = FALSE)
candidate_ids <- unlist(cfg$analysis$selected_candidates, use.names = FALSE)
if (!setequal(unique(trajectories$candidate_id), candidate_ids)) {
  stop("Candidate trajectory identifiers do not match the configured candidates.", call. = FALSE)
}
# Reconcile aggregate counts against the supplied carrier and protein evidence.
concordance <- read_result_table(file.path(mag_root, "10_candidate_annotation_audit", "Table_24_candidate_carrier_MAG_concordance.csv"))
proteins <- read_result_table(file.path(mag_root, "10_candidate_annotation_audit", "Table_23_candidate_representative_protein_clusters.csv"))
if (!all(c("MAG", "candidate_1_carried", "candidate_2_carried", "shared_carrier") %in% names(concordance)) ||
    !all(c("representative_id", "representative_mag", "carries_both_candidate_terms_same_representative") %in% names(proteins)) ||
    anyNA(proteins$carries_both_candidate_terms_same_representative) ||
    anyDuplicated(concordance$MAG) || anyNA(concordance$shared_carrier) ||
    !all(concordance$shared_carrier == (concordance$candidate_1_carried & concordance$candidate_2_carried))) {
  stop("Invalid carrier concordance.", call. = FALSE)
}
expected_counts <- c(sum(concordance$candidate_1_carried & !concordance$candidate_2_carried),
                     sum(concordance$shared_carrier),
                     sum(concordance$candidate_2_carried & !concordance$candidate_1_carried))
classes <- c("KO:K04028 only", "Shared by both candidates", "EC:1.2.1.87 only")
if (nrow(carrier_overlap) != 3L || anyDuplicated(carrier_overlap$display_class) ||
    !setequal(carrier_overlap$display_class, classes) ||
    !all(carrier_overlap$carrier_mags[match(classes, carrier_overlap$display_class)] == expected_counts) ||
    sum(expected_counts) != nrow(concordance) || !nrow(concordance)) {
  stop("Carrier-set totals disagree with supplied carrier membership.", call. = FALSE)
}
both <- proteins[carries_both_candidate_terms_same_representative %in% TRUE]
if (nrow(curation) != length(candidate_ids) || !setequal(curation$candidate_id, candidate_ids) ||
    any(curation$representative_clusters_carrying_both_terms != uniqueN(both$representative_id)) ||
    any(curation$carrier_MAGs_with_both_terms_on_same_representative != uniqueN(both$representative_mag)) ||
    any(curation$carrier_union_MAGs != nrow(concordance)) ||
    any(curation$shared_carrier_MAGs != sum(concordance$shared_carrier))) {
  stop("Curation summary disagrees with supplied protein/carrier evidence.", call. = FALSE)
}
if (!nrow(transitions[layer == cfg$analysis$primary_functional_layer])) {
  stop("No primary-layer adjacent transitions supplied.", call. = FALSE)
}

# GDV colour palettes v0.9.1. The Euler fills use the documented
# colour-blind-safe bivariate family #487DB0 / #9D96BC / #EBADC8.
blue <- "#487DB0"
pink <- "#D21C5E"
soft_pink <- "#EBADC8"
shared_lilac <- "#9D96BC"
orange <- "#FD8D3C"
teal <- "#089099"
dark_purple <- "#443F90"
ink <- "#202124"
grey <- "#A7ADB4"
light_grey <- "#F5DDEB"

# Panel a: descriptive catalogue pattern. Cluster counts are not replicates.
unknown[, size_group := factor(size_group, levels = c("1", "2-5", "6-20", "21-100"), ordered = TRUE)]
unknown[, count_label := ifelse(clusters >= 1e6, sprintf("n=%.1fM", clusters / 1e6), ifelse(clusters >= 1e3, sprintf("n=%.0fk", clusters / 1e3), paste0("n=", clusters)))]
unknown[, axis_label := paste0(as.character(size_group), "\n", sub("^n=", "", count_label))]
panel_a <- ggplot2::ggplot(unknown, ggplot2::aes(size_group, unknown_or_unannotated_fraction)) +
  ggplot2::geom_col(width = 0.66, fill = pink, colour = "#7C1D6F", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = scales::percent(unknown_or_unannotated_fraction, accuracy = 0.1)), vjust = -0.45, size = 2.35, colour = ink) +
  ggplot2::scale_x_discrete(labels = setNames(unknown$axis_label, unknown$size_group)) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, max(0.01, max(unknown$unknown_or_unannotated_fraction, na.rm = TRUE) * 1.15)), expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::labs(tag = "a", title = "Unknown clusters", x = "Member proteins per cluster\n(cluster count below)", y = "Unknown or unannotated clusters") +
  theme_manuscript() +
  ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(size = 5.9, lineheight = 0.88))

# Panel b: descriptive candidate means and standard errors from unnormalized values.
candidate_summary <- trajectories[, .(
  mean_potential = mean(potential),
  standard_error = stats::sd(potential) / sqrt(.N),
  samples = .N
), by = .(candidate_id, Timepoint)]
candidate_summary[, Timepoint := factor(Timepoint, levels = timepoint_levels, ordered = TRUE)]
candidate_summary[, candidate_label := factor(
  candidate_id,
  levels = c("KO:K04028", "EC:1.2.1.87"),
  labels = c("K04028 (EutN-like shell)", "EC 1.2.1.87 (propanal dehydrogenase)")
)]
candidate_screen <- longitudinal_models[
  (ontology == "KO" & term == "K04028") |
    (ontology == "EC" & term == "1.2.1.87"),
  .(candidate_id = paste(ontology, term, sep = ":"), ontology_BH_q = omnibus_p_adj_bh)
]
if (!setequal(candidate_screen$candidate_id, candidate_ids)) {
  stop("Selected candidate screen rows do not match the configured candidates.", call. = FALSE)
}
candidate_q <- stats::setNames(candidate_screen$ontology_BH_q, candidate_screen$candidate_id)
candidate_summary[, legend_label := factor(
  ifelse(candidate_id == "KO:K04028", "K04028", "EC 1.2.1.87"),
  levels = c("K04028", "EC 1.2.1.87")
)]
panel_b <- ggplot2::ggplot(candidate_summary, ggplot2::aes(Timepoint, mean_potential, colour = legend_label, group = legend_label, shape = legend_label)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(0, mean_potential - standard_error), ymax = mean_potential + standard_error), width = 0.12, linewidth = 0.38, alpha = 0.9) +
  ggplot2::geom_line(linewidth = 0.72) +
  ggplot2::geom_point(size = 1.8, fill = "white", stroke = 0.55) +
  ggplot2::scale_colour_manual(values = c("K04028" = blue, "EC 1.2.1.87" = pink)) +
  ggplot2::scale_shape_manual(values = c("K04028" = 21, "EC 1.2.1.87" = 24)) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.13))) +
  ggplot2::annotate("text", x = 1, y = Inf,
    label = sprintf("Ontology-wide BH q: %.3f (K04028); %.3f (EC)", candidate_q[["KO:K04028"]], candidate_q[["EC:1.2.1.87"]]),
    hjust = 0, vjust = 1.3, size = 2.05, colour = "#55585C") +
  ggplot2::labs(tag = "b", title = "Selected functional potentials", x = "Expedition timepoint", y = "Mean unnormalized potential", colour = NULL, shape = NULL) +
  theme_manuscript() +
  ggplot2::theme(legend.position = "top", legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.spacing = grid::unit(0, "pt"))

# Panel c: guide-inspired schematic Venn diagram of carrier-MAG membership.
# Equal circles prioritize readable set relationships; their areas do not
# encode the unequal set totals. Exact counts remain direct labels and data.
carrier_overlap[, fill_class := data.table::fcase(
  display_class == "KO:K04028 only", "K04028 only",
  display_class == "Shared by both candidates", "Shared",
  display_class == "EC:1.2.1.87 only", "EC 1.2.1.87 only"
)]
carrier_overlap[, union_fraction := carrier_mags / sum(carrier_mags)]

ko_only <- carrier_overlap[fill_class == "K04028 only", carrier_mags]
shared_carriers <- carrier_overlap[fill_class == "Shared", carrier_mags]
ec_only <- carrier_overlap[fill_class == "EC 1.2.1.87 only", carrier_mags]
ko_total <- ko_only + shared_carriers
ec_total <- ec_only + shared_carriers
union_total <- ko_only + shared_carriers + ec_only
shared_union_fraction <- shared_carriers / union_total

carrier_overlap[, `:=`(
  segment_label = data.table::fcase(
    fill_class == "K04028 only", as.character(carrier_mags),
    fill_class == "Shared", sprintf("%d MAGs\n%.1f%% of union", carrier_mags, 100 * shared_union_fraction),
    fill_class == "EC 1.2.1.87 only", as.character(carrier_mags)
  ),
  diagram_geometry = "schematic_equal_circles",
  areas_proportional = FALSE,
  union_total = union_total,
  shared_union_fraction = shared_union_fraction
)]

circle_polygon <- function(center_x, radius, points = 720L) {
  angle <- seq(0, 2 * pi, length.out = points)
  data.table::data.table(x = center_x + radius * cos(angle), y = radius * sin(angle))
}

radius_ko <- 1
radius_ec <- 1
center_distance <- 0.72
center_ko <- -center_distance / 2
center_ec <- center_distance / 2

chord_from_ko <- (radius_ko^2 - radius_ec^2 + center_distance^2) / (2 * center_distance)
chord_from_ec <- center_distance - chord_from_ko
theta_ko <- acos(chord_from_ko / radius_ko)
theta_ec <- acos(chord_from_ec / radius_ec)
lens <- data.table::rbindlist(list(
  data.table::data.table(
    x = center_ko + radius_ko * cos(seq(-theta_ko, theta_ko, length.out = 360L)),
    y = radius_ko * sin(seq(-theta_ko, theta_ko, length.out = 360L))
  ),
  data.table::data.table(
    x = center_ec + radius_ec * cos(seq(pi - theta_ec, pi + theta_ec, length.out = 360L)),
    y = radius_ec * sin(seq(pi - theta_ec, pi + theta_ec, length.out = 360L))
  )
))
circle_ko <- circle_polygon(center_ko, radius_ko)
circle_ec <- circle_polygon(center_ec, radius_ec)

panel_c <- ggplot2::ggplot() +
  ggplot2::geom_polygon(data = circle_ko, ggplot2::aes(x, y), fill = scales::alpha(blue, 0.70), colour = NA) +
  ggplot2::geom_polygon(data = circle_ec, ggplot2::aes(x, y), fill = scales::alpha(soft_pink, 0.78), colour = NA) +
  ggplot2::geom_polygon(data = lens, ggplot2::aes(x, y), fill = scales::alpha(shared_lilac, 0.94), colour = NA) +
  ggplot2::geom_path(data = circle_ko, ggplot2::aes(x, y), colour = "#365F95", linewidth = 0.45) +
  ggplot2::geom_path(data = circle_ec, ggplot2::aes(x, y), colour = "#A45C80", linewidth = 0.45) +
  ggplot2::annotate("text", x = -0.86, y = 1.22,
    label = sprintf("K04028\n%d carrier MAGs", ko_total), colour = "#365F95",
    fontface = "bold", size = 2.05, lineheight = 0.9) +
  ggplot2::annotate("text", x = 0.86, y = 1.22,
    label = sprintf("EC 1.2.1.87\n%d carrier MAGs", ec_total), colour = "#8B4167",
    fontface = "bold", size = 2.05, lineheight = 0.9) +
  ggplot2::annotate("text", x = -0.91, y = -0.03,
    label = as.character(ko_only), colour = "white", fontface = "bold", size = 2.75) +
  ggplot2::annotate("text", x = 0, y = 0,
    label = sprintf("%d MAGs\n%.1f%% of union", shared_carriers, 100 * shared_union_fraction),
    colour = "white", fontface = "bold", size = 2.25, lineheight = 0.9) +
  ggplot2::annotate("text", x = 0.91, y = -0.03,
    label = as.character(ec_only), colour = ink, fontface = "bold", size = 2.75) +
  ggplot2::coord_equal(
    xlim = c(-1.58, 1.58),
    ylim = c(-1.25, 1.55),
    expand = FALSE,
    clip = "off"
  ) +
  ggplot2::labs(
    tag = "c", title = "Carrier-MAG overlap",
    subtitle = sprintf("Schematic Venn diagram\nNon-proportional areas; union n = %d", union_total),
    x = NULL, y = NULL
  ) +
  theme_manuscript() +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
    legend.position = "none"
  )

# Panel d: descriptive contributions for the configured shared contrast.
strongest_id <- cfg$analysis$strongest_shared_contrast
selected_contrast <- shared_contrasts[contrast_id == strongest_id]
if (nrow(selected_contrast) != 1L) stop("Configured contrast must have exactly one supplied row.", call. = FALSE)
before_label <- selected_contrast$before_timepoint
after_label <- selected_contrast$after_timepoint
focus_drivers <- copy(genus_drivers[contrast_id == strongest_id])
focus_drivers[, absolute_mean_delta := abs(mean_delta)]
setorder(focus_drivers, -absolute_mean_delta)
top_count <- min(8L, nrow(focus_drivers))
driver <- focus_drivers[seq_len(top_count), .(genus, mean_delta)]
if (nrow(focus_drivers) > top_count) driver <- rbindlist(list(driver, data.table(genus = "Other shared-carrier genera", mean_delta = sum(focus_drivers[(top_count + 1L):.N, mean_delta]))))
driver[, direction := fifelse(mean_delta >= 0, "Higher after", "Lower after")]
if (!nrow(driver)) stop("No genus-driver plot rows found.", call. = FALSE)
driver[, genus_label := sub("^g__", "", genus)]
driver[genus_label == "Other shared-carrier genera", genus_label := "Other genera"]
driver[, genus_label := factor(genus_label, levels = rev(genus_label[order(abs(mean_delta), decreasing = TRUE)]))]
panel_d <- ggplot2::ggplot(driver, ggplot2::aes(mean_delta, genus_label, fill = mean_delta >= 0)) +
  ggplot2::geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, colour = ink) +
  ggplot2::scale_fill_manual(values = c(`TRUE` = orange, `FALSE` = teal), guide = "none") +
  ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.01), expand = ggplot2::expansion(mult = c(0.08, 0.11))) +
  ggplot2::labs(tag = "d", title = paste0(before_label, "–", after_label, " genus drivers"), subtitle = sprintf("Shared %d-MAG component; descriptive", shared_carriers), x = paste0("Mean paired change (", after_label, " − ", before_label, ")"), y = NULL) +
  theme_manuscript() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.y = ggplot2::element_text(size = 6.2),
    plot.margin = ggplot2::margin(5, 6, 5, 14)
  )

# Panel e: primary KO matched adjacent-transition relationship.
ko_transitions <- transitions[layer == cfg$analysis$primary_functional_layer]
ko_sensitivity <- participant_sensitivity[layer == "KO" & sensitivity_family == "participant_spearman_transition_coupling"]
panel_e <- ggplot2::ggplot(ko_transitions, ggplot2::aes(taxonomic_standardized_distance, functional_standardized_distance)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.45, colour = "#6B7075") +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = dark_purple, fill = light_grey, linewidth = 0.65) +
  ggplot2::geom_point(shape = 21, size = 1.55, stroke = 0.35, colour = "#365F95", fill = blue, alpha = 0.66) +
  ggplot2::annotate("label", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.25,
    label = sprintf("Median within-participant rho = %.2f\n%d/%d positive; BH q = %.3g",
      ko_sensitivity$participant_summary, ko_sensitivity$positive_direction_participants,
      ko_sensitivity$participants, ko_sensitivity$p_adj_bh_within_sensitivity_family),
    size = 2.1, fill = "white", colour = ink) +
  ggplot2::labs(tag = "e", title = "Taxonomic–KO movement", subtitle = sprintf("%d adjacent transitions from %d participants", nrow(ko_transitions), uniqueN(ko_transitions$subject_id)), x = "Taxonomic MAG distance (standardized)", y = "KO potential distance (standardized)") +
  theme_manuscript() + ggplot2::theme(legend.position = "none")

layout <- rbind(
  c(1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3),
  c(4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
)
main_figure <- gridExtra::arrangeGrob(
  panel_a, panel_b, panel_c, panel_d, panel_e,
  layout_matrix = layout,
  heights = c(1, 1.05)
)
figure_files <- save_composite(main_figure, "Main_Figure_01_functional_catalogue_and_dynamics", cfg, paths)

# Main Table 1: exact study and catalogue scope.
metric_value <- function(name) catalogue[metric == name, value][[1]]
metric_percent <- function(name) catalogue[metric == name, percent_of_relevant_denominator][[1]]
table1 <- data.table::data.table(
  Domain = c(
    "Longitudinal cohort", "Longitudinal cohort", "Gut abundance layer",
    "Protein catalogue", "Protein catalogue", "Protein catalogue",
    "Protein catalogue", "Annotation", "Annotation", "Annotation"
  ),
  Measure = c(
    "Included samples", "Participants", "Retained gut MAGs",
    "Catalogue MAGs", "Member protein sequences", "Representative clusters",
    "Cross-MAG representative clusters", "Representatives with eggNOG row",
    "Representatives with functional annotation", "Representatives without eggNOG row"
  ),
  Value = c(
    as.character(uniqueN(trajectories$sampleID)), as.character(uniqueN(trajectories$Subject)), format(metric_value("gut_mags_represented_by_members"), big.mark = ",", scientific = FALSE),
    format(metric_value("unique_parsed_mags"), big.mark = ",", scientific = FALSE),
    format(metric_value("unique_member_proteins"), big.mark = ",", scientific = FALSE),
    format(metric_value("unique_representatives"), big.mark = ",", scientific = FALSE),
    format(metric_value("cross_mag_representative_clusters"), big.mark = ",", scientific = FALSE),
    sprintf("%s (%.1f%%)", format(metric_value("representatives_with_annotation_row"), big.mark = ",", scientific = FALSE), metric_percent("representatives_with_annotation_row")),
    sprintf("%s (%.1f%%)", format(metric_value("representatives_with_functional_annotation"), big.mark = ",", scientific = FALSE), metric_percent("representatives_with_functional_annotation")),
    sprintf("%s (%.1f%%)", format(metric_value("representatives_without_annotation_row"), big.mark = ",", scientific = FALSE), metric_percent("representatives_without_annotation_row"))
  ),
  Interpretation = c(
    "T1-T10 only; Ta and Tb excluded", "Repeated longitudinal sampling", "Exact MAG identifiers joined to taxonomy",
    "Catalogue MAG identifiers used for taxonomy matching", "Member sequences are not abundance values", "MMseqs representatives",
    "Representative clusters containing proteins from more than one MAG", "Annotation-row coverage", "At least one retained functional field", "No representative annotation record"
  )
)

tax_global <- global[layer == "taxonomic_MAG" & metric == "hellinger"]
ko_global <- global[layer == "KO" & metric == "hellinger"]
ko_model <- stability_models[layer == "KO" & model_family == "taxonomic_functional_transition_coupling"]
ko_relative <- participant_sensitivity[layer == "KO" & sensitivity_family == "participant_median_relative_movement"]
shared_model <- component_models[component == "shared_carriers" & model_family == "absolute_component_temporal_pattern"]
t2_t6 <- shared_contrasts[contrast_id == strongest_id]
leading_genus <- genus_drivers[contrast_id == strongest_id][order(-abs(mean_delta), genus)][1L]
if (!nrow(leading_genus) || is.na(leading_genus$genus)) stop("No genus contribution for configured contrast.", call. = FALSE)
ko_recovery <- recovery[layer == "KO"]

table2 <- data.table::data.table(
  Analysis = c(
    "Singleton unknown enrichment", "K04028 timepoint screen", "EC 1.2.1.87 timepoint screen",
    sprintf("Shared %d-MAG component", shared_carriers), paste0("Shared component: ", before_label, " vs ", after_label), paste0(before_label, "-", after_label, " ", sub("^g__", "", leading_genus$genus), " contribution"),
    "Taxonomic profile across time", "KO profile across time", "Taxonomic-KO transition coupling",
    "Relative KO movement", "Exploratory T1-T10 return"
  ),
  Effect = c(
    sprintf("OR %.2f (95%% CI %.2f-%.2f)", singleton$odds_ratio, singleton$lower_95, singleton$upper_95),
    "Ontology-wide omnibus", "Ontology-wide omnibus",
    sprintf("Timepoint omnibus across %d MAGs", shared_model$carrier_mags),
    sprintf("%s/%s ratio %.2f", after_label, before_label, t2_t6$ratio_after_vs_before),
    sprintf("%.1f%% of absolute genus-level change", 100 * leading_genus$absolute_change_share),
    sprintf("PERMANOVA R2 = %.3f", tax_global$permanova_r_squared),
    sprintf("PERMANOVA R2 = %.3f", ko_global$permanova_r_squared),
    sprintf("Median participant Spearman rho = %.2f; %d/%d positive", ko_sensitivity$participant_summary, ko_sensitivity$positive_direction_participants, ko_sensitivity$participants),
    sprintf("Median functional/taxonomic ratio = %.2f; %d/%d above 1", ko_relative$ratio, ko_relative$positive_direction_participants, ko_relative$participants),
    sprintf("%d/%d participants showed partial KO return", sum(ko_recovery$recovery_fraction_from_peak > 0), data.table::uniqueN(ko_recovery$subject_id))
  ),
  Statistical_evidence = c(
    sprintf("Descriptive catalogue model; p = %s; status = %s", format(singleton$p_value, scientific = TRUE, digits = 4), singleton$status),
    sprintf("Ontology-specific BH q = %.4f", candidate_q[["KO:K04028"]]),
    sprintf("Ontology-specific BH q = %.4f", candidate_q[["EC:1.2.1.87"]]),
    sprintf("BH q = %.3g", shared_model$p_adj_bh_within_family), sprintf("BH q = %.3g", t2_t6$p_adj_bh_across_45_timepoint_contrasts),
    "Descriptive decomposition; no genus-specific test",
    sprintf("restricted-permutation p = %.4f", tax_global$permanova_p_value), sprintf("BH q = %.3f", ko_global$permanova_p_adj_bh_functional_family),
    sprintf("signed-rank BH q = %.4f", ko_sensitivity$p_adj_bh_within_sensitivity_family), sprintf("signed-rank BH q = %.4f", ko_relative$p_adj_bh_within_sensitivity_family),
    sprintf("Descriptive; paired T1 and T10 n = %d", uniqueN(ko_recovery$subject_id))
  ),
  Interpretation = c(
    sprintf("Singleton/larger unknown-status odds ratio %.2f; clusters are not biological replicates.", singleton$odds_ratio),
    sprintf("Selected predicted potential: ontology-specific q %.4f.", candidate_q[["KO:K04028"]]),
    sprintf("Selected predicted potential: ontology-specific q %.4f.", candidate_q[["EC:1.2.1.87"]]),
    sprintf("Shared-carrier temporal omnibus: reported q %.3g; inspect upstream fit diagnostics.", shared_model$p_adj_bh_within_family),
    "Configured shared-component contrast; selected from the same data, not independent replication.",
    paste0(sub("^g__", "", leading_genus$genus), " has the largest absolute mean contribution; no separate taxon test."),
    sprintf("Taxonomic timepoint R2 %.3f; restricted-permutation p %.4f.", tax_global$permanova_r_squared, tax_global$permanova_p_value),
    sprintf("Functional timepoint R2 %.3f; reported q %.3f.", ko_global$permanova_r_squared, ko_global$permanova_p_adj_bh_functional_family),
    sprintf("Median within-participant association %.2f; %d/%d positive.", ko_sensitivity$participant_summary, ko_sensitivity$positive_direction_participants, ko_sensitivity$participants),
    sprintf("Median relative movement %.2f; interpret with supplied signed-rank evidence.", ko_relative$ratio),
    sprintf("Descriptive return in %d paired participants; no cohort-level recovery inference.", uniqueN(ko_recovery$subject_id))
  )
)
data.table::setnames(table2, "Statistical_evidence", "Statistical evidence")

for (stem in c("Main_Table_01_study_and_catalogue_summary", "Main_Table_02_key_statistical_results")) {
  x <- if (grepl("01_", stem)) table1 else table2
  write_result_table(x, file.path(paths$main, paste0(stem, ".csv")))
  data.table::fwrite(x, file.path(paths$main, paste0(stem, ".tsv")), sep = "\t", na = "")
}
table_preview_files <- c(
  save_table_preview(table1, "Main_Table_01_study_and_catalogue_summary", cfg, paths, c(15, 21, 14, 35), height = 5.5),
  save_table_preview(table2, "Main_Table_02_key_statistical_results", cfg, paths, c(18, 20, 20, 34), height = 7.2)
)

write_result_table(unknown, file.path(paths$reproducibility, "Main_Figure_01_panel_a.csv"))
write_result_table(candidate_summary, file.path(paths$reproducibility, "Main_Figure_01_panel_b.csv"))
write_result_table(carrier_overlap, file.path(paths$reproducibility, "Main_Figure_01_panel_c.csv"))
write_result_table(driver, file.path(paths$reproducibility, "Main_Figure_01_panel_d.csv"))
write_result_table(ko_transitions, file.path(paths$reproducibility, "Main_Figure_01_panel_e.csv"))

index <- c(
  "# Manuscript Functional Synthesis Results Index", "",
  "Generated from completed MAG functional annotation (Stages 0-10) and functional community dynamics (Stages 0-4).", "",
  "- Main Figure 1: `main/Main_Figure_01_functional_catalogue_and_dynamics.{pdf,png,svg}`",
  "- Main Table 1: `main/Main_Table_01_study_and_catalogue_summary.{csv,tsv,pdf,png}`",
  "- Main Table 2: `main/Main_Table_02_key_statistical_results.{csv,tsv,pdf,png}`",
  "- Exact panel data: `reproducibility/Main_Figure_01_panel_[a-e].csv`", "",
  "Panel c is a schematic equal-circle Venn diagram; circle areas are not proportional to set sizes.",
  "The figure reports predicted genomic potential, not expression, enzyme activity, pathway completeness, metabolite production, adaptation, or causality."
)
writeLines(index, file.path(paths$root, "RESULTS_INDEX.md"))

session_path <- file.path(paths$reproducibility, "00_build_main_figure_and_tables_session_info.txt")
capture.output(sessionInfo(), file = session_path)
writeLines(c(
  "Manuscript functional synthesis complete.",
  paste("Main figure files:", length(figure_files)),
  paste("Main table preview files:", length(table_preview_files)),
  "Panel C geometry: schematic equal circles; areas are not proportional.",
  "Source results were read without modification.",
  "Generated results remain uncommitted under repository policy."
), marker)

message("Manuscript synthesis completed: ", file.path(paths$root, "RESULTS_INDEX.md"))
