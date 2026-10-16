#!/usr/bin/env Rscript

# Purpose: Decompose each selected functional-potential candidate into the
# contribution from MAGs shared by both candidates and MAGs specific to one
# candidate carrier set.
# Inputs: Stage-06 carrier and sample-MAG contribution outputs.
# Outputs: Numbered Tables 15-17, Figures 11-12, and an updated follow-up index.
# Dependencies: R, data.table, lme4, ggplot2, scales, yaml, common.R, followup.R.
# Usage: Rscript analysis/mag_functional_annotation/08_decompose_shared_unique_carriers.R --config config/mag_functional_annotation.yml
# Workflow: Exploratory same-data decomposition; carrier overlap is MAG-set
# membership and is not evidence that two annotations occur on the same gene.

source("analysis/mag_functional_annotation/R/common.R")
source("analysis/mag_functional_annotation/R/followup.R")
args <- parse_stage_args()
cfg <- read_magfunc_config(args$config)
assert_packages("lme4")
paths <- output_paths(cfg)
if (!file.exists(file.path(paths$qc, "07_explore_pathway_unknowns.complete"))) {
  stop("Stage 07 completion marker is missing.", call. = FALSE)
}
marker <- file.path(paths$qc, "08_decompose_shared_unique_carriers.complete")
guard_stage(marker, args, cfg)

candidates <- configured_followup_candidates(cfg)
if (nrow(candidates) != 2L) {
  stop("Shared-versus-unique decomposition currently requires exactly two configured candidates.", call. = FALSE)
}
candidate_ids <- candidates$candidate_id
carrier_atlas <- read_table(file.path(
  paths$targeted_validation,
  "Table_08_candidate_carrier_MAG_taxonomy.csv"
))
required_carrier_columns <- c("candidate_id", "MAG", "phylum", "genus", "species")
missing_carrier_columns <- setdiff(required_carrier_columns, names(carrier_atlas))
if (length(missing_carrier_columns) > 0L) {
  stop("Table 08 is missing columns: ", paste(missing_carrier_columns, collapse = ", "), call. = FALSE)
}
carrier_atlas <- carrier_atlas[candidate_id %in% candidate_ids]
if (!setequal(unique(carrier_atlas$candidate_id), candidate_ids)) {
  stop("Table 08 does not contain both configured candidates.", call. = FALSE)
}
if (anyDuplicated(carrier_atlas, by = c("candidate_id", "MAG"))) {
  stop("Table 08 must contain one row per candidate-MAG pair.", call. = FALSE)
}

carrier_sets <- setNames(lapply(candidate_ids, function(id) {
  sort(unique(carrier_atlas[candidate_id == id, MAG]))
}), candidate_ids)
shared_mags <- Reduce(intersect, carrier_sets)
candidate_specific_sets <- setNames(lapply(candidate_ids, function(id) {
  setdiff(carrier_sets[[id]], shared_mags)
}), candidate_ids)
if (length(shared_mags) == 0L || any(lengths(candidate_specific_sets) == 0L)) {
  stop("Decomposition requires shared MAGs and at least one candidate-specific MAG for each candidate.", call. = FALSE)
}

taxonomy_columns <- intersect(
  c("MAG", "domain", "phylum", "class", "order", "ordeer", "family", "genus", "species"),
  names(carrier_atlas)
)
if ("order" %in% taxonomy_columns) {
  taxonomy_columns <- setdiff(taxonomy_columns, "ordeer")
}
taxonomy <- unique(carrier_atlas[, ..taxonomy_columns])
if (anyDuplicated(taxonomy$MAG)) {
  stop("Candidate carrier taxonomy is inconsistent across candidate rows.", call. = FALSE)
}
membership <- carrier_atlas[, .(
  candidate_ids = paste(sort(unique(candidate_id)), collapse = ";"),
  candidate_terms = paste(sort(unique(term)), collapse = ";"),
  candidates_carried = data.table::uniqueN(candidate_id)
), by = MAG]
membership[, `:=`(
  membership_class = data.table::fifelse(
    candidates_carried == 2L,
    "shared_by_both_candidates",
    paste0(candidate_ids, "_only")
  ),
  interpretation = data.table::fifelse(
    candidates_carried == 2L,
    "MAG carries both selected annotation terms; this does not imply the same gene or a complete pathway",
    "MAG carries only one of the two selected annotation terms"
  )
)]
membership <- merge(membership, taxonomy, by = "MAG", all.x = TRUE)
data.table::setorder(membership, -candidates_carried, candidate_ids, MAG)

contribution_path <- file.path(paths$targeted_validation, "candidate_sample_MAG_contributions.csv.gz")
contributions <- read_table(contribution_path)
required_contribution_columns <- c(
  "candidate_id", "ontology", "term", "MAG", "sampleID", "Subject", "Timepoint",
  "analysis_date", "altitude", "days_from_baseline_subject", "mag_value", "candidate_total"
)
missing_contribution_columns <- setdiff(required_contribution_columns, names(contributions))
if (length(missing_contribution_columns) > 0L) {
  stop("Candidate contribution table is missing columns: ", paste(missing_contribution_columns, collapse = ", "), call. = FALSE)
}
contributions <- contributions[candidate_id %in% candidate_ids]
if (anyDuplicated(contributions, by = c("candidate_id", "sampleID", "MAG"))) {
  stop("Candidate contribution rows must be unique by candidate, sample, and MAG.", call. = FALSE)
}
if (any(contributions[, data.table::uniqueN(candidate_total), by = .(candidate_id, sampleID)]$V1 != 1L)) {
  stop("Each candidate-sample must have exactly one original candidate total.", call. = FALSE)
}
contributions[, component := data.table::fifelse(
  MAG %in% shared_mags,
  "shared_carriers",
  "candidate_specific_carriers"
)]
if (any(!contributions[component == "candidate_specific_carriers", MAG] %in%
        unlist(candidate_specific_sets, use.names = FALSE))) {
  stop("A candidate-specific contribution is not present in the expected carrier set.", call. = FALSE)
}

component_samples <- contributions[, .(
  component_potential = sum(mag_value),
  carrier_mags = data.table::uniqueN(MAG),
  original_candidate_total = unique(candidate_total)
), by = .(
  candidate_id, ontology, term, sampleID, Subject, Timepoint,
  analysis_date, altitude, days_from_baseline_subject, component
)]
component_counts <- component_samples[, .N, by = .(candidate_id, sampleID)]
if (any(component_counts$N != 2L)) {
  stop("Every candidate-sample must retain shared and candidate-specific components.", call. = FALSE)
}
component_samples[, reconstructed_candidate_total := sum(component_potential), by = .(candidate_id, sampleID)]
component_samples[, reconstruction_error := reconstructed_candidate_total - original_candidate_total]
maximum_reconstruction_error <- max(abs(component_samples$reconstruction_error), na.rm = TRUE)
if (!is.finite(maximum_reconstruction_error) || maximum_reconstruction_error > 1e-10) {
  stop("Shared plus candidate-specific components do not reconstruct the Stage-06 candidate total.", call. = FALSE)
}
component_samples[, component_share := data.table::fifelse(
  original_candidate_total > 0,
  component_potential / original_candidate_total,
  NA_real_
)]
share_checks <- component_samples[, .(share_sum = sum(component_share)), by = .(candidate_id, sampleID)]
if (any(!is.finite(share_checks$share_sum)) || max(abs(share_checks$share_sum - 1)) > 1e-10) {
  stop("Component shares do not sum to one for every candidate-sample.", call. = FALSE)
}

shared_wide <- data.table::dcast(
  component_samples[component == "shared_carriers"],
  sampleID ~ candidate_id,
  value.var = "component_potential"
)
shared_profile_error <- max(abs(shared_wide[[candidate_ids[[1L]]]] - shared_wide[[candidate_ids[[2L]]]]))
if (!is.finite(shared_profile_error) || shared_profile_error > 1e-12) {
  stop("The shared-carrier component must be identical for both candidate reconstructions.", call. = FALSE)
}

fit_component <- function(d, candidate_id, component, model_family, response, pseudocount, carrier_mags) {
  fit <- fit_timepoint_response_model(d, cfg, response)
  data.table::data.table(
    candidate_id = candidate_id,
    component = component,
    model_family = model_family,
    carrier_mags = as.integer(carrier_mags),
    n_samples = nrow(d),
    n_subjects = data.table::uniqueN(d$Subject),
    n_nonzero_samples = if ("component_potential" %in% names(d)) sum(d$component_potential > 0) else NA_integer_,
    pseudocount = pseudocount,
    omnibus_chisq = fit$chisq %||% NA_real_,
    omnibus_df = fit$df %||% NA_real_,
    omnibus_p_value = fit$p_value %||% NA_real_,
    singular = fit$singular %||% NA,
    converged = fit$converged %||% NA,
    convergence_messages = fit$convergence_messages %||% "",
    status = fit$status,
    inference_scope = "exploratory same-data carrier decomposition; not independent replication"
  )
}

model_rows <- list()
shared_data <- data.table::copy(component_samples[
  candidate_id == candidate_ids[[1L]] & component == "shared_carriers"
])
shared_positive <- shared_data$component_potential[shared_data$component_potential > 0]
shared_pseudocount <- if (length(shared_positive) > 0L) min(shared_positive) / 2 else NA_real_
shared_data[, response_log10 := log10(component_potential + shared_pseudocount)]
model_rows[["shared"]] <- fit_component(
  shared_data, "Both candidates", "shared_carriers", "absolute_component_temporal_pattern",
  "response_log10", shared_pseudocount, length(shared_mags)
)

for (target_candidate_id in candidate_ids) {
  specific <- data.table::copy(component_samples[
    candidate_id == target_candidate_id & component == "candidate_specific_carriers"
  ])
  positive <- specific$component_potential[specific$component_potential > 0]
  pseudocount <- if (length(positive) > 0L) min(positive) / 2 else NA_real_
  specific[, response_log10 := log10(component_potential + pseudocount)]
  model_rows[[paste0(target_candidate_id, "_specific")]] <- fit_component(
    specific, target_candidate_id, "candidate_specific_carriers", "absolute_component_temporal_pattern",
    "response_log10", pseudocount, unique(specific$carrier_mags)
  )

  balance <- data.table::dcast(
    component_samples[candidate_id == target_candidate_id],
    sampleID + Subject + Timepoint + analysis_date + altitude + days_from_baseline_subject ~ component,
    value.var = "component_potential"
  )
  balance_positive <- c(balance$shared_carriers, balance$candidate_specific_carriers)
  balance_positive <- balance_positive[balance_positive > 0]
  balance_pseudocount <- if (length(balance_positive) > 0L) min(balance_positive) / 2 else NA_real_
  balance[, response_balance_log10 := log10(
    (shared_carriers + balance_pseudocount) /
      (candidate_specific_carriers + balance_pseudocount)
  )]
  balance_row <- fit_component(
    balance, target_candidate_id, "shared_vs_candidate_specific_balance", "component_balance_temporal_pattern",
    "response_balance_log10", balance_pseudocount, NA_integer_
  )
  balance_row[, n_nonzero_samples := NA_integer_]
  model_rows[[paste0(target_candidate_id, "_balance")]] <- balance_row
}
component_models <- data.table::rbindlist(model_rows, fill = TRUE)
component_models[, p_adj_bh_within_family := safe_bh(omnibus_p_value), by = model_family]
component_models[, significant_bh := is.finite(p_adj_bh_within_family) &
  p_adj_bh_within_family <= as.numeric(cfg$analysis$fdr_threshold)]
component_models[, adjustment_family := as.character(cfg$followup$shared_unique_adjustment)]
data.table::setorder(component_models, model_family, p_adj_bh_within_family, candidate_id)

component_summary <- component_samples[, .(
  n_samples = .N,
  n_subjects = data.table::uniqueN(Subject),
  mean_component_potential = mean(component_potential),
  standard_error = stats::sd(component_potential) / sqrt(.N),
  median_component_potential = stats::median(component_potential),
  mean_component_share = mean(component_share),
  median_component_share = stats::median(component_share),
  nonzero_fraction = mean(component_potential > 0)
), by = .(candidate_id, ontology, term, Timepoint, component, carrier_mags)]

table15 <- file.path(paths$carrier_decomposition, "Table_15_candidate_carrier_set_membership.csv")
table16 <- file.path(paths$carrier_decomposition, "Table_16_sample_component_decomposition.csv.gz")
table17 <- file.path(paths$carrier_decomposition, "Table_17_component_temporal_models.csv")
write_table(membership, table15)
write_table(component_samples, table16)
write_table(component_models, table17)
write_table(component_summary, file.path(paths$carrier_decomposition, "component_timepoint_summary.csv"))

membership_counts <- membership[, .(carrier_mags = .N), by = membership_class]
membership_counts[, percent_of_union := carrier_mags / sum(carrier_mags)]
display_lookup <- c(
  shared_by_both_candidates = "Shared by both candidates",
  setNames(paste0(candidate_ids, " only"), paste0(candidate_ids, "_only"))
)
membership_counts[, display_class := unname(display_lookup[membership_class])]
membership_counts[, membership_type := data.table::fifelse(
  membership_class == "shared_by_both_candidates", "Shared", "Candidate-specific"
)]
membership_counts[, display_class := factor(
  display_class,
  levels = rev(c("Shared by both candidates", paste0(candidate_ids, " only")))
)]
overlap_colours <- c(
  Shared = magfunc_palette[["orange"]],
  `Candidate-specific` = magfunc_palette[["blue"]]
)
p11 <- ggplot2::ggplot(membership_counts, ggplot2::aes(carrier_mags, display_class, fill = membership_type)) +
  ggplot2::geom_col(width = 0.62, colour = "white", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%d (%.1f%% of union)", carrier_mags, 100 * percent_of_union)),
    hjust = -0.08, size = 3.2
  ) +
  ggplot2::scale_fill_manual(values = overlap_colours) +
  ggplot2::scale_x_continuous(
    limits = c(0, max(membership_counts$carrier_mags) * 1.42),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = "Carrier MAG-set overlap for the two selected candidates",
    subtitle = "Membership is defined at MAG level; shared does not mean the annotations occur on the same gene",
    x = "Number of carrier MAGs", y = NULL, fill = NULL
  ) + theme_nature() + ggplot2::theme(legend.position = "top")
files11 <- save_publication_plot(p11, "Figure_11_candidate_carrier_set_overlap", cfg, width = 7.4, height = 4.2)

component_labels <- c(
  shared_carriers = "Shared-carrier component",
  candidate_specific_carriers = "Candidate-specific component"
)
plot_summary <- data.table::copy(component_summary)
plot_summary[, component_label := unname(component_labels[component])]
plot_summary[, Timepoint := factor(Timepoint, levels = unlist(cfg$analysis$ordered_timepoints))]
component_colours <- c(
  `Shared-carrier component` = magfunc_palette[["orange"]],
  `Candidate-specific component` = magfunc_palette[["blue"]]
)
component_linetypes <- c(
  `Shared-carrier component` = "solid",
  `Candidate-specific component` = "22"
)
p12 <- ggplot2::ggplot(
  plot_summary,
  ggplot2::aes(Timepoint, mean_component_potential, colour = component_label, fill = component_label,
               linetype = component_label, group = component_label)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = pmax(0, mean_component_potential - standard_error),
      ymax = mean_component_potential + standard_error
    ),
    colour = NA, alpha = 0.12
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(shape = 21, size = 1.6, stroke = 0.6) +
  ggplot2::facet_wrap(~candidate_id, ncol = 1) +
  ggplot2::scale_colour_manual(values = component_colours) +
  ggplot2::scale_fill_manual(values = component_colours) +
  ggplot2::scale_linetype_manual(values = component_linetypes) +
  ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.06))) +
  ggplot2::labs(
    title = "Shared and candidate-specific contributions across expedition timepoints",
    subtitle = paste(
      "Arithmetic means ± SE of unnormalized MAG values",
      "The shared-carrier profile is the same source component in both panels",
      sep = "\n"
    ),
    x = "Ordered expedition timepoint", y = "Predicted functional-potential component",
    colour = NULL, fill = NULL, linetype = NULL
  ) + theme_nature() + ggplot2::theme(legend.position = "top")
files12 <- save_publication_plot(p12, "Figure_12_shared_unique_component_trajectories", cfg, width = 7.4, height = 6.8)

write_table(membership_counts, file.path(paths$reproducibility, "Figure_11_plot_data.csv"))
write_table(plot_summary, file.path(paths$reproducibility, "Figure_12_plot_data.csv"))
chart_map_path <- file.path(paths$reproducibility, "followup_publication_chart_map.csv")
chart_map <- if (file.exists(chart_map_path)) read_table(chart_map_path) else data.table::data.table()
new_chart_map <- data.table::data.table(
  figure = c("Figure 11", "Figure 12"),
  analytical_question = c(
    "How many carrier MAGs are shared by or specific to the two selected candidates?",
    "How do shared-carrier and candidate-specific contributions vary across expedition timepoints?"
  ),
  chart_family = c("comparison and composition", "longitudinal decomposition"),
  visual_encoding = c(
    "MAG count and percent of the carrier-MAG union by mutually exclusive membership class",
    "arithmetic mean component potential with between-sample standard error"
  ),
  denominator_or_unit = c(
    "union of MAGs carrying either selected annotation term",
    "sum of unnormalized gut.csv values for MAGs in each carrier component"
  ),
  palette_policy = c(
    "orange shared class; blue candidate-specific classes with direct labels",
    "orange shared component; dashed blue candidate-specific component"
  )
)
chart_map <- data.table::rbindlist(list(chart_map, new_chart_map), fill = TRUE)
chart_map <- unique(chart_map, by = "figure", fromLast = TRUE)
data.table::setorder(chart_map, figure)
write_table(chart_map, chart_map_path)

table_records <- data.table::data.table(
  result_number = sprintf("Table %02d", 15:17),
  category = "shared_unique_carrier_decomposition",
  question = c(
    "Which carrier MAGs are shared by or specific to the two selected candidates?",
    "How much does each carrier component contribute within every candidate-sample?",
    "Do shared, candidate-specific, or component-balance profiles vary across timepoints?"
  ),
  artifact_type = "table",
  path = vapply(c(table15, table16, table17), relative_to_root, character(1)),
  description = c(
    "Mutually exclusive MAG-set membership with taxonomy; overlap does not imply same-gene colocalization.",
    "Exact additive reconstruction of each candidate from shared and candidate-specific carrier MAG values.",
    "Exploratory participant-aware temporal models with BH correction within decomposition model family."
  ),
  stage = "08_decompose_shared_unique_carriers"
)
figure_records <- data.table::rbindlist(list(
  data.table::data.table(
    result_number = "Figure 11", category = "shared_unique_carrier_decomposition",
    question = "How large are the shared and candidate-specific carrier sets?", artifact_type = "figure",
    path = vapply(files11, relative_to_root, character(1)),
    description = "Mutually exclusive carrier-MAG membership counts and percentages of the union.",
    stage = "08_decompose_shared_unique_carriers"
  ),
  data.table::data.table(
    result_number = "Figure 12", category = "shared_unique_carrier_decomposition",
    question = "Which carrier component supplies each candidate potential over time?", artifact_type = "figure",
    path = vapply(files12, relative_to_root, character(1)),
    description = "Descriptive mean unnormalized potentials for shared and candidate-specific MAG components.",
    stage = "08_decompose_shared_unique_carriers"
  )
), fill = TRUE)
register_results(cfg, data.table::rbindlist(list(table_records, figure_records), fill = TRUE))

registry <- data.table::fread(file.path(paths$intermediate, "results_registry.csv"))
followup_registry <- registry[
  result_number %in% c(sprintf("Table %02d", 8:17), sprintf("Figure %02d", 7:12))
]
data.table::setorder(followup_registry, result_number, path)
followup_csv <- file.path(paths$root, "FOLLOWUP_RESULTS_INDEX.csv")
followup_md <- file.path(paths$root, "FOLLOWUP_RESULTS_INDEX.md")
write_table(followup_registry, followup_csv)
escape_md <- function(x) gsub("\\|", "\\\\|", as.character(x))
lines <- c(
  "# MAG Functional Annotation Follow-up Results Index", "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
  "## Analysis contract", "",
  "- The two candidates were selected from the same Stage-03 screen; follow-up p-values are not independent replication.",
  "- The unnormalized gut.csv analysis remains primary; row-sum normalization is sensitivity-only.",
  "- Shared carrier membership is defined at MAG level and does not imply that terms occur on the same gene.",
  "- Shared and candidate-specific components sum MAG values once per carrier and are not protein-copy weighted.",
  "- Co-carriage does not establish pathway completeness, activity, or metabolite production.",
  "- Unknown singleton status does not establish biological novelty.", "",
  "Number | Category | Type | Path | Purpose", "--- | --- | --- | --- | ---"
)
for (i in seq_len(nrow(followup_registry))) {
  lines <- c(lines, paste(
    escape_md(followup_registry$result_number[[i]]),
    escape_md(followup_registry$category[[i]]),
    escape_md(followup_registry$artifact_type[[i]]),
    paste0("`", escape_md(followup_registry$path[[i]]), "`"),
    escape_md(followup_registry$description[[i]]),
    sep = " | "
  ))
}
lines <- c(lines, "", "## Interpretation order", "",
  "1. Read Table 11 before interpreting the robustness of either selected signal.",
  "2. Read Table 09 and Figure 09 to identify participant influence.",
  "3. Use Tables 08 and 10 with Figures 07-08 for carrier and timepoint interpretation.",
  "4. Use Tables 15-17 and Figures 11-12 to decide whether shared or candidate-specific carrier MAGs drive the profiles.",
  "5. Treat Table 12 as exploratory co-carriage, not pathway completeness.",
  "6. Read Table 14 before making any novelty-versus-fragmentation statement from Table 13 or Figure 10."
)
writeLines(lines, followup_md)

writeLines(c(
  "Stage 08 complete",
  paste0("Shared carrier MAGs: ", length(shared_mags)),
  paste0(candidate_ids[[1L]], "-specific carrier MAGs: ", length(candidate_specific_sets[[candidate_ids[[1L]]]])),
  paste0(candidate_ids[[2L]], "-specific carrier MAGs: ", length(candidate_specific_sets[[candidate_ids[[2L]]]])),
  "Shared plus candidate-specific components exactly reconstruct each Stage-06 candidate potential.",
  "Component models are exploratory same-data decompositions, not independent replication."
), marker)
write_session_info(cfg, "08_decompose_shared_unique_carriers")
message("Stage 08 completed: ", followup_md)
