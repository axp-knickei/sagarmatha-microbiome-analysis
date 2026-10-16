#!/usr/bin/env Rscript

# Purpose: Validate the two prespecified BH-supported functional-potential
# candidates with carrier attribution, participant trajectories, leave-one-
# participant-out refits, pairwise timepoint contrasts, and captured-fraction
# sensitivity analyses.
# Inputs: Stage-02 MAG-function intermediates, Stage-03 potential/model outputs,
# gut.csv, sample_info.csv, and GDTB.csv.
# Outputs: Numbered Tables 08-11, Figures 07-09, and reproducibility tables.
# Dependencies: R, data.table, lme4, ggplot2, yaml, common.R, and followup.R.
# Usage: Rscript analysis/mag_functional_annotation/06_validate_selected_candidates.R --config config/mag_functional_annotation.yml
# Workflow: Selection-aware validation; it does not constitute independent replication.

source("analysis/mag_functional_annotation/R/common.R")
source("analysis/mag_functional_annotation/R/followup.R")
args <- parse_stage_args()
cfg <- read_magfunc_config(args$config)
assert_packages("lme4")
paths <- output_paths(cfg)
inputs <- resolve_input_paths(cfg)
if (!file.exists(file.path(paths$qc, "03_longitudinal_analyses.complete"))) {
  stop("Stage 03 completion marker is missing.", call. = FALSE)
}
marker <- file.path(paths$qc, "06_validate_selected_candidates.complete")
guard_stage(marker, args, cfg)

candidates <- configured_followup_candidates(cfg)
original_models <- read_table(file.path(paths$longitudinal, "Table_06_longitudinal_functional_models.csv.gz"))
candidate_models <- merge(candidates, original_models, by = c("ontology", "term"), all.x = TRUE)
if (anyNA(candidate_models$status)) {
  missing <- candidate_models[is.na(status), candidate_id]
  stop("Configured candidates are absent from Table 06: ", paste(missing, collapse = ", "), call. = FALSE)
}
if (any(candidate_models$status != "fitted")) {
  failed <- candidate_models[status != "fitted", paste0(candidate_id, " [", status, "]")]
  stop("Configured candidates do not have fitted Stage-03 models: ", paste(failed, collapse = ", "), call. = FALSE)
}

aligned <- align_gut_metadata(cfg)
taxonomy <- data.table::fread(inputs$taxonomy_file)
tax_id <- cfg$columns$taxonomy_mag_id
if ("ordeer" %in% names(taxonomy) && !"order" %in% names(taxonomy)) taxonomy[, order := ordeer]
if (anyDuplicated(taxonomy[[tax_id]])) stop("Taxonomy MAG identifiers must be unique for follow-up joins.", call. = FALSE)
metadata_join <- aligned$metadata[, .(
  sampleID = get(aligned$sample_col), Subject, Timepoint, analysis_date,
  altitude = altitude_numeric, days_from_baseline_subject
)]
row_sum_table <- data.table::data.table(
  sampleID = aligned$sample_ids,
  abundance_sum_as_provided = aligned$row_sums
)

trajectory_rows <- list()
carrier_rows <- list()
contribution_rows <- list()
lopo_rows <- list()
lopo_prediction_rows <- list()
posthoc_rows <- list()
sensitivity_rows <- list()
timepoint_estimate_rows <- list()

for (candidate_index in seq_len(nrow(candidates))) {
  candidate <- candidates[candidate_index]
  candidate_id <- candidate$candidate_id[[1L]]
  ontology <- candidate$ontology[[1L]]
  term <- candidate$term[[1L]]
  message("Validating selected candidate: ", candidate_id)

  potential_path <- file.path(paths$intermediate, "functional_potential", paste0(ontology, ".rds"))
  potential_all <- read_object(potential_path)
  d <- data.table::as.data.table(potential_all)[term == candidate$term[[1L]]]
  if (nrow(d) == 0L) stop("No potential rows found for ", candidate_id, ".", call. = FALSE)
  d <- merge(d, row_sum_table, by = "sampleID", all.x = TRUE)
  if (anyNA(d$abundance_sum_as_provided) || any(d$abundance_sum_as_provided <= 0)) {
    stop("Candidate rows require positive finite gut.csv row sums.", call. = FALSE)
  }
  d[, `:=`(
    candidate_id = candidate_id,
    selection_basis = candidate$selection_basis[[1L]]
  )]
  trajectory_rows[[candidate_id]] <- d[, .(
    candidate_id, ontology, term, sampleID, Subject, Timepoint, analysis_date,
    altitude, days_from_baseline_subject, potential, abundance_sum_as_provided,
    selection_basis
  )]

  original <- candidate_models[ontology == candidate$ontology[[1L]] & term == candidate$term[[1L]]][1L]
  pseudocount <- as.numeric(original$pseudocount)
  primary_fit <- fit_candidate_timepoint_model(d, cfg, pseudocount)
  if (primary_fit$status != "fitted") stop("Primary candidate refit failed for ", candidate_id, ": ", primary_fit$status, call. = FALSE)
  primary_estimates <- fixed_timepoint_estimates(primary_fit$model, unlist(cfg$analysis$ordered_timepoints))
  primary_estimates[, `:=`(candidate_id = candidate_id, scenario = "primary_unnormalized")]
  timepoint_estimate_rows[[paste0(candidate_id, "_primary")]] <- primary_estimates

  posthoc <- candidate_pairwise_contrasts(primary_fit$model, unlist(cfg$analysis$ordered_timepoints))
  posthoc[, `:=`(
    candidate_id = candidate_id,
    ontology = ontology,
    term = term,
    adjustment_family = as.character(cfg$followup$posthoc_adjustment),
    significant_bh = is.finite(p_adj_bh_within_candidate) &
      p_adj_bh_within_candidate <= as.numeric(cfg$analysis$fdr_threshold),
    inference_scope = "selection-aware post-hoc analysis of a candidate chosen from the same dataset"
  )]
  posthoc_rows[[candidate_id]] <- posthoc

  original_p <- as.numeric(original$omnibus_p_value)
  original_q <- as.numeric(original$omnibus_p_adj_bh)
  sensitivity <- data.table::data.table(
    candidate_id = candidate_id,
    ontology = ontology,
    term = term,
    scenario = "primary_unnormalized_refit",
    response_definition = "log10(sum of provided carrier-MAG values + fixed candidate pseudocount)",
    timepoint_chisq = primary_fit$chisq,
    timepoint_df = primary_fit$df,
    timepoint_p_value = primary_fit$p_value,
    singular = primary_fit$singular,
    converged = primary_fit$converged,
    convergence_messages = primary_fit$convergence_messages,
    original_screen_p_value = original_p,
    original_ontology_bh_q = original_q,
    captured_fraction_coefficient = NA_real_,
    captured_fraction_std_error = NA_real_,
    pooled_spearman_with_row_sum = suppressWarnings(stats::cor(d$potential, d$abundance_sum_as_provided, method = "spearman")),
    within_subject_spearman_with_row_sum = within_subject_correlation(d$potential, d$abundance_sum_as_provided, d$Subject),
    interpretation = "Primary-contract verification; no gut.csv row renormalization"
  )

  adjusted_fit <- fit_candidate_timepoint_model(d, cfg, pseudocount, capture_adjusted = TRUE)
  if (adjusted_fit$status == "fitted") {
    adjusted_coef <- summary(adjusted_fit$model)$coefficients["captured_log10_centered", ]
    sensitivity <- data.table::rbindlist(list(sensitivity, data.table::data.table(
      candidate_id = candidate_id,
      ontology = ontology,
      term = term,
      scenario = "captured_fraction_adjusted",
      response_definition = "primary unnormalized log10 response with centered log10 gut row sum covariate",
      timepoint_chisq = adjusted_fit$chisq,
      timepoint_df = adjusted_fit$df,
      timepoint_p_value = adjusted_fit$p_value,
      singular = adjusted_fit$singular,
      converged = adjusted_fit$converged,
      convergence_messages = adjusted_fit$convergence_messages,
      original_screen_p_value = original_p,
      original_ontology_bh_q = original_q,
      captured_fraction_coefficient = unname(adjusted_coef[["Estimate"]]),
      captured_fraction_std_error = unname(adjusted_coef[["Std. Error"]]),
      pooled_spearman_with_row_sum = suppressWarnings(stats::cor(d$potential, d$abundance_sum_as_provided, method = "spearman")),
      within_subject_spearman_with_row_sum = within_subject_correlation(d$potential, d$abundance_sum_as_provided, d$Subject),
      interpretation = "Sensitivity model only; asks whether the timepoint test persists after captured-fraction adjustment"
    )), fill = TRUE)
    adjusted_estimates <- fixed_timepoint_estimates(adjusted_fit$model, unlist(cfg$analysis$ordered_timepoints))
    adjusted_estimates[, `:=`(candidate_id = candidate_id, scenario = "captured_fraction_adjusted")]
    timepoint_estimate_rows[[paste0(candidate_id, "_adjusted")]] <- adjusted_estimates
  } else {
    sensitivity <- data.table::rbindlist(list(sensitivity, data.table::data.table(
      candidate_id = candidate_id, ontology = ontology, term = term,
      scenario = "captured_fraction_adjusted", response_definition = NA_character_,
      timepoint_chisq = NA_real_, timepoint_df = NA_real_, timepoint_p_value = NA_real_,
      singular = NA, original_screen_p_value = original_p, original_ontology_bh_q = original_q,
      converged = NA, convergence_messages = adjusted_fit$status,
      captured_fraction_coefficient = NA_real_, captured_fraction_std_error = NA_real_,
      pooled_spearman_with_row_sum = NA_real_, within_subject_spearman_with_row_sum = NA_real_,
      interpretation = adjusted_fit$status
    )), fill = TRUE)
  }

  normalized <- data.table::copy(d)
  normalized[, potential := potential / abundance_sum_as_provided]
  normalized_positive <- normalized$potential[normalized$potential > 0]
  normalized_pseudocount <- if (length(normalized_positive) > 0L) min(normalized_positive) / 2 else NA_real_
  normalized_fit <- fit_candidate_timepoint_model(normalized, cfg, normalized_pseudocount)
  if (normalized_fit$status == "fitted") {
    sensitivity <- data.table::rbindlist(list(sensitivity, data.table::data.table(
      candidate_id = candidate_id,
      ontology = ontology,
      term = term,
      scenario = "row_sum_normalized_sensitivity",
      response_definition = "log10(candidate potential divided by the sum of 441 retained MAG values + scale-specific pseudocount)",
      timepoint_chisq = normalized_fit$chisq,
      timepoint_df = normalized_fit$df,
      timepoint_p_value = normalized_fit$p_value,
      singular = normalized_fit$singular,
      converged = normalized_fit$converged,
      convergence_messages = normalized_fit$convergence_messages,
      original_screen_p_value = original_p,
      original_ontology_bh_q = original_q,
      captured_fraction_coefficient = NA_real_,
      captured_fraction_std_error = NA_real_,
      pooled_spearman_with_row_sum = suppressWarnings(stats::cor(normalized$potential, normalized$abundance_sum_as_provided, method = "spearman")),
      within_subject_spearman_with_row_sum = within_subject_correlation(normalized$potential, normalized$abundance_sum_as_provided, normalized$Subject),
      interpretation = "Secondary sensitivity only; does not replace the approved unnormalized primary analysis"
    )), fill = TRUE)
    normalized_estimates <- fixed_timepoint_estimates(normalized_fit$model, unlist(cfg$analysis$ordered_timepoints))
    normalized_estimates[, `:=`(candidate_id = candidate_id, scenario = "row_sum_normalized_sensitivity")]
    timepoint_estimate_rows[[paste0(candidate_id, "_normalized")]] <- normalized_estimates
  } else {
    sensitivity <- data.table::rbindlist(list(sensitivity, data.table::data.table(
      candidate_id = candidate_id, ontology = ontology, term = term,
      scenario = "row_sum_normalized_sensitivity", response_definition = NA_character_,
      timepoint_chisq = NA_real_, timepoint_df = NA_real_, timepoint_p_value = NA_real_,
      singular = NA, original_screen_p_value = original_p, original_ontology_bh_q = original_q,
      converged = NA, convergence_messages = normalized_fit$status,
      captured_fraction_coefficient = NA_real_, captured_fraction_std_error = NA_real_,
      pooled_spearman_with_row_sum = NA_real_, within_subject_spearman_with_row_sum = NA_real_,
      interpretation = normalized_fit$status
    )), fill = TRUE)
  }
  sensitivity_rows[[candidate_id]] <- sensitivity

  base_profile <- primary_estimates[, .(Timepoint, full_estimate_log10 = estimate_log10)]
  subjects <- sort(unique(as.character(d$Subject)))
  for (excluded_subject in subjects) {
    subset <- d[Subject != excluded_subject]
    fit <- fit_candidate_timepoint_model(subset, cfg, pseudocount)
    row <- data.table::data.table(
      candidate_id = candidate_id,
      ontology = ontology,
      term = term,
      excluded_subject = excluded_subject,
      n_subjects_remaining = data.table::uniqueN(subset$Subject),
      n_samples_remaining = nrow(subset),
      omnibus_chisq = if (fit$status == "fitted") fit$chisq else NA_real_,
      omnibus_df = if (fit$status == "fitted") fit$df else NA_real_,
      omnibus_raw_p_value = if (fit$status == "fitted") fit$p_value else NA_real_,
      raw_p_below_diagnostic_threshold = if (fit$status == "fitted")
        fit$p_value < as.numeric(cfg$followup$influence_raw_p_threshold) else NA,
      singular = if (fit$status == "fitted") fit$singular else NA,
      converged = if (fit$status == "fitted") fit$converged else NA,
      convergence_messages = if (fit$status == "fitted") fit$convergence_messages else fit$status,
      maximum_absolute_profile_change_log10 = NA_real_,
      status = fit$status,
      inference_scope = "selection-aware influence diagnostic; no ontology-wide BH refit"
    )
    if (fit$status == "fitted") {
      estimates <- fixed_timepoint_estimates(fit$model, unlist(cfg$analysis$ordered_timepoints))
      estimates <- merge(estimates, base_profile, by = "Timepoint", all.x = TRUE)
      estimates[, `:=`(
        candidate_id = candidate_id,
        excluded_subject = excluded_subject,
        profile_change_log10 = estimate_log10 - full_estimate_log10
      )]
      row[, maximum_absolute_profile_change_log10 := max(abs(estimates$profile_change_log10), na.rm = TRUE)]
      lopo_prediction_rows[[paste(candidate_id, excluded_subject, sep = "_")]] <- estimates
    }
    lopo_rows[[paste(candidate_id, excluded_subject, sep = "_")]] <- row
  }

  mag_function <- read_object(file.path(paths$intermediate, "mag_function_tables", paste0(ontology, ".rds")))
  carriers <- mag_function[term == candidate$term[[1L]], .(
    cluster_count = sum(cluster_count),
    member_gene_count = sum(member_gene_count)
  ), by = MAG]
  if (nrow(carriers) == 0L) stop("No carrier MAGs found for ", candidate_id, ".", call. = FALSE)
  absent_carrier_mags <- setdiff(carriers$MAG, aligned$gut_mags)
  if (length(absent_carrier_mags) > 0L) stop("Carrier MAGs are absent from gut.csv for ", candidate_id, ".", call. = FALSE)
  carrier_values <- data.table::melt(
    aligned$gut,
    id.vars = aligned$sample_col,
    measure.vars = carriers$MAG,
    variable.name = "MAG",
    value.name = "mag_value"
  )
  data.table::setnames(carrier_values, aligned$sample_col, "sampleID")
  carrier_values[, `:=`(candidate_id = candidate_id, ontology = ontology, term = term)]
  carrier_values <- merge(carrier_values, metadata_join, by = "sampleID", all.x = TRUE)
  carrier_values[, candidate_total := sum(mag_value), by = .(candidate_id, sampleID)]
  carrier_values[, share_of_candidate_potential := data.table::fifelse(candidate_total > 0, mag_value / candidate_total, NA_real_)]
  reconciliation <- merge(
    unique(carrier_values[, .(sampleID, candidate_total)]),
    d[, .(sampleID, potential)],
    by = "sampleID", all = TRUE
  )
  maximum_reconciliation_error <- max(abs(reconciliation$candidate_total - reconciliation$potential), na.rm = TRUE)
  if (!is.finite(maximum_reconciliation_error) || maximum_reconciliation_error > 1e-10) {
    stop("Carrier contributions do not reconstruct candidate potential for ", candidate_id, ".", call. = FALSE)
  }
  carrier_values <- merge(carrier_values, taxonomy, by.x = "MAG", by.y = tax_id, all.x = TRUE)
  contribution_rows[[candidate_id]] <- carrier_values

  contribution_summary <- carrier_values[, .(
    mean_mag_value = mean(mag_value),
    total_mag_value_across_samples = sum(mag_value),
    mean_share_of_candidate_potential = mean(share_of_candidate_potential, na.rm = TRUE),
    nonzero_samples = sum(mag_value > 0),
    n_samples = .N
  ), by = MAG]
  carrier_atlas <- merge(carriers, contribution_summary, by = "MAG", all.x = TRUE)
  carrier_atlas <- merge(carrier_atlas, taxonomy, by.x = "MAG", by.y = tax_id, all.x = TRUE)
  carrier_atlas[, `:=`(candidate_id = candidate_id, ontology = ontology, term = term)]
  carrier_rows[[candidate_id]] <- carrier_atlas
}

trajectories <- data.table::rbindlist(trajectory_rows, fill = TRUE)
carrier_atlas <- data.table::rbindlist(carrier_rows, fill = TRUE)
contributions <- data.table::rbindlist(contribution_rows, fill = TRUE)
lopo <- data.table::rbindlist(lopo_rows, fill = TRUE)
lopo_predictions <- data.table::rbindlist(lopo_prediction_rows, fill = TRUE)
posthoc <- data.table::rbindlist(posthoc_rows, fill = TRUE)
sensitivity <- data.table::rbindlist(sensitivity_rows, fill = TRUE)
timepoint_estimates <- data.table::rbindlist(timepoint_estimate_rows, fill = TRUE)

lopo_summary <- lopo[, .(
  exclusions_attempted = .N,
  exclusions_fitted = sum(status == "fitted"),
  raw_p_below_diagnostic_threshold = sum(raw_p_below_diagnostic_threshold %in% TRUE),
  minimum_raw_p = min(omnibus_raw_p_value, na.rm = TRUE),
  maximum_raw_p = max(omnibus_raw_p_value, na.rm = TRUE),
  maximum_profile_change_log10 = max(maximum_absolute_profile_change_log10, na.rm = TRUE),
  singular_refits = sum(singular %in% TRUE),
  nonconverged_refits = sum(converged %in% FALSE),
  interpretation = "Influence diagnostic after same-data candidate selection; raw p-values are not independent confirmation"
), by = .(candidate_id, ontology, term)]

taxon_sample <- contributions[, .(
  taxon_candidate_value = sum(mag_value),
  candidate_total = unique(candidate_total)
), by = .(candidate_id, ontology, term, sampleID, Subject, Timepoint, phylum, genus)]
taxon_sample[, share_of_candidate_potential := data.table::fifelse(candidate_total > 0, taxon_candidate_value / candidate_total, NA_real_)]
taxon_summary <- taxon_sample[, .(
  mean_candidate_value = mean(taxon_candidate_value),
  total_candidate_value_across_samples = sum(taxon_candidate_value),
  mean_share_of_candidate_potential = mean(share_of_candidate_potential, na.rm = TRUE),
  n_samples = data.table::uniqueN(sampleID),
  n_subjects = data.table::uniqueN(Subject)
), by = .(candidate_id, ontology, term, phylum, genus)]

table08 <- file.path(paths$targeted_validation, "Table_08_candidate_carrier_MAG_taxonomy.csv")
table09 <- file.path(paths$targeted_validation, "Table_09_leave_one_participant_out_models.csv")
table10 <- file.path(paths$targeted_validation, "Table_10_candidate_posthoc_timepoint_contrasts.csv")
table11 <- file.path(paths$targeted_validation, "Table_11_captured_MAG_fraction_sensitivity.csv")
write_table(carrier_atlas, table08)
write_table(lopo, table09)
write_table(posthoc, table10)
write_table(sensitivity, table11)
write_table(trajectories, file.path(paths$targeted_validation, "candidate_participant_trajectories.csv.gz"))
write_table(contributions, file.path(paths$targeted_validation, "candidate_sample_MAG_contributions.csv.gz"))
write_table(taxon_summary, file.path(paths$targeted_validation, "candidate_taxon_contribution_summary.csv"))
write_table(lopo_summary, file.path(paths$targeted_validation, "leave_one_participant_out_summary.csv"))
write_table(lopo_predictions, file.path(paths$targeted_validation, "leave_one_participant_out_timepoint_estimates.csv.gz"))
write_table(timepoint_estimates, file.path(paths$targeted_validation, "candidate_sensitivity_timepoint_estimates.csv"))

trajectories[, Timepoint := factor(Timepoint, levels = unlist(cfg$analysis$ordered_timepoints))]
mean_trajectories <- trajectories[, .(
  mean_potential = mean(potential),
  standard_error = stats::sd(potential) / sqrt(.N)
), by = .(candidate_id, Timepoint)]
p7 <- ggplot2::ggplot(trajectories, ggplot2::aes(Timepoint, potential, group = Subject)) +
  ggplot2::geom_line(colour = magfunc_palette[["pale"]], linewidth = 0.45, alpha = 0.9) +
  ggplot2::geom_point(colour = magfunc_palette[["muted"]], size = 0.8, alpha = 0.65) +
  ggplot2::geom_line(
    data = mean_trajectories,
    ggplot2::aes(Timepoint, mean_potential, group = 1),
    colour = magfunc_palette[["orange"]], linewidth = 0.9, inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = mean_trajectories,
    ggplot2::aes(Timepoint, mean_potential),
    shape = 21, fill = "white", colour = magfunc_palette[["orange"]], size = 1.8,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_wrap(~candidate_id, scales = "free_y", ncol = 1) +
  ggplot2::labs(
    title = "Participant trajectories for selected functional-potential candidates",
    subtitle = "Grey lines are participants; orange lines are descriptive means; primary values are not row-renormalized",
    x = "Ordered expedition timepoint", y = "Predicted functional potential"
  ) + theme_nature()
files7 <- save_publication_plot(p7, "Figure_07_candidate_participant_trajectories", cfg, width = 7.4, height = 6.2)

top_genera <- taxon_summary[order(-total_candidate_value_across_samples), head(.SD, 3L), by = candidate_id][, .(candidate_id, genus)]
plot_contributions <- merge(taxon_sample, top_genera[, .(candidate_id, genus, top_genus = TRUE)], by = c("candidate_id", "genus"), all.x = TRUE)
plot_contributions[, genus_display := data.table::fifelse(top_genus %in% TRUE & !is.na(genus), genus, "Other genera")]
plot_contributions <- plot_contributions[, .(
  mean_share = mean(share_of_candidate_potential, na.rm = TRUE)
), by = .(candidate_id, Timepoint, genus_display)]
plot_contributions[, Timepoint := factor(Timepoint, levels = unlist(cfg$analysis$ordered_timepoints))]
plot_contributions[, panel := paste(candidate_id, genus_display, sep = " — ")]
p8 <- ggplot2::ggplot(plot_contributions, ggplot2::aes(Timepoint, mean_share, group = 1)) +
  ggplot2::geom_line(colour = magfunc_palette[["blue"]], linewidth = 0.7) +
  ggplot2::geom_point(shape = 21, fill = "white", colour = magfunc_palette[["blue"]], size = 1.5) +
  ggplot2::facet_wrap(~panel, ncol = 2) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(), limits = c(0, 1),
    expand = ggplot2::expansion(mult = c(0, 0.04))
  ) +
  ggplot2::labs(
    title = "Taxonomic contributions to selected predicted functional potentials",
    subtitle = "Three largest genera plus Other; every panel uses the same 0-100% scale\nExact values are in candidate_taxon_contribution_summary.csv; values are not protein-copy weighted",
    x = "Ordered expedition timepoint", y = "Mean within-sample share"
  ) + theme_nature()
files8 <- save_publication_plot(p8, "Figure_08_candidate_taxonomic_contributions", cfg, width = 7.4, height = 9.0)

lopo[, excluded_subject := factor(excluded_subject, levels = sort(unique(excluded_subject)))]
lopo[, plot_p_value := pmax(omnibus_raw_p_value, .Machine$double.xmin)]
lopo[, model_diagnostic := data.table::fcase(
  converged %in% FALSE, "Convergence warning",
  singular %in% TRUE, "Singular fit",
  default = "No flagged issue"
)]
diagnostic_shapes <- c("No flagged issue" = 21, "Singular fit" = 24, "Convergence warning" = 4)
p9 <- ggplot2::ggplot(lopo, ggplot2::aes(excluded_subject, plot_p_value, shape = model_diagnostic)) +
  ggplot2::geom_hline(yintercept = as.numeric(cfg$followup$influence_raw_p_threshold), linetype = 2, colour = magfunc_palette[["muted"]], linewidth = 0.45) +
  ggplot2::geom_point(colour = magfunc_palette[["blue"]], size = 2) +
  ggplot2::facet_wrap(~candidate_id, ncol = 1) +
  ggplot2::scale_y_log10() +
  ggplot2::scale_shape_manual(values = diagnostic_shapes) +
  ggplot2::labs(
    title = "Leave-one-participant-out influence diagnostics",
    subtitle = "Dashed line: raw p=0.05; candidates were selected from these data\nThese p-values diagnose influence and are not independent confirmation",
    x = "Excluded participant", y = "Omnibus raw p-value", shape = "Model diagnostic"
  ) + theme_nature() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
files9 <- save_publication_plot(p9, "Figure_09_candidate_leave_one_participant_out", cfg, width = 7.4, height = 6.2)

figure_records <- data.table::rbindlist(list(
  data.table::data.table(
    result_number = "Figure 07", category = "targeted_candidate_validation",
    question = "Are the two selected trajectories distributed across participants?", artifact_type = "figure",
    path = vapply(files7, relative_to_root, character(1)),
    description = "Participant-level trajectories with descriptive means for the two same-data selected candidates.",
    stage = "06_validate_selected_candidates"
  ),
  data.table::data.table(
    result_number = "Figure 08", category = "targeted_candidate_validation",
    question = "Which carrier taxa contribute to each selected functional potential?", artifact_type = "figure",
    path = vapply(files8, relative_to_root, character(1)),
    description = "Mean within-sample genus contributions to candidate potential; values are not copy-number weighted.",
    stage = "06_validate_selected_candidates"
  ),
  data.table::data.table(
    result_number = "Figure 09", category = "targeted_candidate_validation",
    question = "Is either candidate strongly influenced by one participant?", artifact_type = "figure",
    path = vapply(files9, relative_to_root, character(1)),
    description = "Selection-aware leave-one-participant-out raw-p influence diagnostic.",
    stage = "06_validate_selected_candidates"
  )
), fill = TRUE)
table_records <- data.table::data.table(
  result_number = sprintf("Table %02d", 8:11),
  category = "targeted_candidate_validation",
  question = c(
    "Which MAGs and taxa carry the two selected functions?",
    "Does excluding one participant materially change each selected timepoint test?",
    "Which T1-T10 pairs differ for each selected candidate after within-candidate multiplicity correction?",
    "Are selected timepoint tests sensitive to the varying fraction represented by the 441 retained MAGs?"
  ),
  artifact_type = "table",
  path = vapply(c(table08, table09, table10, table11), relative_to_root, character(1)),
  description = c(
    "Exact carrier MAGs, catalogue support, taxonomy, and descriptive contribution summaries.",
    "Leave-one-participant-out mixed-model results; raw p-values are influence diagnostics, not independent confirmation.",
    "All pairwise T1-T10 fixed-effect contrasts with BH correction within each selected candidate.",
    "Primary-contract verification, captured-fraction adjustment, and explicitly secondary row-sum-normalized sensitivity."
  ),
  stage = "06_validate_selected_candidates"
)
register_results(cfg, data.table::rbindlist(list(table_records, figure_records), fill = TRUE))
writeLines(c(
  "Stage 06 complete",
  "Candidates were selected from the same Stage-03 screen; follow-up p-values are not independent replication.",
  "The primary analysis remains unnormalized. Row-sum normalization is reported only as a sensitivity analysis.",
  "Carrier contributions sum MAG values once per carrier and are not weighted by protein-copy count."
), marker)
write_session_info(cfg, "06_validate_selected_candidates")
message("Stage 06 completed.")
