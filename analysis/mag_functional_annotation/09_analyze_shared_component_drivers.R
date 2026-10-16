#!/usr/bin/env Rscript

# Purpose: Identify which expedition timepoints differ for the Stage-08 shared
# carrier component, evaluate participant influence, and attribute paired
# changes exactly to shared carrier MAGs and genera.
# Inputs: Stage-06 sample-MAG contributions and Stage-08 Tables 15-17.
# Outputs: Numbered Tables 18-22, Figures 13-15, and updated follow-up indexes.
# Dependencies: R, data.table, lme4, ggplot2, scales, yaml, common.R, followup.R.
# Usage: Rscript analysis/mag_functional_annotation/09_analyze_shared_component_drivers.R --config config/mag_functional_annotation.yml
# Workflow: Selection-aware post-hoc analysis; it is not independent
# replication, expression, activity, adaptation, or causal altitude evidence.

source("analysis/mag_functional_annotation/R/common.R")
source("analysis/mag_functional_annotation/R/followup.R")
args <- parse_stage_args()
cfg <- read_magfunc_config(args$config)
assert_packages("lme4")
paths <- output_paths(cfg)
if (!file.exists(file.path(paths$qc, "08_decompose_shared_unique_carriers.complete"))) {
  stop("Stage 08 completion marker is missing.", call. = FALSE)
}
marker <- file.path(paths$qc, "09_analyze_shared_component_drivers.complete")
guard_stage(marker, args, cfg)

timepoints <- as.character(unlist(cfg$analysis$ordered_timepoints))
fdr_threshold <- as.numeric(cfg$analysis$fdr_threshold)
candidates <- configured_followup_candidates(cfg)
reference_candidate <- candidates$candidate_id[[1L]]

membership <- read_table(file.path(
  paths$carrier_decomposition,
  "Table_15_candidate_carrier_set_membership.csv"
))
required_membership <- c("MAG", "membership_class", "phylum", "genus", "species")
missing_membership <- setdiff(required_membership, names(membership))
if (length(missing_membership) > 0L) {
  stop("Table 15 is missing columns: ", paste(missing_membership, collapse = ", "), call. = FALSE)
}
shared_taxonomy <- unique(membership[
  membership_class == "shared_by_both_candidates",
  .(MAG, phylum, genus, species)
])
if (nrow(shared_taxonomy) == 0L || anyDuplicated(shared_taxonomy$MAG)) {
  stop("Table 15 must contain unique shared carrier MAGs.", call. = FALSE)
}
shared_mags <- shared_taxonomy$MAG

components <- read_table(file.path(
  paths$carrier_decomposition,
  "Table_16_sample_component_decomposition.csv.gz"
))
required_components <- c(
  "candidate_id", "component", "sampleID", "Subject", "Timepoint",
  "analysis_date", "altitude", "days_from_baseline_subject", "component_potential"
)
missing_components <- setdiff(required_components, names(components))
if (length(missing_components) > 0L) {
  stop("Table 16 is missing columns: ", paste(missing_components, collapse = ", "), call. = FALSE)
}
shared_components <- components[component == "shared_carriers"]
if (!setequal(unique(shared_components$candidate_id), candidates$candidate_id)) {
  stop("Table 16 does not contain the shared component for every configured candidate.", call. = FALSE)
}
shared_wide <- data.table::dcast(
  shared_components,
  sampleID ~ candidate_id,
  value.var = "component_potential"
)
shared_profile_error <- max(abs(
  shared_wide[[candidates$candidate_id[[1L]]]] -
    shared_wide[[candidates$candidate_id[[2L]]]]
))
if (!is.finite(shared_profile_error) || shared_profile_error > 1e-12) {
  stop("The Stage-08 shared profiles are not identical between candidates.", call. = FALSE)
}
shared_samples <- data.table::copy(shared_components[candidate_id == reference_candidate])
if (anyDuplicated(shared_samples$sampleID)) {
  stop("The shared component must contain one row per sample.", call. = FALSE)
}

component_models <- read_table(file.path(
  paths$carrier_decomposition,
  "Table_17_component_temporal_models.csv"
))
shared_stage8 <- component_models[
  candidate_id == "Both candidates" & component == "shared_carriers"
]
if (nrow(shared_stage8) != 1L || shared_stage8$status[[1L]] != "fitted") {
  stop("Table 17 must contain one fitted shared-carrier model.", call. = FALSE)
}
pseudocount <- as.numeric(shared_stage8$pseudocount[[1L]])
if (!is.finite(pseudocount) || pseudocount <= 0) {
  positive <- shared_samples$component_potential[shared_samples$component_potential > 0]
  pseudocount <- if (length(positive) > 0L) min(positive) / 2 else NA_real_
}
if (!is.finite(pseudocount) || pseudocount <= 0) {
  stop("A positive shared-component pseudocount could not be established.", call. = FALSE)
}
shared_samples[, response_log10 := log10(component_potential + pseudocount)]
full_fit <- fit_timepoint_response_model(shared_samples, cfg, "response_log10")
if (full_fit$status != "fitted") {
  stop("The shared-component post-hoc model failed: ", full_fit$status, call. = FALSE)
}

pairwise <- candidate_pairwise_contrasts(full_fit$model, timepoints)
data.table::setnames(
  pairwise,
  "p_adj_bh_within_candidate",
  "p_adj_bh_across_45_timepoint_contrasts"
)
pairwise[, `:=`(
  contrast_id = paste0(before_timepoint, "_vs_", after_timepoint),
  significant_bh = is.finite(p_adj_bh_across_45_timepoint_contrasts) &
    p_adj_bh_across_45_timepoint_contrasts <= fdr_threshold,
  shared_omnibus_chisq = as.numeric(shared_stage8$omnibus_chisq[[1L]]),
  shared_omnibus_df = as.numeric(shared_stage8$omnibus_df[[1L]]),
  shared_omnibus_p_value = as.numeric(shared_stage8$omnibus_p_value[[1L]]),
  shared_omnibus_bh_q = as.numeric(shared_stage8$p_adj_bh_within_family[[1L]]),
  shared_omnibus_bh_supported = isTRUE(shared_stage8$significant_bh[[1L]]),
  adjustment_family = as.character(cfg$followup$shared_component_posthoc_adjustment),
  inference_scope = paste(
    "conditional post-hoc analysis of a shared component selected from the same dataset;",
    "not independent replication"
  )
)]
data.table::setcolorder(pairwise, c(
  "contrast_id", "before_timepoint", "after_timepoint",
  setdiff(names(pairwise), c("contrast_id", "before_timepoint", "after_timepoint"))
))
data.table::setorder(pairwise, p_adj_bh_across_45_timepoint_contrasts, p_value)

base_profile <- fixed_timepoint_estimates(full_fit$model, timepoints)[
  , .(Timepoint, full_estimate_log10 = estimate_log10)
]
original_significant <- pairwise[significant_bh %in% TRUE]
original_significant_keys <- original_significant[, paste(before_timepoint, after_timepoint, sep = "|")]
lopo_rows <- list()
lopo_contrast_rows <- list()
lopo_profile_rows <- list()
for (excluded_subject in sort(unique(as.character(shared_samples$Subject)))) {
  subset <- shared_samples[Subject != excluded_subject]
  fit <- fit_timepoint_response_model(subset, cfg, "response_log10")
  row <- data.table::data.table(
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
    significant_contrasts_after_exclusion = NA_integer_,
    original_significant_contrasts = length(original_significant_keys),
    original_significant_contrasts_retained = NA_integer_,
    proportion_original_significant_retained = NA_real_,
    original_significant_directions_preserved = NA,
    status = fit$status,
    inference_scope = "selection-aware participant-influence diagnostic; raw p-values are not independent confirmation"
  )
  if (fit$status == "fitted") {
    estimates <- fixed_timepoint_estimates(fit$model, timepoints)
    estimates <- merge(estimates, base_profile, by = "Timepoint", all.x = TRUE)
    estimates[, `:=`(
      excluded_subject = excluded_subject,
      profile_change_log10 = estimate_log10 - full_estimate_log10
    )]
    row[, maximum_absolute_profile_change_log10 := max(abs(estimates$profile_change_log10))]
    lopo_profile_rows[[excluded_subject]] <- estimates

    contrasts <- candidate_pairwise_contrasts(fit$model, timepoints)
    data.table::setnames(
      contrasts,
      "p_adj_bh_within_candidate",
      "p_adj_bh_across_45_timepoint_contrasts"
    )
    contrasts[, `:=`(
      excluded_subject = excluded_subject,
      contrast_key = paste(before_timepoint, after_timepoint, sep = "|"),
      significant_bh = is.finite(p_adj_bh_across_45_timepoint_contrasts) &
        p_adj_bh_across_45_timepoint_contrasts <= fdr_threshold
    )]
    comparison <- merge(
      contrasts,
      pairwise[, .(
        before_timepoint, after_timepoint,
        full_estimate_log10 = estimate_log10_after_minus_before,
        full_significant_bh = significant_bh
      )],
      by = c("before_timepoint", "after_timepoint"),
      all.x = TRUE
    )
    comparison[, direction_preserved :=
      sign(estimate_log10_after_minus_before) == sign(full_estimate_log10)]
    lopo_contrast_rows[[excluded_subject]] <- comparison
    retained <- comparison[full_significant_bh %in% TRUE & significant_bh %in% TRUE, .N]
    direction_preserved <- if (length(original_significant_keys) > 0L) {
      all(comparison[full_significant_bh %in% TRUE, direction_preserved] %in% TRUE)
    } else {
      NA
    }
    row[, `:=`(
      significant_contrasts_after_exclusion = sum(comparison$significant_bh %in% TRUE),
      original_significant_contrasts_retained = retained,
      proportion_original_significant_retained = if (length(original_significant_keys) > 0L)
        retained / length(original_significant_keys) else NA_real_,
      original_significant_directions_preserved = direction_preserved
    )]
  }
  lopo_rows[[excluded_subject]] <- row
}
lopo <- data.table::rbindlist(lopo_rows, fill = TRUE)
lopo_contrasts <- data.table::rbindlist(lopo_contrast_rows, fill = TRUE)
lopo_profiles <- data.table::rbindlist(lopo_profile_rows, fill = TRUE)

contributions <- read_table(file.path(
  paths$targeted_validation,
  "candidate_sample_MAG_contributions.csv.gz"
))
required_contributions <- c(
  "candidate_id", "MAG", "sampleID", "Subject", "Timepoint", "mag_value"
)
missing_contributions <- setdiff(required_contributions, names(contributions))
if (length(missing_contributions) > 0L) {
  stop("Stage-06 contribution data are missing columns: ", paste(missing_contributions, collapse = ", "), call. = FALSE)
}
shared_values <- contributions[
  candidate_id == reference_candidate & MAG %in% shared_mags,
  .(MAG, sampleID, Subject, Timepoint = as.character(Timepoint), mag_value)
]
if (anyDuplicated(shared_values, by = c("MAG", "sampleID"))) {
  stop("Shared MAG contribution rows must be unique by MAG and sample.", call. = FALSE)
}
if (!setequal(unique(shared_values$MAG), shared_mags)) {
  stop("Stage-06 contribution data do not contain every Stage-08 shared MAG.", call. = FALSE)
}
shared_values <- merge(shared_values, shared_taxonomy, by = "MAG", all.x = TRUE)
shared_reconstruction <- shared_values[, .(
  reconstructed_shared_potential = sum(mag_value)
), by = sampleID]
shared_reconstruction <- merge(
  shared_reconstruction,
  shared_samples[, .(sampleID, component_potential)],
  by = "sampleID", all = TRUE
)
maximum_shared_reconstruction_error <- max(abs(
  shared_reconstruction$reconstructed_shared_potential -
    shared_reconstruction$component_potential
))
if (!is.finite(maximum_shared_reconstruction_error) ||
    maximum_shared_reconstruction_error > 1e-10) {
  stop("Shared MAG values do not reconstruct the Stage-08 shared component.", call. = FALSE)
}

genus_values <- data.table::copy(shared_values)
genus_values[is.na(phylum) | !nzchar(phylum), phylum := "p__Unclassified"]
genus_values[is.na(genus) | !nzchar(genus), genus := "g__Unclassified"]
genus_values <- genus_values[, .(
  mag_value = sum(mag_value)
), by = .(sampleID, Subject, Timepoint, phylum, genus)]

mag_driver_rows <- list()
genus_driver_rows <- list()
for (contrast_index in seq_len(nrow(pairwise))) {
  before_tp <- pairwise$before_timepoint[[contrast_index]]
  after_tp <- pairwise$after_timepoint[[contrast_index]]
  contrast_name <- pairwise$contrast_id[[contrast_index]]
  mag_driver <- paired_change_drivers(
    shared_values, before_tp, after_tp,
    c("MAG", "phylum", "genus", "species")
  )
  genus_driver <- paired_change_drivers(
    genus_values, before_tp, after_tp,
    c("phylum", "genus")
  )
  if (nrow(mag_driver) == 0L || nrow(genus_driver) == 0L) {
    stop("No paired driver rows were available for ", contrast_name, ".", call. = FALSE)
  }
  mag_driver[, contrast_id := contrast_name]
  genus_driver[, contrast_id := contrast_name]
  mag_driver_rows[[contrast_name]] <- mag_driver
  genus_driver_rows[[contrast_name]] <- genus_driver
}
mag_drivers <- data.table::rbindlist(mag_driver_rows, fill = TRUE)
genus_drivers <- data.table::rbindlist(genus_driver_rows, fill = TRUE)
contrast_statistics <- pairwise[, .(
  contrast_id, before_timepoint, after_timepoint,
  estimate_log10_after_minus_before, ratio_after_vs_before, p_value,
  p_adj_bh_across_45_timepoint_contrasts, significant_bh
)]
mag_drivers <- merge(
  mag_drivers,
  contrast_statistics,
  by = c("contrast_id", "before_timepoint", "after_timepoint"),
  all.x = TRUE
)
genus_drivers <- merge(
  genus_drivers,
  contrast_statistics,
  by = c("contrast_id", "before_timepoint", "after_timepoint"),
  all.x = TRUE
)
mag_driver_totals <- mag_drivers[, .(
  mag_net_change = unique(net_change),
  mag_total_absolute_change = unique(total_absolute_change)
), by = contrast_id]
genus_driver_totals <- genus_drivers[, .(
  genus_net_change = unique(net_change),
  genus_total_absolute_change = unique(total_absolute_change)
), by = contrast_id]
driver_reconciliation <- merge(mag_driver_totals, genus_driver_totals, by = "contrast_id")
maximum_driver_reconciliation_error <- max(abs(
  driver_reconciliation$mag_net_change - driver_reconciliation$genus_net_change
))
if (!is.finite(maximum_driver_reconciliation_error) ||
    maximum_driver_reconciliation_error > 1e-10) {
  stop("MAG and genus paired-change decompositions do not reconcile.", call. = FALSE)
}

driver_concentration <- mag_drivers[order(contrast_id, -absolute_change_share), {
  shares <- absolute_change_share[is.finite(absolute_change_share)]
  .(
    n_pairs = unique(n_pairs),
    net_shared_component_change = unique(net_change),
    total_absolute_MAG_change = unique(total_absolute_change),
    cancellation_fraction = if (unique(total_absolute_change) > 0)
      1 - abs(unique(net_change)) / unique(total_absolute_change) else NA_real_,
    top_1_MAG_absolute_change_share = sum(utils::head(shares, 1L)),
    top_3_MAG_absolute_change_share = sum(utils::head(shares, 3L)),
    top_5_MAG_absolute_change_share = sum(utils::head(shares, 5L)),
    top_10_MAG_absolute_change_share = sum(utils::head(shares, 10L)),
    effective_MAG_driver_count = if (sum(shares^2) > 0) 1 / sum(shares^2) else NA_real_
  )
}, by = .(
  contrast_id, before_timepoint, after_timepoint,
  estimate_log10_after_minus_before, ratio_after_vs_before, p_value,
  p_adj_bh_across_45_timepoint_contrasts, significant_bh
)]
data.table::setorder(driver_concentration, p_adj_bh_across_45_timepoint_contrasts, p_value)

table18 <- file.path(paths$shared_component_followup, "Table_18_shared_component_timepoint_contrasts.csv")
table19 <- file.path(paths$shared_component_followup, "Table_19_shared_component_leave_one_participant_out.csv")
table20 <- file.path(paths$shared_component_followup, "Table_20_shared_component_MAG_drivers.csv.gz")
table21 <- file.path(paths$shared_component_followup, "Table_21_shared_component_genus_drivers.csv.gz")
table22 <- file.path(paths$shared_component_followup, "Table_22_shared_component_driver_concentration.csv")
write_table(pairwise, table18)
write_table(lopo, table19)
write_table(mag_drivers, table20)
write_table(genus_drivers, table21)
write_table(driver_concentration, table22)
write_table(lopo_contrasts, file.path(
  paths$shared_component_followup,
  "shared_component_leave_one_participant_out_contrasts.csv.gz"
))
write_table(lopo_profiles, file.path(
  paths$shared_component_followup,
  "shared_component_leave_one_participant_out_timepoint_estimates.csv.gz"
))

heatmap_data <- data.table::copy(pairwise)
heatmap_data[, `:=`(
  before_timepoint = factor(before_timepoint, levels = rev(timepoints)),
  after_timepoint = factor(after_timepoint, levels = timepoints),
  significance_label = data.table::fifelse(significant_bh %in% TRUE, "*", "")
)]
fill_limit <- max(abs(heatmap_data$estimate_log10_after_minus_before), na.rm = TRUE)
p13 <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(after_timepoint, before_timepoint, fill = estimate_log10_after_minus_before)
) +
  ggplot2::geom_tile(colour = "grey85", linewidth = 0.35) +
  ggplot2::geom_text(ggplot2::aes(label = significance_label), size = 4, colour = magfunc_palette[["ink"]]) +
  ggplot2::scale_fill_gradient2(
    low = magfunc_palette[["blue"]], mid = "white", high = magfunc_palette[["orange"]],
    midpoint = 0, limits = c(-fill_limit, fill_limit)
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(
    title = "Shared-component pairwise expedition-timepoint contrasts",
    subtitle = paste(
      "Orange: after timepoint higher; blue: after timepoint lower",
      "* BH q<0.05 across all 45 contrasts; values use the log10 model scale",
      sep = "\n"
    ),
    x = "After timepoint", y = "Before timepoint", fill = "Log10\nchange"
  ) + theme_nature() +
  ggplot2::theme(panel.grid = ggplot2::element_blank())
files13 <- save_publication_plot(
  p13, "Figure_13_shared_component_timepoint_contrasts", cfg,
  width = 7.4, height = 6.2
)

lopo_plot <- data.table::copy(lopo)
lopo_plot[, `:=`(
  excluded_subject = factor(excluded_subject, levels = sort(unique(excluded_subject))),
  plot_p_value = pmax(omnibus_raw_p_value, .Machine$double.xmin),
  retention_label = data.table::fifelse(
    is.finite(proportion_original_significant_retained),
    paste0(original_significant_contrasts_retained, "/", original_significant_contrasts),
    "NA"
  ),
  model_diagnostic = data.table::fcase(
    converged %in% FALSE, "Convergence warning",
    singular %in% TRUE, "Singular fit",
    default = "No flagged issue"
  )
)]
diagnostic_shapes <- c("No flagged issue" = 21, "Singular fit" = 24, "Convergence warning" = 4)
p14 <- ggplot2::ggplot(
  lopo_plot,
  ggplot2::aes(excluded_subject, plot_p_value, shape = model_diagnostic)
) +
  ggplot2::geom_hline(
    yintercept = as.numeric(cfg$followup$influence_raw_p_threshold),
    linetype = 2, colour = magfunc_palette[["muted"]], linewidth = 0.45
  ) +
  ggplot2::geom_point(colour = magfunc_palette[["blue"]], size = 2) +
  ggplot2::geom_text(ggplot2::aes(label = retention_label), vjust = -0.75, size = 2.5, show.legend = FALSE) +
  ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
  ggplot2::scale_shape_manual(values = diagnostic_shapes) +
  ggplot2::labs(
    title = "Shared-component leave-one-participant-out diagnostics",
    subtitle = paste(
      "Labels show retained/full-data BH-significant contrasts; dashed line is raw p=0.05",
      "P-values diagnose influence and are not independent confirmation",
      sep = "\n"
    ),
    x = "Excluded participant", y = "Omnibus raw p-value", shape = "Model diagnostic"
  ) + theme_nature() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
files14 <- save_publication_plot(
  p14, "Figure_14_shared_component_leave_one_participant_out", cfg,
  width = 7.4, height = 5.4
)

focus_supported <- nrow(original_significant) > 0L
focus_contrast <- if (focus_supported) original_significant[1L] else pairwise[1L]
focus_rank_label <- if (focus_supported) "strongest BH-supported" else "top-ranked exploratory"
focus_drivers <- data.table::copy(genus_drivers[contrast_id == focus_contrast$contrast_id[[1L]]])
focus_drivers[, absolute_mean_delta := abs(mean_delta)]
data.table::setorder(focus_drivers, -absolute_mean_delta)
focus_q <- focus_contrast$p_adj_bh_across_45_timepoint_contrasts[[1L]]
focus_q_label <- if (is.finite(focus_q) && focus_q < 0.001) {
  "<0.001"
} else {
  paste0("=", scales::number(focus_q, accuracy = 0.001))
}
top_count <- min(8L, nrow(focus_drivers))
focus_plot <- focus_drivers[seq_len(top_count), .(genus, mean_delta)]
if (nrow(focus_drivers) > top_count) {
  focus_plot <- data.table::rbindlist(list(
    focus_plot,
    data.table::data.table(
      genus = "Other shared-carrier genera",
      mean_delta = sum(focus_drivers[(top_count + 1L):.N, mean_delta])
    )
  ))
}
focus_plot[, `:=`(
  direction = data.table::fifelse(mean_delta >= 0, "Higher after", "Lower after"),
  genus = factor(genus, levels = genus[order(mean_delta)])
)]
direction_colours <- c(
  `Higher after` = magfunc_palette[["orange"]],
  `Lower after` = magfunc_palette[["blue"]]
)
p15 <- ggplot2::ggplot(focus_plot, ggplot2::aes(mean_delta, genus, fill = direction)) +
  ggplot2::geom_vline(xintercept = 0, colour = magfunc_palette[["ink"]], linewidth = 0.4) +
  ggplot2::geom_col(width = 0.68, colour = "white", linewidth = 0.2) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf("%+.4f", mean_delta),
      hjust = data.table::fifelse(mean_delta >= 0, -0.08, 1.08)
    ),
    size = 2.8
  ) +
  ggplot2::scale_fill_manual(values = direction_colours) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.14, 0.18))) +
  ggplot2::labs(
    title = paste("Genus drivers of the", focus_rank_label, "contrast"),
    subtitle = paste0(
      "Shared component; ",
      focus_contrast$before_timepoint[[1L]], " to ", focus_contrast$after_timepoint[[1L]],
      "; paired subjects; BH q", focus_q_label,
      "\nBars sum exactly to the mean paired change in the unnormalized shared component"
    ),
    x = "Mean paired change in provided MAG value (after - before)",
    y = NULL, fill = NULL
  ) + theme_nature() +
  ggplot2::theme(legend.position = "top")
files15 <- save_publication_plot(
  p15, "Figure_15_shared_component_genus_drivers", cfg,
  width = 7.4, height = 5.6
)

write_table(heatmap_data, file.path(paths$reproducibility, "Figure_13_plot_data.csv"))
write_table(lopo_plot, file.path(paths$reproducibility, "Figure_14_plot_data.csv"))
write_table(focus_plot, file.path(paths$reproducibility, "Figure_15_plot_data.csv"))
chart_map_path <- file.path(paths$reproducibility, "followup_publication_chart_map.csv")
chart_map <- if (file.exists(chart_map_path)) read_table(chart_map_path) else data.table::data.table()
new_chart_map <- data.table::data.table(
  figure = sprintf("Figure %02d", 13:15),
  analytical_question = c(
    "Which shared-component timepoint pairs differ after correction across all 45 comparisons?",
    "Is the shared-component timepoint result strongly influenced by one participant?",
    "Which genera supply the paired change for the top-ranked shared-component contrast, with its BH status shown?"
  ),
  chart_family = c("matrix and comparison", "uncertainty and influence", "decomposition and contribution"),
  visual_encoding = c(
    "signed log10 after-minus-before effect with stars for BH-supported contrasts",
    "raw omnibus p-value by excluded participant with retained-contrast labels",
    "signed mean paired genus changes summing to the shared-component change"
  ),
  denominator_or_unit = c(
    "45 fixed-effect T1-T10 contrasts from the participant-random-intercept model",
    "13 leave-one-participant-out model refits; raw p-values are diagnostic",
    "unnormalized gut.csv MAG values among subjects paired for the focus contrast"
  ),
  palette_policy = c(
    "blue negative and orange positive around a white zero midpoint; stars add non-colour evidence",
    "single blue root with diagnostic point shapes and a neutral reference line",
    "blue negative and orange positive with direct signed labels and a zero line"
  )
)
chart_map <- data.table::rbindlist(list(chart_map, new_chart_map), fill = TRUE)
chart_map <- unique(chart_map, by = "figure", fromLast = TRUE)
data.table::setorder(chart_map, figure)
write_table(chart_map, chart_map_path)

table_records <- data.table::data.table(
  result_number = sprintf("Table %02d", 18:22),
  category = "shared_component_followup",
  question = c(
    "Which T1-T10 shared-component pairs differ after correction across all 45 contrasts?",
    "Does excluding one participant materially change the shared-component result?",
    "Which shared carrier MAGs supply each paired timepoint change?",
    "Which shared carrier genera supply each paired timepoint change?",
    "Are paired shared-component changes concentrated in a few MAGs or broadly distributed?"
  ),
  artifact_type = "table",
  path = vapply(c(table18, table19, table20, table21, table22), relative_to_root, character(1)),
  description = c(
    "All 45 shared-component fixed-effect contrasts with one BH correction family.",
    "Selection-aware leave-one-participant-out influence diagnostics with contrast-retention summaries.",
    "Exact participant-paired MAG-level decomposition on the provided unnormalized value scale.",
    "Exact participant-paired genus-level decomposition on the provided unnormalized value scale.",
    "Per-contrast MAG driver concentration, cancellation, and effective-driver summaries."
  ),
  stage = "09_analyze_shared_component_drivers"
)
figure_records <- data.table::rbindlist(list(
  data.table::data.table(
    result_number = "Figure 13", category = "shared_component_followup",
    question = "Which timepoint pairs differ for the shared component?", artifact_type = "figure",
    path = vapply(files13, relative_to_root, character(1)),
    description = "Complete signed 45-contrast heatmap with BH-supported cells marked by stars.",
    stage = "09_analyze_shared_component_drivers"
  ),
  data.table::data.table(
    result_number = "Figure 14", category = "shared_component_followup",
    question = "Is the shared-component evidence participant-sensitive?", artifact_type = "figure",
    path = vapply(files14, relative_to_root, character(1)),
    description = "Leave-one-participant-out raw-p influence diagnostic with retained-contrast counts.",
    stage = "09_analyze_shared_component_drivers"
  ),
  data.table::data.table(
    result_number = "Figure 15", category = "shared_component_followup",
    question = "Which genera drive the strongest supported shared-component contrast?", artifact_type = "figure",
    path = vapply(files15, relative_to_root, character(1)),
    description = "Exact signed paired genus decomposition for the strongest BH-supported contrast.",
    stage = "09_analyze_shared_component_drivers"
  )
), fill = TRUE)
register_results(cfg, data.table::rbindlist(list(table_records, figure_records), fill = TRUE))

followup_csv <- file.path(paths$root, "FOLLOWUP_RESULTS_INDEX.csv")
followup_md <- file.path(paths$root, "FOLLOWUP_RESULTS_INDEX.md")
stage10_marker <- file.path(paths$qc, "10_audit_candidate_annotations.complete")
if (file.exists(stage10_marker)) {
  message("Stage 10 is already complete; preserving its Tables 08-25 follow-up indexes.")
} else {
  registry <- read_table(file.path(paths$intermediate, "results_registry.csv"))
  followup_registry <- registry[
    result_number %in% c(sprintf("Table %02d", 8:22), sprintf("Figure %02d", 7:15))
  ]
  data.table::setorder(followup_registry, result_number, path)
  write_table(followup_registry, followup_csv)
  escape_md <- function(x) gsub("\\|", "\\\\|", as.character(x))
  lines <- c(
    "# MAG Functional Annotation Follow-up Results Index", "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
    "## Analysis contract", "",
    "- The shared component was selected after the same-data candidate screen; Stage-09 p-values are not independent replication.",
    "- All 45 T1-T10 shared-component contrasts form one BH correction family.",
    "- Driver decompositions use paired subjects and reconcile exactly on the provided unnormalized MAG-value scale.",
    "- Leave-one-participant-out raw p-values diagnose influence; they are not new confirmatory tests.",
    "- Results describe predicted genomic potential, not expression, activity, pathway completion, adaptation, or causal altitude effects.", "",
    "Number | Category | Type | Path | Purpose", "--- | --- | --- | --- | ---"
  )
  for (row_index in seq_len(nrow(followup_registry))) {
    lines <- c(lines, paste(
      escape_md(followup_registry$result_number[[row_index]]),
      escape_md(followup_registry$category[[row_index]]),
      escape_md(followup_registry$artifact_type[[row_index]]),
      paste0("`", escape_md(followup_registry$path[[row_index]]), "`"),
      escape_md(followup_registry$description[[row_index]]),
      sep = " | "
    ))
  }
  lines <- c(lines, "", "## Interpretation order", "",
    "1. Read Table 11 before interpreting the robustness of either original selected candidate.",
    "2. Read Table 17 to establish that temporal evidence resides in the shared-carrier component.",
    "3. Use Table 18 and Figure 13 to identify the BH-supported shared-component timepoint pairs.",
    "4. Use Table 19 and Figure 14 to assess participant influence.",
    "5. Use Tables 20-22 and Figure 15 to distinguish broad restructuring from a small number of taxonomic drivers.",
    "6. Treat all Stage-09 findings as selection-aware same-data explanation, not independent replication."
  )
  writeLines(lines, followup_md)
}

writeLines(c(
  "Stage 09 complete",
  paste0("Shared carrier MAGs: ", length(shared_mags)),
  paste0("T1-T10 pairwise contrasts tested: ", nrow(pairwise)),
  paste0("BH-supported pairwise contrasts: ", sum(pairwise$significant_bh %in% TRUE)),
  paste0("Leave-one-participant-out refits fitted: ", sum(lopo$status == "fitted"), "/", nrow(lopo)),
  paste0("Maximum Stage-08 shared reconstruction error: ", format(maximum_shared_reconstruction_error, scientific = TRUE)),
  paste0("Maximum MAG-versus-genus driver reconciliation error: ", format(maximum_driver_reconciliation_error, scientific = TRUE)),
  "All Stage-09 inference is selection-aware and uses unnormalized gut.csv values."
), marker)
write_session_info(cfg, "09_analyze_shared_component_drivers")
message("Stage 09 completed: ", followup_md)
