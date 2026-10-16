#!/usr/bin/env Rscript
source("analysis/functional_community_dynamics/R/common.R")

args <- parse_dynamics_args()
cfg <- read_dynamics_config(args$config)
assert_dynamics_packages()
paths <- dynamics_paths(cfg)
marker <- file.path(paths$qc, "03_analyze_longitudinal_stability.complete")
guard_dynamics_stage(marker, args, cfg)
prepared_file <- file.path(paths$intermediate, "prepared_profiles.rds")
distance_file <- file.path(paths$intermediate, "distance_objects.rds")
if (!file.exists(prepared_file) || !file.exists(distance_file)) stop("Stages 01 and 02 must complete first.", call. = FALSE)
prepared <- readRDS(prepared_file)
distances <- readRDS(distance_file)
metadata <- copy(prepared$metadata)

records <- list()
recovery <- list()
for (key in names(distances)) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  layer <- parts[[1]]
  metric <- parts[[2]]
  distance <- distances[[key]]
  scale <- median(as.numeric(distance))
  if (!is.finite(scale) || scale <= 0) stop("Cannot standardize zero or non-finite median distance for ", key, ".", call. = FALSE)
  extracted <- extract_longitudinal_records(distance, metadata, layer, metric)
  rows <- rbind(extracted$baseline, extracted$adjacent, fill = TRUE)
  rows[, `:=`(distance_scale_median_all_pairs = scale, standardized_distance = distance / scale)]
  records[[key]] <- rows
  rec <- recovery_records(distance, metadata, layer, metric)
  if (nrow(rec)) {
    rec[, `:=`(distance_scale_median_all_pairs = scale, standardized_distance = distance / scale)]
    recovery[[key]] <- rec
  }
}
longitudinal <- rbindlist(c(records, recovery), fill = TRUE)
write_dynamics_table(longitudinal, file.path(paths$longitudinal, "Table_29_participant_distance_reference_and_adjacent_transitions.csv"))

adjacent <- longitudinal[record_type == "adjacent_observed_transition" & metric == cfg$analysis$primary_metric]
tax <- adjacent[layer == "taxonomic_MAG", .(
  subject_id, from_sample, to_sample, from_timepoint, to_timepoint, transition_days,
  taxonomic_distance = distance, taxonomic_standardized_distance = standardized_distance
)]
functional_layers <- setdiff(names(prepared$matrices), "taxonomic_MAG")
model_rows <- list()
matched_rows <- list()
sensitivity_rows <- list()
sensitivity_participant_rows <- list()

for (layer_id in functional_layers) {
  fun <- adjacent[layer == layer_id, .(
    subject_id, from_sample, to_sample, functional_distance = distance,
    functional_standardized_distance = standardized_distance
  )]
  matched <- merge(tax, fun, by = c("subject_id", "from_sample", "to_sample"), all = FALSE)
  matched[, `:=`(
    layer = layer_id,
    log_functional_movement = log1p(functional_standardized_distance),
    log_taxonomic_movement = log1p(taxonomic_standardized_distance),
    log_functional_to_taxonomic_ratio = log((functional_standardized_distance + 1e-8) / (taxonomic_standardized_distance + 1e-8))
  )]
  transition_sd <- stats::sd(matched$transition_days)
  matched[, transition_days_scaled := if (is.finite(transition_sd) && transition_sd > 0) {
    as.numeric(scale(transition_days))
  } else {
    0
  }]
  matched_rows[[layer_id]] <- matched

  layer_role <- prepared$layer_metadata[layer == layer_id, role]
  sensitivity <- participant_transition_sensitivity(matched, layer_id, layer_role)
  sensitivity_rows[[layer_id]] <- sensitivity$summary
  sensitivity_participant_rows[[layer_id]] <- sensitivity$participant

  coupling <- fit_dynamics_lmer(log_functional_movement ~ log_taxonomic_movement + transition_days_scaled + (1 | subject_id), matched)
  if (coupling$status == "fitted") {
    coef <- summary(coupling$fit)$coefficients
    estimate <- unname(coef["log_taxonomic_movement", "Estimate"])
    se <- unname(coef["log_taxonomic_movement", "Std. Error"])
    z <- estimate / se
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  } else estimate <- se <- p <- NA_real_
  model_rows[[length(model_rows) + 1L]] <- data.table(
    layer = layer_id, role = layer_role, model_family = "taxonomic_functional_transition_coupling",
    response = "log1p standardized functional adjacent-transition distance",
    estimate = estimate, standard_error = se, ratio = NA_real_, p_value = p,
    transitions = nrow(matched), subjects = uniqueN(matched$subject_id), singular = coupling$singular,
    converged = coupling$converged, convergence_messages = coupling$messages, status = coupling$status,
    interpretation = "Positive estimates indicate larger functional movement during transitions with larger taxonomic movement."
  )

  smoothing <- fit_dynamics_lmer(log_functional_to_taxonomic_ratio ~ 1 + (1 | subject_id), matched)
  if (smoothing$status == "fitted") {
    coef <- summary(smoothing$fit)$coefficients
    estimate <- unname(coef[1L, "Estimate"])
    se <- unname(coef[1L, "Std. Error"])
    p <- 2 * stats::pnorm(abs(estimate / se), lower.tail = FALSE)
  } else estimate <- se <- p <- NA_real_
  model_rows[[length(model_rows) + 1L]] <- data.table(
    layer = layer_id, role = layer_role, model_family = "functional_smoothing_relative_to_taxonomy",
    response = "log ratio of standardized functional to taxonomic adjacent-transition distance",
    estimate = estimate, standard_error = se, ratio = exp(estimate), p_value = p,
    transitions = nrow(matched), subjects = uniqueN(matched$subject_id), singular = smoothing$singular,
    converged = smoothing$converged, convergence_messages = smoothing$messages, status = smoothing$status,
    interpretation = "Ratios below one indicate smaller standardized functional than taxonomic movement; this supports smoothing, not proof of resilience."
  )
}

models <- rbindlist(model_rows, fill = TRUE)
models[, p_adj_bh_within_model_family := NA_real_]
models[role != "descriptive", p_adj_bh_within_model_family := safe_bh(p_value), by = model_family]
models[, significant_bh := role != "descriptive" & status == "fitted" & p_adj_bh_within_model_family < as.numeric(cfg$analysis$fdr_threshold)]
write_dynamics_table(models, file.path(paths$longitudinal, "Table_30_matched_taxonomic_functional_stability_models.csv"))
sensitivity <- rbindlist(sensitivity_rows, fill = TRUE)
sensitivity[, p_adj_bh_within_sensitivity_family := NA_real_]
sensitivity[role != "descriptive", p_adj_bh_within_sensitivity_family := safe_bh(p_value), by = sensitivity_family]
sensitivity[, significant_bh :=
  role != "descriptive" &
  status == "fitted" &
  p_adj_bh_within_sensitivity_family < as.numeric(cfg$analysis$fdr_threshold)
]
write_dynamics_table(
  sensitivity,
  file.path(paths$longitudinal, "Table_31_participant_level_transition_sensitivity.csv")
)
write_dynamics_table(
  rbindlist(sensitivity_participant_rows, fill = TRUE),
  file.path(paths$reproducibility, "participant_level_transition_sensitivity.csv")
)
write_dynamics_table(rbindlist(matched_rows, fill = TRUE), file.path(paths$reproducibility, "Figure_18_plot_data.csv"))
write_dynamics_table(longitudinal[record_type == "distance_to_earliest_observation" & metric == cfg$analysis$primary_metric], file.path(paths$reproducibility, "Figure_17_plot_data.csv"))
write_dynamics_table(longitudinal[record_type == "exploratory_t1_recovery" & metric == cfg$analysis$primary_metric], file.path(paths$reproducibility, "Figure_19_plot_data.csv"))
write_dynamics_session(cfg, "03_analyze_longitudinal_stability")
writeLines(c(
  "Stage 03 complete", paste("Adjacent matched transitions per functional layer:", nrow(tax)),
  paste("Participant-level sensitivity summaries:", nrow(sensitivity)),
  paste("Participants with paired T1 and T10:", uniqueN(longitudinal[record_type == "exploratory_t1_recovery", subject_id])),
  "Distances are standardized by each layer's median all-sample pairwise distance.",
  "Table 31 treats each participant as one unit using signed-rank tests of participant-level summaries.",
  "T1-to-T10 recovery is exploratory and cannot be generalized beyond paired participants."
), marker)
