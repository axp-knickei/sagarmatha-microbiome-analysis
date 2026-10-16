#!/usr/bin/env Rscript

# Purpose: Quantify internal reproducibility of the exploratory low/initial
# versus high-altitude NetMoss comparison at genus and family levels.
# Inputs: raw gut abundance, metadata, taxonomy, and the dedicated stability YAML.
# Outputs: scenario, edge, module, NetMoss-rank, participant-influence, failure,
# and reproducibility tables under a separate result directory.
# Dependencies: yaml, NetMoss2, rsparcc for production, and base parallel.
# Invocation: run from the repository root; full execution requires an explicit
# exploratory acknowledgement. This workflow does not validate altitude causality.

suppressPackageStartupMessages({
  library(yaml)
  library(dynamicTreeCut)
  library(WGCNA)
})
source("scripts/netmoss/altitude_stability_helpers.R")

usage <- function() {
  cat(
    "Usage: Rscript scripts/netmoss/run_altitude_stability.R [options]\n",
    "  --config PATH       Default: config/netmoss_altitude_stability.yml\n",
    "  --output-dir PATH   Override configured output directory\n",
    "  --threads N         Parallel resampling workers\n",
    "  --bootstrap N       Override bootstrap replicate count\n",
    "  --prepare-only      Validate cohorts and filters without fitting networks\n",
    "  --test-mode         Use repository test data and a small software-only setup\n",
    "  --acknowledge-exploratory  Required for full execution\n",
    "  --force             Replace existing named outputs\n",
    "  --help\n",
    sep = ""
  )
}

parse_args <- function(arguments) {
  parsed <- list(
    config = "config/netmoss_altitude_stability.yml",
    output_dir = NULL,
    threads = NULL,
    bootstrap = NULL,
    prepare_only = FALSE,
    test_mode = FALSE,
    acknowledge = FALSE,
    force = FALSE
  )
  index <- 1L
  while (index <= length(arguments)) {
    argument <- arguments[[index]]
    if (argument == "--help") { usage(); quit(status = 0) }
    if (argument %in% c("--prepare-only", "--test-mode", "--acknowledge-exploratory", "--force")) {
      key <- switch(
        argument,
        `--prepare-only` = "prepare_only",
        `--test-mode` = "test_mode",
        `--acknowledge-exploratory` = "acknowledge",
        `--force` = "force"
      )
      parsed[[key]] <- TRUE
      index <- index + 1L
    } else if (argument %in% c("--config", "--output-dir", "--threads", "--bootstrap")) {
      if (index == length(arguments)) stop("Missing value after ", argument, call. = FALSE)
      key <- switch(
        argument,
        `--config` = "config",
        `--output-dir` = "output_dir",
        `--threads` = "threads",
        `--bootstrap` = "bootstrap"
      )
      parsed[[key]] <- arguments[[index + 1L]]
      index <- index + 2L
    } else stop("Unknown argument: ", argument, call. = FALSE)
  }
  parsed
}

write_matrix_csv <- function(matrix, path) {
  output <- data.frame(taxon = rownames(matrix), matrix, check.names = FALSE)
  utils::write.csv(output, path, row.names = FALSE)
}

safe_frequency <- function(results, field, length_out, requested = length(results)) {
  if (!requested) return(rep(NA_real_, length_out))
  if (!length(results)) return(rep(0, length_out))
  rowSums(do.call(cbind, lapply(results, `[[`, field))) / requested
}

safe_count <- function(results, field, length_out) {
  if (!length(results)) return(rep(0L, length_out))
  rowSums(do.call(cbind, lapply(results, `[[`, field)))
}

arguments <- parse_args(commandArgs(trailingOnly = TRUE))
if (!file.exists(arguments$config)) stop("Configuration not found: ", arguments$config, call. = FALSE)
config <- yaml::read_yaml(arguments$config)
validate_stability_config(config)
if (!arguments$prepare_only && !arguments$acknowledge) {
  stop("Full stability execution requires --acknowledge-exploratory.", call. = FALSE)
}
if (arguments$test_mode) {
  config$inputs$abundance_file <- "data/test_beta_diversity/gut_test.csv"
  config$inputs$metadata_file <- "data/test_beta_diversity/sample_info_test.csv"
  config$inputs$taxonomy_file <- "data/test_beta_diversity/GDTB_test.csv"
  config$netmoss$minimum_module_size <- 2L
  config$netmoss$top_taxa <- 2L
  config$taxonomic_analyses$genus$sensitivity_thresholds <- c(0.5, 0.8)
  config$taxonomic_analyses$family$sensitivity_thresholds <- c(0.5, 0.8)
}

threads <- as.integer(arguments$threads %||% config$execution$default_threads)
bootstrap_replicates <- as.integer(arguments$bootstrap %||% if (arguments$test_mode) {
  config$execution$test_bootstrap_replicates
} else config$stability$bootstrap_replicates)
if (!is.finite(threads) || threads < 1L) stop("Threads must be a positive integer.", call. = FALSE)
if (!is.finite(bootstrap_replicates) || bootstrap_replicates < 1L) {
  stop("Bootstrap replicates must be a positive integer.", call. = FALSE)
}
network_method <- if (arguments$test_mode) config$network$smoke_test_method else config$network$method
output_dir <- arguments$output_dir %||% config$outputs$output_dir
manifest_path <- file.path(output_dir, "reproducibility", "analysis_manifest.yml")
assert_output(manifest_path, arguments$force)
directories <- file.path(output_dir, c("tables", "reference", "logs", "reproducibility"))
invisible(vapply(directories, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))

abundance <- read_csv_checked(config$inputs$abundance_file, "Abundance")
metadata <- read_csv_checked(config$inputs$metadata_file, "Metadata")
taxonomy <- read_csv_checked(config$inputs$taxonomy_file, "Taxonomy")
ranks <- c("genus", "family")
prepared <- stats::setNames(lapply(ranks, function(rank) {
  prepare_rank_profiles(abundance, metadata, taxonomy, config, rank)
}), ranks)

if (!identical(prepared$genus$matched_subjects, prepared$family$matched_subjects)) {
  stop("Matched participant definitions differ between taxonomic ranks.", call. = FALSE)
}
common_subjects <- prepared$genus$matched_subjects
expected_primary_count <- length(common_subjects)
cohort_membership <- do.call(rbind, lapply(sort(union(prepared$genus$subjects_a, prepared$genus$subjects_b)), function(subject) {
  data.frame(
    subject_id = subject,
    in_group_a = subject %in% prepared$genus$subjects_a,
    in_group_b = subject %in% prepared$genus$subjects_b,
    in_primary_matched_cohort = subject %in% common_subjects,
    primary_exclusion_reason = if (subject %in% common_subjects) NA_character_ else "not_observed_in_both_groups",
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(cohort_membership, file.path(output_dir, "tables", "cohort_membership.csv"), row.names = FALSE, na = "NA")

scenario_definitions <- list()
for (rank in ranks) {
  settings <- config$taxonomic_analyses[[rank]]
  thresholds <- sort(unique(as.numeric(unlist(settings$sensitivity_thresholds))))
  primary_threshold <- as.numeric(settings$primary_prevalence_threshold)
  for (threshold in thresholds) {
    scenario_definitions[[length(scenario_definitions) + 1L]] <- list(
      rank = rank,
      cohort = "matched_participants",
      threshold = threshold,
      primary = identical(threshold, primary_threshold)
    )
  }
  scenario_definitions[[length(scenario_definitions) + 1L]] <- list(
    rank = rank,
    cohort = "all_available_participants",
    threshold = primary_threshold,
    primary = FALSE
  )
}

scenario_summaries <- list()
filtering_tables <- list()
replicate_tables <- list()
edge_tables <- list()
taxon_tables <- list()
influence_tables <- list()
reference_scores <- list()
scenario_results <- list()

for (scenario_index in seq_along(scenario_definitions)) {
  definition <- scenario_definitions[[scenario_index]]
  rank_data <- prepared[[definition$rank]]
  group_a <- if (definition$cohort == "matched_participants") rank_data$matched_a else rank_data$all_a
  group_b <- if (definition$cohort == "matched_participants") rank_data$matched_b else rank_data$all_b
  filtered <- filter_rank_profiles(
    group_a,
    group_b,
    threshold = definition$threshold,
    detection = as.numeric(config$taxonomic_analyses$detection_threshold)
  )
  scenario_id <- paste(
    definition$rank,
    if (definition$cohort == "matched_participants") "matched" else "all",
    paste0("p", sprintf("%02d", round(100 * definition$threshold))),
    sep = "_"
  )
  filter_table <- filtered$filtering
  filter_table$scenario_id <- scenario_id
  filter_table$taxonomic_rank <- definition$rank
  filter_table$cohort <- definition$cohort
  filter_table$prevalence_threshold <- definition$threshold
  filtering_tables[[scenario_id]] <- filter_table
  retained_taxa <- nrow(filtered$group_a)
  estimable <- retained_taxa >= as.integer(config$netmoss$minimum_module_size)
  base_summary <- data.frame(
    scenario_id = scenario_id,
    taxonomic_rank = definition$rank,
    cohort = definition$cohort,
    prevalence_threshold = definition$threshold,
    primary_scenario = definition$primary,
    participants_group_a = ncol(filtered$group_a),
    participants_group_b = ncol(filtered$group_b),
    matched_participants = length(intersect(colnames(filtered$group_a), colnames(filtered$group_b))),
    retained_taxa = retained_taxa,
    minimum_module_size = as.integer(config$netmoss$minimum_module_size),
    estimable = estimable,
    status = if (estimable) if (arguments$prepare_only) "prepared" else "pending" else "not_estimable_too_few_taxa",
    reference_error = NA_character_,
    bootstrap_requested = bootstrap_replicates,
    bootstrap_successes = NA_integer_,
    bootstrap_success_fraction = NA_real_,
    leave_one_out_requested = NA_integer_,
    leave_one_out_successes = NA_integer_,
    median_bootstrap_network_similarity_a = NA_real_,
    median_bootstrap_network_similarity_b = NA_real_,
    median_bootstrap_module_ari_a = NA_real_,
    median_bootstrap_module_ari_b = NA_real_,
    reference_module_degenerate_a = NA,
    reference_module_degenerate_b = NA,
    operational_module_stability = FALSE,
    operational_network_stability = FALSE,
    stringsAsFactors = FALSE
  )
  if (!estimable || arguments$prepare_only) {
    scenario_summaries[[scenario_id]] <- base_summary
    next
  }

  seed <- as.integer(config$network$random_seed) + scenario_index * 100000L
  reference <- tryCatch(
    run_netmoss_analysis(filtered$group_a, filtered$group_b, network_method, config, seed),
    error = identity
  )
  if (inherits(reference, "error")) {
    base_summary$status <- "reference_failed"
    base_summary$reference_error <- conditionMessage(reference)
    scenario_summaries[[scenario_id]] <- base_summary
    next
  }
  reference$module_degenerate_a <- module_is_degenerate(
    reference$module_a, config$stability$degenerate_module_fraction
  )
  reference$module_degenerate_b <- module_is_degenerate(
    reference$module_b, config$stability$degenerate_module_fraction
  )
  reference_dir <- file.path(output_dir, "reference", scenario_id)
  dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
  write_matrix_csv(reference$network_a, file.path(reference_dir, "group_a_network.csv"))
  write_matrix_csv(reference$network_b, file.path(reference_dir, "group_b_network.csv"))
  write_matrix_csv(close_columns(filtered$group_a), file.path(reference_dir, "group_a_abundance.csv"))
  write_matrix_csv(close_columns(filtered$group_b), file.path(reference_dir, "group_b_abundance.csv"))

  reference_rank <- stats::setNames(rep(Inf, length(reference$scores)), names(reference$scores))
  reference_rank[reference$scored] <- rank(-reference$scores[reference$scored], ties.method = "min")
  rank_output <- rep(NA_integer_, length(reference_rank))
  rank_output[is.finite(reference_rank)] <- as.integer(reference_rank[is.finite(reference_rank)])
  reference_score_table <- data.frame(
    scenario_id = scenario_id,
    taxonomic_rank = definition$rank,
    cohort = definition$cohort,
    prevalence_threshold = definition$threshold,
    taxon = names(reference$scores),
    netmoss_score = as.numeric(reference$scores),
    netmoss_scored = as.logical(reference$scored),
    reference_rank = rank_output,
    reference_top_taxon = reference$scored & reference_rank <= as.integer(config$netmoss$top_taxa),
    module_group_a = as.integer(reference$module_a),
    module_group_b = as.integer(reference$module_b),
    stringsAsFactors = FALSE
  )
  reference_scores[[scenario_id]] <- reference_score_table

  specifications <- make_resample_specifications(
    filtered$group_a,
    filtered$group_b,
    bootstrap_replicates,
    seed + 1000L
  )
  worker <- function(index) analyze_resample(
    specifications[[index]], filtered$group_a, filtered$group_b,
    reference, config, network_method
  )
  if (.Platform$OS.type == "unix" && threads > 1L) {
    results <- parallel::mclapply(
      seq_along(specifications), worker,
      mc.cores = min(threads, length(specifications)),
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    if (threads > 1L) warning("Parallel workers unavailable; using one worker.", call. = FALSE)
    results <- lapply(seq_along(specifications), worker)
  }
  metric_table <- do.call(rbind, lapply(results, function(result) {
    data.frame(
      scenario_id = scenario_id,
      resample_type = result$type,
      resample_id = result$id,
      omitted_participant = result$omitted,
      success = result$success,
      error = result$error,
      network_similarity_a = result$network_similarity_a %||% NA_real_,
      network_similarity_b = result$network_similarity_b %||% NA_real_,
      score_rank_correlation = result$score_rank_correlation %||% NA_real_,
      module_ari_a = result$module_ari_a %||% NA_real_,
      module_ari_b = result$module_ari_b %||% NA_real_,
      degenerate_module_a = result$degenerate_module_a %||% NA,
      degenerate_module_b = result$degenerate_module_b %||% NA,
      stringsAsFactors = FALSE
    )
  }))
  replicate_tables[[scenario_id]] <- metric_table
  bootstrap_results <- results[vapply(results, function(x) x$success && x$type == "bootstrap", logical(1))]
  loo_results <- results[vapply(results, function(x) x$success && x$type == "leave_one_out", logical(1))]
  requested_bootstrap <- sum(vapply(specifications, function(x) x$type == "bootstrap", logical(1)))
  requested_loo <- sum(vapply(specifications, function(x) x$type == "leave_one_out", logical(1)))
  bootstrap_success_fraction <- length(bootstrap_results) / requested_bootstrap
  median_a <- if (length(bootstrap_results)) stats::median(vapply(bootstrap_results, `[[`, numeric(1), "network_similarity_a"), na.rm = TRUE) else NA_real_
  median_b <- if (length(bootstrap_results)) stats::median(vapply(bootstrap_results, `[[`, numeric(1), "network_similarity_b"), na.rm = TRUE) else NA_real_
  finite_median <- function(results, field) {
    if (!length(results)) return(NA_real_)
    values <- vapply(results, function(x) x[[field]], numeric(1))
    if (any(is.finite(values))) stats::median(values[is.finite(values)]) else NA_real_
  }
  module_median_a <- finite_median(bootstrap_results, "module_ari_a")
  module_median_b <- finite_median(bootstrap_results, "module_ari_b")
  operational_network_stability <- operational_network_pass(
    bootstrap_success_fraction, median_a, median_b, config
  )
  operational_module_stability <- !reference$module_degenerate_a &&
    !reference$module_degenerate_b &&
    is.finite(module_median_a) && module_median_a >= as.numeric(config$stability$module_ari_threshold) &&
    is.finite(module_median_b) && module_median_b >= as.numeric(config$stability$module_ari_threshold)

  base_summary$status <- "completed_internal_stability"
  base_summary$bootstrap_successes <- length(bootstrap_results)
  base_summary$bootstrap_success_fraction <- bootstrap_success_fraction
  base_summary$leave_one_out_requested <- requested_loo
  base_summary$leave_one_out_successes <- length(loo_results)
  base_summary$median_bootstrap_network_similarity_a <- median_a
  base_summary$median_bootstrap_network_similarity_b <- median_b
  base_summary$median_bootstrap_module_ari_a <- module_median_a
  base_summary$median_bootstrap_module_ari_b <- module_median_b
  base_summary$reference_module_degenerate_a <- reference$module_degenerate_a
  base_summary$reference_module_degenerate_b <- reference$module_degenerate_b
  base_summary$operational_module_stability <- operational_module_stability
  base_summary$operational_network_stability <- operational_network_stability
  scenario_summaries[[scenario_id]] <- base_summary

  taxa <- names(reference$scores)
  top_bootstrap_frequency <- safe_frequency(
    bootstrap_results, "top_taxa", length(taxa), requested_bootstrap
  )
  top_loo_count <- safe_count(loo_results, "top_taxa", length(taxa))
  loo_required <- ceiling(as.numeric(config$stability$selection_frequency_threshold) * requested_loo)
  taxon_table <- reference_score_table
  taxon_table$bootstrap_top10_frequency <- top_bootstrap_frequency
  taxon_table$leave_one_out_top10_count <- top_loo_count
  taxon_table$leave_one_out_required <- loo_required
  taxon_table$operationally_stable_taxon <- operational_selection(
    taxon_table$reference_top_taxon,
    top_bootstrap_frequency,
    top_loo_count,
    requested_loo,
    operational_network_stability,
    config
  )
  taxon_tables[[scenario_id]] <- taxon_table

  edge_a <- edge_index(taxa, config$design$group_a_label, reference$network_a)
  edge_b <- edge_index(taxa, config$design$group_b_label, reference$network_b)
  edge_a$bootstrap_selection_frequency <- safe_frequency(
    bootstrap_results, "edge_same_a", nrow(edge_a), requested_bootstrap
  )
  edge_b$bootstrap_selection_frequency <- safe_frequency(
    bootstrap_results, "edge_same_b", nrow(edge_b), requested_bootstrap
  )
  edge_a$leave_one_out_selection_count <- safe_count(loo_results, "edge_same_a", nrow(edge_a))
  edge_b$leave_one_out_selection_count <- safe_count(loo_results, "edge_same_b", nrow(edge_b))
  edges <- rbind(edge_a, edge_b)
  edges$scenario_id <- scenario_id
  edges$reference_strong <- abs(edges$reference_correlation) >= as.numeric(config$network$edge_magnitude_threshold)
  edges$leave_one_out_required <- loo_required
  edges$operationally_stable_edge <- operational_selection(
    edges$reference_strong,
    edges$bootstrap_selection_frequency,
    edges$leave_one_out_selection_count,
    requested_loo,
    operational_network_stability,
    config
  )
  edge_tables[[scenario_id]] <- edges

  influence <- metric_table[metric_table$resample_type == "leave_one_out", , drop = FALSE]
  influence$mean_network_similarity <- rowMeans(
    influence[c("network_similarity_a", "network_similarity_b")], na.rm = TRUE
  )
  influence$influence_rank <- rank(influence$mean_network_similarity, ties.method = "min", na.last = "keep")
  influence_tables[[scenario_id]] <- influence
  scenario_results[[scenario_id]] <- list(reference = reference, definition = definition)
}

scenario_summary <- do.call(rbind, scenario_summaries)
filtering <- do.call(rbind, filtering_tables)
utils::write.csv(scenario_summary, file.path(output_dir, "tables", "scenario_summary.csv"), row.names = FALSE, na = "NA")
utils::write.csv(filtering, file.path(output_dir, "tables", "taxon_filtering.csv"), row.names = FALSE, na = "NA")

status <- if (arguments$prepare_only) "PREPARED_ONLY" else "COMPLETED_EXPLORATORY_STABILITY"
if (!arguments$prepare_only) {
  combined <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
  replicate_summary <- combined(replicate_tables)
  edge_stability <- combined(edge_tables)
  taxon_stability <- combined(taxon_tables)
  participant_influence <- combined(influence_tables)
  reference_score_table <- combined(reference_scores)
  utils::write.csv(replicate_summary, file.path(output_dir, "tables", "replicate_summary.csv"), row.names = FALSE, na = "NA")
  utils::write.csv(
    replicate_summary[c(
      "scenario_id", "resample_type", "resample_id", "omitted_participant",
      "success", "module_ari_a", "module_ari_b", "degenerate_module_a",
      "degenerate_module_b"
    )],
    file.path(output_dir, "tables", "module_stability.csv"),
    row.names = FALSE,
    na = "NA"
  )
  utils::write.csv(
    replicate_summary[!replicate_summary$success, , drop = FALSE],
    file.path(output_dir, "tables", "failed_resamples.csv"),
    row.names = FALSE,
    na = "NA"
  )
  utils::write.csv(edge_stability, file.path(output_dir, "tables", "edge_stability.csv"), row.names = FALSE, na = "NA")
  utils::write.csv(taxon_stability, file.path(output_dir, "tables", "taxon_rank_stability.csv"), row.names = FALSE, na = "NA")
  utils::write.csv(participant_influence, file.path(output_dir, "tables", "participant_influence.csv"), row.names = FALSE, na = "NA")
  utils::write.csv(reference_score_table, file.path(output_dir, "tables", "reference_netmoss_scores.csv"), row.names = FALSE, na = "NA")
  utils::write.csv(edge_stability[edge_stability$operationally_stable_edge %in% TRUE, , drop = FALSE], file.path(output_dir, "tables", "operationally_stable_edges.csv"), row.names = FALSE, na = "NA")
  utils::write.csv(taxon_stability[taxon_stability$operationally_stable_taxon %in% TRUE, , drop = FALSE], file.path(output_dir, "tables", "operationally_stable_taxa.csv"), row.names = FALSE, na = "NA")

  genus_primary <- taxon_stability[
    taxon_stability$scenario_id == "genus_matched_p80", , drop = FALSE
  ]
  family_primary <- taxon_stability[
    taxon_stability$scenario_id == "family_matched_p50", , drop = FALSE
  ]
  if (nrow(genus_primary) && nrow(family_primary)) {
    raw_mapping <- data.frame(
      genus = as.character(taxonomy$genus),
      family = as.character(taxonomy$family),
      stringsAsFactors = FALSE
    )
    raw_mapping <- raw_mapping[
      !is.na(raw_mapping$genus) & nzchar(raw_mapping$genus) &
        !is.na(raw_mapping$family) & nzchar(raw_mapping$family),
      ,
      drop = FALSE
    ]
    family_by_genus <- split(raw_mapping$family, raw_mapping$genus)
    mapping <- do.call(rbind, lapply(names(family_by_genus), function(genus) {
      families <- sort(unique(family_by_genus[[genus]]))
      data.frame(
        genus = genus,
        family = if (length(families) == 1L) families else NA_character_,
        taxonomy_mapping_ambiguous = length(families) != 1L,
        stringsAsFactors = FALSE
      )
    }))
    cross_rank <- merge(genus_primary, mapping, by.x = "taxon", by.y = "genus", all.x = TRUE)
    family_fields <- family_primary[c("taxon", "reference_rank", "bootstrap_top10_frequency", "leave_one_out_top10_count", "operationally_stable_taxon")]
    names(family_fields) <- c("family", "family_reference_rank", "family_bootstrap_top10_frequency", "family_leave_one_out_top10_count", "family_operationally_stable")
    cross_rank <- merge(cross_rank, family_fields, by = "family", all.x = TRUE)
    cross_rank$cross_rank_status <- ifelse(
      cross_rank$operationally_stable_taxon & cross_rank$family_operationally_stable %in% TRUE,
      "stable_at_genus_and_parent_family",
      ifelse(cross_rank$operationally_stable_taxon, "stable_at_genus_only", ifelse(
        cross_rank$family_operationally_stable %in% TRUE, "stable_at_parent_family_only", "not_stable_at_either_level"
      ))
    )
  } else {
    cross_rank <- data.frame(
      taxon = character(0), family = character(0), cross_rank_status = character(0),
      stringsAsFactors = FALSE
    )
  }
  utils::write.csv(cross_rank, file.path(output_dir, "tables", "cross_rank_concordance.csv"), row.names = FALSE, na = "NA")
}

manifest <- list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  status = status,
  analysis_status = config$project$analysis_status,
  network_method = network_method,
  bootstrap_replicates = bootstrap_replicates,
  threads = threads,
  matched_participants = common_subjects,
  matched_participant_count = expected_primary_count,
  group_a_all_participants = prepared$genus$subjects_a,
  group_b_all_participants = prepared$genus$subjects_b,
  group_a_raw_samples = prepared$genus$raw_sample_count_a,
  group_b_raw_samples = prepared$genus$raw_sample_count_b,
  scenario_count = nrow(scenario_summary),
  estimable_scenarios = sum(scenario_summary$estimable),
  operationally_stable_scenarios = sum(scenario_summary$operational_network_stability, na.rm = TRUE),
  input_md5 = as.list(tools::md5sum(unlist(config$inputs))),
  netmoss2_commit_expected = config$netmoss$package_commit,
  netmoss2_remote_sha_observed = if (requireNamespace("NetMoss2", quietly = TRUE)) utils::packageDescription("NetMoss2")$RemoteSha %||% "not_recorded" else "not_installed",
  git_commit = tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(error) "unknown"),
  scientific_boundary = paste(
    "Internal reproducibility assessment only; not external validation, not causal altitude inference,",
    "and not permission to reinterpret unstable networks as biomarkers."
  )
)
yaml::write_yaml(manifest, manifest_path)

old_timezone <- Sys.getenv("TZ", unset = "")
Sys.setenv(TZ = "UTC")
suppressWarnings(utils::capture.output(utils::sessionInfo(), file = file.path(output_dir, "reproducibility", "sessionInfo.txt")))
if (nzchar(old_timezone)) Sys.setenv(TZ = old_timezone) else Sys.unsetenv("TZ")
writeLines(c(
  "NetMoss altitude network-stability workflow",
  paste0("Status: ", status),
  paste0("Matched primary participants: ", expected_primary_count),
  paste0("Scenarios: ", nrow(scenario_summary)),
  paste0("Estimable scenarios: ", sum(scenario_summary$estimable)),
  paste0("Operationally stable scenarios: ", sum(scenario_summary$operational_network_stability, na.rm = TRUE)),
  "Interpretation: internal reproducibility only; failure is evidence of network fragility."
), file.path(output_dir, "logs", "analysis_summary.txt"))

cat("NetMoss altitude stability status: ", status, "\n", sep = "")
cat("Output directory: ", normalizePath(output_dir, mustWork = TRUE), "\n", sep = "")
