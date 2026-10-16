#!/usr/bin/env Rscript

# Purpose: Run a complete, matched-participant genus-abundance comparison for
# every ordered pair of observed expedition timepoints.
#
# Inputs: gut MAG abundance, sample metadata, exact MAG-to-taxonomy mapping, and
# config/netmoss_matched_timepoint_companion.yml.
# Outputs: pair support, full genus results, participant-level changes,
# candidates passing predeclared safeguards, logs, and reproducibility records.
# Dependencies: base R, yaml; parallel execution uses base parallel on Unix.
# Invocation: run from the repository root; see
# docs/netmoss-methods.md and the dedicated OpenPBS job.
# Relationship: this is an abundance-based companion to NetMoss, not a NetMoss
# network analysis. It preserves temporal pairing without treating repeated
# samples as independent network observations.

suppressPackageStartupMessages(library(yaml))
source("scripts/netmoss/matched_timepoint_helpers.R")

usage <- function() {
  cat(
    "Usage: Rscript scripts/netmoss/run_matched_timepoint_companion.R [options]\n",
    "  --config PATH       Default: config/netmoss_matched_timepoint_companion.yml\n",
    "  --output-dir PATH   Override configured output directory\n",
    "  --mode MODE         standard or advanced; default from configuration\n",
    "  --threads N         Parallel timepoint-pair workers\n",
    "  --test-mode         Use repository test data and small bootstrap count\n",
    "  --force             Replace existing named outputs\n",
    "  --help\n",
    sep = ""
  )
}

parse_args <- function(arguments) {
  parsed <- list(
    config = "config/netmoss_matched_timepoint_companion.yml",
    output_dir = NULL,
    mode = NULL,
    threads = NULL,
    test_mode = FALSE,
    force = FALSE
  )
  index <- 1L
  while (index <= length(arguments)) {
    argument <- arguments[[index]]
    if (argument == "--help") {
      usage()
      quit(status = 0)
    } else if (argument %in% c("--test-mode", "--force")) {
      parsed[[if (argument == "--test-mode") "test_mode" else "force"]] <- TRUE
      index <- index + 1L
    } else if (argument %in% c("--config", "--output-dir", "--mode", "--threads")) {
      if (index == length(arguments)) stop("Missing value after ", argument, call. = FALSE)
      key <- switch(
        argument,
        `--config` = "config",
        `--output-dir` = "output_dir",
        `--mode` = "mode",
        `--threads` = "threads"
      )
      parsed[[key]] <- arguments[[index + 1L]]
      index <- index + 2L
    } else {
      stop("Unknown argument: ", argument, call. = FALSE)
    }
  }
  parsed
}

arguments <- parse_args(commandArgs(trailingOnly = TRUE))
if (!file.exists(arguments$config)) {
  stop("Configuration not found: ", arguments$config, call. = FALSE)
}
config <- yaml::read_yaml(arguments$config)

if (arguments$test_mode) {
  config$inputs$abundance_file <- "data/test_beta_diversity/gut_test.csv"
  config$inputs$metadata_file <- "data/test_beta_diversity/sample_info_test.csv"
  config$inputs$taxonomy_file <- "data/test_beta_diversity/GDTB_test.csv"
  config$design$timepoint_order <- c("T2", "T3", "T7", "T10")
  config$design$minimum_pairs_for_inference <- 2L
  config$design$minimum_detected_participants <- 1L
}

mode <- arguments$mode %||% config$execution$default_mode %||% "standard"
if (!mode %in% c("standard", "advanced")) {
  stop("Mode must be standard or advanced.", call. = FALSE)
}
bootstrap_replicates <- if (arguments$test_mode) {
  100L
} else if (mode == "advanced") {
  as.integer(config$inference$advanced_bootstrap_replicates)
} else {
  as.integer(config$inference$standard_bootstrap_replicates)
}
threads <- as.integer(arguments$threads %||% config$execution$default_threads %||% 1L)
if (!is.finite(threads) || threads < 1L) stop("Threads must be a positive integer.", call. = FALSE)

output_dir <- arguments$output_dir %||% config$outputs$output_dir
manifest_path <- file.path(output_dir, "reproducibility", "analysis_manifest.yml")
assert_output(manifest_path, arguments$force)
directories <- file.path(output_dir, c("tables", "logs", "reproducibility"))
invisible(vapply(directories, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))

abundance <- read_csv_checked(config$inputs$abundance_file, "Abundance")
metadata <- read_csv_checked(config$inputs$metadata_file, "Metadata")
taxonomy <- read_csv_checked(config$inputs$taxonomy_file, "Taxonomy")
excluded_timepoints <- unlist(config$design$excluded_timepoints %||% character(0))
excluded_sample_ids <- character(0)
excluded_abundance_samples <- 0L
if (length(excluded_timepoints)) {
  assert_columns(
    metadata,
    c(config$columns$sample_id, config$columns$timepoint),
    "Metadata"
  )
  excluded_rows <- as.character(metadata[[config$columns$timepoint]]) %in%
    excluded_timepoints
  excluded_sample_ids <- as.character(metadata[[config$columns$sample_id]][excluded_rows])
  excluded_abundance_samples <- sum(
    as.character(abundance[[config$columns$sample_id]]) %in% excluded_sample_ids
  )
  metadata <- metadata[!excluded_rows, , drop = FALSE]
  abundance <- abundance[
    !as.character(abundance[[config$columns$sample_id]]) %in% excluded_sample_ids,
    ,
    drop = FALSE
  ]
}
prepared <- prepare_genus_matrices(abundance, metadata, taxonomy, config)

configured_timepoints <- unlist(config$design$timepoint_order)
observed_timepoints <- unique(prepared$metadata$timepoint_internal)
missing_timepoints <- setdiff(configured_timepoints, observed_timepoints)
if (length(missing_timepoints)) {
  stop("Configured timepoints are absent from matched abundance data: ",
       paste(missing_timepoints, collapse = ", "), call. = FALSE)
}
timepoint_pairs <- make_timepoint_pairs(configured_timepoints)
expected_pair_count <- choose(length(configured_timepoints), 2L)
if (nrow(timepoint_pairs) != expected_pair_count) {
  stop("Timepoint-pair enumeration failed.", call. = FALSE)
}

worker <- function(pair_index) {
  set.seed(as.integer(config$inference$random_seed) + pair_index)
  pair <- timepoint_pairs[pair_index, , drop = FALSE]
  analyze_timepoint_pair(pair, prepared, config, bootstrap_replicates)
}

if (.Platform$OS.type == "unix" && threads > 1L) {
  analyses <- parallel::mclapply(
    seq_len(nrow(timepoint_pairs)),
    worker,
    mc.cores = min(threads, nrow(timepoint_pairs)),
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  if (threads > 1L) warning("Parallel workers are unavailable; using one worker.", call. = FALSE)
  analyses <- lapply(seq_len(nrow(timepoint_pairs)), worker)
}

pair_summary <- do.call(rbind, lapply(analyses, `[[`, "pair_summary"))
genus_results <- do.call(rbind, lapply(analyses, `[[`, "genus_results"))
participant_results <- do.call(rbind, lapply(analyses, `[[`, "participant_results"))

genus_results$global_bh_q_value <- NA_real_
tested <- genus_results$testable & is.finite(genus_results$exact_p_value)
genus_results$global_bh_q_value[tested] <- stats::p.adjust(
  genus_results$exact_p_value[tested],
  method = "BH"
)
ci_excludes_zero <- genus_results$bootstrap_ci_lower > 0 |
  genus_results$bootstrap_ci_upper < 0
genus_results$predeclared_candidate <- tested &
  genus_results$global_bh_q_value <= as.numeric(config$inference$fdr_threshold) &
  ci_excludes_zero &
  genus_results$direction_concordance >= as.numeric(
    config$inference$direction_concordance_threshold
  ) &
  genus_results$leave_one_out_sign_stability >= as.numeric(
    config$inference$leave_one_out_sign_stability_threshold
  )
genus_results$interpretation_status <- ifelse(
  !genus_results$inferential_pair,
  "descriptive_low_pair_support",
  ifelse(
    !genus_results$testable,
    "descriptive_low_detection_support",
    ifelse(genus_results$predeclared_candidate, "candidate", "tested_not_candidate")
  )
)

candidate_results <- genus_results[genus_results$predeclared_candidate, , drop = FALSE]
candidate_participants <- participant_results[
  paste(participant_results$pair_id, participant_results$genus) %in%
    paste(candidate_results$pair_id, candidate_results$genus),
  ,
  drop = FALSE
]

timepoint_support <- as.data.frame(table(
  factor(prepared$metadata$timepoint_internal, levels = configured_timepoints)
))
names(timepoint_support) <- c("timepoint", "observed_participants")

tables <- list(
  timepoint_support = timepoint_support,
  timepoint_pair_support = pair_summary,
  all_genus_pair_results = genus_results,
  all_participant_genus_changes = participant_results,
  predeclared_candidates = candidate_results,
  candidate_participant_changes = candidate_participants
)
for (name in names(tables)) {
  utils::write.csv(
    tables[[name]],
    file.path(output_dir, "tables", paste0(name, ".csv")),
    row.names = FALSE,
    na = "NA"
  )
}

manifest <- list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  status = "COMPLETED_EXPLORATORY_COMPANION",
  analysis_status = config$project$analysis_status,
  mode = mode,
  threads = threads,
  bootstrap_replicates = bootstrap_replicates,
  timepoints = configured_timepoints,
  excluded_timepoints = excluded_timepoints,
  excluded_metadata_records = length(excluded_sample_ids),
  excluded_abundance_samples = excluded_abundance_samples,
  enumerated_timepoint_pairs = nrow(timepoint_pairs),
  sample_count = ncol(prepared$raw),
  mag_count = prepared$feature_count,
  genus_count = prepared$genus_count,
  zero_replacement = prepared$zero_replacement,
  tested_genus_pair_hypotheses = sum(tested),
  predeclared_candidates = nrow(candidate_results),
  netmoss_feasible_timepoint_pairs = sum(pair_summary$netmoss_network_feasible),
  input_md5 = as.list(tools::md5sum(unlist(config$inputs))),
  git_commit = tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(error) "unknown"),
  scientific_boundary = paste(
    "Paired genus-abundance companion analysis; not NetMoss, not participant-level",
    "network inference, and not validated biomarker discovery."
  )
)
yaml::write_yaml(manifest, manifest_path)

old_timezone <- Sys.getenv("TZ", unset = "")
Sys.setenv(TZ = "UTC")
suppressWarnings(utils::capture.output(
  utils::sessionInfo(),
  file = file.path(output_dir, "reproducibility", "sessionInfo.txt")
))
if (nzchar(old_timezone)) Sys.setenv(TZ = old_timezone) else Sys.unsetenv("TZ")

writeLines(c(
  "Matched-timepoint NetMoss companion analysis",
  paste0("Status: COMPLETED_EXPLORATORY_COMPANION"),
  paste0("Mode: ", mode),
  paste0("Timepoint pairs enumerated: ", nrow(timepoint_pairs)),
  paste0("Matched participants per pair: ", min(pair_summary$matched_participants),
         "-", max(pair_summary$matched_participants)),
  paste0("Genus-pair hypotheses tested: ", sum(tested)),
  paste0("Predeclared candidates: ", nrow(candidate_results)),
  paste0("Timepoint pairs feasible for NetMoss networks: ",
         sum(pair_summary$netmoss_network_feasible)),
  "Interpretation: abundance-based paired exploration only; this is not NetMoss."
), file.path(output_dir, "logs", "analysis_summary.txt"))

cat("Matched-timepoint companion analysis completed.\n")
cat("Output directory: ", normalizePath(output_dir, mustWork = TRUE), "\n", sep = "")
