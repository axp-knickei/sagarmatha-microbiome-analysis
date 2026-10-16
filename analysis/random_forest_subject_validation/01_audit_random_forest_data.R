#!/usr/bin/env Rscript
# Purpose: Audit the long-format Random Forest input, existing split, class
# balance, subject coverage, and altitude/time confounding before validation.
# Inputs: Data and existing-run paths from random_forest_subject_validation.yml.
# Outputs: Audit CSVs and reports/analysis_audit.md in RF_VALIDATION_RUN_DIR.
# Dependencies: R/common.R and packages declared in environment.yml.
# Usage: Normally invoked by run_random_forest_subject_validation.R; set
# RF_VALIDATION_RUN_DIR when running the stage directly.
# Workflow: Stage 01; establishes data-integrity evidence for later stages.

source("analysis/random_forest_subject_validation/R/common.R")
args <- commandArgs(trailingOnly = TRUE)
cfg <- read_validation_config(if (length(args)) args[[1]] else "config/random_forest_subject_validation.yml")
run_dir <- Sys.getenv("RF_VALIDATION_RUN_DIR", unset = "")
if (!nzchar(run_dir)) run_dir <- make_run_directory(cfg)
out <- file.path(run_dir, "audit"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "reports"), recursive = TRUE, showWarnings = FALSE)
d <- read_long_data(cfg); m <- sample_metadata(d, cfg)
sid <- cfg$data$sample_id; sub <- cfg$data$subject_id; fid <- cfg$data$feature_id; y <- cfg$data$abundance
alt <- cfg$data$outcome; tp <- cfg$data$timepoint; date <- cfg$data$sampling_date
days <- cfg$data$days_from_baseline; global_days <- cfg$data$days_from_global_baseline
dataset_summary <- tibble(metric = c("samples", "subjects", "altitude_classes", "timepoints", "MAG_predictors", "long_rows", "missing_values", "zero_abundance_proportion"),
 value = as.character(c(n_distinct(m[[sid]]), n_distinct(m[[sub]]), n_distinct(m[[alt]]), n_distinct(m[[tp]]), n_distinct(d[[fid]]), nrow(d), sum(is.na(d)), mean(d[[y]] == 0, na.rm = TRUE))))
write.csv(dataset_summary, file.path(out, "dataset_summary.csv"), row.names = FALSE)
sample_feature_counts <- d |> count(.data[[sid]], name = "n_feature_rows")
data_integrity <- tibble(
  check = c(
    "missing_sample_ids", "missing_subject_ids", "missing_outcomes",
    "duplicate_sample_feature_pairs", "non_numeric_abundance", "missing_abundance",
    "negative_abundance", "infinite_abundance", "minimum_features_per_sample",
    "maximum_features_per_sample", "expected_features_per_sample",
    "all_samples_have_complete_feature_grid"
  ),
  value = as.character(c(
    sum(is.na(d[[sid]]) | !nzchar(d[[sid]])),
    sum(is.na(d[[sub]]) | !nzchar(as.character(d[[sub]]))),
    sum(is.na(d[[alt]])), sum(duplicated(d[c(sid, fid)])),
    !is.numeric(d[[y]]), sum(is.na(d[[y]])),
    sum(d[[y]] < 0, na.rm = TRUE), sum(is.infinite(d[[y]])),
    min(sample_feature_counts$n_feature_rows), max(sample_feature_counts$n_feature_rows),
    n_distinct(d[[fid]]),
    all(sample_feature_counts$n_feature_rows == n_distinct(d[[fid]]))
  ))
)
write.csv(data_integrity, file.path(out, "data_integrity_checks.csv"), row.names = FALSE)
write.csv(sample_feature_counts, file.path(out, "sample_feature_counts.csv"), row.names = FALSE)
samples_by_subject <- m |> count(.data[[sub]], name = "n_samples")
samples_by_altitude <- m |> count(.data[[alt]], name = "n_samples") |> mutate(proportion = n_samples / sum(n_samples))
subjects_by_altitude <- m |> group_by(.data[[alt]]) |> summarise(n_subjects = n_distinct(.data[[sub]]), subject_ids = paste(sort(unique(.data[[sub]])), collapse = ";"), .groups = "drop")
coverage <- m |> count(.data[[sub]], .data[[alt]], .data[[tp]], name = "n_samples")
feature_prev <- d |> group_by(.data[[fid]]) |> summarise(sample_prevalence = mean(.data[[y]] > 0), subject_prevalence = n_distinct(.data[[sub]][.data[[y]] > 0]), zero_proportion = mean(.data[[y]] == 0), min = min(.data[[y]]), median = median(.data[[y]]), mean = mean(.data[[y]]), max = max(.data[[y]]), .groups = "drop")
write.csv(samples_by_subject, file.path(out, "samples_by_subject.csv"), row.names = FALSE)
write.csv(samples_by_altitude, file.path(out, "samples_by_altitude.csv"), row.names = FALSE)
write.csv(subjects_by_altitude, file.path(out, "subjects_by_altitude.csv"), row.names = FALSE)
write.csv(coverage, file.path(out, "subject_altitude_coverage.csv"), row.names = FALSE)
write.csv(feature_prev, file.path(out, "feature_prevalence_summary.csv"), row.names = FALSE)
split_path <- file.path(cfg$data$existing_run, "sample_split.csv")
if (!file.exists(split_path)) stop("Configured legacy split file is unavailable: ", split_path, call. = FALSE)
split <- read.csv(split_path, stringsAsFactors = FALSE) |> rename(!!sid := sample_id, !!sub := subject_id, !!alt := response)
split_audit <- split |> group_by(split, .data[[alt]]) |> summarise(n_samples = n(), n_subjects = n_distinct(.data[[sub]]), subject_ids = paste(sort(unique(.data[[sub]])), collapse = ";"), .groups = "drop") |> group_by(split) |> mutate(proportion = n_samples / sum(n_samples)) |> ungroup()
write.csv(split_audit, file.path(out, "existing_split_audit.csv"), row.names = FALSE)
overlap <- intersect(split[[sub]][split$split == "train"], split[[sub]][split$split == "test"])
altitude_by_time <- all((m |> group_by(.data[[tp]]) |> summarise(n = n_distinct(.data[[alt]]), .groups = "drop"))$n == 1)
altitude_by_date <- all((m |> group_by(.data[[date]]) |> summarise(n = n_distinct(.data[[alt]]), .groups = "drop"))$n == 1)
conf <- tibble(
  test = c(
    "altitude_determined_by_timepoint", "altitude_determined_by_date",
    "subject_overlap_existing_split", "pearson_altitude_vs_global_days",
    "pearson_altitude_vs_subject_days"
  ),
  value = c(
    altitude_by_time, altitude_by_date, length(overlap) > 0,
    cor(m[[alt]], m[[global_days]]), cor(m[[alt]], m[[days]])
  )
)
write.csv(conf, file.path(out, "altitude_time_confounding_summary.csv"), row.names = FALSE)
sample_sums <- tapply(d[[y]], d[[sid]], sum)
writeLines(c(
  "# Analysis audit", "",
  paste0("Input: `", attr(d, "input_path"), "`"), "",
  paste0("The data contain ", nrow(m), " samples from ", n_distinct(m[[sub]]), " subjects and ", n_distinct(d[[fid]]), " MAGs."),
  paste0("Zero abundance proportion: ", sprintf("%.4f", mean(d[[y]] == 0)), "."),
  paste0("Every sample has exactly ", n_distinct(d[[fid]]), " explicit MAG rows; absent sample-MAG combinations are rejected rather than zero-filled."),
  paste0("Sample abundance sums range from ", sprintf("%.4f", min(sample_sums)), " to ", sprintf("%.4f", max(sample_sums)),
    ", consistent with relative-abundance-like rather than raw-count data; the upstream normalization is not recorded in this RDS."),
  paste0("Existing split subject overlap: ", ifelse(length(overlap), paste(overlap, collapse = ", "), "none"), "."),
  "",
  "Altitude is perfectly determined by timepoint and date. Prediction cannot separate altitude from expedition chronology in this design.",
  "No expedition-phase or taxonomy field is available; those analyses are not defined."
), file.path(run_dir, "reports", "analysis_audit.md"))
cat(run_dir, "\n")
