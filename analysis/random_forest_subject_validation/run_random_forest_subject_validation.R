#!/usr/bin/env Rscript
# Purpose: Orchestrate the complete subject-aware Random Forest validation and
# record stage completion in a run manifest.
# Inputs: config/random_forest_subject_validation.yml, its configured data, and
# an optional mode (smoke, development, or production) plus run name.
# Outputs: A named run tree under results/random_forest_subject_validation/ with
# audits, predictions, comparisons, features, figures, reports, and manifest.
# Dependencies: R/common.R, stages 01-12, and environment.yml packages.
# Usage: Rscript analysis/random_forest_subject_validation/run_random_forest_subject_validation.R smoke
# Workflow: Canonical local orchestrator; the HPC wrapper calls this entrypoint.

source("analysis/random_forest_subject_validation/R/common.R")
args <- commandArgs(trailingOnly = TRUE)
requested_mode <- if (length(args)) args[[1]] else ""
if (nzchar(requested_mode) && !requested_mode %in% c("smoke", "development", "production")) {
  stop("First argument must be smoke, development, or production.", call. = FALSE)
}
if (nzchar(requested_mode)) Sys.setenv(RF_VALIDATION_MODE = requested_mode)
cfg <- read_validation_config()
run_dir <- make_run_directory(cfg, if (length(args) > 1) args[[2]] else NULL)
run_dir <- normalizePath(run_dir)
env <- c(
  paste0("RF_VALIDATION_RUN_DIR=", run_dir),
  paste0("RF_VALIDATION_MODE=", cfg$analysis$mode)
)
scripts <- c(
  "01_audit_random_forest_data.R",
  "03_run_baseline_models.R",
  "04_run_loso_random_forest.R",
  "05_run_repeated_grouped_cv.R",
  "06_compare_class_balancing.R",
  "07_compare_feature_strategies.R",
  "11_run_tree_stability.R",
  "08_compare_outcome_formulations.R",
  "09_summarize_subject_aware_results.R",
  "10_create_validation_figures.R",
  "12_validate_outputs.R"
)
input_path <- resolve_input(cfg)
resolved_cfg <- cfg
resolved_cfg$.config_path <- NULL
yaml::write_yaml(resolved_cfg, file.path(run_dir, "reports", "resolved_config.yml"))
git_commit <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[[1]], error = function(e) NA_character_)
git_status <- tryCatch(system2("git", c("status", "--short"), stdout = TRUE, stderr = FALSE), error = function(e) "unavailable")
packages <- c("dplyr", "tidyr", "yaml", "randomForest", "ggplot2", "testthat")
package_versions <- setNames(vapply(packages, function(package) {
  if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else NA_character_
}, character(1)), packages)
manifest_path <- file.path(run_dir, "reports", "run_manifest.yml")
started <- Sys.time()
manifest <- list(
  run_id = basename(run_dir), mode = cfg$analysis$mode, status = "running",
  started = format(started, "%Y-%m-%d %H:%M:%S %Z"), completed = NULL,
  seed = cfg$analysis$seed, input = input_path,
  input_md5 = unname(tools::md5sum(input_path)),
  config = cfg$.config_path, config_md5 = unname(tools::md5sum(cfg$.config_path)),
  git_commit = git_commit, git_dirty = length(git_status) > 0L,
  package_versions = as.list(package_versions),
  completed_stages = character(), failed_stage = NULL,
  caution = paste(
    "All results are exploratory prediction among altitude/time-confounded expedition categories.",
    "Production mode is not independent confirmation and fixed hyperparameters are not nested-tuned."
  )
)
write_manifest <- function(x) {
  temporary <- paste0(manifest_path, ".tmp")
  yaml::write_yaml(x, temporary)
  if (!file.rename(temporary, manifest_path)) stop("Could not atomically update run manifest.", call. = FALSE)
}
write_manifest(manifest)
for (s in scripts) {
  status <- system2("Rscript", file.path("analysis/random_forest_subject_validation", s), env = env)
  if (status != 0) {
    manifest$status <- "failed"
    manifest$failed_stage <- s
    manifest$completed <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    write_manifest(manifest)
    stop("Validation stage failed: ", s, call. = FALSE)
  }
  manifest$completed_stages <- c(manifest$completed_stages, s)
  write_manifest(manifest)
}
manifest$status <- "complete"
manifest$completed <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
manifest$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
write_manifest(manifest)
cat(run_dir, "\n")
