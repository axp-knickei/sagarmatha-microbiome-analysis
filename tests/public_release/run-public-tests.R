#!/usr/bin/env Rscript
# Purpose: Run the approved public software checks with artificial inputs only.
# Inputs: Released code/configuration and synthetic beta fixtures; no study data.
# Outputs: Per-suite logs outside the source tree; beta-generated test results.
# Dependencies: environment.yml test/runtime packages, Rscript, Bash.
# Usage: Rscript tests/public_release/run-public-tests.R [external-log-directory]
# Scope: Software smoke/unit checks, not study-result or historical reproduction.
args <- commandArgs(trailingOnly = TRUE)
log_dir <- if (length(args)) args[[1]] else tempfile("public-test-logs-")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_dir <- normalizePath(log_dir)
if (!file.exists("config/beta_diversity_test.yml")) stop("Run from the repository root.")
rscript <- file.path(R.home("bin"), "Rscript")
r_expr <- function(expression) list(command = rscript, args = c("-e", shQuote(expression)))
r_file <- function(path) list(command = rscript, args = shQuote(path))
suites <- list(
  random_forest = r_expr('testthat::test_dir("tests/testthat", reporter="summary", stop_on_failure=TRUE)'),
  public_renderers = r_file("tests/public_release/test-public-renderers.R"),
  public_synthesis = r_file("tests/public_release/test-public-synthesis.R"),
  functional_helpers = r_file("analysis/mag_functional_annotation/run_tests.R"),
  community_helpers = r_file("analysis/functional_community_dynamics/run_tests.R"),
  netmoss_hpc_publication = r_file("scripts/netmoss/tests/testthat.R"),
  netmoss_stability = r_file("scripts/netmoss/tests/testthat-altitude-stability.R"),
  netmoss_matched_timepoint = r_file("scripts/netmoss/tests/testthat-matched-timepoint.R"),
  beta_synthetic = list(command = "bash", args = c("scripts/run_beta_diversity_local.sh",
                          "--config", "config/beta_diversity_test.yml", "--force"))
)
for (name in names(suites)) {
  suite <- suites[[name]]
  log <- file.path(log_dir, paste0(name, ".log"))
  status <- system2(suite$command, suite$args, stdout = log, stderr = log)
  cat(if (status == 0L) "PASS" else "FAIL", name, "log:", log, "\n")
  if (status != 0L) stop("Public test suite failed: ", name, ". Inspect ", log, call. = FALSE)
}
cat("All public suites passed. Generated test artifacts are not release-source files.\n")
