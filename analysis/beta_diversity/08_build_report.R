#!/usr/bin/env Rscript
# Purpose: Assemble a run inventory and per-directory manifests after the
# beta-diversity calculations complete.
# Inputs: Existing audit, QC, distance, ordination, statistics, figure, log, and
# provenance outputs in the configured results directory.
# Outputs: reports/beta_diversity_run_report.md and YAML manifests; this report
# inventories execution and is not a biological significance statement.
# Dependencies: Shared beta-diversity I/O and provenance helpers.
# Usage: Rscript analysis/beta_diversity/08_build_report.R --config config/beta_diversity.yml
# Workflow: Final reporting stage of scripts/run_beta_diversity_local.sh.

source("analysis/beta_diversity/utils/io_helpers.R")
source("analysis/beta_diversity/utils/provenance_helpers.R")
args <- parse_beta_args(); config <- load_beta_config(args$config); ensure_beta_dirs(config)
with_stage_log(config, "08_build_report", {
  distance_manifests <- list.files(beta_output_dir(config, "distances"), pattern = "_manifest.yml$", full.names = FALSE)
  figures <- list.files(beta_output_dir(config, "figures"), pattern = "\\.(pdf|png)$")
  stats <- list.files(beta_output_dir(config, "statistics"), pattern = "\\.tsv$")
  lines <- c("# Beta-diversity run report", "", paste("Generated:", timestamp_beta()), "",
             "This report is an execution inventory, not a biological significance statement.", "",
             "## Outputs", "", paste0("- Distance manifests: ", length(distance_manifests), "."),
             paste0("- Statistical tables: ", length(stats), "."), paste0("- Figure files: ", length(figures), "."), "",
             "## Scientific cautions", "", "- Altitude and timepoint are linked by the expedition schedule.",
             "- Restricted permutation accounts for subject exchangeability but is not a random-intercept PERMANOVA.",
             "- Significant PERMANOVA results must be interpreted alongside dispersion tests.",
             "- UniFrac is unavailable until a compatible branch-length MAG tree is supplied.", "",
             "## Files", "", paste0("- `", c(distance_manifests, stats, figures), "`"))
  writeLines(lines, beta_output_dir(config, "reports", "beta_diversity_run_report.md"))
  write_yaml_manifest(config, beta_output_dir(config, "reports", "manifest.yml"),
                      list(status = "COMPLETED", distance_manifests = distance_manifests, statistics = stats, figures = figures))
  for (directory in c("audit", "qc", "prepared_data", "distances", "ordination", "statistics", "figures", "tables", "reports", "logs", "provenance")) {
    path <- beta_output_dir(config, directory)
    files <- setdiff(list.files(path, recursive = TRUE), "manifest.yml")
    write_yaml_manifest(config, file.path(path, "manifest.yml"), list(directory = directory, status = "COMPLETED", files = files))
  }
})
