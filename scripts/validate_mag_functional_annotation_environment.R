#!/usr/bin/env Rscript

# Purpose: Verify the R dependencies required by the MAG functional workflow.
# Inputs: Active R/Conda environment.
# Outputs: Console-only version summary; no scientific files.
# Dependencies: R.
# Usage: Rscript scripts/validate_mag_functional_annotation_environment.R
# Workflow: Run before the OpenPBS analysis stages.

packages <- c("data.table", "yaml", "ggplot2", "ggrepel", "gridExtra", "scales", "svglite", "Matrix", "lme4")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
cat("MAG functional annotation environment validated\n")
cat("R:", R.version.string, "\n")
for (package in packages) cat(package, as.character(utils::packageVersion(package)), "\n")
