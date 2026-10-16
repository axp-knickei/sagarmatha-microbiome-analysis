#!/usr/bin/env Rscript

# Validate packages needed for NetMoss manuscript-only post-processing.

required <- c(
  "dplyr", "tidyr", "ggplot2", "ggrepel", "igraph", "scales",
  "gridExtra", "grid", "svglite", "yaml"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing R packages: ", paste(missing, collapse = ", "),
    ". Update the existing glmm-hpc environment from environment.yml.",
    call. = FALSE
  )
}
cat("NetMoss manuscript synthesis environment validated.\n")
