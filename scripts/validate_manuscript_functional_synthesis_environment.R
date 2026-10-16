#!/usr/bin/env Rscript
# Validate the lightweight runtime used by the manuscript synthesis workflow.

packages <- c("data.table", "ggplot2", "gridExtra", "grid", "yaml", "svglite", "scales")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
cat("Manuscript functional synthesis environment validated\n")
cat("R:", R.version.string, "\n")
for (package in setdiff(packages, "grid")) cat(package, as.character(utils::packageVersion(package)), "\n")
