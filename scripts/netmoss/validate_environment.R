#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
require_netmoss <- "--require-netmoss" %in% args

base_packages <- c("yaml", "dplyr", "tidyr", "ggplot2", "svglite", "scales", "igraph", "gridExtra")
netmoss_packages <- c(
  "NetMoss2", "rsparcc", "WGCNA", "coin", "psych", "reshape2", "doBy",
  "WeightedROC", "pROC", "plotROC", "Rmisc", "ggforce", "ggpubr", "randomForest"
)
required <- if (require_netmoss) c(base_packages, netmoss_packages) else base_packages
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)

cat("NetMoss workflow environment validation passed.\n")
for (package in required) {
  cat(sprintf("%-14s %s\n", package, as.character(utils::packageVersion(package))))
}
