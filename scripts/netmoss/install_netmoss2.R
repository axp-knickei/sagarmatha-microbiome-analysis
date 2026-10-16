#!/usr/bin/env Rscript

# Install immutable NetMoss2 and r-sparcc revisions after the Conda packages in
# environment.yml have been resolved. Run this only during environment setup.

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check-only" %in% args

netmoss_repo <- "xiaolw95/NetMoss2"
netmoss_sha <- "79da905efccfc8bec9c08a4f447bf1ddb12223f6"
rsparcc_repo <- "MPBA/r-sparcc"
rsparcc_sha <- "4ac4c6c6e83219d0e970b2b8b28a9fef8c0eb30b"

cran_dependencies <- c("WeightedROC", "plotROC")
runtime_dependencies <- c(
  "yaml", "coin", "psych", "WGCNA", "pROC", "ggplot2", "reshape2",
  "Rmisc", "ggforce", "igraph", "ggpubr", "randomForest", "doBy"
)

verify_package <- function(package, expected_sha = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) return(FALSE)
  if (is.null(expected_sha)) return(TRUE)
  description <- utils::packageDescription(package)
  installed_sha <- description$RemoteSha %||% ""
  identical(tolower(installed_sha), tolower(expected_sha))
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) y else x

if (check_only) {
  required <- c(runtime_dependencies, cran_dependencies, "rsparcc", "NetMoss2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing NetMoss2 environment packages: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!verify_package("rsparcc", rsparcc_sha)) stop("rsparcc is installed but not at pinned commit ", rsparcc_sha, call. = FALSE)
  if (!verify_package("NetMoss2", netmoss_sha)) stop("NetMoss2 is installed but not at pinned commit ", netmoss_sha, call. = FALSE)
  cat("NetMoss2 environment check passed.\n")
  cat("NetMoss2: ", as.character(utils::packageVersion("NetMoss2")), " @ ", netmoss_sha, "\n", sep = "")
  cat("rsparcc: ", as.character(utils::packageVersion("rsparcc")), " @ ", rsparcc_sha, "\n", sep = "")
  quit(status = 0)
}

repos <- getOption("repos")
if (is.null(repos) || identical(unname(repos[["CRAN"]]), "@CRAN@")) {
  repos["CRAN"] <- "https://cloud.r-project.org"
  options(repos = repos)
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("Package 'remotes' is required. Update environment.yml before running this installer.", call. = FALSE)
}

missing_cran <- cran_dependencies[!vapply(cran_dependencies, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, dependencies = TRUE)
}

missing_runtime <- runtime_dependencies[!vapply(runtime_dependencies, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_runtime)) {
  stop(
    "Conda environment is missing compiled/common dependencies: ",
    paste(missing_runtime, collapse = ", "),
    ". Re-run conda env update from environment.yml.", call. = FALSE
  )
}

if (!verify_package("rsparcc", rsparcc_sha)) {
  remotes::install_github(paste0(rsparcc_repo, "@", rsparcc_sha), dependencies = FALSE, upgrade = "never")
}
if (!verify_package("NetMoss2", netmoss_sha)) {
  remotes::install_github(paste0(netmoss_repo, "@", netmoss_sha), dependencies = FALSE, upgrade = "never")
}

if (!verify_package("rsparcc", rsparcc_sha) || !verify_package("NetMoss2", netmoss_sha)) {
  stop("Pinned package verification failed after installation.", call. = FALSE)
}

cat("Pinned NetMoss2 environment installation completed successfully.\n")
