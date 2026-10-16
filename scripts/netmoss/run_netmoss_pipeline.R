#!/usr/bin/env Rscript

# Prepare a conservative, subject-aggregated two-group dataset, construct two
# association networks, and run pinned NetMoss2. This is an exploratory wrapper;
# it does not make the current study design suitable for primary biomarker claims.

suppressPackageStartupMessages({
  library(dplyr)
  library(yaml)
  library(dynamicTreeCut)
  library(WGCNA)
})

`%||%` <- function(x, y) if (is.null(x) || !length(x) || all(is.na(x)) || identical(x, "")) y else x

usage <- function() {
  cat(
    "Usage: Rscript scripts/netmoss/run_netmoss_pipeline.R [options]\n",
    "  --config PATH                 Default: config/netmoss_analysis.yml\n",
    "  --output-dir PATH             Override configured output directory\n",
    "  --prepare-only                Validate and write inputs/networks; do not call NetMoss2\n",
    "  --test-mode                   Use repository test data and Spearman smoke network\n",
    "  --acknowledge-exploratory     Required for a full NetMoss2 run\n",
    "  --force                       Replace existing named outputs\n",
    "  --help\n",
    sep = ""
  )
}

parse_args <- function(args) {
  out <- list(config = "config/netmoss_analysis.yml", output_dir = NULL, prepare_only = FALSE,
              test_mode = FALSE, acknowledge = FALSE, force = FALSE)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") { usage(); quit(status = 0) }
    if (arg %in% c("--prepare-only", "--test-mode", "--acknowledge-exploratory", "--force")) {
      key <- switch(arg, `--prepare-only` = "prepare_only", `--test-mode` = "test_mode",
                    `--acknowledge-exploratory` = "acknowledge", `--force` = "force")
      out[[key]] <- TRUE
      i <- i + 1L
    } else if (arg %in% c("--config", "--output-dir")) {
      if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
      key <- if (arg == "--config") "config" else "output_dir"
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    } else stop("Unknown argument: ", arg, call. = FALSE)
  }
  out
}

read_csv <- function(path, label) {
  if (!file.exists(path)) stop(label, " file not found: ", path, call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

assert_columns <- function(df, columns, label) {
  missing <- setdiff(columns, names(df))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

assert_output <- function(path, force) {
  if (file.exists(path) && !force) stop("Refusing to overwrite ", path, ". Use --force.", call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

write_matrix_csv <- function(matrix, path, id_name = "taxon_names") {
  out <- data.frame(taxon_names = rownames(matrix), matrix, check.names = FALSE)
  names(out)[[1]] <- id_name
  utils::write.csv(out, path, row.names = FALSE, na = "NA")
}

close_columns <- function(matrix) {
  totals <- colSums(matrix, na.rm = TRUE)
  if (any(!is.finite(totals) | totals <= 0)) stop("One or more profiles have non-positive totals.", call. = FALSE)
  sweep(matrix, 2, totals, "/")
}

aggregate_subject_profiles <- function(matrix, metadata, subject_col, method) {
  subjects <- unique(as.character(metadata[[subject_col]]))
  result <- vapply(subjects, function(subject) {
    ids <- as.character(metadata$sample_id[as.character(metadata[[subject_col]]) == subject])
    values <- matrix[, ids, drop = FALSE]
    if (method == "mean") rowMeans(values, na.rm = TRUE)
    else if (method == "median") apply(values, 1, stats::median, na.rm = TRUE)
    else stop("Unsupported subject aggregation: ", method, call. = FALSE)
  }, numeric(nrow(matrix)))
  if (is.null(dim(result))) result <- matrix(result, ncol = 1, dimnames = list(rownames(matrix), subjects))
  rownames(result) <- rownames(matrix)
  colnames(result) <- subjects
  result
}

validate_network <- function(network, taxa, label) {
  if (!is.matrix(network) || nrow(network) != ncol(network)) stop(label, " is not square.", call. = FALSE)
  if (!identical(rownames(network), taxa) || !identical(colnames(network), taxa)) stop(label, " taxon ordering mismatch.", call. = FALSE)
  if (any(!is.finite(network))) stop(label, " contains non-finite values.", call. = FALSE)
  if (max(abs(network - t(network))) > 1e-8) stop(label, " is not symmetric.", call. = FALSE)
  if (max(abs(diag(network) - 1)) > 1e-8) stop(label, " diagonal is not one.", call. = FALSE)
  if (max(abs(network)) > 1 + 1e-8) stop(label, " contains correlations outside [-1, 1].", call. = FALSE)
  invisible(TRUE)
}

build_network <- function(abundance, method) {
  if (method == "sparcc") {
    if (!requireNamespace("rsparcc", quietly = TRUE)) stop("rsparcc is required for SparCC networks.", call. = FALSE)
    network <- rsparcc::sparcc(t(abundance))$CORR
  } else if (method == "spearman") {
    network <- stats::cor(t(abundance), method = "spearman", use = "pairwise.complete.obs")
  } else {
    stop("Unsupported network method: ", method, call. = FALSE)
  }
  network[!is.finite(network)] <- 0
  network <- (network + t(network)) / 2
  diag(network) <- 1
  rownames(network) <- rownames(abundance)
  colnames(network) <- rownames(abundance)
  network
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!file.exists(args$config)) stop("Configuration not found: ", args$config, call. = FALSE)
config <- yaml::read_yaml(args$config)

if (!args$prepare_only && !args$acknowledge) {
  stop(
    "Full NetMoss2 execution requires --acknowledge-exploratory. The readiness assessment does not support primary biomarker claims.",
    call. = FALSE
  )
}

if (args$test_mode) {
  config$inputs$abundance_file <- "data/test_beta_diversity/gut_test.csv"
  config$inputs$metadata_file <- "data/test_beta_diversity/sample_info_test.csv"
  config$inputs$taxonomy_file <- "data/test_beta_diversity/GDTB_test.csv"
  config$scenario$group_a_values <- c("T2", "T3")
  config$scenario$group_b_values <- c("T7", "T10")
  config$preprocessing$minimum_samples_per_group <- 3
  config$preprocessing$minimum_features <- 5
  config$preprocessing$prevalence_threshold <- 0.20
}

output_dir <- args$output_dir %||% config$outputs$output_dir
network_method <- if (args$test_mode) config$network$smoke_test_method else config$network$method
dirs <- file.path(output_dir, c("prepared", "networks", "tables", "objects", "logs", "reproducibility"))
invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
manifest_path <- file.path(output_dir, "reproducibility", "netmoss_manifest.yml")
assert_output(manifest_path, args$force)

abundance <- read_csv(config$inputs$abundance_file, "Abundance")
metadata <- read_csv(config$inputs$metadata_file, "Metadata")
taxonomy <- read_csv(config$inputs$taxonomy_file, "Taxonomy")
cols <- config$columns
assert_columns(abundance, cols$sample_id, "Abundance")
assert_columns(metadata, c(cols$sample_id, cols$subject_id, cols$timepoint), "Metadata")
assert_columns(taxonomy, c(cols$taxonomy_feature_id, cols$taxonomy_rank), "Taxonomy")

if (anyDuplicated(abundance[[cols$sample_id]])) stop("Duplicated abundance sample IDs.", call. = FALSE)
if (anyDuplicated(metadata[[cols$sample_id]])) stop("Duplicated metadata sample IDs.", call. = FALSE)

feature_ids <- setdiff(names(abundance), cols$sample_id)
numeric_ok <- vapply(abundance[feature_ids], is.numeric, logical(1))
if (!all(numeric_ok)) stop("Non-numeric abundance columns: ", paste(feature_ids[!numeric_ok], collapse = ", "), call. = FALSE)
sample_ids <- as.character(abundance[[cols$sample_id]])
metadata_ids <- as.character(metadata[[cols$sample_id]])
if (length(setdiff(sample_ids, metadata_ids))) stop("Abundance samples lack metadata.", call. = FALSE)

sample_matrix <- as.matrix(abundance[feature_ids])
storage.mode(sample_matrix) <- "double"
rownames(sample_matrix) <- sample_ids
if (any(!is.finite(sample_matrix)) || any(sample_matrix < 0)) stop("Abundance must be finite and non-negative.", call. = FALSE)

taxonomy_index <- match(feature_ids, as.character(taxonomy[[cols$taxonomy_feature_id]]))
if (anyNA(taxonomy_index)) stop("Some abundance features lack taxonomy.", call. = FALSE)
taxon_labels <- trimws(as.character(taxonomy[[cols$taxonomy_rank]][taxonomy_index]))
missing_taxon <- is.na(taxon_labels) | !nzchar(taxon_labels)
taxon_labels[missing_taxon] <- feature_ids[missing_taxon]

# rowsum on the transposed matrix aggregates MAGs to the configured rank.
taxon_by_sample <- rowsum(t(sample_matrix), group = taxon_labels, reorder = TRUE)
metadata$sample_id <- as.character(metadata[[cols$sample_id]])
metadata <- metadata[match(sample_ids, metadata$sample_id), , drop = FALSE]

excluded_timepoints <- unlist(config$scenario$excluded_values %||% character(0))
if (length(intersect(
  excluded_timepoints,
  c(unlist(config$scenario$group_a_values), unlist(config$scenario$group_b_values))
))) {
  stop("An excluded timepoint is also configured in a NetMoss group.", call. = FALSE)
}
if (length(excluded_timepoints)) {
  keep_samples <- !as.character(metadata[[cols$timepoint]]) %in% excluded_timepoints
  metadata <- metadata[keep_samples, , drop = FALSE]
  taxon_by_sample <- taxon_by_sample[, metadata$sample_id, drop = FALSE]
}

group_a_meta <- metadata[as.character(metadata[[cols$timepoint]]) %in% unlist(config$scenario$group_a_values), , drop = FALSE]
group_b_meta <- metadata[as.character(metadata[[cols$timepoint]]) %in% unlist(config$scenario$group_b_values), , drop = FALSE]
if (!nrow(group_a_meta) || !nrow(group_b_meta)) stop("Configured scenario produces an empty group.", call. = FALSE)

aggregation <- config$preprocessing$subject_aggregation
group_a <- aggregate_subject_profiles(taxon_by_sample, group_a_meta, cols$subject_id, aggregation)
group_b <- aggregate_subject_profiles(taxon_by_sample, group_b_meta, cols$subject_id, aggregation)
colnames(group_a) <- paste0(config$scenario$group_a_label, "__", colnames(group_a))
colnames(group_b) <- paste0(config$scenario$group_b_label, "__", colnames(group_b))

minimum_samples <- as.integer(config$preprocessing$minimum_samples_per_group)
if (ncol(group_a) < minimum_samples || ncol(group_b) < minimum_samples) {
  stop("Too few subject-aggregated profiles: ", ncol(group_a), " and ", ncol(group_b), call. = FALSE)
}
if (min(ncol(group_a), ncol(group_b)) < 20) {
  warning("Each group has fewer than 20 independent profiles; network results are highly unstable and exploratory.", call. = FALSE)
}

detection <- as.numeric(config$preprocessing$detection_threshold)
prevalence <- pmax(rowMeans(group_a > detection), rowMeans(group_b > detection))
keep <- prevalence >= as.numeric(config$preprocessing$prevalence_threshold) &
  apply(cbind(group_a, group_b), 1, stats::var) > 0
group_a <- group_a[keep, , drop = FALSE]
group_b <- group_b[keep, , drop = FALSE]
if (nrow(group_a) < as.integer(config$preprocessing$minimum_features)) {
  stop("Filtering retains only ", nrow(group_a), " taxa.", call. = FALSE)
}

if (isTRUE(config$preprocessing$close_to_relative_abundance)) {
  group_a <- close_columns(group_a)
  group_b <- close_columns(group_b)
}

set.seed(as.integer(config$network$random_seed))
network_a <- build_network(group_a, network_method)
network_b <- build_network(group_b, network_method)
validate_network(network_a, rownames(group_a), "Group A network")
validate_network(network_b, rownames(group_b), "Group B network")

group_a_file <- file.path(output_dir, "prepared", "group_a_abundance.csv")
group_b_file <- file.path(output_dir, "prepared", "group_b_abundance.csv")
network_a_file <- file.path(output_dir, "networks", "group_a_network.csv")
network_b_file <- file.path(output_dir, "networks", "group_b_network.csv")
write_matrix_csv(group_a, group_a_file)
write_matrix_csv(group_b, group_b_file)
write_matrix_csv(network_a, network_a_file)
write_matrix_csv(network_b, network_b_file)

sample_manifest <- bind_rows(
  group_a_meta %>% transmute(sample_id, subject_id = as.character(.data[[cols$subject_id]]), group = config$scenario$group_a_label),
  group_b_meta %>% transmute(sample_id, subject_id = as.character(.data[[cols$subject_id]]), group = config$scenario$group_b_label)
)
utils::write.csv(sample_manifest, file.path(output_dir, "prepared", "sample_group_manifest.csv"), row.names = FALSE)
filter_table <- data.frame(taxon_names = names(prevalence), maximum_group_prevalence = prevalence,
                           retained = names(prevalence) %in% rownames(group_a))
utils::write.csv(filter_table, file.path(output_dir, "prepared", "taxon_filtering.csv"), row.names = FALSE)

status <- "PREPARED_ONLY"
if (!args$prepare_only) {
  if (!requireNamespace("NetMoss2", quietly = TRUE)) stop("NetMoss2 is not installed; run hpc/update_netmoss_env.pbs.", call. = FALSE)
  case_data <- data.frame(X = rownames(group_b), group_b, check.names = FALSE)
  control_data <- data.frame(X = rownames(group_a), group_a, check.names = FALSE)
  result <- NetMoss2::NetMoss(
    case_dir = case_data,
    control_dir = control_data,
    net_case_dir = network_b,
    net_control_dir = network_a,
    scaled = isTRUE(config$netmoss$scaled),
    deepSplit = as.integer(config$netmoss$deep_split),
    minModuleSize = as.integer(config$netmoss$minimum_module_size),
    soft_pow = as.numeric(config$netmoss$soft_power)
  )
  saveRDS(result, file.path(output_dir, "objects", "netmoss_result.rds"))
  scores <- as.data.frame(result[[1]], stringsAsFactors = FALSE)
  if (!"taxon_names" %in% names(scores)) scores$taxon_names <- rownames(scores)
  utils::write.csv(scores, file.path(output_dir, "tables", "netmoss_scores.csv"), row.names = FALSE)
  status <- "COMPLETED_EXPLORATORY"
}

manifest <- list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  status = status,
  analysis_status = config$project$analysis_status,
  scenario = config$scenario$name,
  group_a = list(label = config$scenario$group_a_label, timepoints = unlist(config$scenario$group_a_values),
                 raw_samples = nrow(group_a_meta), independent_profiles = ncol(group_a)),
  group_b = list(label = config$scenario$group_b_label, timepoints = unlist(config$scenario$group_b_values),
                 raw_samples = nrow(group_b_meta), independent_profiles = ncol(group_b)),
  taxonomy_rank = cols$taxonomy_rank,
  subject_aggregation = aggregation,
  excluded_timepoints = excluded_timepoints,
  prevalence_threshold = config$preprocessing$prevalence_threshold,
  retained_taxa = nrow(group_a),
  network_method = network_method,
  random_seed = config$network$random_seed,
  netmoss2_commit = config$netmoss$package_commit,
  git_commit = tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "unknown"),
  caution = "Exploratory/supporting only; the readiness assessment found insufficient independent participants for robust primary NetMoss inference."
)
yaml::write_yaml(manifest, manifest_path)
old_tz <- Sys.getenv("TZ", unset = "")
Sys.setenv(TZ = "UTC")
suppressWarnings(utils::capture.output(utils::sessionInfo(), file = file.path(output_dir, "reproducibility", "sessionInfo.txt")))
if (nzchar(old_tz)) Sys.setenv(TZ = old_tz) else Sys.unsetenv("TZ")
writeLines(c(
  paste0("Status: ", status),
  paste0("Scenario: ", config$scenario$name),
  paste0("Independent profiles: ", ncol(group_a), " versus ", ncol(group_b)),
  paste0("Retained ", cols$taxonomy_rank, " taxa: ", nrow(group_a)),
  paste0("Network method: ", network_method),
  "Interpretation: exploratory/supporting only; do not report as validated biomarkers."
), file.path(output_dir, "logs", "netmoss_run_summary.txt"))

cat("NetMoss pipeline status: ", status, "\n", sep = "")
cat("Output directory: ", normalizePath(output_dir, mustWork = TRUE), "\n", sep = "")
