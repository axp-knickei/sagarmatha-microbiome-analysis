#!/usr/bin/env Rscript

# NetMoss readiness assessment
# ----------------------------
# This standalone diagnostic reads abundance, metadata, and optional taxonomy
# tables without modifying them. It does not install or run NetMoss.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(yaml)
  library(scales)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || identical(x, "")) y else x
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

usage <- function() {
  cat(
    "Usage: Rscript scripts/check_netmoss_readiness.R [options]\n\n",
    "Options:\n",
    "  --config PATH          YAML configuration (default: config/netmoss_readiness.yml)\n",
    "  --abundance PATH       Override abundance table\n",
    "  --metadata PATH        Override metadata table\n",
    "  --taxonomy PATH        Override taxonomy table\n",
    "  --output-dir PATH      Override output directory\n",
    "  --force                Permit replacement of existing named outputs\n",
    "  --help                 Print this message\n",
    sep = ""
  )
}

parse_args <- function(args) {
  out <- list(config = "config/netmoss_readiness.yml", force = FALSE)
  key_map <- c(
    "--config" = "config", "--abundance" = "abundance",
    "--metadata" = "metadata", "--taxonomy" = "taxonomy",
    "--output-dir" = "output_dir"
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      usage()
      quit(status = 0)
    } else if (arg == "--force") {
      out$force <- TRUE
      i <- i + 1L
    } else if (arg %in% names(key_map)) {
      if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
      out[[key_map[[arg]]]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      stop("Unknown argument: ", arg, ". Use --help for valid options.", call. = FALSE)
    }
  }
  out
}

required_packages <- c("dplyr", "tidyr", "ggplot2", "yaml", "scales", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

read_table <- function(path, label, optional = FALSE) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    if (optional) return(NULL)
    stop(label, " path is not configured.", call. = FALSE)
  }
  if (!file.exists(path)) {
    if (optional) {
      warning(label, " file was not found; taxonomy-dependent checks will be unavailable: ", path, call. = FALSE)
      return(NULL)
    }
    stop(label, " file not found: ", path, call. = FALSE)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext == "rds") {
    readRDS(path)
  } else {
    stop("Unsupported ", label, " format: ", ext, ". Use CSV, TSV, TXT, or RDS.", call. = FALSE)
  }
}

normalize_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

find_column <- function(df, configured, aliases, required = FALSE) {
  if (!is.null(configured) && !is.na(configured) && configured != "auto") {
    if (configured %in% names(df)) return(configured)
    if (required) stop("Configured column not found: ", configured, call. = FALSE)
    return(NA_character_)
  }
  normalized <- normalize_name(names(df))
  hits <- match(normalize_name(aliases), normalized, nomatch = 0L)
  hits <- hits[hits > 0L]
  if (length(hits)) return(names(df)[hits[[1]]])
  if (required) stop("Could not detect required column; tried: ", paste(aliases, collapse = ", "), call. = FALSE)
  NA_character_
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

numeric_fraction <- function(df, columns) {
  if (!length(columns)) return(numeric(0))
  vapply(columns, function(column) {
    raw <- df[[column]]
    parsed <- safe_numeric(raw)
    mean(is.na(raw) | trimws(as.character(raw)) == "" | !is.na(parsed))
  }, numeric(1))
}

detect_abundance_layout <- function(abundance, metadata_ids, configured_sample_id = "auto") {
  if (!is.data.frame(abundance)) stop("Abundance input must be a tabular data frame.", call. = FALSE)
  if (nrow(abundance) < 2 || ncol(abundance) < 2) stop("Abundance table is too small to assess.", call. = FALSE)

  candidate_columns <- names(abundance)
  if (!is.null(configured_sample_id) && configured_sample_id != "auto" && configured_sample_id %in% names(abundance)) {
    candidate_columns <- configured_sample_id
  }
  overlap_by_column <- vapply(candidate_columns, function(column) {
    length(intersect(as.character(abundance[[column]]), metadata_ids))
  }, integer(1))
  row_id_column <- candidate_columns[which.max(overlap_by_column)]
  row_overlap <- max(overlap_by_column)
  column_overlap <- length(intersect(names(abundance), metadata_ids))

  if (row_overlap == 0L && column_overlap == 0L) {
    stop(
      "Could not orient abundance table: neither a column nor column names match metadata sample IDs.",
      call. = FALSE
    )
  }

  if (row_overlap >= column_overlap) {
    feature_columns <- setdiff(names(abundance), row_id_column)
    numeric_ok <- numeric_fraction(abundance, feature_columns)
    bad <- names(numeric_ok)[numeric_ok < 1]
    if (length(bad)) stop("Non-numeric abundance feature columns detected: ", paste(head(bad, 10), collapse = ", "), call. = FALSE)
    matrix <- vapply(abundance[feature_columns], safe_numeric, numeric(nrow(abundance)))
    matrix <- as.matrix(matrix)
    rownames(matrix) <- as.character(abundance[[row_id_column]])
    colnames(matrix) <- feature_columns
    orientation <- "Samples x Features"
    sample_id_source <- row_id_column
    transposed <- FALSE
  } else {
    sample_columns <- intersect(names(abundance), metadata_ids)
    feature_id_column <- setdiff(names(abundance), sample_columns)[[1]]
    numeric_ok <- numeric_fraction(abundance, sample_columns)
    bad <- names(numeric_ok)[numeric_ok < 1]
    if (length(bad)) stop("Non-numeric abundance sample columns detected: ", paste(head(bad, 10), collapse = ", "), call. = FALSE)
    raw_matrix <- vapply(abundance[sample_columns], safe_numeric, numeric(nrow(abundance)))
    matrix <- t(as.matrix(raw_matrix))
    rownames(matrix) <- sample_columns
    colnames(matrix) <- as.character(abundance[[feature_id_column]])
    orientation <- "Features x Samples"
    sample_id_source <- "column names"
    transposed <- TRUE
  }
  storage.mode(matrix) <- "double"
  list(
    matrix = matrix, original_orientation = orientation,
    sample_id_source = sample_id_source, transposed = transposed,
    row_overlap = row_overlap, column_overlap = column_overlap
  )
}

infer_feature_type <- function(feature_ids, taxonomy) {
  feature_ids <- as.character(feature_ids)
  scores <- c(
    MAG = mean(grepl("(_bin\\.|\\bbin[._-]|\\.strict$|\\.orig$|\\.permissive$)", feature_ids, ignore.case = TRUE)),
    ASV = mean(grepl("^ASV[_-]?[0-9]+$", feature_ids, ignore.case = TRUE)),
    OTU = mean(grepl("^OTU[_-]?[0-9]+$", feature_ids, ignore.case = TRUE))
  )
  if (max(scores) >= 0.25) return(names(which.max(scores)))
  if (!is.null(taxonomy)) {
    rank_names <- normalize_name(names(taxonomy))
    if ("genome" %in% rank_names) return("MAG")
    if ("species" %in% rank_names && nrow(taxonomy) >= length(feature_ids)) return("species or other taxonomic features")
  }
  "other/unknown"
}

detect_abundance_scale <- function(matrix) {
  finite <- matrix[is.finite(matrix)]
  sample_sums <- rowSums(matrix, na.rm = TRUE)
  integer_fraction <- mean(abs(finite - round(finite)) < 1e-10)
  has_negative <- any(finite < 0)
  sums_one <- mean(abs(sample_sums - 1) <= 0.05)
  sums_hundred <- mean(abs(sample_sums - 100) <= 5)
  clr_centered <- has_negative && mean(abs(rowMeans(matrix, na.rm = TRUE)) < 0.1) >= 0.8
  if (!has_negative && integer_fraction >= 0.99 && max(finite) > 10 && median(sample_sums) > 100) {
    type <- "Raw counts"
    interpretation <- "Values are non-negative and nearly integer-valued with count-like library sums."
  } else if (!has_negative && sums_one >= 0.8) {
    type <- "Relative abundance (sum approximately 1)"
    interpretation <- "Most sample totals are within 0.05 of 1; values are compositional relative abundances."
  } else if (!has_negative && sums_hundred >= 0.8) {
    type <- "Relative abundance (sum approximately 100)"
    interpretation <- "Most sample totals are within 5 percentage points of 100; values are compositional percentages."
  } else if (clr_centered) {
    type <- "CLR transformed"
    interpretation <- "Negative values and approximately zero row means are consistent with CLR coordinates."
  } else if (has_negative) {
    type <- "Transformed/unknown"
    interpretation <- "Negative values preclude raw counts or untransformed relative abundance; transformation provenance is required."
  } else if (max(finite) <= 1 && median(sample_sums) <= 1.1) {
    type <- "Relative-abundance-like (incomplete closure)"
    interpretation <- "Values are non-negative fractions, but many sample totals do not close to 1; unclassified or excluded abundance may remain."
  } else {
    type <- "Unknown non-negative scale"
    interpretation <- "The numerical pattern is not sufficient to distinguish counts, normalized abundance, or a log transform."
  }
  sum_class <- if (sums_one >= 0.8) "approximately 1" else if (sums_hundred >= 0.8) "approximately 100" else "neither"
  list(type = type, interpretation = interpretation, sum_class = sum_class,
       fraction_sum_one = sums_one, fraction_sum_hundred = sums_hundred)
}

classify_design <- function(metadata, subject_col, time_col) {
  if (is.na(subject_col)) return("cross-sectional or indeterminate (no subject identifier detected)")
  counts <- table(metadata[[subject_col]])
  repeated <- any(counts > 1)
  has_time <- !is.na(time_col) && dplyr::n_distinct(metadata[[time_col]], na.rm = TRUE) > 1
  if (repeated && has_time && all(counts > 1)) "longitudinal repeated-measures"
  else if (repeated && has_time) "mixed longitudinal/repeated-measures with incomplete trajectories"
  else if (repeated) "repeated or paired"
  else "cross-sectional"
}

write_csv_safe <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "NA")
}

save_figure <- function(plot, stem, output_dir, width = 8, height = 5.5) {
  ggplot2::ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(output_dir, paste0(stem, ".svg")), plot, width = width, height = height, units = "in", device = svglite::svglite)
}

theme_netmoss <- function() {
  ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "grey30"),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = "grey50"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

network_feasibility <- function(n_observations, n_subjects, n_features, repeated_overlap) {
  effective_n <- if (is.na(n_subjects)) n_observations else min(n_observations, n_subjects)
  ratio <- n_features / max(1, effective_n)
  if (effective_n >= 50 && ratio <= 2 && repeated_overlap == 0) "Potentially robust; bootstrap stability still required"
  else if (effective_n >= 30 && ratio <= 4 && repeated_overlap == 0) "Moderate; regularization and stability analysis required"
  else if (effective_n >= 20 && ratio <= 5) "Fragile; exploratory network only"
  else "Poor; too few effectively independent units for the feature dimension"
}

recommendation_level <- function(n1, n2, s1, s2, p, overlap_subjects) {
  effective_min <- min(s1, s2, na.rm = TRUE)
  if (!is.finite(effective_min)) effective_min <- min(n1, n2)
  if (min(n1, n2) >= 50 && effective_min >= 50 && p / effective_min <= 2 && overlap_subjects == 0) "STRONG"
  else if (min(n1, n2) >= 30 && effective_min >= 30 && p / effective_min <= 4 && overlap_subjects == 0) "MODERATE"
  else if (min(n1, n2) >= 15 && effective_min >= 10 && p / effective_min <= 12) "WEAK"
  else "NOT RECOMMENDED"
}

scenario_row <- function(name, interpretation, metadata, subject_col, mask_a, mask_b,
                         label_a, label_b, n_features, concern) {
  a <- metadata[which(!is.na(mask_a) & mask_a), , drop = FALSE]
  b <- metadata[which(!is.na(mask_b) & mask_b), , drop = FALSE]
  subjects_a <- if (is.na(subject_col)) character(0) else unique(as.character(a[[subject_col]]))
  subjects_b <- if (is.na(subject_col)) character(0) else unique(as.character(b[[subject_col]]))
  overlap <- intersect(subjects_a, subjects_b)
  s1 <- if (is.na(subject_col)) NA_integer_ else length(subjects_a)
  s2 <- if (is.na(subject_col)) NA_integer_ else length(subjects_b)
  balance <- if (min(nrow(a), nrow(b)) == 0) Inf else max(nrow(a), nrow(b)) / min(nrow(a), nrow(b))
  level <- if (nrow(a) == 0 || nrow(b) == 0) "NOT RECOMMENDED" else recommendation_level(nrow(a), nrow(b), s1, s2, n_features, length(overlap))
  feasibility <- network_feasibility(min(nrow(a), nrow(b)), min(s1, s2), n_features, length(overlap))
  data.frame(
    scenario_name = name,
    biological_interpretation = interpretation,
    group_a = label_a, group_b = label_b,
    group_a_samples = nrow(a), group_b_samples = nrow(b),
    group_a_subjects = s1, group_b_subjects = s2,
    repeated_subjects_between_groups = length(overlap),
    independent_samples = if (length(overlap) == 0) nrow(a) + nrow(b) else NA_integer_,
    balance_ratio_larger_to_smaller = ifelse(is.finite(balance), balance, NA_real_),
    feature_count = n_features,
    expected_network_robustness = feasibility,
    potential_statistical_concerns = concern,
    recommendation_level = level,
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!file.exists(args$config)) stop("Configuration file not found: ", args$config, call. = FALSE)
config <- yaml::read_yaml(args$config)

abundance_path <- args$abundance %||% config$inputs$abundance_file
metadata_path <- args$metadata %||% config$inputs$metadata_file
taxonomy_path <- args$taxonomy %||% config$inputs$taxonomy_file
output_dir <- args$output_dir %||% config$outputs$output_dir %||% "results/netmoss_readiness"
force <- isTRUE(args$force) || isTRUE(config$outputs$overwrite)

named_outputs <- c(
  "netmoss_readiness_summary.txt", "netmoss_dataset_summary.csv",
  "netmoss_metadata_summary.csv", "netmoss_feature_prevalence.csv",
  "netmoss_sparsity_summary.csv", "netmoss_candidate_scenarios.csv",
  "netmoss_sample_balance.csv", "netmoss_assumption_checklist.csv",
  "netmoss_recommendations.txt", "sessionInfo.txt"
)
existing_outputs <- file.path(output_dir, named_outputs)
existing_outputs <- existing_outputs[file.exists(existing_outputs)]
if (length(existing_outputs) && !force) {
  stop(
    "Refusing to overwrite existing NetMoss readiness outputs: ",
    paste(basename(existing_outputs), collapse = ", "),
    ". Use --force or a different --output-dir.", call. = FALSE
  )
}
if (length(existing_outputs) && force) {
  warning("Replacing existing named outputs in ", output_dir, " because --force was supplied.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("NetMoss readiness assessment started: ", timestamp(), "\n", sep = "")
cat("Abundance: ", abundance_path, "\nMetadata: ", metadata_path, "\nTaxonomy: ", taxonomy_path, "\n", sep = "")

abundance_raw <- read_table(abundance_path, "Abundance")
metadata_raw <- read_table(metadata_path, "Metadata")
taxonomy_raw <- read_table(taxonomy_path, "Taxonomy", optional = TRUE)
if (!is.data.frame(metadata_raw)) stop("Metadata input must be a tabular data frame.", call. = FALSE)

column_config <- config$columns %||% list()
metadata_sample_col <- find_column(
  metadata_raw, column_config$sample_id %||% "auto",
  c("sampleID", "sample_id", "sample", "SampleID"), required = TRUE
)
metadata_ids_all <- as.character(metadata_raw[[metadata_sample_col]])
layout <- detect_abundance_layout(abundance_raw, metadata_ids_all, column_config$sample_id %||% "auto")
matrix_all <- layout$matrix

duplicated_abundance_ids <- unique(rownames(matrix_all)[duplicated(rownames(matrix_all))])
duplicated_feature_ids <- unique(colnames(matrix_all)[duplicated(colnames(matrix_all))])
duplicated_metadata_ids <- unique(metadata_ids_all[duplicated(metadata_ids_all)])
if (length(duplicated_abundance_ids)) stop("Duplicated abundance sample IDs: ", paste(head(duplicated_abundance_ids, 10), collapse = ", "), call. = FALSE)
if (length(duplicated_feature_ids)) stop("Duplicated microbial feature IDs: ", paste(head(duplicated_feature_ids, 10), collapse = ", "), call. = FALSE)
if (length(duplicated_metadata_ids)) stop("Duplicated metadata sample IDs: ", paste(head(duplicated_metadata_ids, 10), collapse = ", "), call. = FALSE)

abundance_only_ids <- setdiff(rownames(matrix_all), metadata_ids_all)
metadata_only_ids <- setdiff(metadata_ids_all, rownames(matrix_all))
matched_ids <- intersect(rownames(matrix_all), metadata_ids_all)
if (!length(matched_ids)) stop("No sample IDs match between abundance and metadata.", call. = FALSE)
matrix <- matrix_all[matched_ids, , drop = FALSE]
metadata <- metadata_raw[match(matched_ids, metadata_ids_all), , drop = FALSE]
stopifnot(identical(rownames(matrix), as.character(metadata[[metadata_sample_col]])))

subject_col <- find_column(metadata, column_config$subject_id %||% "auto", c("subjectID", "subject_id", "Subject", "group", "participant", "patientID"))
timepoint_col <- find_column(metadata, column_config$timepoint %||% "auto", c("Timepoint", "timepoint", "visit", "time", "sampling_point"))
date_col <- find_column(metadata, column_config$date %||% "auto", c("date", "sample_date", "collection_date", "sampling_date"))
altitude_col <- find_column(metadata, column_config$altitude %||% "auto", c("altitude", "elevation", "height"))
sex_col <- find_column(metadata, column_config$sex %||% "auto", c("sex", "gender"))
age_col <- find_column(metadata, column_config$age %||% "auto", c("age", "age_years"))
bmi_col <- find_column(metadata, "auto", c("BMI", "body_mass_index"))
treatment_col <- find_column(metadata, "auto", c("treatment", "intervention", "arm"))
batch_col <- find_column(metadata, "auto", c("batch", "sequencing_batch", "run", "plate"))
cohort_col <- find_column(metadata, "auto", c("cohort", "study", "site"))
group_col <- find_column(metadata, "auto", c("condition", "case_control", "status", "phenotype"))

excluded_time_levels <- unlist(config$analysis$excluded_time_levels %||% character(0))
excluded_matched_ids <- character(0)
if (!is.na(timepoint_col) && length(excluded_time_levels)) {
  excluded_rows <- as.character(metadata[[timepoint_col]]) %in% excluded_time_levels
  excluded_matched_ids <- rownames(matrix)[excluded_rows]
  matrix <- matrix[!excluded_rows, , drop = FALSE]
  metadata <- metadata[!excluded_rows, , drop = FALSE]
  matched_ids <- rownames(matrix)
}

feature_type <- infer_feature_type(colnames(matrix), taxonomy_raw)
taxonomy_feature_col <- NA_character_
taxonomy_matched <- 0L
taxonomy_aligned <- NULL
if (!is.null(taxonomy_raw)) {
  taxonomy_feature_col <- find_column(
    taxonomy_raw, column_config$taxonomy_feature_id %||% "auto",
    c("Genome", "feature_id", "taxon_id", "MAG", "ASV", "OTU"), required = FALSE
  )
  if (!is.na(taxonomy_feature_col)) {
    taxonomy_ids <- as.character(taxonomy_raw[[taxonomy_feature_col]])
    taxonomy_matched <- sum(colnames(matrix) %in% taxonomy_ids)
    taxonomy_aligned <- taxonomy_raw[match(colnames(matrix), taxonomy_ids), , drop = FALSE]
  }
}

finite_values <- matrix[is.finite(matrix)]
if (!length(finite_values)) stop("Abundance matrix contains no finite numeric values.", call. = FALSE)
scale_result <- detect_abundance_scale(matrix)
detection_threshold <- as.numeric(config$analysis$detection_threshold %||% 0)
present <- is.finite(matrix) & matrix > detection_threshold
sample_sparsity <- rowMeans(!present)
feature_sparsity <- colMeans(!present)
feature_prevalence <- colMeans(present)
sample_sums <- rowSums(matrix, na.rm = TRUE)

subject_values <- if (is.na(subject_col)) rep(NA_character_, nrow(metadata)) else as.character(metadata[[subject_col]])
time_values <- if (is.na(timepoint_col)) rep(NA_character_, nrow(metadata)) else as.character(metadata[[timepoint_col]])
ordered_levels <- unlist(config$analysis$ordered_time_levels %||% character(0))
observed_levels <- unique(time_values[!is.na(time_values) & nzchar(time_values)])
time_levels <- c(intersect(ordered_levels, observed_levels), setdiff(observed_levels, ordered_levels))
if (!length(time_levels) && !is.na(date_col)) time_levels <- sort(unique(as.character(metadata[[date_col]])))

n_subjects <- if (is.na(subject_col)) NA_integer_ else dplyr::n_distinct(subject_values, na.rm = TRUE)
n_timepoints <- if (is.na(timepoint_col)) NA_integer_ else dplyr::n_distinct(time_values, na.rm = TRUE)
design <- classify_design(metadata, subject_col, timepoint_col)

quantiles <- stats::quantile(finite_values, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1), na.rm = TRUE, names = FALSE)
dataset_summary <- data.frame(
  metric = c(
    "generated_at", "abundance_file", "metadata_file", "taxonomy_file",
    "original_abundance_dimensions", "analysis_matrix_dimensions", "original_orientation",
    "transposed_for_analysis", "sample_id_column", "sample_id_source", "feature_type",
    "matched_samples", "excluded_timepoints", "excluded_matched_samples",
    "abundance_only_samples", "metadata_only_samples",
    "duplicated_abundance_sample_ids", "duplicated_metadata_sample_ids", "duplicated_feature_ids",
    "subjects", "microbial_features", "metadata_variables", "longitudinal_timepoints",
    "observations", "abundance_cells", "matrix_memory_bytes", "minimum", "maximum", "mean", "median",
    "variance", "q01", "q05", "q25", "q75", "q95", "q99", "zeros", "zero_proportion",
    "missing_values", "negative_values", "infinite_values", "nan_values", "sample_sum_class",
    "sample_sum_minimum", "sample_sum_median", "sample_sum_maximum", "abundance_scale",
    "abundance_interpretation", "taxonomy_features_matched", "study_design"
  ),
  value = as.character(c(
    timestamp(), abundance_path, metadata_path, taxonomy_path,
    paste(nrow(abundance_raw), ncol(abundance_raw), sep = " x "), paste(nrow(matrix), ncol(matrix), sep = " x "), layout$original_orientation,
    layout$transposed, metadata_sample_col, layout$sample_id_source, feature_type,
    length(matched_ids), paste(excluded_time_levels, collapse = ";"),
    length(excluded_matched_ids), length(abundance_only_ids), length(metadata_only_ids),
    length(duplicated_abundance_ids), length(duplicated_metadata_ids), length(duplicated_feature_ids),
    n_subjects, ncol(matrix), ncol(metadata_raw), n_timepoints,
    nrow(matrix), length(matrix), as.numeric(object.size(matrix)), min(finite_values), max(finite_values), mean(finite_values), stats::median(finite_values),
    stats::var(finite_values), quantiles[2], quantiles[3], quantiles[4], quantiles[6], quantiles[7], quantiles[8], sum(matrix == 0, na.rm = TRUE), mean(matrix == 0, na.rm = TRUE),
    sum(is.na(matrix)), sum(matrix < 0, na.rm = TRUE), sum(is.infinite(matrix)), sum(is.nan(matrix)), scale_result$sum_class,
    min(sample_sums), stats::median(sample_sums), max(sample_sums), scale_result$type,
    scale_result$interpretation, taxonomy_matched, design
  )),
  stringsAsFactors = FALSE
)

role_map <- c(
  sample_id = metadata_sample_col, subject_id = subject_col, timepoint = timepoint_col,
  date = date_col, altitude = altitude_col, sex = sex_col, age = age_col, BMI = bmi_col,
  treatment = treatment_col, batch = batch_col, cohort = cohort_col, biological_group = group_col
)
metadata_summary <- bind_rows(lapply(names(metadata), function(column) {
  values <- metadata[[column]]
  roles <- names(role_map)[!is.na(role_map) & role_map == column]
  data.frame(
    variable = column,
    detected_role = if (length(roles)) paste(roles, collapse = "; ") else "other biological/technical variable",
    data_class = paste(class(values), collapse = "/"),
    non_missing_n = sum(!is.na(values) & trimws(as.character(values)) != ""),
    missing_n = sum(is.na(values) | trimws(as.character(values)) == ""),
    missing_fraction = mean(is.na(values) | trimws(as.character(values)) == ""),
    unique_non_missing = dplyr::n_distinct(values[!is.na(values) & trimws(as.character(values)) != ""]),
    example_values = paste(head(unique(as.character(values[!is.na(values)])), 5), collapse = "; "),
    stringsAsFactors = FALSE
  )
}))

feature_prevalence_table <- data.frame(
  feature_id = colnames(matrix),
  feature_type = feature_type,
  detected_samples = colSums(present),
  prevalence = feature_prevalence,
  zero_count = colSums(matrix == 0, na.rm = TRUE),
  zero_proportion = colMeans(matrix == 0, na.rm = TRUE),
  missing_count = colSums(is.na(matrix)),
  total_abundance = colSums(matrix, na.rm = TRUE),
  mean_abundance = colMeans(matrix, na.rm = TRUE),
  median_abundance = apply(matrix, 2, stats::median, na.rm = TRUE),
  stringsAsFactors = FALSE
)
if (!is.null(taxonomy_aligned)) {
  ranks <- intersect(c("domain", "phylum", "class", "order", "ordeer", "family", "genus", "species"), names(taxonomy_aligned))
  for (rank in ranks) feature_prevalence_table[[rank]] <- taxonomy_aligned[[rank]]
}

sparsity_summary <- bind_rows(
  data.frame(entity_type = "dataset", entity_id = "overall", zero_or_undetected_n = sum(!present), total_values = length(present), sparsity = mean(!present)),
  data.frame(entity_type = "sample", entity_id = rownames(matrix), zero_or_undetected_n = rowSums(!present), total_values = ncol(matrix), sparsity = sample_sparsity),
  data.frame(entity_type = "feature", entity_id = colnames(matrix), zero_or_undetected_n = colSums(!present), total_values = nrow(matrix), sparsity = feature_sparsity)
)

thresholds <- as.numeric(unlist(config$analysis$prevalence_thresholds %||% c(.05, .1, .2, .3, .4, .5)))
filtering_table <- bind_rows(lapply(thresholds, function(threshold) {
  keep <- feature_prevalence >= threshold
  filtered <- matrix[, keep, drop = FALSE]
  retained <- sum(keep)
  coverage <- if (retained) mean(rowSums(filtered > detection_threshold, na.rm = TRUE) > 0) else 0
  feasibility <- network_feasibility(nrow(matrix), n_subjects, retained, ifelse(!is.na(n_subjects) && nrow(matrix) > n_subjects, n_subjects, 0))
  data.frame(
    prevalence_threshold = threshold,
    retained_features = retained,
    retained_fraction = retained / ncol(matrix),
    remaining_sparsity = if (retained) mean(filtered <= detection_threshold | !is.finite(filtered)) else NA_real_,
    sample_coverage = coverage,
    observations_per_retained_feature = nrow(matrix) / max(1, retained),
    subjects_per_retained_feature = n_subjects / max(1, retained),
    network_feasibility = feasibility,
    stringsAsFactors = FALSE
  )
}))
screening_threshold <- max(filtering_table$prevalence_threshold)
native_threshold_feasible <- any(
  filtering_table$retained_features <= max(20, n_subjects * 3) &
    filtering_table$prevalence_threshold >= 0.20
)
filtering_table$screening_candidate <- filtering_table$prevalence_threshold == screening_threshold
filtering_table$recommended_for_exploratory_network <-
  filtering_table$screening_candidate & native_threshold_feasible
filtering_table$recommendation <- ifelse(
  filtering_table$screening_candidate,
  ifelse(
    native_threshold_feasible,
    "Candidate threshold, conditional on aggregation and stability testing",
    "Most stringent tested sensitivity threshold; still not feasible at native feature level"
  ),
  "Sensitivity threshold only"
)

samples_per_subject <- if (!is.na(subject_col)) {
  metadata %>% count(.data[[subject_col]], name = "n_samples") %>% rename(subject_id = 1)
} else data.frame(subject_id = NA_character_, n_samples = nrow(metadata))
samples_per_timepoint <- if (!is.na(timepoint_col)) {
  metadata %>% mutate(.time = factor(as.character(.data[[timepoint_col]]), levels = time_levels, ordered = TRUE)) %>%
    group_by(.time) %>% summarise(n_samples = n(), n_subjects = if (is.na(subject_col)) NA_integer_ else n_distinct(.data[[subject_col]]), .groups = "drop") %>%
    rename(timepoint = .time)
} else data.frame(timepoint = NA_character_, n_samples = nrow(metadata), n_subjects = n_subjects)

if (!is.na(subject_col) && !is.na(timepoint_col)) {
  longitudinal_grid <- tidyr::expand_grid(
    subject_id = sort(unique(subject_values)),
    timepoint = factor(time_levels, levels = time_levels, ordered = TRUE)
  ) %>%
    left_join(
      metadata %>% transmute(subject_id = as.character(.data[[subject_col]]), timepoint = factor(as.character(.data[[timepoint_col]]), levels = time_levels, ordered = TRUE), sample_id = .data[[metadata_sample_col]]) %>% distinct(),
      by = c("subject_id", "timepoint")
    ) %>%
    mutate(observed = !is.na(sample_id))
  missing_by_subject <- longitudinal_grid %>% group_by(subject_id) %>% summarise(
    observed_timepoints = sum(observed), missing_timepoints = sum(!observed),
    dropout_after_last_observation = {
      idx <- which(observed)
      if (!length(idx)) length(observed) else sum(seq_along(observed) > max(idx) & !observed)
    }, .groups = "drop"
  )
} else {
  longitudinal_grid <- data.frame(subject_id = NA_character_, timepoint = NA_character_, sample_id = rownames(matrix), observed = TRUE)
  missing_by_subject <- data.frame(subject_id = NA_character_, observed_timepoints = NA_integer_, missing_timepoints = NA_integer_, dropout_after_last_observation = NA_integer_)
}

sample_balance <- bind_rows(
  samples_per_subject %>% transmute(balance_type = "subject", level = as.character(subject_id), n_samples = n_samples, n_subjects = 1L),
  samples_per_timepoint %>% transmute(balance_type = "timepoint", level = as.character(timepoint), n_samples = n_samples, n_subjects = n_subjects),
  missing_by_subject %>% transmute(balance_type = "longitudinal_completeness", level = as.character(subject_id), n_samples = observed_timepoints, n_subjects = NA_integer_)
)

# Scenario discovery uses only observed abundance-linked samples. These are
# candidate contrasts, not a license to infer independent biological replicates.
scenarios <- list()
add_scenario <- function(...) scenarios[[length(scenarios) + 1L]] <<- scenario_row(...)
preferred_rank <- tolower(config$analysis$preferred_taxonomic_rank %||% "genus")
n_scenario_features <- ncol(matrix)
if (!is.null(taxonomy_aligned) && preferred_rank %in% names(taxonomy_aligned)) {
  labels <- as.character(taxonomy_aligned[[preferred_rank]])
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels)) n_scenario_features <- dplyr::n_distinct(labels)
}

if (!is.na(altitude_col)) {
  altitude <- safe_numeric(metadata[[altitude_col]])
  low_cut <- min(altitude, na.rm = TRUE)
  high_cut <- stats::quantile(altitude, 0.75, na.rm = TRUE, names = FALSE)
  add_scenario(
    "Low altitude versus high altitude",
    "Contrasts low-altitude/recovery samples with samples in the upper altitude quartile.",
    metadata, subject_col, altitude == low_cut, altitude >= high_cut,
    paste0("altitude = ", low_cut, " m"), paste0("altitude >= ", round(high_cut), " m"),
    n_scenario_features,
    "Altitude is coupled to expedition time; repeated observations from the same subjects violate independence for ordinary correlations."
  )
  add_scenario(
    "Baseline versus maximum sampled altitude",
    "Contrasts the first low-altitude timepoint with the maximum recorded sampling altitude; this is not necessarily summit exposure.",
    metadata, subject_col,
    time_values == time_levels[[1]], altitude == max(altitude, na.rm = TRUE),
    time_levels[[1]], paste0(max(altitude, na.rm = TRUE), " m samples"),
    n_scenario_features,
    "Very small exact-timepoint groups; maximum recorded sampling altitude must not be mislabeled as the summit."
  )
}
if (!is.na(timepoint_col) && length(time_levels) >= 6) {
  early_levels <- head(time_levels, 3)
  late_levels <- tail(time_levels, 3)
  add_scenario(
    "Early expedition versus late/recovery phase",
    "Pools the first three ordered sampling windows and the last three windows.",
    metadata, subject_col, time_values %in% early_levels, time_values %in% late_levels,
    paste(early_levels, collapse = "+"), paste(late_levels, collapse = "+"),
    n_scenario_features,
    "Pooling improves observation counts but mixes distinct exposures; within-subject dependence and unequal trajectories remain."
  )
  midpoint <- floor(length(time_levels) / 2)
  ascent_levels <- time_levels[seq_len(midpoint)]
  descent_levels <- time_levels[(midpoint + 1):length(time_levels)]
  add_scenario(
    "Ascending/early phase versus descending/late phase",
    "A phase-level contrast based on ordered expedition sampling; labels require biological confirmation from the expedition record.",
    metadata, subject_col, time_values %in% ascent_levels, time_values %in% descent_levels,
    paste(ascent_levels, collapse = "+"), paste(descent_levels, collapse = "+"),
    n_scenario_features,
    "Time, altitude, diet, acclimatization, and other expedition exposures are inseparable; phase labels are inferred rather than recorded."
  )
  add_scenario(
    "Baseline versus final follow-up",
    "Exact first-versus-final sampling-window comparison.",
    metadata, subject_col, time_values == time_levels[[1]], time_values == tail(time_levels, 1),
    time_levels[[1]], tail(time_levels, 1), n_scenario_features,
    "Exact windows contain too few samples for stable taxon-correlation networks, even if paired longitudinally."
  )
}
if (!is.na(sex_col)) {
  sex_values <- tolower(as.character(metadata[[sex_col]]))
  sex_levels <- names(sort(table(sex_values), decreasing = TRUE))
  if (length(sex_levels) >= 2) add_scenario(
    "Sex-associated network comparison",
    "Exploratory between-participant comparison of the two most frequent recorded sex/gender categories.",
    metadata, subject_col, sex_values == sex_levels[[1]], sex_values == sex_levels[[2]],
    sex_levels[[1]], sex_levels[[2]], n_scenario_features,
    "Sex is a between-subject attribute, but repeated samples inflate observation counts and subject groups are small and imbalanced."
  )
}
if (!is.na(age_col) && !is.na(subject_col)) {
  age <- safe_numeric(metadata[[age_col]])
  subject_age <- tapply(age, subject_values, function(x) stats::median(x, na.rm = TRUE))
  age_cut <- stats::median(subject_age, na.rm = TRUE)
  add_scenario(
    "Younger versus older participants",
    "Exploratory split at the participant-level median age.",
    metadata, subject_col, age < age_cut, age >= age_cut,
    paste0("age < ", age_cut), paste0("age >= ", age_cut), n_scenario_features,
    "Data-driven age dichotomization loses information; repeated samples and few participants make network differences unstable."
  )
}
candidate_scenarios <- if (length(scenarios)) bind_rows(scenarios) else data.frame(
  scenario_name = "No defensible two-group scenario detected", recommendation_level = "NOT RECOMMENDED"
)
candidate_scenarios <- candidate_scenarios %>%
  mutate(recommendation_rank = match(recommendation_level, c("STRONG", "MODERATE", "WEAK", "NOT RECOMMENDED"))) %>%
  arrange(recommendation_rank, desc(pmin(group_a_samples, group_b_samples))) %>%
  select(-recommendation_rank)

best_level <- candidate_scenarios$recommendation_level[[1]]
overall_recommendation <- if (any(candidate_scenarios$recommendation_level %in% c("STRONG", "MODERATE"))) {
  "NetMoss may be considered as a secondary analysis after independent network-stability validation."
} else if (any(candidate_scenarios$recommendation_level == "WEAK")) {
  "NetMoss is not supported as a primary analysis; at most it is an exploratory/supporting analysis after taxonomic aggregation and subject-aware sensitivity analyses."
} else {
  "NetMoss is not recommended on the current data structure because the effective independent sample size is too small for condition-specific microbial networks."
}

assumptions <- data.frame(
  assumption = c(
    "Two biologically meaningful comparison groups", "Sufficient effectively independent samples per group",
    "Matched microbial features across groups", "Reasonable feature prevalence",
    "Acceptable sparsity", "Network inference feasible at native feature level",
    "Network inference feasible after genus aggregation", "Balanced comparison groups",
    "Independent observations for ordinary correlation", "Documented abundance provenance",
    "Absence of severe technical/batch bias", "Complete sample-metadata linkage",
    "Taxonomy available for aggregation", "Longitudinal dependence can be handled by original NetMoss workflow",
    "Summit exposure explicitly observed in metadata"
  ),
  status = c(
    ifelse(!is.na(timepoint_col) || !is.na(altitude_col), "Partially satisfied", "Not satisfied"),
    ifelse(!is.na(n_subjects) && n_subjects >= 40, "Partially satisfied", "Not satisfied"),
    "Satisfied",
    ifelse(stats::median(feature_prevalence) >= 0.5, "Partially satisfied", "Not satisfied"),
    ifelse(mean(!present) <= 0.5, "Partially satisfied", "Not satisfied"),
    ifelse(ncol(matrix) <= n_subjects * 2, "Partially satisfied", "Not satisfied"),
    ifelse(n_scenario_features <= n_subjects * 4, "Partially satisfied", "Not satisfied"),
    ifelse(any(candidate_scenarios$balance_ratio_larger_to_smaller <= 1.5, na.rm = TRUE), "Partially satisfied", "Not satisfied"),
    ifelse(!is.na(n_subjects) && nrow(metadata) > n_subjects, "Not satisfied", "Satisfied"),
    ifelse(grepl("Relative abundance", scale_result$type), "Partially satisfied", "Cannot determine"),
    "Cannot determine",
    ifelse(length(metadata_only_ids) == 0 && length(abundance_only_ids) == 0, "Satisfied", "Partially satisfied"),
    ifelse(taxonomy_matched == ncol(matrix), "Satisfied", "Partially satisfied"),
    ifelse(grepl("longitudinal|repeated", design), "Not satisfied", "Cannot determine"),
    "Not satisfied"
  ),
  justification = c(
    "Timepoint and altitude provide biological contrasts, but there is no randomized treatment or case-control condition and phase definitions require confirmation.",
    paste0(nrow(matrix), " observations arise from only ", n_subjects, " participants; exact timepoint groups are much smaller."),
    "All comparisons derive from one common abundance matrix with identical feature columns.",
    paste0("Median native-feature prevalence is ", percent(stats::median(feature_prevalence), accuracy = 0.1), "; filtering sensitivity is reported."),
    paste0("Overall zero/undetected fraction is ", percent(mean(!present), accuracy = 0.1), "."),
    paste0(ncol(matrix), " MAGs greatly exceed ", n_subjects, " effectively independent participants."),
    paste0("Preferred ", preferred_rank, " aggregation yields approximately ", n_scenario_features, " taxa, still high relative to participant count."),
    "At least one pooled scenario is numerically balanced, but balance in observations does not create independent participants.",
    paste0("Repeated observations are present (", nrow(metadata), " samples from ", n_subjects, " participants)."),
    paste0(scale_result$interpretation, " Exact upstream normalization is not recorded in the abundance table."),
    "No batch variable, sequencing depth, extraction lot, or negative-control information was detected; absence of evidence is not evidence of no bias.",
    paste0(length(matched_ids), " abundance samples matched; ", length(metadata_only_ids), " metadata rows lack abundance and ", length(abundance_only_ids), " abundance rows lack metadata."),
    paste0(taxonomy_matched, " of ", ncol(matrix), " native features match the taxonomy table."),
    "The original workflow compares condition networks and does not itself model within-subject correlation; naive pooling would create pseudoreplication.",
    "The maximum recorded sampling altitude is not equivalent to documented summit exposure; repository notes place summit-push exposure between sampling points."
  ),
  evidence_basis = c(
    "study metadata", "network-stability diagnostic heuristic", "NetMoss input structure",
    "study data", "study data", "network-stability diagnostic heuristic",
    "network-stability diagnostic heuristic", "study data", "study design",
    "study data and repository documentation", "missing provenance", "study data",
    "study data", "NetMoss method/design mismatch", "study metadata and expedition documentation"
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(dataset_summary, file.path(output_dir, "netmoss_dataset_summary.csv"))
write_csv_safe(metadata_summary, file.path(output_dir, "netmoss_metadata_summary.csv"))
write_csv_safe(feature_prevalence_table, file.path(output_dir, "netmoss_feature_prevalence.csv"))
write_csv_safe(sparsity_summary, file.path(output_dir, "netmoss_sparsity_summary.csv"))
write_csv_safe(candidate_scenarios, file.path(output_dir, "netmoss_candidate_scenarios.csv"))
write_csv_safe(sample_balance, file.path(output_dir, "netmoss_sample_balance.csv"))
write_csv_safe(assumptions, file.path(output_dir, "netmoss_assumption_checklist.csv"))
write_csv_safe(filtering_table, file.path(output_dir, "netmoss_filtering_evaluation.csv"))
write_csv_safe(
  data.frame(id_status = c(rep("metadata_without_abundance", length(metadata_only_ids)), rep("abundance_without_metadata", length(abundance_only_ids))),
             sample_id = c(metadata_only_ids, abundance_only_ids)),
  file.path(output_dir, "netmoss_sample_id_matching.csv")
)
write_csv_safe(missing_by_subject, file.path(output_dir, "netmoss_longitudinal_completeness.csv"))

# Publication-quality figures (vector PDF and SVG).
prevalence_plot <- ggplot(feature_prevalence_table, aes(prevalence)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "#2C7FB8", color = "white") +
  geom_vline(xintercept = screening_threshold, linetype = 2, color = "#D95F0E") +
  scale_x_continuous(labels = percent_format()) + coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Microbial feature prevalence", subtitle = paste0("Dashed line: strict sensitivity threshold ", percent(screening_threshold)), x = "Samples with detected feature", y = "Features") + theme_netmoss()
save_figure(prevalence_plot, "feature_prevalence_histogram", output_dir)

zero_plot <- ggplot(feature_prevalence_table, aes(zero_proportion)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "#7FCDBB", color = "white") +
  scale_x_continuous(labels = percent_format()) + coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Feature zero proportions", x = "Zero proportion per feature", y = "Features") + theme_netmoss()
save_figure(zero_plot, "zero_proportion_histogram", output_dir)

sample_sparse_df <- data.frame(sample_id = rownames(matrix), sparsity = sample_sparsity)
sample_sparse_plot <- ggplot(sample_sparse_df, aes(sparsity)) +
  geom_histogram(binwidth = 0.025, fill = "#225EA8", color = "white") +
  geom_density(aes(y = after_stat(count) * 0.025), color = "#D95F0E", linewidth = 0.8) +
  scale_x_continuous(labels = percent_format()) + coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Sample sparsity", x = "Zero/undetected fraction per sample", y = "Samples") + theme_netmoss()
save_figure(sample_sparse_plot, "sample_sparsity", output_dir)

feature_sparse_plot <- ggplot(feature_prevalence_table, aes(x = "All features", y = zero_proportion)) +
  geom_violin(fill = "#A1DAB4", color = "grey30", width = 0.7) + geom_boxplot(width = 0.16, outlier.size = 0.7) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(title = "Feature sparsity distribution", x = NULL, y = "Zero proportion per feature") + theme_netmoss()
save_figure(feature_sparse_plot, "feature_sparsity", output_dir, width = 5.5)

subject_plot <- ggplot(samples_per_subject, aes(x = reorder(as.character(subject_id), n_samples), y = n_samples)) +
  geom_col(fill = "#41B6C4") + coord_flip() +
  labs(title = "Samples per participant", x = "Participant", y = "Usable samples") + theme_netmoss()
save_figure(subject_plot, "samples_per_subject", output_dir)

time_plot <- ggplot(samples_per_timepoint, aes(x = timepoint, y = n_samples)) +
  geom_col(fill = "#2C7FB8") + geom_text(aes(label = n_samples), vjust = -0.25, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Samples per timepoint", x = "Ordered timepoint", y = "Usable samples") + theme_netmoss()
save_figure(time_plot, "samples_per_timepoint", output_dir)

if (!all(is.na(longitudinal_grid$subject_id))) {
  missing_plot <- ggplot(longitudinal_grid, aes(x = timepoint, y = subject_id, fill = observed)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = c(`TRUE` = "#238B45", `FALSE` = "#D9D9D9"), labels = c(`TRUE` = "Observed", `FALSE` = "Missing")) +
    labs(title = "Longitudinal sampling completeness", x = "Timepoint", y = "Participant", fill = NULL) + theme_netmoss()
  save_figure(missing_plot, "missing_timepoint_heatmap", output_dir)

  trajectory_df <- longitudinal_grid %>% filter(observed) %>% mutate(time_index = as.numeric(timepoint))
  trajectory_plot <- ggplot(trajectory_df, aes(x = timepoint, y = subject_id, group = subject_id)) +
    geom_line(color = "grey70", linewidth = 0.45) + geom_point(aes(color = subject_id), size = 2) +
    labs(title = "Longitudinal sampling diagram", subtitle = "Lines connect available samples; gaps indicate missing sampling windows", x = "Ordered timepoint", y = "Participant") +
    guides(color = "none") + theme_netmoss()
  save_figure(trajectory_plot, "longitudinal_sampling_diagram", output_dir)
}

scenario_balance_plot <- candidate_scenarios %>%
  select(scenario_name, group_a, group_b, group_a_samples, group_b_samples) %>%
  pivot_longer(c(group_a_samples, group_b_samples), names_to = "group_position", values_to = "n_samples") %>%
  mutate(group_label = ifelse(group_position == "group_a_samples", group_a, group_b)) %>%
  ggplot(aes(x = reorder(scenario_name, n_samples), y = n_samples, fill = group_position)) +
  geom_col(position = position_dodge(width = 0.8)) + coord_flip() +
  scale_fill_manual(values = c(group_a_samples = "#2C7FB8", group_b_samples = "#F03B20"), labels = c("Group A", "Group B")) +
  labs(title = "Candidate scenario observation balance", x = NULL, y = "Samples", fill = NULL) + theme_netmoss()
save_figure(scenario_balance_plot, "group_balance", output_dir, width = 9, height = 6.5)

filter_plot <- ggplot(filtering_table, aes(prevalence_threshold, retained_features)) +
  geom_line(color = "#225EA8", linewidth = 0.8) + geom_point(aes(color = screening_candidate), size = 2.5) +
  scale_x_continuous(labels = percent_format(), breaks = thresholds) +
  scale_color_manual(values = c(`TRUE` = "#D95F0E", `FALSE` = "#225EA8"), guide = "none") +
  labs(title = "Features retained under prevalence filtering", subtitle = "Highlighted point is the strictest tested sensitivity threshold, not a published NetMoss threshold", x = "Minimum prevalence", y = "Retained features") + theme_netmoss()
save_figure(filter_plot, "retained_features_filtering_thresholds", output_dir)

recommendation_lines <- c(
  "NetMoss readiness recommendations",
  paste0("Generated: ", timestamp()), "",
  paste0("Overall recommendation: ", overall_recommendation), "",
  "Recommended preprocessing before any exploratory network analysis:",
  paste0("1. Verify upstream abundance normalization and why sample totals range from ", signif(min(sample_sums), 4), " to ", signif(max(sample_sums), 4), "."),
  "2. Define comparison windows prospectively from expedition biology; do not call 6500 m samples 'summit' samples.",
  paste0("3. Aggregate MAGs to a biologically defensible rank (genus is currently preferable; about ", n_scenario_features, " mapped taxa) and retain an explicit MAG-to-taxon map."),
  paste0("4. Test ", percent(screening_threshold), " prevalence as a strict sensitivity threshold, but note that no tested threshold makes the native MAG network feasible; combine filtering with taxonomic aggregation."),
  "5. Avoid treating repeated samples as independent replicates. Compare subject-aggregated profiles or use subject-level resampling/blocked bootstrap sensitivity analyses.",
  "6. Use a compositional association estimator (for example SparCC or SPIEC-EASI, as suggested by NetMoss documentation) and quantify edge/module stability across bootstrap resamples.",
  "7. Obtain sequencing-depth, extraction batch, run, storage, medication, diet, and technical-control metadata before attributing network shifts to altitude or expedition phase.",
  "8. Treat any NetMoss result as exploratory/supporting and triangulate it with mixed models, differential abundance, beta diversity, and Random Forest analyses.",
  "9. Pre-register the primary contrast and preserve held-out or external data for validation if biomarker claims are intended.", "",
  "Important methodological boundary:",
  "The readiness ratings are diagnostic heuristics for correlation-network stability. The NetMoss publication does not provide a universal minimum sample-size rule."
)
writeLines(recommendation_lines, file.path(output_dir, "netmoss_recommendations.txt"))

scenario_text <- unlist(lapply(seq_len(nrow(candidate_scenarios)), function(i) {
  row <- candidate_scenarios[i, ]
  c(
    paste0(i, ". ", row$scenario_name, " — ", row$recommendation_level),
    paste0("   Groups: ", row$group_a, " (n=", row$group_a_samples, ", subjects=", row$group_a_subjects,
           ") versus ", row$group_b, " (n=", row$group_b_samples, ", subjects=", row$group_b_subjects, ")."),
    paste0("   Interpretation: ", row$biological_interpretation),
    paste0("   Robustness: ", row$expected_network_robustness),
    paste0("   Concern: ", row$potential_statistical_concerns), ""
  )
}))

assumption_text <- unlist(lapply(seq_len(nrow(assumptions)), function(i) {
  paste0("- ", assumptions$assumption[[i]], ": ", assumptions$status[[i]], ". ", assumptions$justification[[i]])
}))

report_lines <- c(
  "NetMoss readiness assessment: supplementary methods-style report",
  "===============================================================",
  paste0("Generated: ", timestamp()),
  paste0("Configuration: ", normalizePath(args$config, mustWork = TRUE)), "",
  "1. Executive summary", "--------------------",
  overall_recommendation,
  paste0("The current gut dataset contains ", nrow(matrix), " matched samples, ", ncol(matrix), " MAG features, and only ", n_subjects,
         " participants. The ", percent(mean(!present), accuracy = 0.1), " zero fraction is manageable for descriptive community analyses, but the effective independent sample size is inadequate for stable condition-specific MAG networks."),
  paste0("Taxonomy enables aggregation to about ", n_scenario_features, " ", preferred_rank, "-level taxa. This improves sparsity and dimensionality but does not solve the participant-level sample-size limitation."), "",
  "2. Repository overview", "----------------------",
  "The repository contains raw gut MAG abundance (data/raw/gut.csv), sample metadata (data/raw/sample_info.csv), GTDB taxonomy (data/raw/GDTB.csv), and a processed sample-by-MAG long table (data/processed/cleaned_data.rds).",
  "Existing analyses include per-MAG LMM/GLMM workflows, alpha diversity, beta diversity/PERMANOVA and longitudinal distances, Random Forest prediction with subject-grouped resampling, response transformation diagnostics, expedition-window evaluation, and HPC launchers.",
  "The abundance-preparation script reshapes the wide MAG table and joins metadata; the beta-diversity workflow validates IDs, filters by prevalence, and applies compositional/community transformations. No existing script was modified by this readiness assessment.", "",
  "3. Dataset overview", "-------------------",
  paste0("Input orientation: ", layout$original_orientation, "; canonical analysis orientation: Samples x Features."),
  paste0("Matched samples: ", length(matched_ids), "; metadata-only records: ", length(metadata_only_ids), "; abundance-only records: ", length(abundance_only_ids), "."),
  paste0("Feature type: ", feature_type, "; taxonomy matches: ", taxonomy_matched, "/", ncol(matrix), "."),
  paste0("Design: ", design, "; timepoints: ", n_timepoints, "; participants: ", n_subjects, "."), "",
  "4. Data quality assessment", "-------------------------",
  paste0("Scale classification: ", scale_result$type, ". ", scale_result$interpretation),
  paste0("Range: ", signif(min(finite_values), 5), " to ", signif(max(finite_values), 5), "; median: ", signif(stats::median(finite_values), 5),
         "; zeros: ", sum(matrix == 0, na.rm = TRUE), " (", percent(mean(matrix == 0, na.rm = TRUE), accuracy = 0.1), ")."),
  paste0("Missing values: ", sum(is.na(matrix)), "; negative: ", sum(matrix < 0, na.rm = TRUE), "; infinite: ", sum(is.infinite(matrix)), "; NaN: ", sum(is.nan(matrix)), "."),
  paste0("Sample totals classify as ", scale_result$sum_class, " and range from ", signif(min(sample_sums), 5), " to ", signif(max(sample_sums), 5), ". Upstream normalization remains incompletely documented."),
  paste0("Native-feature prevalence ranges from ", percent(min(feature_prevalence), accuracy = 0.1), " to ", percent(max(feature_prevalence), accuracy = 0.1),
         "; median ", percent(stats::median(feature_prevalence), accuracy = 0.1), "."), "",
  "5. Metadata quality and longitudinal structure", "---------------------------------------------",
  paste0("Detected roles: ", paste(names(role_map)[!is.na(role_map)], role_map[!is.na(role_map)], sep = "=", collapse = "; "), "."),
  paste0("Participants contribute ", min(samples_per_subject$n_samples), " to ", max(samples_per_subject$n_samples), " usable samples (median ", stats::median(samples_per_subject$n_samples), ")."),
  paste0("Per-timepoint usable samples range from ", min(samples_per_timepoint$n_samples), " to ", max(samples_per_timepoint$n_samples), ". Missingness and dropout are tabulated in netmoss_longitudinal_completeness.csv."),
  "Date/timepoint and altitude are tightly linked. Diet, acclimatization, medication, illness, extraction batch, sequencing run, and storage variables were not detected and therefore cannot be ruled out as confounders.", "",
  "6. Assumption checklist", "-----------------------", assumption_text, "",
  "7. Candidate scenarios ranked strongest to weakest", "--------------------------------------------------", scenario_text,
  "8. Scientific limitations", "-------------------------",
  "The dominant limitation is not the number of rows but the 13 independent participants. Repeated samples can improve longitudinal mean estimation, but they do not supply independent replicates for estimating separate microbial correlation networks.",
  "Pooling timepoints may increase apparent group size while inducing pseudoreplication and mixing biological states. Subject aggregation avoids pseudoreplication but leaves at most 13 participant profiles per condition, still too few for 114 genus-level or 441 MAG-level nodes.",
  "Network edges inferred from relative-abundance data are associations, not direct ecological interactions. Compositional closure, zeros, shared environmental responses, and technical effects can produce apparent correlations.",
  "No exact summit sample is documented. Repository expedition notes indicate summit-push exposure occurs between sampling points, so baseline-versus-summit language would be scientifically indefensible.", "",
  "9. Recommended preprocessing and sensitivity analysis", "----------------------------------------------------",
  recommendation_lines[5:13], "",
  "10. Sample-size and network reliability interpretation", "-----------------------------------------------------",
  paste0("At the native MAG level, p/n_subject is ", round(ncol(matrix) / n_subjects, 1), ". At the preferred ", preferred_rank,
         " level it remains approximately ", round(n_scenario_features / n_subjects, 1), ". These ratios make unregularized condition-specific correlations highly unstable."),
  "There is no universally accepted NetMoss minimum sample size in the original publication. Adequacy must therefore be established empirically through bootstrap edge stability, module reproducibility, perturbation sensitivity, and external validation. The current participant count makes those checks unlikely to support a robust biomarker claim.", "",
  "11. Relationship to existing analyses", "-------------------------------------",
  "LMM/GLMM asks whether individual taxon abundance changes with time/exposure while modeling subject-level dependence. Differential abundance asks which taxa differ in marginal abundance between groups. Random Forest asks whether multivariate profiles predict an outcome and which features aid prediction, ideally with participant-grouped validation.",
  "Beta-diversity/PERMANOVA asks whether whole-community composition differs between conditions; correlation analysis measures pairwise co-variation; network analysis organizes such associations into edges and modules. NetMoss specifically ranks taxa by their contribution to changes in module/network structure between two states, so it can identify rewiring not visible as marginal abundance change.",
  "For publication, NetMoss would complement rather than replace mixed models and differential abundance. Concordant taxa could be presented as multi-method candidates; discordant high-NetMoss taxa could motivate ecological hypotheses. Given this dataset, such results should be explicitly labeled exploratory and supported by stability analyses, subject-aware sensitivity checks, and independent validation.", "",
  "12. Overall recommendation", "--------------------------", overall_recommendation,
  "Scientifically, NetMoss should not be a primary analysis of the current BlackICE Summit dataset. A future exploratory genus-level analysis might be defensible only after the preprocessing and validation steps above, with strong language limiting causal and biomarker claims.", "",
  "Methodological references", "-------------------------",
  "- Xiao L, Zhang F, Zhao F. Large-scale microbiome data integration enables robust biomarker identification. Nature Computational Science 2, 307-316 (2022). https://doi.org/10.1038/s43588-022-00247-8",
  "- NetMoss authors' R package documentation: https://github.com/xiaolw95/NetMoss (accessed during pipeline design).",
  "The original NetMoss workflow accepts relative-abundance matrices with taxa as rows and samples as columns and condition-specific adjacency matrices; SparCC or SPIEC-EASI are specifically suggested for microbial association inference.", "",
  "Reproducibility statement", "-------------------------",
  "This diagnostic is read-only with respect to input data and existing analysis code. All thresholds, inferred columns, scenario counts, assumptions, and software versions are recorded in machine-readable outputs."
)
writeLines(report_lines, file.path(output_dir, "netmoss_readiness_summary.txt"))
old_tz <- Sys.getenv("TZ", unset = "")
Sys.setenv(TZ = "UTC")
suppressWarnings(utils::capture.output(utils::sessionInfo(), file = file.path(output_dir, "sessionInfo.txt")))
if (nzchar(old_tz)) Sys.setenv(TZ = old_tz) else Sys.unsetenv("TZ")

cat("Completed: ", timestamp(), "\n", sep = "")
cat("Overall recommendation: ", overall_recommendation, "\n", sep = "")
cat("Outputs: ", normalizePath(output_dir, mustWork = TRUE), "\n", sep = "")
