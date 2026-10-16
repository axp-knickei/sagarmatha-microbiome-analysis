# Prepare validated alpha-diversity inputs from the gut MAG abundance table.

required_packages <- c("dplyr", "tidyr", "lubridate", "yaml", "here")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

library(dplyr)
library(tidyr)
library(lubridate)
library(yaml)
library(here)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

env_value <- function(...) {
  for (key in c(...)) {
    value <- Sys.getenv(key, unset = NA_character_)
    if (!is.na(value) && nzchar(value)) return(value)
  }
  NULL
}

is_absolute_path <- function(path) {
  grepl("^/", path) || grepl("^[A-Za-z]:[/\\\\]", path)
}

translate_windows_path <- function(path) {
  if (.Platform$OS.type != "unix" || !grepl("^[A-Za-z]:[/\\\\]", path)) {
    return(path)
  }

  drive <- tolower(substr(path, 1, 1))
  rest <- substring(path, 3)
  rest <- gsub("\\\\", "/", rest)
  file.path("/mnt", drive, rest)
}

resolve_path <- function(path) {
  path <- translate_windows_path(path)
  if (is_absolute_path(path)) path else here::here(path)
}

read_csv_file <- function(path) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN"))
}

write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
}

is_count_like_matrix <- function(x) {
  values <- as.numeric(x)
  values <- values[!is.na(values)]
  if (length(values) == 0 || any(values < 0)) return(FALSE)
  all(abs(values - round(values)) < 1e-8) && max(values) > 1
}

config <- yaml::read_yaml(here::here("config", "config.yml"))
raw_dir_override <- env_value("BLACKICE_RAW_DIR", "RAW_DATA_DIR")
raw_dir <- raw_dir_override %||% config$data$raw_dir %||% "data/raw"
abundance_file <- env_value("BLACKICE_ABUNDANCE_FILE", "ABUNDANCE_FILE") %||% config$data$abundance %||% "gut.csv"
metadata_file <- env_value("BLACKICE_METADATA_FILE", "METADATA_FILE") %||% config$data$metadata %||% "sample_info.csv"
taxonomy_file <- env_value("BLACKICE_TAXONOMY_FILE", "TAXONOMY_FILE") %||% "GDTB.csv"

raw_dir_candidates <- raw_dir
if (is.null(raw_dir_override)) {
  raw_dir_candidates <- unique(c(raw_dir_candidates, "data/raw"))
}

for (candidate_raw_dir in raw_dir_candidates) {
  candidate_path <- resolve_path(candidate_raw_dir)
  if (
    file.exists(file.path(candidate_path, abundance_file)) &&
      file.exists(file.path(candidate_path, metadata_file))
  ) {
    if (!identical(candidate_raw_dir, raw_dir)) {
      message("Configured raw data directory was not usable; using project data/raw instead: ", candidate_path)
    }
    raw_dir <- candidate_raw_dir
    break
  }
}

raw_dir_path <- resolve_path(raw_dir)
abundance_path <- file.path(raw_dir_path, abundance_file)
metadata_path <- file.path(raw_dir_path, metadata_file)
taxonomy_path <- file.path(raw_dir_path, taxonomy_file)
abundance_relpath <- file.path(raw_dir, abundance_file)
metadata_relpath <- file.path(raw_dir, metadata_file)
taxonomy_relpath <- file.path(raw_dir, taxonomy_file)

output_dir <- env_value("BLACKICE_ALPHA_OUTPUT_DIR") %||% "results/alpha_diversity"
tables_dir <- file.path(output_dir, "tables")
logs_dir <- file.path(output_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("docs/reports", recursive = TRUE, showWarnings = FALSE)

if (!file.exists(abundance_path)) stop("Abundance file not found: ", abundance_path, call. = FALSE)
if (!file.exists(metadata_path)) stop("Metadata file not found: ", metadata_path, call. = FALSE)

abundance <- read_csv_file(abundance_path)
metadata <- read_csv_file(metadata_path)
taxonomy <- if (file.exists(taxonomy_path)) read_csv_file(taxonomy_path) else NULL

if (!"sampleID" %in% names(abundance)) stop("Abundance table must contain sampleID.", call. = FALSE)
if (!"sampleID" %in% names(metadata)) stop("Metadata table must contain sampleID.", call. = FALSE)

abundance$sampleID <- as.character(abundance$sampleID)
metadata$sampleID <- as.character(metadata$sampleID)
feature_cols <- setdiff(names(abundance), "sampleID")
feature_data <- abundance[, feature_cols, drop = FALSE]

converted <- lapply(feature_data, function(x) suppressWarnings(as.numeric(as.character(x))))
non_numeric_cols <- names(feature_data)[vapply(
  seq_along(converted),
  function(i) {
    original <- trimws(as.character(feature_data[[i]]))
    any(is.na(converted[[i]]) & nzchar(original) & !is.na(original))
  },
  logical(1)
)]
feature_data <- as.data.frame(converted, check.names = FALSE)
value_matrix <- as.matrix(feature_data)
storage.mode(value_matrix) <- "numeric"
rownames(value_matrix) <- abundance$sampleID

missing_sample_ids <- sum(is.na(abundance$sampleID) | !nzchar(abundance$sampleID))
duplicate_sample_ids <- unique(abundance$sampleID[duplicated(abundance$sampleID)])
duplicate_metadata_ids <- unique(metadata$sampleID[duplicated(metadata$sampleID)])
duplicate_mag_ids <- unique(feature_cols[duplicated(feature_cols)])
negative_values <- sum(value_matrix < 0, na.rm = TRUE)
all_zero_samples <- rownames(value_matrix)[rowSums(value_matrix, na.rm = TRUE) == 0]
all_zero_features <- colnames(value_matrix)[colSums(value_matrix, na.rm = TRUE) == 0]

if (length(duplicate_sample_ids) > 0) stop("Duplicate sampleID values in abundance table: ", paste(head(duplicate_sample_ids, 10), collapse = ", "), call. = FALSE)
if (length(duplicate_metadata_ids) > 0) stop("Duplicate sampleID values in metadata table: ", paste(head(duplicate_metadata_ids, 10), collapse = ", "), call. = FALSE)
if (length(duplicate_mag_ids) > 0) stop("Duplicate MAG IDs in abundance columns: ", paste(head(duplicate_mag_ids, 10), collapse = ", "), call. = FALSE)
if (length(non_numeric_cols) > 0) stop("Non-numeric abundance columns detected: ", paste(head(non_numeric_cols, 10), collapse = ", "), call. = FALSE)
if (negative_values > 0) stop("Negative abundance values detected: ", negative_values, call. = FALSE)

metadata_only <- setdiff(metadata$sampleID, abundance$sampleID)
abundance_only <- setdiff(abundance$sampleID, metadata$sampleID)
matched_ids <- intersect(abundance$sampleID, metadata$sampleID)

metadata_prepared <- metadata %>%
  mutate(
    subjectID = if ("subjectID" %in% names(.)) as.character(subjectID) else as.character(group),
    date_converted = if ("date" %in% names(.)) lubridate::ymd(as.character(date)) else as.Date(NA),
    Timepoint = if ("Timepoint" %in% names(.)) as.character(Timepoint) else NA_character_,
    altitude = if ("altitude" %in% names(.)) suppressWarnings(as.numeric(altitude)) else NA_real_,
    age = if ("age" %in% names(.)) suppressWarnings(as.numeric(age)) else NA_real_,
    gender = if ("gender" %in% names(.)) as.character(gender) else NA_character_
  ) %>%
  group_by(subjectID) %>%
  mutate(days_from_baseline = as.numeric(date_converted - min(date_converted, na.rm = TRUE))) %>%
  ungroup()

if (all(is.infinite(metadata_prepared$days_from_baseline) | is.na(metadata_prepared$days_from_baseline))) {
  metadata_prepared$days_from_baseline <- NA_real_
}

keep_samples <- matched_ids
if (length(all_zero_samples) > 0) keep_samples <- setdiff(keep_samples, all_zero_samples)
if (length(all_zero_features) > 0) value_matrix <- value_matrix[, setdiff(colnames(value_matrix), all_zero_features), drop = FALSE]
value_matrix <- value_matrix[keep_samples, , drop = FALSE]

nonmissing_values <- as.numeric(value_matrix)
nonmissing_values <- nonmissing_values[!is.na(nonmissing_values)]
row_sums <- rowSums(value_matrix, na.rm = TRUE)
count_like <- is_count_like_matrix(value_matrix)
relative_like <- !count_like &&
  length(nonmissing_values) > 0 &&
  min(nonmissing_values, na.rm = TRUE) >= 0 &&
  max(nonmissing_values, na.rm = TRUE) <= 1 &&
  stats::median(row_sums, na.rm = TRUE) > 0.5 &&
  stats::median(row_sums, na.rm = TRUE) < 1.5
scale_decision <- if (count_like) "count_like" else if (relative_like) "relative_abundance_like" else "continuous_non_count"

prepared <- list(
  abundance_matrix = value_matrix,
  metadata = metadata_prepared %>% filter(sampleID %in% rownames(value_matrix)),
  diagnostics = list(
    abundance_file = abundance_path,
    metadata_file = metadata_path,
    taxonomy_file = if (file.exists(taxonomy_path)) taxonomy_path else NA_character_,
    abundance_relpath = abundance_relpath,
    metadata_relpath = metadata_relpath,
    taxonomy_relpath = if (file.exists(taxonomy_path)) taxonomy_relpath else NA_character_,
    abundance_orientation = "samples_as_rows_MAGs_as_columns",
    n_samples_raw = nrow(abundance),
    n_mag_features_raw = length(feature_cols),
    n_samples_prepared = nrow(value_matrix),
    n_mag_features_prepared = ncol(value_matrix),
    n_metadata_rows = nrow(metadata),
    n_samples_matched = length(matched_ids),
    abundance_only_sample_ids = abundance_only,
    metadata_only_sample_ids = metadata_only,
    all_zero_samples = all_zero_samples,
    all_zero_features = all_zero_features,
    abundance_scale_decision = scale_decision,
    chao1_ace_decision = if (count_like) "calculated" else "not_calculated_relative_or_non_count_input"
  )
)

saveRDS(prepared, file.path(tables_dir, "alpha_diversity_prepared_inputs.rds"))
utils::write.csv(cbind(sampleID = rownames(value_matrix), value_matrix), file.path(tables_dir, "alpha_diversity_abundance_matrix.csv"), row.names = FALSE)
utils::write.csv(prepared$metadata, file.path(tables_dir, "alpha_diversity_metadata_prepared.csv"), row.names = FALSE)

diagnostic_table <- data.frame(
  item = c(
    "abundance_file", "metadata_file", "taxonomy_file", "abundance_orientation",
    "samples_are_rows", "MAGs_are_columns", "abundance_scale_decision",
    "raw_samples", "raw_MAG_features", "prepared_samples", "prepared_MAG_features",
    "metadata_rows", "matched_samples", "abundance_only_samples", "metadata_only_samples",
    "missing_sample_ids", "duplicated_abundance_sample_ids", "duplicated_metadata_sample_ids",
    "duplicated_MAG_ids", "negative_values", "all_zero_samples", "all_zero_MAGs",
    "taxonomy_rows", "taxonomy_MAG_matches"
  ),
  value = c(
    abundance_relpath, metadata_relpath, if (file.exists(taxonomy_path)) taxonomy_relpath else "not_found",
    "samples_as_rows_MAGs_as_columns", TRUE, TRUE, scale_decision,
    nrow(abundance), length(feature_cols), nrow(value_matrix), ncol(value_matrix),
    nrow(metadata), length(matched_ids), length(abundance_only), length(metadata_only),
    missing_sample_ids, length(duplicate_sample_ids), length(duplicate_metadata_ids),
    length(duplicate_mag_ids), negative_values, length(all_zero_samples), length(all_zero_features),
    if (is.null(taxonomy)) NA_integer_ else nrow(taxonomy),
    if (is.null(taxonomy) || !"Genome" %in% names(taxonomy)) NA_integer_ else sum(colnames(value_matrix) %in% taxonomy$Genome)
  )
)
utils::write.csv(diagnostic_table, file.path(tables_dir, "alpha_diversity_diagnostic_summary.csv"), row.names = FALSE)

warning_lines <- if (!count_like) {
  c(
    "",
    "WARNING: The abundance matrix is not count-like. Values are continuous and relative-abundance-like.",
    "Chao1 and ACE are not calculated in the primary workflow because those estimators assume count/incidence information."
  )
} else character(0)

summary_lines <- c(
  "# Alpha-Diversity Input Diagnostic Summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("- Gut abundance table: `", abundance_relpath, "`"),
  paste0("- Metadata table: `", metadata_relpath, "`"),
  paste0("- Taxonomy table: `", if (file.exists(taxonomy_path)) taxonomy_relpath else "not found", "`"),
  "- Detected orientation: samples are rows; MAGs/features are columns.",
  paste0("- Raw samples: ", nrow(abundance)),
  paste0("- Raw MAG/features: ", length(feature_cols)),
  paste0("- Prepared samples: ", nrow(value_matrix)),
  paste0("- Prepared MAG/features: ", ncol(value_matrix)),
  paste0("- Matched sample IDs: ", length(matched_ids)),
  paste0("- Abundance-only sample IDs: ", length(abundance_only)),
  paste0("- Metadata-only sample IDs: ", length(metadata_only)),
  paste0("- All-zero samples excluded: ", length(all_zero_samples)),
  paste0("- All-zero MAG/features excluded: ", length(all_zero_features)),
  paste0("- Detected abundance scale: ", scale_decision),
  paste0("- Chao1/ACE decision: ", if (count_like) "calculated" else "not calculated for primary output"),
  warning_lines
)
write_lines(summary_lines, file.path(logs_dir, "alpha_diversity_diagnostic_summary.md"))

cat(paste(summary_lines, collapse = "\n"), "\n")
