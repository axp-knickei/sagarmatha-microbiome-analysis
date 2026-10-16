# ============================================
# 01_data_prep.R - Data Preparation
# ============================================
# Reads metadata plus one wide MAG abundance table, reshapes to model-ready
# long format, validates critical columns, and saves processed data.
#
# Inputs: config/config.yml plus sample_info.csv and the configured abundance
# CSV from data/raw/ or the BLACKICE_RAW_DIR override.
# Outputs: data/processed/cleaned_data.rds, unless config$data$processed changes
# the path.
# Dependencies: dplyr, tidyr, lubridate, yaml, here, and optionally readr.
# Usage: Rscript scripts/01_data_prep.R
# Workflow: Shared preparation stage for the primary and legacy model paths.

library(dplyr)
library(tidyr)
library(lubridate)
library(yaml)
library(here)

has_readr <- requireNamespace("readr", quietly = TRUE)
if (has_readr) {
  library(readr)
} else {
  warning("Package 'readr' is not installed; falling back to utils::read.csv().", call. = FALSE)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

env_value <- function(...) {
  keys <- c(...)
  for (key in keys) {
    value <- Sys.getenv(key, unset = NA_character_)
    if (!is.na(value) && nzchar(value)) {
      return(value)
    }
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
  if (is_absolute_path(path)) path else here(path)
}

count_missing <- function(df, columns) {
  available_columns <- intersect(columns, names(df))
  missing_columns <- setdiff(columns, names(df))

  counts <- vapply(
    available_columns,
    function(column) sum(is.na(df[[column]])),
    integer(1)
  )

  if (length(missing_columns) > 0) {
    counts <- c(counts, stats::setNames(rep(NA_integer_, length(missing_columns)), missing_columns))
  }

  counts
}

print_missing_counts <- function(label, df, columns) {
  cat(label, "\n", sep = "")
  counts <- count_missing(df, columns)
  for (column in names(counts)) {
    cat(sprintf("  %s: %s\n", column, counts[[column]]))
  }
}

read_csv_file <- function(path) {
  if (has_readr) {
    return(readr::read_csv(path, show_col_types = FALSE))
  }

  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

stop_with_examples <- function(message, values, max_examples = 10) {
  examples <- utils::head(values, max_examples)
  stop(
    message,
    paste(examples, collapse = ", "),
    if (length(values) > max_examples) paste0(" ... plus ", length(values) - max_examples, " more") else "",
    call. = FALSE
  )
}

check_duplicate_sample_ids <- function(values, label) {
  values <- as.character(values)
  duplicated_values <- unique(values[!is.na(values) & duplicated(values)])
  if (length(duplicated_values) > 0) {
    stop_with_examples(
      paste0(label, " contains duplicate sampleID values: "),
      duplicated_values
    )
  }
}

config <- yaml::read_yaml(here("config", "config.yml"))

# The raw dataset is stored outside the repository. Prefer config$data$raw_dir.
# Older config files may still use config$data$input, so keep that as a fallback.
configured_raw_dir <- config$data$raw_dir %||% dirname(config$data$input %||% "")
raw_dir_override <- env_value("BLACKICE_RAW_DIR", "RAW_DATA_DIR")
raw_dir <- raw_dir_override %||% configured_raw_dir
if (raw_dir %in% c("", ".", "data/raw")) {
  raw_dir <- "data/raw"
}

# Select one abundance table for this preparation run. Use "gut.csv" first,
# then repeat with "oral.csv", "skin.csv", or diff_* files as needed.
abundance_filename <- env_value("BLACKICE_ABUNDANCE_FILE", "ABUNDANCE_FILE") %||%
  config$data$abundance %||%
  config$data$abundance_file %||%
  "gut.csv"
metadata_filename <- env_value("BLACKICE_METADATA_FILE", "METADATA_FILE") %||%
  config$data$metadata %||%
  config$data$metadata_file %||%
  "sample_info.csv"

processed_path <- resolve_path(env_value("BLACKICE_PROCESSED_PATH", "PROCESSED_DATA_PATH") %||% config$data$processed)

raw_dir_candidates <- raw_dir
if (is.null(raw_dir_override)) {
  raw_dir_candidates <- unique(c(raw_dir_candidates, "data/raw"))
}

for (candidate_raw_dir in raw_dir_candidates) {
  candidate_path <- resolve_path(candidate_raw_dir)
  if (
    file.exists(file.path(candidate_path, metadata_filename)) &&
      file.exists(file.path(candidate_path, abundance_filename))
  ) {
    if (!identical(candidate_raw_dir, raw_dir)) {
      cat(
        "Configured raw data directory was not usable; using project data/raw instead:",
        candidate_path,
        "\n"
      )
    }
    raw_dir <- candidate_raw_dir
    break
  }
}

metadata_path <- file.path(resolve_path(raw_dir), metadata_filename)
abundance_path <- file.path(resolve_path(raw_dir), abundance_filename)

cat("Selected abundance file:", abundance_filename, "\n")
cat("Raw data directory:", resolve_path(raw_dir), "\n")
cat("Metadata file:", metadata_path, "\n")
cat("Abundance file:", abundance_path, "\n\n")

if (!file.exists(metadata_path)) {
  stop("Metadata file not found: ", metadata_path, call. = FALSE)
}

if (!file.exists(abundance_path)) {
  stop("Abundance file not found: ", abundance_path, call. = FALSE)
}

# --- Load raw data ---
sample_info <- read_csv_file(metadata_path)
abundance_wide <- read_csv_file(abundance_path)

cat(sprintf("Raw metadata dimensions: %d rows x %d columns\n", nrow(sample_info), ncol(sample_info)))
cat(sprintf("Raw abundance dimensions: %d rows x %d columns\n", nrow(abundance_wide), ncol(abundance_wide)))

required_metadata_columns <- c("sampleID", "group", "Timepoint", "date", "gender", "age", "altitude")
missing_metadata_columns <- setdiff(required_metadata_columns, names(sample_info))
if (length(missing_metadata_columns) > 0) {
  stop(
    "sample_info.csv is missing required columns: ",
    paste(missing_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

if (!"sampleID" %in% names(abundance_wide)) {
  stop("Selected abundance file is missing required column: sampleID", call. = FALSE)
}

sample_info$sampleID <- as.character(sample_info$sampleID)
abundance_wide$sampleID <- as.character(abundance_wide$sampleID)

check_duplicate_sample_ids(sample_info$sampleID, "Metadata file")
check_duplicate_sample_ids(abundance_wide$sampleID, "Abundance file")

metadata_ids <- unique(sample_info$sampleID)
abundance_ids <- unique(abundance_wide$sampleID)
abundance_without_metadata <- setdiff(abundance_ids, metadata_ids)
metadata_without_abundance <- setdiff(metadata_ids, abundance_ids)

if (length(abundance_without_metadata) > 0) {
  stop_with_examples(
    "Abundance file contains sampleID values that are missing from metadata: ",
    abundance_without_metadata
  )
}

if (length(metadata_without_abundance) > 0) {
  warning(
    "Metadata contains sampleID values that are not present in the abundance file: ",
    paste(utils::head(metadata_without_abundance, 10), collapse = ", "),
    if (length(metadata_without_abundance) > 10) paste0(" ... plus ", length(metadata_without_abundance) - 10, " more") else "",
    call. = FALSE
  )
}

# Abundance tables are wide: one row per sampleID, one column per MAG/bin.
# Mixed-model input needs one response value per row, so pivot MAG/bin columns
# into two columns: MAG identifier and numeric abundance response Y.
abundance_long <- abundance_wide %>%
  tidyr::pivot_longer(
    cols = -sampleID,
    names_to = "MAG",
    values_to = "Y"
  )

raw_y_text <- as.character(abundance_long$Y)
parsed_y <- suppressWarnings(as.numeric(raw_y_text))
non_numeric_y <- !is.na(abundance_long$Y) & nzchar(trimws(raw_y_text)) & is.na(parsed_y)
if (any(non_numeric_y)) {
  bad_values <- unique(raw_y_text[non_numeric_y])
  stop_with_examples(
    "Abundance values must be numeric; found non-numeric values such as: ",
    bad_values
  )
}
abundance_long$Y <- parsed_y

cat(sprintf("Long abundance dimensions: %d rows x %d columns\n", nrow(abundance_long), ncol(abundance_long)))

# Join metadata after pivoting so each sampleID x MAG observation carries the
# subject, time, date, and candidate covariates needed for modeling.
df_joined <- abundance_long %>%
  dplyr::left_join(sample_info, by = "sampleID")

critical_columns <- c("sampleID", "MAG", "Y", "Subject", "Timepoint", "date_converted")

df_typed <- df_joined %>%
  dplyr::mutate(
    Y = as.numeric(Y),
    Subject = factor(group),
    Timepoint = factor(Timepoint, levels = paste0("T", 1:10), ordered = TRUE),
    date_converted = lubridate::ymd(as.character(date)),
    gender = factor(gender),
    age = as.numeric(age),
    altitude = as.numeric(altitude)
  ) %>%
  dplyr::mutate(
    days_from_baseline_global = as.numeric(date_converted - min(date_converted, na.rm = TRUE))
  ) %>%
  dplyr::group_by(Subject) %>%
  dplyr::mutate(
    days_from_baseline_subject = as.numeric(date_converted - min(date_converted, na.rm = TRUE))
  ) %>%
  dplyr::ungroup()

print_missing_counts(
  "\nMissing values in model-critical columns before filtering:",
  df_typed,
  critical_columns
)

initial_rows <- nrow(df_typed)

df <- df_typed %>%
  dplyr::filter(
    !is.na(sampleID),
    !is.na(MAG),
    !is.na(Y),
    !is.na(Subject),
    !is.na(Timepoint),
    !is.na(date_converted)
  )

rows_removed <- initial_rows - nrow(df)

print_missing_counts(
  "\nMissing values in model-critical columns after filtering:",
  df,
  critical_columns
)

# --- Validation and summaries ---
cat("\nCleaning summary:\n")
cat(sprintf("Cleaned data dimensions: %d rows x %d columns\n", nrow(df), ncol(df)))
cat(sprintf("Rows removed for missing model-critical values: %d\n", rows_removed))
cat(sprintf("Unique sampleID: %d\n", dplyr::n_distinct(df$sampleID)))
cat(sprintf("Unique Subject: %d\n", dplyr::n_distinct(df$Subject)))
cat(sprintf("Unique MAG: %d\n", dplyr::n_distinct(df$MAG)))
cat("Observed Timepoint levels:", paste(levels(droplevels(df$Timepoint)), collapse = ", "), "\n")
cat(
  "Date range:",
  format(min(df$date_converted, na.rm = TRUE)),
  "to",
  format(max(df$date_converted, na.rm = TRUE)),
  "\n"
)

if (!is.null(config$model$formula)) {
  model_formula <- as.formula(config$model$formula)
  formula_columns <- all.vars(model_formula)
  missing_formula_columns <- setdiff(formula_columns, names(df))

  if (length(missing_formula_columns) > 0) {
    warning(
      "Processed data does not contain all variables in config$model$formula: ",
      paste(missing_formula_columns, collapse = ", "),
      ". This is expected until config/config.yml is updated for the inspected dataset.",
      call. = FALSE
    )
  }
}

# --- Save processed ---
dir.create(dirname(processed_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(df, processed_path)
cat("Saved cleaned data to ", processed_path, "\n", sep = "")
