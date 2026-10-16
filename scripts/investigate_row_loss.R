# ============================================
# investigate_row_loss.R - Sample Row Loss Trace
# ============================================
# Compares raw abundance sampleID values with processed sampleID values and
# reports why any raw samples were removed during data preparation.

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

read_csv_file <- function(path) {
  if (has_readr) {
    return(readr::read_csv(path, show_col_types = FALSE))
  }

  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

classify_removed_sample <- function(row) {
  if (is.na(row$in_sample_info) || !row$in_sample_info) {
    return("missing from sample_info.csv")
  }

  if (is.na(row$Timepoint_raw) || identical(row$Timepoint_raw, "")) {
    return("Timepoint missing")
  }

  if (is.na(row$Timepoint_model)) {
    return("Timepoint outside T1-T10")
  }

  if (is.na(row$date_converted)) {
    return("date conversion failed")
  }

  "another reason"
}

config <- yaml::read_yaml(here("config", "config.yml"))

raw_dir <- config$data$raw_dir %||% dirname(config$data$input %||% "")
metadata_filename <- config$data$metadata %||% config$data$metadata_file %||% "sample_info.csv"
abundance_filename <- config$data$abundance %||% config$data$abundance_file %||% "gut.csv"
processed_path <- resolve_path(config$data$processed %||% "data/processed/cleaned_data.rds")

metadata_path <- file.path(resolve_path(raw_dir), metadata_filename)
abundance_path <- file.path(resolve_path(raw_dir), abundance_filename)
report_path <- here("data", "processed", "row_loss_investigation.txt")

dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

sink(report_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Row Loss Investigation\n")
cat("======================\n\n")
cat("Abundance file:", abundance_path, "\n")
cat("Metadata file:", metadata_path, "\n")
cat("Processed file:", processed_path, "\n")
cat("Report file:", report_path, "\n\n")

if (!file.exists(abundance_path)) {
  stop("Abundance file not found: ", abundance_path, call. = FALSE)
}
if (!file.exists(metadata_path)) {
  stop("Metadata file not found: ", metadata_path, call. = FALSE)
}
if (!file.exists(processed_path)) {
  stop("Processed data not found: ", processed_path, call. = FALSE)
}

gut_raw <- read_csv_file(abundance_path)
sample_info <- read_csv_file(metadata_path)
cleaned <- readRDS(processed_path)

if (!"sampleID" %in% names(gut_raw)) {
  stop("Raw abundance file is missing sampleID.", call. = FALSE)
}
if (!"sampleID" %in% names(sample_info)) {
  stop("sample_info.csv is missing sampleID.", call. = FALSE)
}
if (!"sampleID" %in% names(cleaned)) {
  stop("cleaned_data.rds is missing sampleID.", call. = FALSE)
}

raw_sample_ids <- unique(gut_raw$sampleID)
cleaned_sample_ids <- unique(cleaned$sampleID)
removed_sample_ids <- setdiff(raw_sample_ids, cleaned_sample_ids)

cat("Raw gut.csv sampleID count:", length(raw_sample_ids), "\n")
cat("Cleaned sampleID count:", length(cleaned_sample_ids), "\n")
cat("Removed sampleID count:", length(removed_sample_ids), "\n")
cat("MAG columns in raw gut.csv:", ncol(gut_raw) - 1, "\n")
cat("Expected removed long rows:", length(removed_sample_ids) * (ncol(gut_raw) - 1), "\n\n")

cat("Removed sampleID values:\n")
print(removed_sample_ids)
cat("\n")

sample_info_lookup <- sample_info
sample_info_lookup$in_sample_info <- TRUE

removed <- data.frame(sampleID = removed_sample_ids, stringsAsFactors = FALSE)
removed_joined <- merge(
  removed,
  sample_info_lookup,
  by = "sampleID",
  all.x = TRUE,
  sort = FALSE
)

removed_joined$in_sample_info[is.na(removed_joined$in_sample_info)] <- FALSE

if ("Timepoint" %in% names(removed_joined)) {
  removed_joined$Timepoint_raw <- removed_joined$Timepoint
  removed_joined$Timepoint_model <- factor(
    removed_joined$Timepoint,
    levels = paste0("T", 1:10),
    ordered = TRUE
  )
} else {
  removed_joined$Timepoint_raw <- NA_character_
  removed_joined$Timepoint_model <- factor(NA_character_, levels = paste0("T", 1:10), ordered = TRUE)
}

if ("date" %in% names(removed_joined) && requireNamespace("lubridate", quietly = TRUE)) {
  removed_joined$date_converted <- lubridate::ymd(as.character(removed_joined$date))
} else if ("date" %in% names(removed_joined)) {
  removed_joined$date_converted <- as.Date(as.character(removed_joined$date), format = "%Y%m%d")
} else {
  removed_joined$date_converted <- as.Date(NA)
}

removed_joined$removal_reason <- vapply(
  seq_len(nrow(removed_joined)),
  function(i) classify_removed_sample(removed_joined[i, , drop = FALSE]),
  character(1)
)

cat("Removed sampleID joined to sample_info.csv:\n")
print(removed_joined)
cat("\n")

cat("Removal reason counts:\n")
print(table(removed_joined$removal_reason, useNA = "ifany"))
cat("\n")

cat("Compact removal table:\n")
compact_columns <- intersect(
  c("sampleID", "group", "Timepoint_raw", "date", "date_converted", "gender", "age", "altitude", "in_sample_info", "removal_reason"),
  names(removed_joined)
)
print(removed_joined[, compact_columns, drop = FALSE], row.names = FALSE)
cat("\n")

all_outside_timepoint <- all(removed_joined$removal_reason == "Timepoint outside T1-T10")
all_have_metadata <- all(removed_joined$in_sample_info)
all_dates_ok <- all(!is.na(removed_joined$date_converted))

cat("Interpretation\n")
cat("--------------\n")
if (all_outside_timepoint && all_have_metadata && all_dates_ok) {
  cat(
    "All removed sampleID values were present in sample_info.csv and had valid dates, ",
    "but their Timepoint values were outside the configured T1-T10 levels. ",
    "This row loss is expected under the current filtering rule, but it is analytically concerning ",
    "because those sample/time points are being excluded by design and should be confirmed before modeling.\n",
    sep = ""
  )
} else {
  cat(
    "The removed samples have mixed removal reasons. Review the compact removal table above before changing filters or modeling.\n"
  )
}

cat("\nSaved row-loss investigation to data/processed/row_loss_investigation.txt\n")
