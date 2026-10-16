# ============================================
# validate_processed_data.R - Processed Data Validation
# ============================================
# Reads data/processed/cleaned_data.rds, prints inspection summaries, validates
# model-ready columns and types, and writes a readable validation report.

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

check_item <- function(ok, label, detail = NULL) {
  box <- if (isTRUE(ok)) "[x]" else "[ ]"
  line <- paste(box, label)
  if (!is.null(detail) && nzchar(detail)) {
    line <- paste(line, "-", detail)
  }
  line
}

count_missing <- function(df, columns) {
  stats::setNames(
    vapply(columns, function(column) {
      if (!column %in% names(df)) return(NA_integer_)
      sum(is.na(df[[column]]))
    }, integer(1)),
    columns
  )
}

print_named_counts <- function(title, counts) {
  cat(title, "\n", sep = "")
  for (name in names(counts)) {
    cat(sprintf("  %s: %s\n", name, counts[[name]]))
  }
}

read_csv_file <- function(path) {
  if (has_readr) {
    return(readr::read_csv(path, show_col_types = FALSE))
  }

  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

config <- yaml::read_yaml(here("config", "config.yml"))
processed_path <- resolve_path(config$data$processed %||% "data/processed/cleaned_data.rds")
report_path <- here("data", "processed", "validation_report.txt")

raw_dir <- config$data$raw_dir %||% dirname(config$data$input %||% "")
abundance_filename <- config$data$abundance %||% config$data$abundance_file %||% "gut.csv"
abundance_path <- file.path(resolve_path(raw_dir), abundance_filename)

dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

sink(report_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Processed Data Validation Report\n")
cat("================================\n\n")
cat("Processed file:", processed_path, "\n")
cat("Selected abundance file:", abundance_filename, "\n")
cat("Validation report:", report_path, "\n\n")

if (!file.exists(processed_path)) {
  stop("Processed data file not found: ", processed_path, call. = FALSE)
}

df <- readRDS(processed_path)

cat("Basic Inspection\n")
cat("----------------\n")
cat("Dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
cat("Column names:\n")
print(names(df))
cat("\nStructure:\n")
str(df)
cat("\nSummary:\n")
print(summary(df))
cat("\nHead:\n")
print(utils::head(df))

cat("\nKey Counts\n")
cat("----------\n")
cat("Unique sampleID:", length(unique(df$sampleID)), "\n")
cat("Unique Subject:", length(unique(df$Subject)), "\n")
cat("Unique MAG:", length(unique(df$MAG)), "\n")
cat("Observed Timepoint values:", paste(sort(unique(as.character(df$Timepoint))), collapse = ", "), "\n")
cat(
  "days_from_baseline_subject range:",
  paste(range(df$days_from_baseline_subject, na.rm = TRUE), collapse = " to "),
  "\n"
)
cat(
  "Date range:",
  paste(format(range(df$date_converted, na.rm = TRUE)), collapse = " to "),
  "\n"
)

required_cols <- c(
  "sampleID",
  "MAG",
  "Y",
  "Subject",
  "Timepoint",
  "date_converted",
  "days_from_baseline_subject",
  "gender",
  "age",
  "altitude"
)

missing_required_cols <- setdiff(required_cols, names(df))

cat("\nRequired Columns\n")
cat("----------------\n")
if (length(missing_required_cols) == 0) {
  cat("All required columns are present.\n")
} else {
  cat("Missing required columns:", paste(missing_required_cols, collapse = ", "), "\n")
}

type_checks <- c(
  Y_numeric = "Y" %in% names(df) && is.numeric(df$Y),
  Subject_factor = "Subject" %in% names(df) && is.factor(df$Subject),
  MAG_character_or_factor = "MAG" %in% names(df) && (is.character(df$MAG) || is.factor(df$MAG)),
  Timepoint_factor = "Timepoint" %in% names(df) && is.factor(df$Timepoint),
  date_converted_Date = "date_converted" %in% names(df) && inherits(df$date_converted, "Date"),
  days_from_baseline_subject_numeric = "days_from_baseline_subject" %in% names(df) &&
    is.numeric(df$days_from_baseline_subject),
  gender_factor = "gender" %in% names(df) && is.factor(df$gender),
  age_numeric = "age" %in% names(df) && is.numeric(df$age),
  altitude_numeric = "altitude" %in% names(df) && is.numeric(df$altitude)
)

cat("\nVariable Types\n")
cat("--------------\n")
for (name in names(type_checks)) {
  cat(check_item(type_checks[[name]], name), "\n")
}

cat("\nRepeated Observations\n")
cat("---------------------\n")
rows_per_subject <- sort(table(df$Subject), decreasing = TRUE)
timepoints_per_subject <- stats::aggregate(Timepoint ~ Subject, df, function(x) length(unique(x)))
sampleids_per_subject <- stats::aggregate(sampleID ~ Subject, df, function(x) length(unique(x)))

cat("Rows per Subject:\n")
print(rows_per_subject)
cat("\nUnique Timepoint count per Subject:\n")
print(timepoints_per_subject)
cat("\nUnique sampleID count per Subject:\n")
print(sampleids_per_subject)

subjects_with_one_sample <- sampleids_per_subject$Subject[sampleids_per_subject$sampleID <= 1]
if (length(subjects_with_one_sample) > 0) {
  cat(
    "Warning: some subjects have only one sample/timepoint:",
    paste(subjects_with_one_sample, collapse = ", "),
    "\n"
  )
  warning(
    "Some subjects have only one sample/timepoint: ",
    paste(subjects_with_one_sample, collapse = ", "),
    call. = FALSE
  )
} else {
  cat("All subjects have repeated sampleID observations.\n")
}

cat("\nTime Validation\n")
cat("---------------\n")
cat("Timepoint table:\n")
print(table(df$Timepoint))
cat("\nUnique days_from_baseline_subject values:\n")
print(sort(unique(df$days_from_baseline_subject)))

baseline_by_subject <- stats::aggregate(days_from_baseline_subject ~ Subject, df, min, na.rm = TRUE)
subjects_without_baseline_zero <- baseline_by_subject$Subject[
  is.na(baseline_by_subject$days_from_baseline_subject) |
    baseline_by_subject$days_from_baseline_subject != 0
]

if (length(subjects_without_baseline_zero) > 0) {
  cat(
    "Warning: some subjects do not have days_from_baseline_subject starting at 0:",
    paste(subjects_without_baseline_zero, collapse = ", "),
    "\n"
  )
  warning(
    "Some subjects do not have days_from_baseline_subject starting at 0: ",
    paste(subjects_without_baseline_zero, collapse = ", "),
    call. = FALSE
  )
} else {
  cat("days_from_baseline_subject starts at 0 for every Subject.\n")
}

model_critical_cols <- c(
  "sampleID",
  "MAG",
  "Y",
  "Subject",
  "Timepoint",
  "date_converted",
  "days_from_baseline_subject"
)

cat("\nMissing Values\n")
cat("--------------\n")
critical_missing <- count_missing(df, model_critical_cols)
print_named_counts("Model-critical missing values:", critical_missing)

optional_covariates <- c("gender", "age", "altitude")
optional_missing <- count_missing(df, optional_covariates)
print_named_counts("\nOptional covariate missing values:", optional_missing)

cat("\nRow Loss\n")
cat("--------\n")
expected_long_rows <- NA_integer_
if (file.exists(abundance_path)) {
  abundance_wide <- read_csv_file(abundance_path)
  expected_long_rows <- nrow(abundance_wide) * (ncol(abundance_wide) - 1)
  cat("Expected long rows:", expected_long_rows, "\n")
  cat("Observed processed rows:", nrow(df), "\n")
  cat("Rows removed:", expected_long_rows - nrow(df), "\n")

  row_loss_fraction <- if (expected_long_rows > 0) (expected_long_rows - nrow(df)) / expected_long_rows else NA_real_
  if (!is.na(row_loss_fraction) && row_loss_fraction > 0.05) {
    cat(sprintf("Warning: more than 5%% of expected long rows were removed: %.2f%%\n", row_loss_fraction * 100))
    warning(
      sprintf("More than 5%% of expected long rows were removed: %.2f%%", row_loss_fraction * 100),
      call. = FALSE
    )
  } else {
    cat("No unexpected massive row loss detected.\n")
  }
} else {
  cat("Could not calculate expected long rows because raw abundance file was not found:\n")
  cat("  ", abundance_path, "\n", sep = "")
}

cat("\nGroup Column\n")
cat("------------\n")
if ("Group" %in% names(df)) {
  print(table(df$Group))
} else {
  cat("No Group column found. This is expected for the current dataset because no treatment/condition Group has been identified.\n")
}

subject_has_repeats <- exists("sampleids_per_subject") &&
  nrow(sampleids_per_subject) > 0 &&
  all(sampleids_per_subject$sampleID > 1)
time_has_expected_values <- "Timepoint" %in% names(df) &&
  all(na.omit(as.character(df$Timepoint)) %in% paste0("T", 1:10)) &&
  length(subjects_without_baseline_zero) == 0
no_critical_missing <- all(critical_missing == 0, na.rm = FALSE)
no_massive_row_loss <- is.na(expected_long_rows) ||
  expected_long_rows == 0 ||
  ((expected_long_rows - nrow(df)) / expected_long_rows) <= 0.05

cat("\nChecklist Summary\n")
cat("-----------------\n")
cat(check_item(type_checks[["Y_numeric"]], "Y is numeric"), "\n")
cat(check_item(subject_has_repeats, "Subject has repeated observations"), "\n")
cat("[ ] Group has expected levels - not applicable yet / no true Group identified\n")
cat(check_item(time_has_expected_values, "Time has expected values"), "\n")
cat(check_item(no_massive_row_loss, "No unexpected massive row loss occurred"), "\n")
cat(check_item(length(missing_required_cols) == 0 && no_critical_missing, "No important columns are missing"), "\n")

cat("\nValidation report saved to data/processed/validation_report.txt\n")
