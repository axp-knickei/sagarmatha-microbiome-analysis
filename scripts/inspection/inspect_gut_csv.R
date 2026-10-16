# Inspect the gut.csv MAG abundance table and write reproducible summaries.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) y else x
}

parse_args <- function(args) {
  parsed <- list(
    input = Sys.getenv("BLACKICE_GUT_CSV", unset = ""),
    output_dir = "results/data_inspection/gut_csv"
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      cat(
        "Usage: Rscript scripts/inspection/inspect_gut_csv.R [options]\n\n",
        "Options:\n",
        "  --input PATH       Path to gut.csv. Default: BLACKICE_GUT_CSV, data/gut.csv, then data/raw/gut.csv.\n",
        "  --output-dir PATH  Output directory. Default: results/data_inspection/gut_csv\n",
        "  --help             Show this help.\n",
        sep = ""
      )
      quit(status = 0)
    } else if (arg %in% c("--input", "--output-dir")) {
      if (i == length(args)) {
        stop("Missing value for argument: ", arg, call. = FALSE)
      }
      value <- args[[i + 1]]
      if (arg == "--input") parsed$input <- value
      if (arg == "--output-dir") parsed$output_dir <- value
      i <- i + 2
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }

  if (!nzchar(parsed$input)) {
    candidates <- c("data/gut.csv", "data/raw/gut.csv")
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0) {
      parsed$input <- existing[[1]]
    }
  }

  parsed
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (!nzchar(args$input) || !file.exists(args$input)) {
  stop(
    "Could not find gut.csv. Use --input PATH or set BLACKICE_GUT_CSV. Tried data/gut.csv and data/raw/gut.csv.",
    call. = FALSE
  )
}

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading: ", args$input)
gut <- utils::read.csv(
  args$input,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)

if (!"sampleID" %in% names(gut)) {
  stop("Expected a sampleID column in gut.csv.", call. = FALSE)
}

sample_id <- as.character(gut$sampleID)
feature_cols <- setdiff(names(gut), "sampleID")
feature_data <- gut[, feature_cols, drop = FALSE]

numeric_flags <- vapply(feature_data, is.numeric, logical(1))
if (!all(numeric_flags)) {
  converted <- lapply(feature_data, function(x) suppressWarnings(as.numeric(as.character(x))))
  nonnumeric_after_conversion <- vapply(converted, function(x) all(is.na(x)), logical(1))
  if (any(nonnumeric_after_conversion)) {
    stop(
      "Some feature columns could not be converted to numeric: ",
      paste(names(feature_data)[nonnumeric_after_conversion], collapse = ", "),
      call. = FALSE
    )
  }
  feature_data <- as.data.frame(converted, check.names = FALSE)
}

value_matrix <- as.matrix(feature_data)
missing_by_col <- colSums(is.na(value_matrix))
missing_by_row <- rowSums(is.na(value_matrix))
zero_by_col <- colSums(value_matrix == 0, na.rm = TRUE)
zero_by_row <- rowSums(value_matrix == 0, na.rm = TRUE)
nonzero_by_col <- colSums(value_matrix != 0 & !is.na(value_matrix))
nonzero_by_row <- rowSums(value_matrix != 0 & !is.na(value_matrix))

all_values <- as.numeric(value_matrix)
nonmissing_values <- all_values[!is.na(all_values)]
nonzero_values <- nonmissing_values[nonmissing_values != 0]

sample_duplicates <- sample_id[duplicated(sample_id)]
duplicate_row_count <- sum(duplicated(gut))
duplicate_feature_names <- feature_cols[duplicated(feature_cols)]

overview <- data.frame(
  metric = c(
    "input_file",
    "n_rows",
    "n_columns",
    "n_sample_id_columns",
    "n_feature_columns",
    "n_duplicate_sample_ids",
    "n_duplicate_rows",
    "n_duplicate_feature_names",
    "n_missing_values_total",
    "n_missing_sample_ids",
    "n_missing_feature_values",
    "n_zero_values",
    "n_nonzero_values",
    "zero_fraction",
    "min_value",
    "max_value",
    "min_nonzero_value",
    "max_nonzero_value",
    "mean_value",
    "median_value",
    "mean_nonzero_value",
    "median_nonzero_value",
    "row_unit_inference",
    "column_unit_inference",
    "value_scale_inference"
  ),
  value = c(
    args$input,
    nrow(gut),
    ncol(gut),
    sum(names(gut) == "sampleID"),
    length(feature_cols),
    length(unique(sample_duplicates)),
    duplicate_row_count,
    length(unique(duplicate_feature_names)),
    sum(is.na(gut)),
    sum(is.na(sample_id) | !nzchar(sample_id)),
    sum(is.na(value_matrix)),
    sum(value_matrix == 0, na.rm = TRUE),
    sum(value_matrix != 0 & !is.na(value_matrix)),
    if (length(nonmissing_values) > 0) mean(nonmissing_values == 0) else NA_real_,
    if (length(nonmissing_values) > 0) min(nonmissing_values) else NA_real_,
    if (length(nonmissing_values) > 0) max(nonmissing_values) else NA_real_,
    if (length(nonzero_values) > 0) min(nonzero_values) else NA_real_,
    if (length(nonzero_values) > 0) max(nonzero_values) else NA_real_,
    if (length(nonmissing_values) > 0) mean(nonmissing_values) else NA_real_,
    if (length(nonmissing_values) > 0) stats::median(nonmissing_values) else NA_real_,
    if (length(nonzero_values) > 0) mean(nonzero_values) else NA_real_,
    if (length(nonzero_values) > 0) stats::median(nonzero_values) else NA_real_,
    "rows are microbiome samples identified by sampleID",
    "columns after sampleID are MAG abundance features",
    "non-negative sparse continuous relative-abundance-like values"
  ),
  stringsAsFactors = FALSE
)

column_summary <- data.frame(
  column_name = names(gut),
  inferred_role = ifelse(names(gut) == "sampleID", "sample_identifier", "MAG_abundance_feature"),
  r_type = vapply(gut, function(x) paste(class(x), collapse = "|"), character(1)),
  missing_count = vapply(gut, function(x) sum(is.na(x)), integer(1)),
  unique_count = vapply(gut, function(x) length(unique(x)), integer(1)),
  stringsAsFactors = FALSE
)

mag_summary <- data.frame(
  MAG = feature_cols,
  missing_count = missing_by_col,
  zero_count = zero_by_col,
  nonzero_count = nonzero_by_col,
  zero_fraction = zero_by_col / nrow(gut),
  min_value = apply(value_matrix, 2, min, na.rm = TRUE),
  max_value = apply(value_matrix, 2, max, na.rm = TRUE),
  mean_value = colMeans(value_matrix, na.rm = TRUE),
  median_value = apply(value_matrix, 2, stats::median, na.rm = TRUE),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_summary <- data.frame(
  sampleID = sample_id,
  missing_feature_count = missing_by_row,
  zero_feature_count = zero_by_row,
  nonzero_feature_count = nonzero_by_row,
  zero_fraction = zero_by_row / length(feature_cols),
  row_sum = rowSums(value_matrix, na.rm = TRUE),
  min_value = apply(value_matrix, 1, min, na.rm = TRUE),
  max_value = apply(value_matrix, 1, max, na.rm = TRUE),
  stringsAsFactors = FALSE
)

duplicate_samples <- data.frame(
  sampleID = unique(sample_duplicates),
  occurrence_count = as.integer(table(sample_id)[unique(sample_duplicates)]),
  stringsAsFactors = FALSE
)

utils::write.csv(overview, file.path(args$output_dir, "gut_csv_overview.csv"), row.names = FALSE)
utils::write.csv(column_summary, file.path(args$output_dir, "gut_csv_column_summary.csv"), row.names = FALSE)
utils::write.csv(mag_summary, file.path(args$output_dir, "gut_csv_mag_summary.csv"), row.names = FALSE)
utils::write.csv(sample_summary, file.path(args$output_dir, "gut_csv_sample_summary.csv"), row.names = FALSE)
utils::write.csv(duplicate_samples, file.path(args$output_dir, "gut_csv_duplicate_sample_ids.csv"), row.names = FALSE)

summary_md <- file.path(args$output_dir, "gut_csv_inspection_summary.md")
writeLines(
  c(
    "# gut.csv Inspection Summary",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    paste0("- Input file: `", args$input, "`"),
    paste0("- Rows: ", nrow(gut)),
    paste0("- Columns: ", ncol(gut)),
    paste0("- MAG abundance feature columns: ", length(feature_cols)),
    paste0("- Duplicate sample IDs: ", length(unique(sample_duplicates))),
    paste0("- Duplicate full rows: ", duplicate_row_count),
    paste0("- Missing values in feature matrix: ", sum(is.na(value_matrix))),
    paste0("- Zero values in feature matrix: ", sum(value_matrix == 0, na.rm = TRUE)),
    paste0("- Non-zero values in feature matrix: ", sum(value_matrix != 0 & !is.na(value_matrix))),
    paste0("- Zero fraction: ", format(round(mean(nonmissing_values == 0), 6), nsmall = 6)),
    paste0("- Value range: ", format(min(nonmissing_values), scientific = TRUE), " to ", format(max(nonmissing_values), scientific = TRUE)),
    paste0("- Non-zero value range: ", format(min(nonzero_values), scientific = TRUE), " to ", format(max(nonzero_values), scientific = TRUE)),
    "",
    "## Inferred Organization",
    "",
    "Rows appear to represent gut microbiome samples identified by `sampleID`.",
    "Columns after `sampleID` appear to represent MAG abundance features.",
    "Values are non-negative, sparse, continuous, and mostly less than one, consistent with relative-abundance-like MAG abundance data rather than raw integer counts.",
    "",
    "## Output Files",
    "",
    "- `gut_csv_overview.csv`",
    "- `gut_csv_column_summary.csv`",
    "- `gut_csv_mag_summary.csv`",
    "- `gut_csv_sample_summary.csv`",
    "- `gut_csv_duplicate_sample_ids.csv`"
  ),
  summary_md
)

message("Wrote inspection outputs to: ", args$output_dir)
