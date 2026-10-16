# Combine per-workflow reproducibility execution records into a supplementary
# table-friendly project-level CSV and Markdown file.

required_columns <- c(
  "analysis_name",
  "script_path",
  "scheduler",
  "job_id",
  "start_time",
  "completion_time",
  "runtime_seconds",
  "input_data_file",
  "window_definition_file",
  "config_files",
  "conda_environment",
  "conda_prefix",
  "loaded_modules",
  "number_of_tests",
  "fdr_family",
  "result_directory",
  "main_results_file",
  "summary_log_file",
  "session_info_file",
  "git_branch",
  "git_commit",
  "git_status_short",
  "notes"
)

markdown_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  gsub("\\|", "\\\\|", x)
}

write_markdown_table <- function(df, path) {
  lines <- c(
    "# Combined HPC Execution Records",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    paste(required_columns, collapse = " | "),
    paste(rep("---", length(required_columns)), collapse = " | ")
  )
  if (nrow(df) > 0) {
    for (i in seq_len(nrow(df))) {
      lines <- c(lines, paste(markdown_escape(df[i, required_columns, drop = TRUE]), collapse = " | "))
    }
  }
  writeLines(lines, path)
}

record_files <- list.files(
  "results",
  pattern = "^hpc_execution_record\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
record_files <- record_files[grepl("/reproducibility/hpc_execution_record\\.csv$", record_files)]

records <- lapply(record_files, function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(required_columns, names(df))
  for (column in missing) {
    df[[column]] <- NA_character_
  }
  df[, required_columns, drop = FALSE]
})

combined <- if (length(records) == 0) {
  as.data.frame(stats::setNames(rep(list(character()), length(required_columns)), required_columns))
} else {
  do.call(rbind, records)
}

if (nrow(combined) > 0) {
  combined <- combined[order(combined$analysis_name, combined$completion_time), , drop = FALSE]
}

csv_path <- file.path("results", "hpc_execution_records_combined.csv")
md_path <- file.path("results", "hpc_execution_records_combined.md")
utils::write.csv(combined, csv_path, row.names = FALSE, na = "")
write_markdown_table(combined, md_path)

cat("Combined HPC execution records:", nrow(combined), "\n")
cat("CSV:", csv_path, "\n")
cat("Markdown:", md_path, "\n")
