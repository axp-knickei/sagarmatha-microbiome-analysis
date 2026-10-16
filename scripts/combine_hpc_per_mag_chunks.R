# ============================================
# combine_hpc_per_mag_chunks.R
# ============================================
# Merges chunk outputs from scripts/hpc_run_per_mag_lmertest.R and recomputes
# final Benjamini-Hochberg FDR across the combined MAG set.

parse_args <- function(args) {
  parsed <- list(
    input_dir = "results/hpc_per_mag_lmertest",
    output_dir = "results/hpc_per_mag_lmertest"
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      cat(
        "Usage: Rscript scripts/combine_hpc_per_mag_chunks.R [options]\n\n",
        "Options:\n",
        "  --input-dir PATH    Directory containing chunks/. Default: results/hpc_per_mag_lmertest\n",
        "  --output-dir PATH   Directory for combined outputs. Default: results/hpc_per_mag_lmertest\n",
        "  --help              Show this help.\n",
        sep = ""
      )
      quit(status = 0)
    } else if (arg %in% c("--input-dir", "--output-dir")) {
      if (i == length(args)) {
        stop("Missing value for argument: ", arg, call. = FALSE)
      }
      value <- args[[i + 1]]
      if (arg == "--input-dir") {
        parsed$input_dir <- sub("/+$", "", value)
      } else if (arg == "--output-dir") {
        parsed$output_dir <- sub("/+$", "", value)
      }
      i <- i + 2
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }

  parsed
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

args <- parse_args(commandArgs(trailingOnly = TRUE))

chunk_dir <- file.path(args$input_dir, "chunks")
chunk_files <- list.files(
  chunk_dir,
  pattern = "^per_mag_lmertest_hpc_results\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(chunk_files) == 0) {
  stop("No chunk result files found under: ", chunk_dir, call. = FALSE)
}

tables_dir <- file.path(args$output_dir, "tables")
logs_dir <- file.path(args$output_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

combined_csv <- file.path(tables_dir, "per_mag_lmertest_hpc_results_combined.csv")
summary_path <- file.path(logs_dir, "per_mag_lmertest_hpc_combined_summary.txt")

chunk_results <- lapply(chunk_files, function(path) {
  result <- utils::read.csv(path, stringsAsFactors = FALSE)
  result$chunk_source_file <- path
  result
})
results <- do.call(rbind, chunk_results)

required_cols <- c("MAG", "status", "time_p_value")
missing_cols <- setdiff(required_cols, names(results))
if (length(missing_cols) > 0) {
  stop("Chunk outputs are missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

duplicated_mags <- unique(results$MAG[duplicated(results$MAG)])
if (length(duplicated_mags) > 0) {
  stop_with_examples(
    "MAG IDs appeared in more than one chunk output; cannot compute final FDR safely: ",
    duplicated_mags
  )
}

results$time_p_adj_bh_all <- stats::p.adjust(results$time_p_value, method = "BH")
results$primary_analysis <- results$status == "lmm_ok"
results$time_p_adj_bh_primary <- NA_real_
results$time_p_adj_bh_primary[results$primary_analysis] <- stats::p.adjust(
  results$time_p_value[results$primary_analysis],
  method = "BH"
)
results$primary_significant <- results$primary_analysis & results$time_p_adj_bh_primary < 0.05
results$nominal_p_lt_0_05 <- results$time_p_value < 0.05
results$fdr_scope <- "combined_across_chunks"

utils::write.csv(results, combined_csv, row.names = FALSE)

summary_connection <- file(summary_path, open = "wt")
sink(summary_connection, split = TRUE)

cat("Combined HPC Per-MAG lmerTest LMM Summary\n")
cat("=========================================\n\n")
cat("Input directory:", args$input_dir, "\n")
cat("Chunk result files:", length(chunk_files), "\n")
cat("MAGs combined:", nrow(results), "\n")
cat("Output CSV:", combined_csv, "\n")
cat("Summary file:", summary_path, "\n\n")

cat("Status counts:\n")
print(table(results$status))
cat("\n")

cat("Significance counts:\n")
cat("Raw p < 0.05:", sum(results$nominal_p_lt_0_05 %in% TRUE, na.rm = TRUE), "\n")
cat("BH q < 0.05 across all combined MAGs:", sum(results$time_p_adj_bh_all < 0.05, na.rm = TRUE), "\n")
cat("BH q < 0.05 in primary-analysis set:", sum(results$time_p_adj_bh_primary < 0.05, na.rm = TRUE), "\n")
cat("Primary significant MAGs:", sum(results$primary_significant %in% TRUE, na.rm = TRUE), "\n\n")

cat("Top 20 MAGs ranked by primary BH q-value, then raw p-value:\n")
ranked <- results[order(
  is.na(results$time_p_adj_bh_primary),
  results$time_p_adj_bh_primary,
  results$time_p_value
), ]
top_cols <- intersect(
  c(
    "MAG", "status", "primary_analysis", "primary_significant",
    "time_estimate", "time_std_error", "time_df", "time_statistic",
    "time_p_value", "time_p_adj_bh_all", "time_p_adj_bh_primary",
    "zero_fraction", "nonzero_sample_count", "subjects_with_nonzero",
    "response_variance", "singular", "chunk_source_file"
  ),
  names(ranked)
)
print(utils::head(ranked[, top_cols], 20), row.names = FALSE)

sink()
close(summary_connection)

cat("Saved combined chunk results to ", combined_csv, "\n", sep = "")
cat("Saved combined summary to ", summary_path, "\n", sep = "")
