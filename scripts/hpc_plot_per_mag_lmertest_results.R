#!/usr/bin/env Rscript

# Orchestrate post-processing plots for the full T1-T10 per-MAG LMM workflow.
#
# This script ensures the canonical result files exist, optionally combines
# chunked outputs, and then runs the per-MAG plotting scripts in a stable order.

parse_args <- function(args) {
  parsed <- list(
    input_dir = "results/hpc_per_mag_lmertest",
    results_csv = "",
    output_dir = "results/hpc_per_mag_lmertest",
    combine_chunks = TRUE,
    write_root_metadata = TRUE
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      cat(
        "Usage: Rscript scripts/hpc_plot_per_mag_lmertest_results.R [options]\n\n",
        "Options:\n",
        "  --input-dir PATH         Directory containing HPC outputs. Default: results/hpc_per_mag_lmertest\n",
        "  --results-csv PATH       Canonical or combined result CSV to use.\n",
        "  --output-dir PATH        Output directory. Default: results/hpc_per_mag_lmertest\n",
        "  --no-combine-chunks      Do not attempt to combine chunk outputs automatically.\n",
        "  --no-root-metadata       Do not rewrite the canonical summary/session/record files.\n",
        "  --help                   Show this help.\n",
        sep = ""
      )
      quit(status = 0)
    } else if (arg %in% c("--input-dir", "--results-csv", "--output-dir")) {
      if (i == length(args)) stop("Missing value for argument: ", arg, call. = FALSE)
      value <- args[[i + 1]]
      if (arg == "--input-dir") {
        parsed$input_dir <- sub("/+$", "", value)
      } else if (arg == "--results-csv") {
        parsed$results_csv <- value
      } else if (arg == "--output-dir") {
        parsed$output_dir <- sub("/+$", "", value)
      }
      i <- i + 2
    } else if (arg == "--no-combine-chunks") {
      parsed$combine_chunks <- FALSE
      i <- i + 1
    } else if (arg == "--no-root-metadata") {
      parsed$write_root_metadata <- FALSE
      i <- i + 1
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }

  parsed
}

run_rscript <- function(script, args = character()) {
  message("Running: Rscript ", script, if (length(args) > 0) paste0(" ", paste(args, collapse = " ")) else "")
  status <- system2("Rscript", c(script, args))
  exit_status <- attr(status, "status")
  if (!is.null(exit_status) && !identical(exit_status, 0L)) {
    stop("Failed to run: ", script, call. = FALSE)
  }
  invisible(TRUE)
}

combine_chunk_outputs <- function(input_dir, output_dir) {
  run_rscript(
    file.path("scripts", "combine_hpc_per_mag_chunks.R"),
    c("--input-dir", input_dir, "--output-dir", output_dir)
  )
}

copy_if_needed <- function(src, dest) {
  if (!file.exists(src)) {
    stop("Missing required file: ", src, call. = FALSE)
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(src, dest, overwrite = TRUE)) {
    stop("Failed to copy ", src, " to ", dest, call. = FALSE)
  }
  invisible(dest)
}

write_session_info <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = path)
  invisible(path)
}

ensure_root_metadata <- function(result_dir, results_csv, summary_log_file) {
  canonical_results <- file.path(result_dir, "tables", "per_mag_lmertest_hpc_results.csv")
  canonical_summary <- file.path(result_dir, "logs", "per_mag_lmertest_hpc_summary.txt")
  canonical_session <- file.path(result_dir, "session", "session_info.txt")
  if (file.exists(results_csv) && normalizePath(results_csv, winslash = "/", mustWork = FALSE) != normalizePath(canonical_results, winslash = "/", mustWork = FALSE)) {
    copy_if_needed(results_csv, canonical_results)
  }
  if (file.exists(summary_log_file) && normalizePath(summary_log_file, winslash = "/", mustWork = FALSE) != normalizePath(canonical_summary, winslash = "/", mustWork = FALSE)) {
    copy_if_needed(summary_log_file, canonical_summary)
  }

  write_session_info(canonical_session)

  execution_record_helper <- file.path("scripts", "utils", "write_hpc_execution_record.R")
  if (file.exists(execution_record_helper)) {
    tryCatch(source(execution_record_helper), error = function(e) warning("Could not source execution recorder: ", conditionMessage(e), call. = FALSE))
  }
  if (exists("write_hpc_execution_record", mode = "function")) {
    write_hpc_execution_record(
      analysis_name = "full_T1-T10_per_MAG_lmertest_plot_postprocess",
      script_path = "scripts/hpc_plot_per_mag_lmertest_results.R",
      result_dir = result_dir,
      input_data_file = results_csv,
      config_files = c("config/model_decision_rules.yml"),
      number_of_tests = NA,
      fdr_family = "Benjamini-Hochberg across all modeled MAGs in the full per-MAG LMM run",
      main_results_file = canonical_results,
      summary_log_file = canonical_summary,
      session_info_file = canonical_session,
      notes = c("Post-processing plotting stage for the per-MAG LMM workflow"),
      start_time = Sys.time()
    )
  }

  invisible(list(
    results_csv = canonical_results,
    summary_log = canonical_summary,
    session_info = canonical_session
  ))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

input_dir <- args$input_dir
output_dir <- args$output_dir
tables_dir <- file.path(output_dir, "tables")
logs_dir <- file.path(output_dir, "logs")
figures_dir <- file.path(output_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

combined_csv <- file.path(tables_dir, "per_mag_lmertest_hpc_results_combined.csv")
canonical_csv <- file.path(tables_dir, "per_mag_lmertest_hpc_results.csv")
combined_summary <- file.path(logs_dir, "per_mag_lmertest_hpc_combined_summary.txt")

resolve_results_csv <- function(requested_csv) {
  candidates <- unique(stats::na.omit(c(requested_csv, canonical_csv, combined_csv)))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }

  if (!args$combine_chunks) {
    return(character())
  }

  chunk_dir <- file.path(input_dir, "chunks")
  chunk_files <- list.files(chunk_dir, pattern = "^per_mag_lmertest_hpc_results\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(chunk_files) == 0) {
    return(character())
  }

  combine_chunk_outputs(input_dir, output_dir)

  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }

  character()
}

if (!nzchar(args$results_csv)) {
  args$results_csv <- resolve_results_csv(character())
} else {
  resolved_csv <- resolve_results_csv(args$results_csv)
  if (!nzchar(resolved_csv)) {
    stop("Requested per-MAG result CSV does not exist and could not be rebuilt: ", args$results_csv, call. = FALSE)
  }
  args$results_csv <- resolved_csv
}

if (!nzchar(args$results_csv)) {
  stop("No per-MAG result CSV found in ", input_dir, call. = FALSE)
}

if (args$write_root_metadata) {
  ensure_root_metadata(output_dir, args$results_csv, if (file.exists(combined_summary)) combined_summary else file.path(logs_dir, "per_mag_lmertest_hpc_summary.txt"))
}

plot_scripts <- c(
  "scripts/plot_figure_1_hpc_per_mag_lmertest_results.R",
  "scripts/plot_figure_2_hpc_per_mag_lmertest_results.R",
  "scripts/plot_figure_3_hpc_per_mag_lmertest_results.R",
  "scripts/plot_figure_4_hpc_per_mag_lmertest_results_public.R",
  "scripts/plot_figure_5_hpc_per_mag_lmertest_results_public.R",
  "scripts/plot_figure_6_hpc_per_mag_lmertest_results_public.R",
  "scripts/plot_figure_7_hpc_per_mag_lmertest_results_public.R",
  "scripts/plot_figure_8_hpc_per_mag_lmertest_results_public.R",
  "scripts/plot_figure_9_hpc_per_mag_lmertest_results.R",
  "scripts/plot_figure_10_hpc_per_mag_lmertest_results_public.R",
  "scripts/summarize_hpc_candidate_signals.R"
)

for (script in plot_scripts) {
  if (!file.exists(script)) {
    stop("Missing plotting script: ", script, call. = FALSE)
  }
  run_rscript(script)
}

message("Per-MAG plotting post-process complete.")
