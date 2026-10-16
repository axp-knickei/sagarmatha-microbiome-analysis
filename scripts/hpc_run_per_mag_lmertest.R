# ============================================
# hpc_run_per_mag_lmertest.R
# ============================================
# HPC-ready per-MAG LMM workflow using the primary working response scale:
# log10(Y + pseudocount), where pseudocount is half of the smallest non-zero
# abundance in the processed dataset.
#
# This script writes to results/hpc_per_mag_lmertest/ by default and does not
# overwrite the exploratory per-MAG result files.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    return(y)
  }
  if (is.atomic(x) && length(x) == 1 && (is.na(x) || identical(x, ""))) {
    return(y)
  }
  x
}

get_nested <- function(x, path, default = NULL) {
  current <- x
  for (key in path) {
    if (is.null(current) || is.null(current[[key]])) {
      return(default)
    }
    current <- current[[key]]
  }
  current
}

parse_args <- function(args) {
  parsed <- list(
    input = "data/processed/cleaned_data.rds",
    output_dir = "results/hpc_per_mag_lmertest",
    config = "config/model_decision_rules.yml",
    max_mags = NA_integer_,
    chunk_count = 1L,
    chunk_index = 1L
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      cat(
        "Usage: Rscript scripts/hpc_run_per_mag_lmertest.R [options]\n\n",
        "Options:\n",
        "  --input PATH          Processed RDS input. Default: data/processed/cleaned_data.rds\n",
        "  --output-dir PATH     Output directory. Default: results/hpc_per_mag_lmertest\n",
        "  --config PATH         Decision-rule YAML. Default: config/model_decision_rules.yml\n",
        "  --max-mags N          Limit to first N MAGs for dry runs.\n",
        "  --chunk-count N       Split selected MAGs into N chunks.\n",
        "  --chunk-index N       Run 1-based chunk N. Use with --chunk-count.\n",
        "  --help                Show this help.\n",
        sep = ""
      )
      quit(status = 0)
    } else if (arg %in% c("--input", "--output-dir", "--config", "--max-mags", "--chunk-count", "--chunk-index")) {
      if (i == length(args)) {
        stop("Missing value for argument: ", arg, call. = FALSE)
      }
      value <- args[[i + 1]]
      if (arg == "--input") {
        parsed$input <- value
      } else if (arg == "--output-dir") {
        parsed$output_dir <- value
      } else if (arg == "--config") {
        parsed$config <- value
      } else if (arg == "--max-mags") {
        parsed$max_mags <- as.integer(value)
        if (is.na(parsed$max_mags) || parsed$max_mags < 1) {
          stop("--max-mags must be a positive integer.", call. = FALSE)
        }
      } else if (arg == "--chunk-count") {
        parsed$chunk_count <- as.integer(value)
        if (is.na(parsed$chunk_count) || parsed$chunk_count < 1) {
          stop("--chunk-count must be a positive integer.", call. = FALSE)
        }
      } else if (arg == "--chunk-index") {
        parsed$chunk_index <- as.integer(value)
        if (is.na(parsed$chunk_index) || parsed$chunk_index < 1) {
          stop("--chunk-index must be a positive integer.", call. = FALSE)
        }
      }
      i <- i + 2
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }

  parsed$output_dir <- sub("/+$", "", parsed$output_dir)
  if (parsed$chunk_index > parsed$chunk_count) {
    stop("--chunk-index must be less than or equal to --chunk-count.", call. = FALSE)
  }
  parsed
}

execution_start_time <- Sys.time()
args <- parse_args(commandArgs(trailingOnly = TRUE))

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(yaml)
})

execution_record_helper <- file.path("scripts", "utils", "write_hpc_execution_record.R")
if (file.exists(execution_record_helper)) {
  tryCatch(
    source(execution_record_helper),
    error = function(e) warning("Could not source HPC execution recorder; continuing analysis: ", conditionMessage(e), call. = FALSE)
  )
}

decision_rules <- if (file.exists(args$config)) yaml::read_yaml(args$config) else list()
decision_rules <- decision_rules$model_decision_rules %||% decision_rules

chunk_mode <- args$chunk_count > 1
chunk_label <- sprintf("chunk_%03d_of_%03d", args$chunk_index, args$chunk_count)
run_output_dir <- if (chunk_mode) file.path(args$output_dir, "chunks", chunk_label) else args$output_dir

tables_dir <- file.path(run_output_dir, "tables")
logs_dir <- file.path(run_output_dir, "logs")
session_dir <- file.path(run_output_dir, "session")

output_csv <- file.path(tables_dir, "per_mag_lmertest_hpc_results.csv")
log_path <- file.path(logs_dir, "per_mag_lmertest_hpc_log.txt")
summary_path <- file.path(logs_dir, "per_mag_lmertest_hpc_summary.txt")
session_path <- file.path(session_dir, "session_info.txt")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)

too_sparse_nonzero_samples <- get_nested(decision_rules, c("thresholds", "too_sparse_nonzero_samples"), 10)
too_sparse_nonzero_subjects <- get_nested(decision_rules, c("thresholds", "too_sparse_nonzero_subjects"), 3)
low_variance_threshold <- get_nested(decision_rules, c("thresholds", "low_variance_threshold"), 1e-10)
too_sparse_nonzero_samples <- as.integer(too_sparse_nonzero_samples)
too_sparse_nonzero_subjects <- as.integer(too_sparse_nonzero_subjects)
low_variance_threshold <- as.numeric(low_variance_threshold)
p_value_method <- get_nested(decision_rules, c("p_value_method"), "lmerTest_Satterthwaite")
fdr_method <- get_nested(decision_rules, c("multiple_testing"), "Benjamini-Hochberg FDR")
formula_text <- get_nested(
  decision_rules,
  c("primary_formula"),
  "log10(Y + pseudocount) ~ days_from_baseline_subject + (1 | Subject)"
)
supported_formula_text <- "log10(Y + pseudocount) ~ days_from_baseline_subject + (1 | Subject)"
if (!identical(formula_text, supported_formula_text)) {
  stop(
    "Unsupported primary_formula in decision config: ",
    formula_text,
    ". This HPC entrypoint currently supports: ",
    supported_formula_text,
    call. = FALSE
  )
}
pseudocount_rule <- get_nested(decision_rules, c("pseudocount", "rule"), "half_smallest_nonzero_Y")
pseudocount_rule <- as.character(pseudocount_rule)

status_for_result <- function(result) {
  if (isTRUE(result$too_sparse)) {
    return("too_sparse")
  }
  if (isTRUE(result$low_variance)) {
    return("low_variance")
  }
  if (!result$fit_succeeded || !result$required_output_extracted) {
    return("model_failed")
  }
  if (isTRUE(result$singular)) {
    return("lmm_singular")
  }
  "lmm_ok"
}

extract_time_effect <- function(fit) {
  coef_table <- summary(fit)$coefficients
  if (!"days_from_baseline_subject" %in% rownames(coef_table)) {
    return(NULL)
  }

  statistic_column <- intersect(c("t value", "z value"), colnames(coef_table))
  p_column <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(coef_table))

  list(
    estimate = coef_table["days_from_baseline_subject", "Estimate"],
    std_error = coef_table["days_from_baseline_subject", "Std. Error"],
    df = if ("df" %in% colnames(coef_table)) coef_table["days_from_baseline_subject", "df"] else NA_real_,
    statistic = if (length(statistic_column) > 0) coef_table["days_from_baseline_subject", statistic_column[[1]]] else NA_real_,
    p_value = if (length(p_column) > 0) coef_table["days_from_baseline_subject", p_column[[1]]] else NA_real_
  )
}

fit_mag <- function(df, mag_id, index, total, pseudocount) {
  mag_df <- df[df$MAG == mag_id, , drop = FALSE]
  mag_df <- mag_df[stats::complete.cases(mag_df[, c("Y", "days_from_baseline_subject", "Subject")]), , drop = FALSE]
  mag_df$Subject <- droplevels(factor(mag_df$Subject))
  mag_df$response_value <- log10(mag_df$Y + pseudocount)

  nonzero <- mag_df$Y != 0 & !is.na(mag_df$Y)
  subject_nonzero <- tapply(nonzero, mag_df$Subject, any)
  response_variance <- stats::var(mag_df$response_value, na.rm = TRUE)

  result <- data.frame(
    MAG = mag_id,
    n_rows = nrow(mag_df),
    n_subjects = length(unique(mag_df$Subject)),
    n_sampleID = length(unique(mag_df$sampleID)),
    n_timepoints = length(unique(mag_df$Timepoint)),
    zero_fraction = mean(mag_df$Y == 0, na.rm = TRUE),
    nonzero_sample_count = sum(nonzero),
    subjects_with_nonzero = sum(subject_nonzero, na.rm = TRUE),
    response_variance = response_variance,
    fit_succeeded = FALSE,
    required_output_extracted = FALSE,
    singular = NA,
    random_intercept_variance = NA_real_,
    time_estimate = NA_real_,
    time_std_error = NA_real_,
    time_df = NA_real_,
    time_statistic = NA_real_,
    time_p_value = NA_real_,
    p_value_method = p_value_method,
    warnings = "",
    errors = "",
    too_sparse = sum(nonzero) < too_sparse_nonzero_samples ||
      sum(subject_nonzero, na.rm = TRUE) < too_sparse_nonzero_subjects,
    low_variance = is.na(response_variance) || response_variance < low_variance_threshold,
    status = "",
    primary_analysis = FALSE,
    time_p_adj_bh_all = NA_real_,
    time_p_adj_bh_primary = NA_real_,
    primary_significant = FALSE,
    nominal_p_lt_0_05 = FALSE,
    stringsAsFactors = FALSE
  )

  warning_messages <- character(0)
  error_messages <- character(0)

  cat(sprintf("[%d/%d] Modeling MAG: %s\n", index, total, mag_id))

  fit <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(
        response_value ~ days_from_baseline_subject + (1 | Subject),
        data = mag_df,
        control = lme4::lmerControl(optimizer = "bobyqa")
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_messages <<- c(error_messages, conditionMessage(e))
      NULL
    }
  )

  if (!is.null(fit)) {
    result$fit_succeeded <- TRUE
    result$singular <- lme4::isSingular(fit)

    vc <- as.data.frame(lme4::VarCorr(fit))
    subject_intercept <- subset(vc, grp == "Subject" & var1 == "(Intercept)" & is.na(var2))
    if (nrow(subject_intercept) > 0) {
      result$random_intercept_variance <- subject_intercept$vcov[[1]]
    }

    time_effect <- extract_time_effect(fit)
    if (!is.null(time_effect)) {
      result$time_estimate <- time_effect$estimate
      result$time_std_error <- time_effect$std_error
      result$time_df <- time_effect$df
      result$time_statistic <- time_effect$statistic
      result$time_p_value <- time_effect$p_value
      result$required_output_extracted <- is.finite(result$time_estimate) &&
        is.finite(result$time_std_error) &&
        !is.na(result$time_p_value)
    }
  }

  result$warnings <- paste(unique(warning_messages), collapse = " | ")
  result$errors <- paste(unique(error_messages), collapse = " | ")
  result$status <- status_for_result(result)
  result$primary_analysis <- result$status == "lmm_ok"

  cat("  Status:", result$status, "\n")
  if (nzchar(result$warnings)) cat("  Warnings:", result$warnings, "\n")
  if (nzchar(result$errors)) cat("  Errors:", result$errors, "\n")

  result
}

log_connection <- file(log_path, open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")

cat("HPC Per-MAG lmerTest LMM Workflow\n")
cat("=================================\n\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Input file:", args$input, "\n")
cat("Output directory:", run_output_dir, "\n")
cat("Decision-rule config:", args$config, "\n")
cat("Formula:", formula_text, "\n")
cat("P-value method:", p_value_method, "\n")
cat("FDR method:", fdr_method, "\n")
cat("Max MAGs:", ifelse(is.na(args$max_mags), "all", args$max_mags), "\n\n")
cat("Chunk count:", args$chunk_count, "\n")
cat("Chunk index:", args$chunk_index, "\n")
if (chunk_mode) {
  cat("Chunk mode: FDR columns are left NA in chunk outputs; run scripts/combine_hpc_per_mag_chunks.R after all chunks finish.\n")
}
cat("\n")

if (!file.exists(args$input)) {
  stop("Input file not found: ", args$input, call. = FALSE)
}

df <- readRDS(args$input)
required_cols <- c("MAG", "Y", "days_from_baseline_subject", "Subject", "sampleID", "Timepoint")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Input data is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

if (!identical(pseudocount_rule, "half_smallest_nonzero_Y")) {
  stop(
    "Unsupported pseudocount rule in decision config: ",
    pseudocount_rule,
    ". Supported rule: half_smallest_nonzero_Y.",
    call. = FALSE
  )
}

smallest_nonzero_y <- min(df$Y[df$Y > 0], na.rm = TRUE)
if (!is.finite(smallest_nonzero_y)) {
  stop("Cannot define pseudocount because no positive Y values were found.", call. = FALSE)
}
pseudocount <- smallest_nonzero_y / 2

cat("Pseudocount rule: half of smallest non-zero Y\n")
cat("Smallest non-zero Y:", format(smallest_nonzero_y, scientific = TRUE), "\n")
cat("Pseudocount:", format(pseudocount, scientific = TRUE), "\n\n")
cat("Thresholds:\n")
cat("  too_sparse_nonzero_samples:", too_sparse_nonzero_samples, "\n")
cat("  too_sparse_nonzero_subjects:", too_sparse_nonzero_subjects, "\n")
cat("  low_variance_threshold:", format(low_variance_threshold, scientific = TRUE), "\n\n")

mag_ids <- sort(unique(df$MAG))
total_available_mags <- length(mag_ids)
if (!is.na(args$max_mags)) {
  mag_ids <- utils::head(mag_ids, args$max_mags)
}
total_selected_before_chunk <- length(mag_ids)
if (chunk_mode) {
  mag_positions <- seq_along(mag_ids)
  mag_ids <- mag_ids[((mag_positions - 1L) %% args$chunk_count) + 1L == args$chunk_index]
}

cat("Total available MAGs:", total_available_mags, "\n")
cat("MAGs selected before chunking:", total_selected_before_chunk, "\n")
cat("MAGs selected for this run:", length(mag_ids), "\n\n")

if (length(mag_ids) == 0) {
  stop("No MAGs selected for this run. Check --max-mags, --chunk-count, and --chunk-index.", call. = FALSE)
}

results <- do.call(
  rbind,
  Map(fit_mag, mag_ids, seq_along(mag_ids), MoreArgs = list(
    df = df,
    total = length(mag_ids),
    pseudocount = pseudocount
  ))
)

results$time_p_adj_bh_all <- NA_real_
results$primary_analysis <- results$status == "lmm_ok"
results$time_p_adj_bh_primary <- NA_real_
results$primary_significant <- FALSE
results$nominal_p_lt_0_05 <- results$time_p_value < 0.05
results$fdr_scope <- if (chunk_mode) "not_computed_in_chunk_mode" else "single_run_all_modeled_mags"

if (!chunk_mode) {
  results$time_p_adj_bh_all <- stats::p.adjust(results$time_p_value, method = "BH")
  results$time_p_adj_bh_primary[results$primary_analysis] <- stats::p.adjust(
    results$time_p_value[results$primary_analysis],
    method = "BH"
  )
  results$primary_significant <- results$primary_analysis & results$time_p_adj_bh_primary < 0.05
}

utils::write.csv(results, output_csv, row.names = FALSE)

cat("\nFinished modeling:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Saved results to:", output_csv, "\n")

sink(type = "message")
sink()
close(log_connection)

summary_connection <- file(summary_path, open = "wt")
sink(summary_connection, split = TRUE)

cat("HPC Per-MAG lmerTest LMM Summary\n")
cat("================================\n\n")
cat("Input file:", args$input, "\n")
cat("Output directory:", run_output_dir, "\n")
cat("Decision-rule config:", args$config, "\n")
cat("Date/time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Total available MAGs:", total_available_mags, "\n")
cat("MAGs selected before chunking:", total_selected_before_chunk, "\n")
cat("MAGs modeled:", nrow(results), "\n")
cat("Run mode:", ifelse(is.na(args$max_mags), "full", paste0("max-mags ", args$max_mags)), "\n\n")
cat("Chunk mode:", chunk_mode, "\n")
cat("Chunk count:", args$chunk_count, "\n")
cat("Chunk index:", args$chunk_index, "\n\n")

cat("Pseudocount rule: half of smallest non-zero Y\n")
cat("Smallest non-zero Y:", format(smallest_nonzero_y, scientific = TRUE), "\n")
cat("Pseudocount value:", format(pseudocount, scientific = TRUE), "\n")
cat("Formula:", formula_text, "\n")
cat("P-value method:", p_value_method, "\n")
cat("FDR method:", fdr_method, "\n\n")

cat("Status counts:\n")
cat("Successful models:", sum(results$fit_succeeded %in% TRUE, na.rm = TRUE), "\n")
cat("Failed models:", sum(results$status == "model_failed", na.rm = TRUE), "\n")
cat("lmm_ok:", sum(results$status == "lmm_ok", na.rm = TRUE), "\n")
cat("lmm_singular:", sum(results$status == "lmm_singular", na.rm = TRUE), "\n")
cat("too_sparse:", sum(results$status == "too_sparse", na.rm = TRUE), "\n")
cat("low_variance:", sum(results$status == "low_variance", na.rm = TRUE), "\n\n")

cat("Significance counts:\n")
cat("Raw p < 0.05:", sum(results$nominal_p_lt_0_05 %in% TRUE, na.rm = TRUE), "\n")
if (chunk_mode) {
  cat("BH q < 0.05 across all modeled MAGs: not computed for chunk outputs\n")
  cat("BH q < 0.05 in primary-analysis set: not computed for chunk outputs\n")
  cat("Primary significant MAGs: not computed for chunk outputs\n")
  cat("Chunk note: run scripts/combine_hpc_per_mag_chunks.R after all chunks finish to compute final FDR.\n\n")
} else {
  cat("BH q < 0.05 across all modeled MAGs:", sum(results$time_p_adj_bh_all < 0.05, na.rm = TRUE), "\n")
  cat("BH q < 0.05 in primary-analysis set:", sum(results$time_p_adj_bh_primary < 0.05, na.rm = TRUE), "\n")
  cat("Primary significant MAGs:", sum(results$primary_significant %in% TRUE, na.rm = TRUE), "\n\n")
}

if (chunk_mode) {
  cat("Top 20 MAGs in this chunk ranked by raw p-value because final FDR is computed after chunk merge:\n")
  ranked <- results[order(is.na(results$time_p_value), results$time_p_value), ]
} else {
  cat("Top 20 MAGs ranked by primary BH q-value, then raw p-value:\n")
  ranked <- results[order(
    is.na(results$time_p_adj_bh_primary),
    results$time_p_adj_bh_primary,
    results$time_p_value
  ), ]
}
top_cols <- c(
  "MAG", "status", "primary_analysis", "primary_significant",
  "time_estimate", "time_std_error", "time_df", "time_statistic",
  "time_p_value", "time_p_adj_bh_all", "time_p_adj_bh_primary",
  "zero_fraction", "nonzero_sample_count", "subjects_with_nonzero",
  "response_variance", "singular"
)
print(utils::head(ranked[, top_cols], 20), row.names = FALSE)

cat("\nImportant warning:\n")
cat("This script implements the working primary response scale and model decision rules, but biological interpretation still requires completion of manual plot review, response-scale caveat review, and biological plausibility checks. Flagged statuses are retained in output but excluded from primary significance claims.\n")

cat("\nSaved summary to:", summary_path, "\n")

sink()
close(summary_connection)

session_connection <- file(session_path, open = "wt")
sink(session_connection)
cat("R Session Information\n")
cat("=====================\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
print(sessionInfo())
sink()
close(session_connection)

cat("Saved HPC per-MAG results to ", output_csv, "\n", sep = "")
cat("Saved summary to ", summary_path, "\n", sep = "")
cat("Saved session info to ", session_path, "\n", sep = "")

tryCatch(
  {
    if (exists("write_hpc_execution_record", mode = "function")) {
      total_modeled_tests <- nrow(results)
      primary_eligible_tests <- sum(results$primary_analysis %in% TRUE, na.rm = TRUE)
      config_files <- c(args$config)
      config_files <- c(config_files, c("config/response_scale_decision.yml", "config/config.yml")[file.exists(c("config/response_scale_decision.yml", "config/config.yml"))])

      write_hpc_execution_record(
        analysis_name = "full_T1_T10_per_MAG_lmertest_LMM",
        script_path = "scripts/hpc_run_per_mag_lmertest.R",
        result_dir = run_output_dir,
        input_data_file = args$input,
        config_files = config_files,
        number_of_tests = total_modeled_tests,
        fdr_family = "Benjamini-Hochberg across all modeled MAGs in the full per-MAG LMM run",
        main_results_file = output_csv,
        summary_log_file = summary_path,
        session_info_file = session_path,
        notes = c(
          paste0("Total modeled MAG tests: ", total_modeled_tests),
          paste0("Primary eligible lmm_ok tests: ", primary_eligible_tests),
          paste0("Chunk mode: ", chunk_mode),
          if (chunk_mode) "FDR columns are not computed in chunk outputs; combine chunks before final inference." else NA_character_
        ),
        start_time = execution_start_time
      )
    }
  },
  error = function(e) warning("Failed to write HPC execution record; continuing analysis: ", conditionMessage(e), call. = FALSE)
)
