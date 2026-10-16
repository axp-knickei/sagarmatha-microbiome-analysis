# ============================================
# hpc_run_class1_window_lmertest.R
# ============================================
# HPC-ready and locally testable Class 1 predefined-window per-MAG LMM
# workflow. This is a secondary stage-local analysis using windows defined in
# results/metadata_expedition_windows/proposed_biological_windows.csv.
#
# This script does not redefine windows and does not run exploratory all-window
# scanning.

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
    windows_csv = "results/metadata_expedition_windows/proposed_biological_windows.csv",
    output_dir = "results/class1_window_lmertest",
    config = "config/model_decision_rules.yml",
    test = FALSE,
    max_mags = NA_integer_,
    max_windows = NA_integer_
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      cat(
        "Usage: Rscript scripts/hpc_run_class1_window_lmertest.R [options]\n\n",
        "Options:\n",
        "  --test                Run local validation mode: first 2 windows and first 5 MAGs.\n",
        "  --input PATH          Processed RDS input. Default: data/processed/cleaned_data.rds\n",
        "  --windows-csv PATH    Class 1 window definition CSV. Default: results/metadata_expedition_windows/proposed_biological_windows.csv\n",
        "  --output-dir PATH     Output directory. Default: results/class1_window_lmertest\n",
        "  --config PATH         Decision-rule YAML. Default: config/model_decision_rules.yml\n",
        "  --max-mags N          Limit to first N MAGs.\n",
        "  --max-windows N       Limit to first N predefined windows.\n",
        "  --help                Show this help.\n",
        sep = ""
      )
      quit(status = 0)
    } else if (arg == "--test") {
      parsed$test <- TRUE
      i <- i + 1
    } else if (arg %in% c("--input", "--windows-csv", "--output-dir", "--config", "--max-mags", "--max-windows")) {
      if (i == length(args)) {
        stop("Missing value for argument: ", arg, call. = FALSE)
      }
      value <- args[[i + 1]]
      if (arg == "--input") {
        parsed$input <- value
      } else if (arg == "--windows-csv") {
        parsed$windows_csv <- value
      } else if (arg == "--output-dir") {
        parsed$output_dir <- value
      } else if (arg == "--config") {
        parsed$config <- value
      } else if (arg == "--max-mags") {
        parsed$max_mags <- as.integer(value)
        if (is.na(parsed$max_mags) || parsed$max_mags < 1) {
          stop("--max-mags must be a positive integer.", call. = FALSE)
        }
      } else if (arg == "--max-windows") {
        parsed$max_windows <- as.integer(value)
        if (is.na(parsed$max_windows) || parsed$max_windows < 1) {
          stop("--max-windows must be a positive integer.", call. = FALSE)
        }
      }
      i <- i + 2
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }

  if (isTRUE(parsed$test)) {
    parsed$max_mags <- if (is.na(parsed$max_mags)) 5L else min(parsed$max_mags, 5L)
    parsed$max_windows <- if (is.na(parsed$max_windows)) 2L else min(parsed$max_windows, 2L)
  }

  parsed$output_dir <- sub("/+$", "", parsed$output_dir)
  parsed
}

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
  if (isTRUE(result$singular_fit)) {
    return("lmm_singular")
  }
  "lmm_ok"
}

direction_for_estimate <- function(estimate) {
  if (!is.finite(estimate) || is.na(estimate)) {
    return("flat_or_zero")
  }
  if (estimate > 0) {
    return("increasing")
  }
  if (estimate < 0) {
    return("decreasing")
  }
  "flat_or_zero"
}

split_timepoints <- function(x) {
  parts <- unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE)
  trimws(parts)
}

extract_window_order_effect <- function(fit) {
  coef_table <- summary(fit)$coefficients
  if (!"window_order" %in% rownames(coef_table)) {
    return(NULL)
  }

  statistic_column <- intersect(c("t value", "z value"), colnames(coef_table))
  p_column <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(coef_table))

  list(
    estimate = coef_table["window_order", "Estimate"],
    std_error = coef_table["window_order", "Std. Error"],
    df = if ("df" %in% colnames(coef_table)) coef_table["window_order", "df"] else NA_real_,
    statistic = if (length(statistic_column) > 0) coef_table["window_order", statistic_column[[1]]] else NA_real_,
    p_value = if (length(p_column) > 0) coef_table["window_order", p_column[[1]]] else NA_real_
  )
}

complete_subject_count <- function(df, timepoints) {
  if (nrow(df) == 0) {
    return(0L)
  }
  sample_level <- unique(df[, c("Subject", "Timepoint"), drop = FALSE])
  sample_level$Timepoint <- as.character(sample_level$Timepoint)
  tab <- table(sample_level$Subject, sample_level$Timepoint)
  if (!all(timepoints %in% colnames(tab))) {
    return(0L)
  }
  sum(rowSums(tab[, timepoints, drop = FALSE] > 0) == length(timepoints))
}

fit_mag_window <- function(df, window_row, mag_id, index, total, pseudocount, thresholds, p_value_method) {
  timepoints <- split_timepoints(window_row$timepoints)
  if (length(timepoints) != 3) {
    stop("Window does not contain exactly three timepoints: ", window_row$window_label, call. = FALSE)
  }

  mag_window_df <- df[df$MAG == mag_id & as.character(df$Timepoint) %in% timepoints, , drop = FALSE]
  mag_window_df <- mag_window_df[stats::complete.cases(mag_window_df[, c("Y", "Subject", "Timepoint")]), , drop = FALSE]
  mag_window_df$Timepoint <- as.character(mag_window_df$Timepoint)
  mag_window_df$window_order <- match(mag_window_df$Timepoint, timepoints) - 1L
  mag_window_df$Subject <- droplevels(factor(mag_window_df$Subject))
  mag_window_df$response_value <- log10(mag_window_df$Y + pseudocount)

  nonzero <- mag_window_df$Y != 0 & !is.na(mag_window_df$Y)
  subject_nonzero <- if (nrow(mag_window_df) > 0) tapply(nonzero, mag_window_df$Subject, any) else logical(0)
  response_variance <- stats::var(mag_window_df$response_value, na.rm = TRUE)
  complete_subjects <- complete_subject_count(mag_window_df, timepoints)

  too_sparse <- sum(nonzero, na.rm = TRUE) < thresholds$too_sparse_nonzero_samples ||
    sum(subject_nonzero, na.rm = TRUE) < thresholds$too_sparse_nonzero_subjects
  low_variance <- is.na(response_variance) || response_variance < thresholds$low_variance_threshold

  result <- data.frame(
    window_label = window_row$window_label,
    forward_timepoint_path = window_row$forward_timepoint_path,
    timepoints = paste(timepoints, collapse = ","),
    evidence_type = window_row$evidence_type %||% NA_character_,
    evidence_type_explanation = window_row$evidence_type_explanation %||% NA_character_,
    MAG = mag_id,
    n_samples_window = length(unique(mag_window_df$sampleID)),
    n_subjects_window = length(unique(mag_window_df$Subject)),
    n_complete_subjects_all_three_timepoints = complete_subjects,
    zero_fraction_window = mean(mag_window_df$Y == 0, na.rm = TRUE),
    nonzero_samples_window = sum(nonzero, na.rm = TRUE),
    subjects_with_nonzero_window = sum(subject_nonzero, na.rm = TRUE),
    pseudocount = pseudocount,
    response_variance = response_variance,
    fit_succeeded = FALSE,
    required_output_extracted = FALSE,
    singular_fit = NA,
    random_intercept_variance = NA_real_,
    estimate_window_order = NA_real_,
    std_error = NA_real_,
    df = NA_real_,
    statistic = NA_real_,
    raw_p_value = NA_real_,
    q_global_class1_bh = NA_real_,
    q_descriptive_within_window_bh = NA_real_,
    p_value_method = p_value_method,
    warnings = "",
    errors = "",
    too_sparse = too_sparse,
    low_variance = low_variance,
    status = "",
    primary_analysis = FALSE,
    primary_significant = FALSE,
    nominal_p_lt_0_05 = FALSE,
    direction = "flat_or_zero",
    stringsAsFactors = FALSE
  )

  warning_messages <- character(0)
  error_messages <- character(0)

  cat(sprintf(
    "[%d/%d] Modeling window: %s | MAG: %s\n",
    index, total, window_row$window_label, mag_id
  ))

  if (!too_sparse && !low_variance) {
    fit <- tryCatch(
      withCallingHandlers(
        lmerTest::lmer(
          response_value ~ window_order + (1 | Subject),
          data = mag_window_df,
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
      result$singular_fit <- lme4::isSingular(fit)

      vc <- as.data.frame(lme4::VarCorr(fit))
      subject_intercept <- subset(vc, grp == "Subject" & var1 == "(Intercept)" & is.na(var2))
      if (nrow(subject_intercept) > 0) {
        result$random_intercept_variance <- subject_intercept$vcov[[1]]
      }

      window_effect <- extract_window_order_effect(fit)
      if (!is.null(window_effect)) {
        result$estimate_window_order <- window_effect$estimate
        result$std_error <- window_effect$std_error
        result$df <- window_effect$df
        result$statistic <- window_effect$statistic
        result$raw_p_value <- window_effect$p_value
        result$required_output_extracted <- is.finite(result$estimate_window_order) &&
          is.finite(result$std_error) &&
          !is.na(result$raw_p_value)
      }
    }
  }

  result$warnings <- paste(unique(warning_messages), collapse = " | ")
  result$errors <- paste(unique(error_messages), collapse = " | ")
  result$status <- status_for_result(result)
  result$primary_analysis <- result$status == "lmm_ok"
  result$direction <- direction_for_estimate(result$estimate_window_order)

  cat("  Status:", result$status, "\n")
  if (nzchar(result$warnings)) cat("  Warnings:", result$warnings, "\n")
  if (nzchar(result$errors)) cat("  Errors:", result$errors, "\n")

  result
}

safe_write_method_note <- function(path) {
  note <- c(
    "# Class 1 Window lmerTest Method Note",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "## Purpose",
    "",
    "This secondary workflow tests predefined Class 1 expedition windows for stage-local MAG abundance associations. The windows are imported from `results/metadata_expedition_windows/proposed_biological_windows.csv` and are not redefined by this script.",
    "",
    "## Overlapping Windows",
    "",
    "Class 1 windows partially overlap by design because expedition stages are biologically related and adjacent. Overlap is allowed, but the resulting tests are not independent.",
    "",
    "Because predefined windows partially overlap, window-level tests are interpreted as related secondary stage-local analyses rather than independent biological discoveries. Benjamini-Hochberg FDR correction is applied globally across all predefined Class 1 MAG-window tests.",
    "",
    "## Window Order",
    "",
    "`window_order` is coded within each three-timepoint window as 0, 1, and 2. The fixed effect estimates whether a MAG tends to increase or decrease across that predefined biological window.",
    "",
    "## Why Not Primary days_from_baseline Inside Windows",
    "",
    "The primary full-run model uses expedition time across T1-T10. In the Class 1 analysis, the goal is a stage-local contrast across a predefined three-sample biological window. `window_order` emphasizes the ordered stage-local progression and avoids treating irregular date gaps as the only primary exposure within these short windows.",
    "",
    "## FDR",
    "",
    "The main FDR family is all eligible MAG-window tests across all predefined Class 1 windows. Per-window q-values, if present, are descriptive only and are not the primary significance claim.",
    "",
    "## Claims",
    "",
    "Allowed language includes stage-local association, predefined biological window, secondary analysis, exploratory candidate, MAG-window signal, and hypothesis-generating result. Do not claim causality, mediation, or a summit MAG signature from this workflow alone.",
    "",
    "## Limitations",
    "",
    "This is a secondary stage-local analysis, not final biological proof. Interpretation requires trajectory review, diagnostics, taxonomy/biology review, response-scale sensitivity checks, and careful handling of the distinction between sampling altitude and exposure altitude trajectory."
  )
  writeLines(note, path)
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

too_sparse_nonzero_samples <- as.integer(get_nested(decision_rules, c("thresholds", "too_sparse_nonzero_samples"), 10))
too_sparse_nonzero_subjects <- as.integer(get_nested(decision_rules, c("thresholds", "too_sparse_nonzero_subjects"), 3))
low_variance_threshold <- as.numeric(get_nested(decision_rules, c("thresholds", "low_variance_threshold"), 1e-10))
p_value_method <- get_nested(decision_rules, c("p_value_method"), "lmerTest_Satterthwaite")
fdr_method <- get_nested(decision_rules, c("multiple_testing"), "Benjamini-Hochberg FDR")
pseudocount_rule <- as.character(get_nested(decision_rules, c("pseudocount", "rule"), "half_smallest_nonzero_Y"))

thresholds <- list(
  too_sparse_nonzero_samples = too_sparse_nonzero_samples,
  too_sparse_nonzero_subjects = too_sparse_nonzero_subjects,
  low_variance_threshold = low_variance_threshold
)

tables_dir <- file.path(args$output_dir, "tables")
logs_dir <- file.path(args$output_dir, "logs")
candidates_dir <- file.path(args$output_dir, "candidates")
figures_dir <- file.path(args$output_dir, "figures")
session_dir <- file.path(args$output_dir, "session")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(candidates_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)

results_csv <- file.path(tables_dir, "class1_window_lmertest_results.csv")
primary_csv <- file.path(tables_dir, "class1_window_lmertest_results_primary_only.csv")
nominal_csv <- file.path(tables_dir, "class1_window_nominal_candidates_raw_p_lt_0_05.csv")
fdr_csv <- file.path(tables_dir, "class1_window_fdr_significant_candidates.csv")
direction_csv <- file.path(tables_dir, "class1_window_direction_summary.csv")
summary_path <- file.path(logs_dir, "class1_window_lmertest_summary.txt")
log_path <- file.path(logs_dir, "class1_window_lmertest_log.txt")
session_path <- file.path(session_dir, "session_info.txt")
method_note_path <- file.path(args$output_dir, "class1_window_lmertest_method_note.md")

log_connection <- file(log_path, open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")

cat("Class 1 Predefined-Window lmerTest LMM Workflow\n")
cat("===============================================\n\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Run mode:", ifelse(args$test, "test", "full"), "\n")
cat("Input file:", args$input, "\n")
cat("Window definition file:", args$windows_csv, "\n")
cat("Output directory:", args$output_dir, "\n")
cat("Decision-rule config:", args$config, "\n")
cat("Formula: log10(Y + pseudocount) ~ window_order + (1 | Subject)\n")
cat("P-value method:", p_value_method, "\n")
cat("FDR method:", fdr_method, "\n\n")

if (!file.exists(args$input)) {
  stop("Input file not found: ", args$input, call. = FALSE)
}
if (!file.exists(args$windows_csv)) {
  stop("Window definition CSV not found: ", args$windows_csv, call. = FALSE)
}
if (!identical(pseudocount_rule, "half_smallest_nonzero_Y")) {
  stop("Unsupported pseudocount rule: ", pseudocount_rule, call. = FALSE)
}

df <- readRDS(args$input)
windows <- utils::read.csv(args$windows_csv, stringsAsFactors = FALSE, check.names = FALSE)

required_df_cols <- c("MAG", "Y", "Subject", "sampleID", "Timepoint")
missing_df_cols <- setdiff(required_df_cols, names(df))
if (length(missing_df_cols) > 0) {
  stop("Input data is missing required columns: ", paste(missing_df_cols, collapse = ", "), call. = FALSE)
}

required_window_cols <- c("window_label", "timepoints", "forward_timepoint_path", "evidence_type")
missing_window_cols <- setdiff(required_window_cols, names(windows))
if (length(missing_window_cols) > 0) {
  stop("Window CSV is missing required columns: ", paste(missing_window_cols, collapse = ", "), call. = FALSE)
}
if (!"evidence_type_explanation" %in% names(windows)) {
  windows$evidence_type_explanation <- ""
}

window_lengths <- vapply(windows$timepoints, function(x) length(split_timepoints(x)), integer(1))
if (any(window_lengths != 3)) {
  stop("All Class 1 windows must contain exactly three timepoints.", call. = FALSE)
}

total_predefined_windows <- nrow(windows)
selected_windows <- windows
if (!is.na(args$max_windows)) {
  selected_windows <- utils::head(selected_windows, args$max_windows)
}

mag_ids <- sort(unique(df$MAG))
total_available_mags <- length(mag_ids)
if (!is.na(args$max_mags)) {
  mag_ids <- utils::head(mag_ids, args$max_mags)
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
cat("Predefined windows available:", total_predefined_windows, "\n")
cat("Windows selected:", nrow(selected_windows), "\n")
cat("Total MAGs available:", total_available_mags, "\n")
cat("MAGs selected:", length(mag_ids), "\n\n")

if (nrow(selected_windows) == 0 || length(mag_ids) == 0) {
  stop("No windows or MAGs selected for modeling.", call. = FALSE)
}

task_grid <- expand.grid(
  window_index = seq_len(nrow(selected_windows)),
  MAG = mag_ids,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_tests <- nrow(task_grid)

results_list <- vector("list", total_tests)
for (i in seq_len(total_tests)) {
  window_row <- selected_windows[task_grid$window_index[[i]], , drop = FALSE]
  results_list[[i]] <- fit_mag_window(
    df = df,
    window_row = window_row,
    mag_id = task_grid$MAG[[i]],
    index = i,
    total = total_tests,
    pseudocount = pseudocount,
    thresholds = thresholds,
    p_value_method = p_value_method
  )
}

results <- do.call(rbind, results_list)
results$primary_analysis <- results$status == "lmm_ok"
results$nominal_p_lt_0_05 <- results$raw_p_value < 0.05
results$q_global_class1_bh <- NA_real_
results$q_descriptive_within_window_bh <- NA_real_
results$primary_significant <- FALSE

if (any(results$primary_analysis, na.rm = TRUE)) {
  results$q_global_class1_bh[results$primary_analysis] <- stats::p.adjust(
    results$raw_p_value[results$primary_analysis],
    method = "BH"
  )
  results$primary_significant <- results$primary_analysis & results$q_global_class1_bh < 0.05

  for (window_label in unique(results$window_label)) {
    idx <- results$window_label == window_label & results$primary_analysis
    if (any(idx, na.rm = TRUE)) {
      results$q_descriptive_within_window_bh[idx] <- stats::p.adjust(results$raw_p_value[idx], method = "BH")
    }
  }
}

primary_results <- results[results$primary_analysis, , drop = FALSE]
nominal_candidates <- results[results$primary_analysis & results$raw_p_value < 0.05, , drop = FALSE]
nominal_candidates <- nominal_candidates[order(nominal_candidates$raw_p_value), , drop = FALSE]
fdr_candidates <- results[results$primary_significant, , drop = FALSE]
fdr_candidates <- fdr_candidates[order(fdr_candidates$q_global_class1_bh, fdr_candidates$raw_p_value), , drop = FALSE]

direction_summary <- as.data.frame.matrix(table(
  factor(results$direction, levels = c("increasing", "decreasing", "flat_or_zero")),
  results$status
))
direction_summary$direction <- rownames(direction_summary)
rownames(direction_summary) <- NULL
direction_summary <- direction_summary[, c("direction", setdiff(names(direction_summary), "direction")), drop = FALSE]

nominal_direction_summary <- as.data.frame.matrix(table(
  factor(nominal_candidates$direction, levels = c("increasing", "decreasing", "flat_or_zero")),
  nominal_candidates$window_label
))
nominal_direction_summary$direction <- rownames(nominal_direction_summary)
rownames(nominal_direction_summary) <- NULL
direction_summary$nominal_raw_p_lt_0_05_total <- vapply(
  direction_summary$direction,
  function(x) sum(nominal_candidates$direction == x, na.rm = TRUE),
  integer(1)
)

utils::write.csv(results, results_csv, row.names = FALSE)
utils::write.csv(primary_results, primary_csv, row.names = FALSE)
utils::write.csv(nominal_candidates, nominal_csv, row.names = FALSE)
utils::write.csv(fdr_candidates, fdr_csv, row.names = FALSE)
utils::write.csv(direction_summary, direction_csv, row.names = FALSE)
safe_write_method_note(method_note_path)

plot_warning <- NULL
tryCatch({
  if (requireNamespace("ggplot2", quietly = TRUE) && nrow(nominal_candidates) > 0) {
    top_candidates <- utils::head(nominal_candidates, 20)
    top_keys <- paste(top_candidates$window_label, top_candidates$MAG, sep = "___")
    plot_df <- data.frame()
    for (j in seq_len(nrow(top_candidates))) {
      tps <- split_timepoints(top_candidates$timepoints[[j]])
      subset_df <- df[df$MAG == top_candidates$MAG[[j]] & as.character(df$Timepoint) %in% tps, , drop = FALSE]
      subset_df$window_label <- top_candidates$window_label[[j]]
      subset_df$window_order <- match(as.character(subset_df$Timepoint), tps) - 1L
      subset_df$response_value <- log10(subset_df$Y + pseudocount)
      subset_df$MAG_window <- paste(top_candidates$window_label[[j]], top_candidates$MAG[[j]], sep = "\n")
      plot_df <- rbind(plot_df, subset_df)
    }
    if (nrow(plot_df) > 0) {
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = window_order, y = response_value, group = Subject)) +
        ggplot2::geom_line(alpha = 0.35) +
        ggplot2::geom_point(alpha = 0.65, size = 0.9) +
        ggplot2::stat_summary(ggplot2::aes(group = 1), fun = mean, geom = "line", linewidth = 0.8, color = "red") +
        ggplot2::facet_wrap(~MAG_window, scales = "free_y") +
        ggplot2::scale_x_continuous(breaks = c(0, 1, 2)) +
        ggplot2::labs(
          title = "Top nominal Class 1 MAG-window trajectories",
          x = "Window order",
          y = "log10(Y + pseudocount)"
        ) +
        ggplot2::theme_bw(base_size = 9)
      ggplot2::ggsave(
        file.path(figures_dir, "top_class1_window_candidate_trajectories.png"),
        p,
        width = 12,
        height = 10,
        dpi = 300
      )
      ggplot2::ggsave(
        file.path(figures_dir, "top_class1_window_candidate_trajectories.svg"),
        p,
        width = 12,
        height = 10,
        device = svglite::svglite
      )
    }
  }
}, error = function(e) {
  plot_warning <<- conditionMessage(e)
})

cat("\nFinished modeling:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Saved results to:", results_csv, "\n")
if (!is.null(plot_warning)) {
  cat("Plot warning:", plot_warning, "\n")
}

sink(type = "message")
sink()
close(log_connection)

summary_connection <- file(summary_path, open = "wt")
sink(summary_connection, split = TRUE)

cat("Class 1 Predefined-Window lmerTest LMM Summary\n")
cat("==============================================\n\n")
cat("Date/time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Run mode:", ifelse(args$test, "test", "full"), "\n")
cat("Input file:", args$input, "\n")
cat("Window definition file:", args$windows_csv, "\n")
cat("Output directory:", args$output_dir, "\n")
cat("Decision-rule config:", args$config, "\n")
cat("Number of predefined windows available:", total_predefined_windows, "\n")
cat("Number of predefined windows modeled:", nrow(selected_windows), "\n\n")

cat("Window labels and paths:\n")
for (i in seq_len(nrow(selected_windows))) {
  cat(" - ", selected_windows$window_label[[i]], ": ", selected_windows$forward_timepoint_path[[i]], " [", selected_windows$evidence_type[[i]], "]\n", sep = "")
}
cat("\n")

cat("Total MAGs available:", total_available_mags, "\n")
cat("Total MAGs modeled:", length(mag_ids), "\n")
cat("Total MAG-window tests attempted:", nrow(results), "\n")
cat("Successful models:", sum(results$fit_succeeded %in% TRUE, na.rm = TRUE), "\n")
cat("Failed models:", sum(results$status == "model_failed", na.rm = TRUE), "\n\n")

cat("Status counts:\n")
status_counts <- table(results$status)
status_count <- function(status_name) {
  if (status_name %in% names(status_counts)) {
    return(unname(status_counts[[status_name]]))
  }
  0L
}
for (status_name in c("lmm_ok", "lmm_singular", "too_sparse", "low_variance", "model_failed")) {
  cat(status_name, ":", status_count(status_name), "\n", sep = "")
}
cat("\n")

cat("Significance counts:\n")
cat("Raw p < 0.05 among primary-analysis tests:", sum(results$primary_analysis & results$raw_p_value < 0.05, na.rm = TRUE), "\n")
cat("q_global_class1_bh < 0.05:", sum(results$q_global_class1_bh < 0.05, na.rm = TRUE), "\n")
cat("Primary significant MAG-window tests:", sum(results$primary_significant %in% TRUE, na.rm = TRUE), "\n")
best_q <- suppressWarnings(min(results$q_global_class1_bh, na.rm = TRUE))
cat("Best global Class 1 q-value:", ifelse(is.finite(best_q), format(best_q, digits = 6), "NA"), "\n\n")

cat("Direction summary among nominal candidates:\n")
if (nrow(nominal_candidates) == 0) {
  cat("No nominal raw p < 0.05 primary-analysis candidates.\n")
} else {
  print(as.data.frame(table(nominal_candidates$direction)), row.names = FALSE)
}
cat("\n")

cat("Top 20 primary-analysis MAG-window tests ranked by global Class 1 BH q-value, then raw p-value:\n")
ranked <- primary_results[order(
  is.na(primary_results$q_global_class1_bh),
  primary_results$q_global_class1_bh,
  primary_results$raw_p_value
), ]
top_cols <- c(
  "window_label", "forward_timepoint_path", "evidence_type", "MAG",
  "status", "primary_significant", "direction", "estimate_window_order",
  "std_error", "df", "statistic", "raw_p_value", "q_global_class1_bh",
  "q_descriptive_within_window_bh", "zero_fraction_window",
  "nonzero_samples_window", "subjects_with_nonzero_window",
  "n_complete_subjects_all_three_timepoints"
)
if (nrow(ranked) > 0) {
  print(utils::head(ranked[, top_cols], 20), row.names = FALSE)
} else {
  cat("No primary-analysis rows available.\n")
}

cat("\nInterpretation principle:\n")
cat("Because predefined windows partially overlap, window-level tests are interpreted as related secondary stage-local analyses rather than independent biological discoveries. Benjamini-Hochberg FDR correction is applied globally across all predefined Class 1 MAG-window tests.\n\n")

cat("Warnings:\n")
cat("- Partially overlapping windows are related secondary stage-local tests, not independent biological discoveries.\n")
cat("- Sampling altitude and exposure altitude trajectory differ; T7 -> T8 summit-push exposure is figure-informed and not directly encoded as a metadata sampling altitude.\n")
cat("- Missing subjects and incomplete three-timepoint windows can reduce LMM stability; review n_complete_subjects_all_three_timepoints before interpretation.\n")
cat("- Raw p < 0.05 rows are exploratory candidates only unless q_global_class1_bh < 0.05 and diagnostics support interpretation.\n")
if (!is.null(plot_warning)) {
  cat("- Plotting warning:", plot_warning, "\n")
}

cat("\nSaved outputs:\n")
cat("- ", results_csv, "\n", sep = "")
cat("- ", primary_csv, "\n", sep = "")
cat("- ", nominal_csv, "\n", sep = "")
cat("- ", fdr_csv, "\n", sep = "")
cat("- ", direction_csv, "\n", sep = "")
cat("- ", method_note_path, "\n", sep = "")
cat("- ", session_path, "\n", sep = "")

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

cat("Class 1 predefined-window LMM complete\n")
cat("======================================\n")
cat("Run mode:", ifelse(args$test, "test", "full"), "\n")
cat("Files created:\n")
cat(" - ", results_csv, "\n", sep = "")
cat(" - ", primary_csv, "\n", sep = "")
cat(" - ", nominal_csv, "\n", sep = "")
cat(" - ", fdr_csv, "\n", sep = "")
cat(" - ", direction_csv, "\n", sep = "")
cat(" - ", summary_path, "\n", sep = "")
cat(" - ", method_note_path, "\n", sep = "")
cat(" - ", session_path, "\n", sep = "")
cat("Test windows:", nrow(selected_windows), "\n")
cat("Test MAGs:", length(mag_ids), "\n")
cat("Model fitting succeeded:", sum(results$fit_succeeded %in% TRUE, na.rm = TRUE), "of", nrow(results), "attempted tests\n")
cat("FDR columns exist:", all(c("q_global_class1_bh", "q_descriptive_within_window_bh") %in% names(results)), "\n")
cat("Warnings for review:\n")
cat(" - Predefined windows partially overlap and are not independent.\n")
cat(" - Sampling altitude is not the full exposure altitude trajectory.\n")
cat(" - Review sparse/low-variance/singular statuses before interpretation.\n")

tryCatch(
  {
    if (exists("write_hpc_execution_record", mode = "function")) {
      total_attempted_tests <- nrow(results)
      primary_eligible_tests <- sum(results$primary_analysis %in% TRUE, na.rm = TRUE)
      config_files <- c(args$config)
      config_files <- c(config_files, c("config/response_scale_decision.yml", "config/config.yml")[file.exists(c("config/response_scale_decision.yml", "config/config.yml"))])

      write_hpc_execution_record(
        analysis_name = "class1_predefined_window_per_MAG_lmertest_LMM",
        script_path = "scripts/hpc_run_class1_window_lmertest.R",
        result_dir = args$output_dir,
        input_data_file = args$input,
        window_definition_file = args$windows_csv,
        config_files = config_files,
        number_of_tests = primary_eligible_tests,
        fdr_family = "Benjamini-Hochberg globally across all eligible predefined Class 1 MAG-window tests",
        main_results_file = results_csv,
        summary_log_file = summary_path,
        session_info_file = session_path,
        notes = c(
          paste0("Total attempted MAG-window rows: ", total_attempted_tests),
          paste0("Eligible primary Class 1 MAG-window tests: ", primary_eligible_tests),
          paste0("Run mode: ", ifelse(args$test, "test", "full"))
        ),
        start_time = execution_start_time
      )
    }
  },
  error = function(e) warning("Failed to write HPC execution record; continuing analysis: ", conditionMessage(e), call. = FALSE)
)
