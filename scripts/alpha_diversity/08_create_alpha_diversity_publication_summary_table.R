#!/usr/bin/env Rscript

# Create a publication-ready alpha-diversity summary table from the prepared
# alpha-diversity results and sample metadata.

options(stringsAsFactors = FALSE)
Sys.setenv(TZ = Sys.getenv("TZ", unset = "Asia/Jakarta"))

required_packages <- c("dplyr", "tidyr", "lme4", "knitr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them in the project/HPC environment before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  if (requireNamespace("lmerTest", quietly = TRUE)) {
    library(lmerTest)
  }
})

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

path_or_env <- function(env_name, default) {
  value <- Sys.getenv(env_name, unset = default)
  normalizePath(file.path(repo_root, value), winslash = "/", mustWork = FALSE)
}

alpha_input <- path_or_env(
  "BLACKICE_ALPHA_SUMMARY_INPUT",
  "results/alpha_diversity/tables/alpha_diversity_with_metadata.csv"
)
metadata_input <- path_or_env("BLACKICE_METADATA_PATH", "data/raw/sample_info.csv")
output_dir <- path_or_env("BLACKICE_ALPHA_OUTPUT_DIR", "results/alpha_diversity")
tables_dir <- file.path(output_dir, "tables")
logs_dir <- file.path(output_dir, "logs")
reports_dir <- file.path(repo_root, "docs", "reports")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(logs_dir, "alpha_diversity_publication_summary_table.log")
log_message <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  message(line)
  cat(line, "\n", file = log_path, append = TRUE)
}

read_csv_base <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file does not exist: ", path, call. = FALSE)
  }
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

find_first_column <- function(data, candidates) {
  matches <- candidates[candidates %in% names(data)]
  if (length(matches) == 0) NA_character_ else matches[[1]]
}

standardize_metrics <- function(data) {
  metric_map <- list(
    Observed_richness = c("Observed_richness", "observed_richness", "Observed", "observed", "Sobs", "ObservedRichness"),
    Chao1 = c("Chao1", "chao1", "Chao_1", "chao_1"),
    ACE = c("ACE", "ace"),
    Shannon = c("Shannon", "shannon", "Shannon_diversity", "shannon_diversity"),
    Simpson = c("Simpson", "simpson", "Simpson_diversity", "simpson_diversity")
  )

  for (standard_name in names(metric_map)) {
    found <- find_first_column(data, metric_map[[standard_name]])
    if (!is.na(found) && found != standard_name) {
      names(data)[names(data) == found] <- standard_name
    }
  }
  data
}

coerce_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
}

format_number <- function(x, digits = 3) {
  if (is.na(x)) return(NA_character_)
  formatC(x, digits = digits, format = "f")
}

metric_digits <- function(metric, values) {
  if (metric %in% c("Shannon", "Simpson")) return(3)
  values <- values[is.finite(values)]
  if (length(values) > 0 && all(abs(values - round(values)) < 1e-8)) 1 else 2
}

format_median_iqr <- function(values, metric) {
  values <- values[is.finite(values)]
  if (length(values) == 0) return(NA_character_)
  digits <- metric_digits(metric, values)
  q <- as.numeric(stats::quantile(values, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE))
  paste0(
    format_number(q[[2]], digits),
    " [",
    format_number(q[[1]], digits),
    ", ",
    format_number(q[[3]], digits),
    "]"
  )
}

format_min_max <- function(values, metric) {
  values <- values[is.finite(values)]
  if (length(values) == 0) return(NA_character_)
  digits <- metric_digits(metric, values)
  paste0(format_number(min(values), digits), "-", format_number(max(values), digits))
}

format_p <- function(x) {
  if (is.na(x)) return(NA_character_)
  if (x < 0.001) return("<0.001")
  formatC(signif(x, 3), digits = 3, format = "fg", flag = "#")
}

format_effect <- function(x) {
  if (is.na(x)) return(NA_character_)
  formatC(signif(x, 4), digits = 4, format = "fg", flag = "#")
}

clean_timepoint_order <- function(timepoints, data) {
  tp <- unique(as.character(timepoints[!is.na(timepoints)]))
  if (length(tp) == 0) return(character())

  canonical <- c("T1", "T2", "Ta", "T3", "Tb", paste0("T", 4:10))
  ordered <- canonical[canonical %in% tp]
  remainder <- setdiff(tp, ordered)

  if ("date_converted" %in% names(data) || "date" %in% names(data) || "days_from_baseline" %in% names(data)) {
    date_col <- find_first_column(data, c("date_converted", "date"))
    parse_dates <- function(x) {
      x_chr <- as.character(x)
      parsed <- suppressWarnings(as.Date(x_chr))
      if (all(is.na(parsed))) {
        parsed <- suppressWarnings(as.Date(x_chr, format = "%Y%m%d"))
      }
      parsed
    }
    summarise_min_date <- function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0) return(as.Date(NA))
      min(x)
    }
    summarise_median_numeric <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) return(NA_real_)
      median(x)
    }
    order_data <- data %>%
      mutate(
        .timepoint = as.character(.data$Timepoint),
        .date_value = if (!is.na(date_col)) parse_dates(.data[[date_col]]) else as.Date(NA),
        .day_value = if ("days_from_baseline" %in% names(data)) coerce_numeric(.data$days_from_baseline) else NA_real_
      ) %>%
      group_by(.data$.timepoint) %>%
      summarise(
        .date_value = summarise_min_date(.data$.date_value),
        .day_value = summarise_median_numeric(.data$.day_value),
        .groups = "drop"
      ) %>%
      arrange(.data$.date_value, .data$.day_value, .data$.timepoint)

    inferred <- order_data$.timepoint[order_data$.timepoint %in% remainder]
    remainder <- c(inferred, setdiff(remainder, inferred))
  }

  c(ordered, sort(remainder))
}

write_csv_base <- function(data, path) {
  utils::write.csv(data, path, row.names = FALSE, na = "")
}

write_tsv_base <- function(data, path) {
  utils::write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

session_info_path <- file.path(logs_dir, "alpha_diversity_publication_summary_table_session_info.txt")
cat("Alpha-diversity publication summary table session\n", file = session_info_path)
cat("Run time: ", as.character(Sys.time()), "\n\n", file = session_info_path, append = TRUE, sep = "")
suppressWarnings(capture.output(sessionInfo(), file = session_info_path, append = TRUE))

log_message("Starting alpha-diversity publication summary table workflow")
log_message("Project root: ", repo_root)
log_message("Primary alpha-diversity input: ", alpha_input)

alpha_data <- read_csv_base(alpha_input)
metadata_used <- NA_character_

if (!all(c("sampleID", "subjectID") %in% names(alpha_data)) && file.exists(metadata_input)) {
  log_message("Alpha table is missing metadata columns; attempting metadata join with: ", metadata_input)
  metadata <- read_csv_base(metadata_input)
  sample_col_alpha <- find_first_column(alpha_data, c("sampleID", "SampleID", "sample_id", "Sample_ID"))
  sample_col_meta <- find_first_column(metadata, c("sampleID", "SampleID", "sample_id", "Sample_ID"))
  if (is.na(sample_col_alpha) || is.na(sample_col_meta)) {
    stop("Cannot join metadata because a sample ID column was not detected in both tables.", call. = FALSE)
  }
  if (sample_col_alpha != "sampleID") names(alpha_data)[names(alpha_data) == sample_col_alpha] <- "sampleID"
  if (sample_col_meta != "sampleID") names(metadata)[names(metadata) == sample_col_meta] <- "sampleID"
  alpha_data <- left_join(alpha_data, metadata, by = "sampleID", suffix = c("", "_metadata"))
  metadata_used <- metadata_input
} else if (file.exists(metadata_input)) {
  metadata_used <- "metadata already present in alpha-diversity table; fallback metadata file available but not joined"
}

alpha_data <- standardize_metrics(alpha_data)

sample_col <- find_first_column(alpha_data, c("sampleID", "SampleID", "sample_id", "Sample_ID"))
subject_col <- find_first_column(alpha_data, c("subjectID", "SubjectID", "subject_id", "Subject_ID", "group"))
timepoint_col <- find_first_column(alpha_data, c("Timepoint", "timepoint", "TimePoint", "time_point"))
days_col <- find_first_column(alpha_data, c("days_from_baseline", "Days_from_baseline", "day_from_baseline", "days"))
altitude_col <- find_first_column(alpha_data, c("altitude", "Altitude", "elevation", "Elevation"))

if (is.na(sample_col)) stop("No sample ID column detected in alpha-diversity table.", call. = FALSE)
if (sample_col != "sampleID") names(alpha_data)[names(alpha_data) == sample_col] <- "sampleID"
if (!is.na(subject_col) && subject_col != "subjectID") names(alpha_data)[names(alpha_data) == subject_col] <- "subjectID"
if (!is.na(timepoint_col) && timepoint_col != "Timepoint") names(alpha_data)[names(alpha_data) == timepoint_col] <- "Timepoint"
if (!is.na(days_col) && days_col != "days_from_baseline") names(alpha_data)[names(alpha_data) == days_col] <- "days_from_baseline"
if (!is.na(altitude_col) && altitude_col != "altitude") names(alpha_data)[names(alpha_data) == altitude_col] <- "altitude"

if (any(is.na(alpha_data$sampleID) | alpha_data$sampleID == "")) {
  stop("Missing sample IDs detected; cannot create a manuscript summary table safely.", call. = FALSE)
}
if (anyDuplicated(alpha_data$sampleID) > 0) {
  stop("Duplicated sample IDs detected in alpha-diversity table.", call. = FALSE)
}
if (!"subjectID" %in% names(alpha_data)) {
  stop("No subjectID column detected. Longitudinal model inference requires subject-level repeated measures.", call. = FALSE)
}
if (!("days_from_baseline" %in% names(alpha_data)) && !("Timepoint" %in% names(alpha_data)) && !("altitude" %in% names(alpha_data))) {
  stop("No usable time/phase variable detected: expected days_from_baseline, Timepoint, or altitude.", call. = FALSE)
}

metric_order <- c("Observed_richness", "Chao1", "ACE", "Shannon", "Simpson")
available_metrics <- metric_order[metric_order %in% names(alpha_data)]
if (length(available_metrics) == 0) {
  stop("No supported alpha-diversity metric columns detected.", call. = FALSE)
}

for (metric in available_metrics) {
  alpha_data[[metric]] <- coerce_numeric(alpha_data[[metric]])
}

if ("days_from_baseline" %in% names(alpha_data)) {
  alpha_data$days_from_baseline <- coerce_numeric(alpha_data$days_from_baseline)
}
if ("altitude" %in% names(alpha_data)) {
  alpha_data$altitude <- coerce_numeric(alpha_data$altitude)
}

baseline_label <- NA_character_
comparison_label <- NA_character_
baseline_definition <- NA_character_
comparison_definition <- NA_character_
descriptive_variable <- NA_character_

if ("Timepoint" %in% names(alpha_data)) {
  alpha_data$Timepoint <- as.character(alpha_data$Timepoint)
  timepoint_order <- clean_timepoint_order(alpha_data$Timepoint, alpha_data)
  alpha_data$Timepoint <- factor(alpha_data$Timepoint, levels = timepoint_order, ordered = TRUE)
  descriptive_variable <- "Timepoint"

  if ("T1" %in% timepoint_order) {
    baseline_label <- "T1"
  } else {
    baseline_label <- timepoint_order[[1]]
  }
  baseline_definition <- paste0("Baseline = ", baseline_label)

  if ("altitude" %in% names(alpha_data) && any(is.finite(alpha_data$altitude))) {
    altitude_by_timepoint <- alpha_data %>%
      filter(!is.na(.data$Timepoint), is.finite(.data$altitude)) %>%
      group_by(.data$Timepoint) %>%
      summarise(
        median_altitude = median(.data$altitude, na.rm = TRUE),
        first_day = if ("days_from_baseline" %in% names(alpha_data)) median(.data$days_from_baseline, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      ) %>%
      arrange(desc(.data$median_altitude), .data$first_day, .data$Timepoint)

    if (nrow(altitude_by_timepoint) > 0 && max(altitude_by_timepoint$median_altitude, na.rm = TRUE) > min(altitude_by_timepoint$median_altitude, na.rm = TRUE)) {
      comparison_label <- as.character(altitude_by_timepoint$Timepoint[[1]])
      comparison_definition <- paste0(
        "Comparison = highest-altitude timepoint ",
        comparison_label,
        " (median altitude ",
        format_number(altitude_by_timepoint$median_altitude[[1]], 0),
        " m)"
      )
    }
  }

  if (is.na(comparison_label) && "T10" %in% timepoint_order) {
    comparison_label <- "T10"
    comparison_definition <- "Comparison = final timepoint T10"
  }
  if (is.na(comparison_label)) {
    comparison_label <- tail(timepoint_order, 1)
    comparison_definition <- paste0("Comparison = latest available timepoint ", comparison_label)
  }
} else if ("altitude" %in% names(alpha_data) && any(is.finite(alpha_data$altitude))) {
  descriptive_variable <- "altitude"
  baseline_label <- min(alpha_data$altitude, na.rm = TRUE)
  comparison_label <- max(alpha_data$altitude, na.rm = TRUE)
  baseline_definition <- paste0("Baseline/low altitude = ", baseline_label, " m")
  comparison_definition <- paste0("Comparison/high altitude = ", comparison_label, " m")
} else {
  descriptive_variable <- "days_from_baseline"
  baseline_label <- min(alpha_data$days_from_baseline, na.rm = TRUE)
  comparison_label <- max(alpha_data$days_from_baseline, na.rm = TRUE)
  baseline_definition <- paste0("Baseline = minimum days from baseline (", baseline_label, ")")
  comparison_definition <- paste0("Comparison = maximum days from baseline (", comparison_label, ")")
}

if ("days_from_baseline" %in% names(alpha_data) && any(is.finite(alpha_data$days_from_baseline))) {
  model_variable <- "days_from_baseline"
  model_variable_label <- "days_from_baseline"
} else if ("Timepoint" %in% names(alpha_data)) {
  alpha_data$timepoint_numeric <- as.numeric(alpha_data$Timepoint)
  model_variable <- "timepoint_numeric"
  model_variable_label <- "ordered Timepoint"
} else {
  model_variable <- NA_character_
  model_variable_label <- descriptive_variable
}

extract_lmm_result <- function(model, term, formula_text, n_samples, n_subjects, warnings) {
  summary_coef <- as.data.frame(coef(summary(model)))
  if (!term %in% rownames(summary_coef)) {
    stop("Model term not found in coefficient table: ", term, call. = FALSE)
  }

  estimate <- summary_coef[term, "Estimate"]
  std_error <- summary_coef[term, "Std. Error"]
  statistic_col <- intersect(c("t value", "z value"), colnames(summary_coef))
  statistic <- if (length(statistic_col) > 0) summary_coef[term, statistic_col[[1]]] else estimate / std_error
  df <- if ("df" %in% colnames(summary_coef)) summary_coef[term, "df"] else NA_real_

  p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(summary_coef))
  p_value <- if (length(p_col) > 0) {
    summary_coef[term, p_col[[1]]]
  } else {
    2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  }

  conv_messages <- model@optinfo$conv$lme4$messages
  converged <- length(conv_messages) == 0
  all_warnings <- unique(c(warnings, conv_messages))
  if (length(p_col) == 0) {
    all_warnings <- unique(c(all_warnings, "lmerTest unavailable; p-value uses normal approximation from lme4 t statistic"))
  }

  data.frame(
    Metric = NA_character_,
    Model_or_test = "Linear mixed model",
    Formula = formula_text,
    Term = term,
    Estimate = estimate,
    Std_error = std_error,
    Statistic = statistic,
    DF = df,
    P_value = p_value,
    FDR = NA_real_,
    n_samples = n_samples,
    n_subjects = n_subjects,
    Converged = converged,
    Warnings = paste(all_warnings, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

fit_metric_model <- function(data, metric) {
  metric_values <- data[[metric]]
  n_total_missing <- sum(!is.finite(metric_values))
  model_data <- data %>%
    mutate(metric_value = .data[[metric]]) %>%
    filter(is.finite(.data$metric_value), !is.na(.data$subjectID), .data$subjectID != "")

  if (!is.na(model_variable)) {
    model_data <- model_data %>% filter(is.finite(.data[[model_variable]]))
  }

  n_samples <- nrow(model_data)
  n_subjects <- dplyr::n_distinct(model_data$subjectID)
  if (n_samples < 10 || n_subjects < 3 || length(unique(model_data$metric_value)) < 2 || is.na(model_variable)) {
    reason <- paste0(
      "Insufficient data for LMM or no model time variable; n_samples=",
      n_samples,
      ", n_subjects=",
      n_subjects,
      ", missing_values=",
      n_total_missing
    )
    return(data.frame(
      Metric = metric,
      Model_or_test = "Not modeled",
      Formula = NA_character_,
      Term = NA_character_,
      Estimate = NA_real_,
      Std_error = NA_real_,
      Statistic = NA_real_,
      DF = NA_real_,
      P_value = NA_real_,
      FDR = NA_real_,
      n_samples = n_samples,
      n_subjects = n_subjects,
      Converged = FALSE,
      Warnings = reason,
      stringsAsFactors = FALSE
    ))
  }

  formula_text <- paste0("metric_value ~ ", model_variable, " + (1 | subjectID)")
  formula_obj <- stats::as.formula(formula_text)
  warnings <- character()

  lmer_fun <- if (requireNamespace("lmerTest", quietly = TRUE)) lmerTest::lmer else lme4::lmer
  fit <- tryCatch(
    withCallingHandlers(
      lmer_fun(formula_obj, data = model_data, REML = TRUE, control = lme4::lmerControl(check.conv.singular = "ignore")),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      Metric = metric,
      Model_or_test = "Linear mixed model failed",
      Formula = formula_text,
      Term = model_variable,
      Estimate = NA_real_,
      Std_error = NA_real_,
      Statistic = NA_real_,
      DF = NA_real_,
      P_value = NA_real_,
      FDR = NA_real_,
      n_samples = n_samples,
      n_subjects = n_subjects,
      Converged = FALSE,
      Warnings = conditionMessage(fit),
      stringsAsFactors = FALSE
    ))
  }

  result <- extract_lmm_result(fit, model_variable, formula_text, n_samples, n_subjects, warnings)
  result$Metric <- metric
  if (lme4::isSingular(fit, tol = 1e-4)) {
    result$Warnings <- paste(unique(c(result$Warnings, "singular fit")), collapse = "; ")
  }
  result
}

summary_rows <- lapply(available_metrics, function(metric) {
  values <- alpha_data[[metric]]

  if (descriptive_variable == "Timepoint") {
    baseline_values <- values[as.character(alpha_data$Timepoint) == baseline_label]
    comparison_values <- values[as.character(alpha_data$Timepoint) == comparison_label]
  } else if (descriptive_variable == "altitude") {
    baseline_values <- values[alpha_data$altitude == baseline_label]
    comparison_values <- values[alpha_data$altitude == comparison_label]
  } else {
    baseline_values <- values[alpha_data$days_from_baseline == baseline_label]
    comparison_values <- values[alpha_data$days_from_baseline == comparison_label]
  }

  data.frame(
    Metric = metric,
    n_samples = sum(is.finite(values)),
    n_subjects = dplyr::n_distinct(alpha_data$subjectID[is.finite(values) & !is.na(alpha_data$subjectID)]),
    Baseline_median_IQR = format_median_iqr(baseline_values, metric),
    High_altitude_or_late_phase_median_IQR = format_median_iqr(comparison_values, metric),
    Overall_min_max = format_min_max(values, metric),
    Missing_values = sum(!is.finite(values)),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

statistical_results <- lapply(available_metrics, function(metric) fit_metric_model(alpha_data, metric)) %>%
  bind_rows()

modeled <- !is.na(statistical_results$P_value)
statistical_results$FDR[modeled] <- stats::p.adjust(statistical_results$P_value[modeled], method = "BH")

effect_direction <- function(estimate) {
  if (is.na(estimate) || abs(estimate) < 1e-12) return("No clear direction")
  if (estimate > 0) "Increase" else "Decrease"
}

interpret_result <- function(metric, estimate, p_value, fdr, model_or_test) {
  if (model_or_test == "Not modeled" || is.na(estimate)) {
    return("Metric unavailable or insufficient data for formal modeling.")
  }
  direction <- tolower(effect_direction(estimate))
  if (!is.na(fdr) && fdr < 0.05) {
    paste0("FDR-significant ", direction, " over ", model_variable_label, ".")
  } else if (!is.na(p_value) && p_value < 0.05) {
    paste0("Nominal ", direction, " over ", model_variable_label, "; not FDR-significant.")
  } else {
    paste0("No clear evidence of directional change over ", model_variable_label, ".")
  }
}

caution_for_metric <- function(metric) {
  if (metric %in% c("Chao1", "ACE")) {
    "Not interpreted here because Chao1/ACE were not calculated from the relative-abundance-like input."
  } else if (metric == "Observed_richness") {
    "Richness is sensitive to filtering and detection thresholds; interpret with sequencing-depth context."
  } else {
    "Interpret alongside longitudinal trajectories and model diagnostics."
  }
}

main_table <- summary_rows %>%
  left_join(
    statistical_results %>%
      select(
        "Metric",
        Model_or_test = "Model_or_test",
        Effect_estimate = "Estimate",
        P_value = "P_value",
        FDR = "FDR"
      ),
    by = "Metric"
  ) %>%
  mutate(
    Comparison_variable = paste0(model_variable_label, "; descriptive ", baseline_label, " vs ", comparison_label),
    Effect_direction = vapply(.data$Effect_estimate, effect_direction, character(1)),
    Interpretation = mapply(
      interpret_result,
      .data$Metric,
      .data$Effect_estimate,
      .data$P_value,
      .data$FDR,
      .data$Model_or_test,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    ),
    Caution = vapply(.data$Metric, caution_for_metric, character(1))
  ) %>%
  select(
    "Metric",
    "Comparison_variable",
    "n_samples",
    "n_subjects",
    "Baseline_median_IQR",
    "High_altitude_or_late_phase_median_IQR",
    "Overall_min_max",
    "Missing_values",
    "Effect_direction",
    "Effect_estimate",
    "P_value",
    "FDR",
    "Interpretation",
    "Caution"
  )

main_table_export <- main_table %>%
  mutate(
    Effect_estimate = vapply(.data$Effect_estimate, format_effect, character(1)),
    P_value = vapply(.data$P_value, format_p, character(1)),
    FDR = vapply(.data$FDR, format_p, character(1))
  )

statistical_results_export <- statistical_results %>%
  mutate(
    Estimate = vapply(.data$Estimate, format_effect, character(1)),
    Std_error = vapply(.data$Std_error, format_effect, character(1)),
    Statistic = vapply(.data$Statistic, format_effect, character(1)),
    DF = vapply(.data$DF, format_effect, character(1)),
    P_value = vapply(.data$P_value, format_p, character(1)),
    FDR = vapply(.data$FDR, format_p, character(1))
  )

main_csv <- file.path(tables_dir, "alpha_diversity_publication_summary_table.csv")
main_tsv <- file.path(tables_dir, "alpha_diversity_publication_summary_table.tsv")
full_csv <- file.path(tables_dir, "alpha_diversity_statistical_results_full.csv")
markdown_table_path <- file.path(reports_dir, "alpha_diversity_publication_summary_table.md")
report_path <- file.path(reports_dir, "alpha_diversity_publication_summary_table_report.md")

write_csv_base(main_table_export, main_csv)
write_tsv_base(main_table_export, main_tsv)
write_csv_base(statistical_results_export, full_csv)

markdown_display <- main_table_export %>%
  transmute(
    Metric = .data$Metric,
    `Time variable` = .data$Comparison_variable,
    `n samples` = .data$n_samples,
    `n subjects` = .data$n_subjects,
    `Baseline median [IQR]` = .data$Baseline_median_IQR,
    `Comparison median [IQR]` = .data$High_altitude_or_late_phase_median_IQR,
    Direction = .data$Effect_direction,
    `Effect estimate` = .data$Effect_estimate,
    `P value` = .data$P_value,
    FDR = .data$FDR,
    Interpretation = .data$Interpretation
  )

markdown_text <- c(
  "# Alpha-Diversity Publication Summary Table",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("Input table: `", sub(paste0("^", repo_root, "/?"), "", alpha_input), "`"),
  "",
  paste0("Baseline definition: ", baseline_definition),
  "",
  paste0("Comparison definition: ", comparison_definition),
  "",
  paste0("Primary model/test variable: ", model_variable_label),
  "",
  knitr::kable(markdown_display, format = "markdown", align = "l"),
  "",
  "Notes: P values are from the primary longitudinal model when available. FDR values use Benjamini-Hochberg correction across tested alpha-diversity metrics. Chao1 and ACE are retained for transparency but are not modeled when unavailable."
)
writeLines(markdown_text, markdown_table_path)

output_paths <- c(main_csv, main_tsv, full_csv, markdown_table_path, report_path, session_info_path, log_path)

report_lines <- c(
  "# Alpha-Diversity Publication Summary Table Report",
  "",
  paste0("Date and time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Inputs",
  "",
  paste0("- Alpha-diversity table: `", sub(paste0("^", repo_root, "/?"), "", alpha_input), "`"),
  paste0("- Metadata: ", ifelse(is.na(metadata_used), "not used", paste0("`", metadata_used, "`"))),
  paste0("- Detected metadata columns: ", paste(intersect(c("sampleID", "subjectID", "Timepoint", "days_from_baseline", "altitude", "date", "date_converted", "age", "gender"), names(alpha_data)), collapse = ", ")),
  "",
  "## Summary",
  "",
  paste0("- Metrics summarized: ", paste(available_metrics, collapse = ", ")),
  paste0("- Samples in input table: ", nrow(alpha_data)),
  paste0("- Subjects in input table: ", dplyr::n_distinct(alpha_data$subjectID)),
  paste0("- Baseline definition: ", baseline_definition),
  paste0("- Comparison phase definition: ", comparison_definition),
  paste0("- Statistical model/test used where possible: metric ~ ", model_variable_label, " + (1 | subjectID)"),
  "- FDR correction: Benjamini-Hochberg across metrics with valid model/test P values.",
  "",
  "## Outputs",
  "",
  paste0("- `", sub(paste0("^", repo_root, "/?"), "", main_csv), "`"),
  paste0("- `", sub(paste0("^", repo_root, "/?"), "", main_tsv), "`"),
  paste0("- `", sub(paste0("^", repo_root, "/?"), "", full_csv), "`"),
  paste0("- `", sub(paste0("^", repo_root, "/?"), "", markdown_table_path), "`"),
  paste0("- `", sub(paste0("^", repo_root, "/?"), "", session_info_path), "`"),
  "",
  "## Model Notes and Limitations",
  "",
  paste0("- ", paste(unique(statistical_results$Warnings[statistical_results$Warnings != ""]), collapse = "\n- ")),
  "- Alpha diversity summarizes within-sample diversity; it does not identify which MAGs changed.",
  "- Richness metrics such as Observed richness, Chao1, and ACE are sensitive to sequencing depth, filtering, and feature detection thresholds.",
  "- Chao1 and ACE should be interpreted cautiously when calculated from relative-abundance rather than raw count-like data. In this run, they were unavailable because the alpha-diversity workflow flagged the input as relative-abundance-like.",
  "- Because this is a longitudinal dataset, subject-level repeated sampling is handled with a random intercept where data are sufficient.",
  "- The publication summary table is a compact overview; use the full statistical table for detailed model interpretation.",
  "",
  "## Run Instructions",
  "",
  "Local:",
  "",
  "```bash",
  "Rscript scripts/alpha_diversity/08_create_alpha_diversity_publication_summary_table.R",
  "```",
  "",
  "HPC:",
  "",
  "```bash",
  "qsub hpc/job_alpha_diversity_publication_summary_table.pbs",
  "```"
)
writeLines(report_lines, report_path)

log_message("Wrote main summary table: ", main_csv)
log_message("Wrote detailed statistical table: ", full_csv)
log_message("Wrote Markdown summary table: ", markdown_table_path)
log_message("Wrote report: ", report_path)
log_message("Finished alpha-diversity publication summary table workflow")
