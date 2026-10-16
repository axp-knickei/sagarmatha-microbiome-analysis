# Prepare validated alpha-diversity input for sample-level longitudinal LMMs.

required_packages <- c("dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages(library(dplyr))

repo_root <- getwd()
alpha_dir <- Sys.getenv("BLACKICE_ALPHA_OUTPUT_DIR", unset = "results/alpha_diversity")
output_dir <- Sys.getenv("BLACKICE_ALPHA_LMM_OUTPUT_DIR", unset = "results/alpha_diversity_lmm")
tables_dir <- file.path(output_dir, "tables")
logs_dir <- file.path(output_dir, "logs")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(alpha_dir, "tables", "alpha_diversity_with_metadata.csv")
vegan_path <- file.path(alpha_dir, "tables", "alpha_diversity_vegan.csv")
diagnostic_path <- file.path(alpha_dir, "tables", "alpha_diversity_diagnostic_summary.csv")

log_path <- file.path(logs_dir, "alpha_lmm_prepare.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Alpha-diversity LMM input preparation started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Working directory:", repo_root, "\n")

if (file.exists(input_path)) {
  alpha <- utils::read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
  input_source <- input_path
  alpha_reused <- TRUE
} else if (file.exists(vegan_path)) {
  stop(
    "Found alpha-diversity metrics but not the metadata-joined table: ",
    vegan_path,
    ". Run scripts/alpha_diversity/04_compare_alpha_diversity_tools.R before alpha LMM modeling.",
    call. = FALSE
  )
} else {
  stop(
    "No alpha-diversity table was found. Run the alpha-diversity workflow first, for example: ",
    "Rscript scripts/alpha_diversity/01_prepare_alpha_diversity_inputs.R; ",
    "Rscript scripts/alpha_diversity/02_calculate_alpha_diversity_vegan.R; ",
    "python scripts/alpha_diversity/03_calculate_alpha_diversity_scikitbio.py; ",
    "Rscript scripts/alpha_diversity/04_compare_alpha_diversity_tools.R",
    call. = FALSE
  )
}

diagnostic <- if (file.exists(diagnostic_path)) {
  utils::read.csv(diagnostic_path, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  data.frame(item = character(), value = character(), stringsAsFactors = FALSE)
}

metric_columns <- c("Shannon", "Simpson", "Observed_richness", "Chao1", "ACE")
metadata_columns <- c("sampleID", "subjectID", "Timepoint", "date", "days_from_baseline", "altitude", "age", "gender")
required_columns <- c("sampleID", "subjectID")
missing_required <- setdiff(required_columns, names(alpha))
if (length(missing_required) > 0) {
  stop("Alpha-diversity input is missing required column(s): ", paste(missing_required, collapse = ", "), call. = FALSE)
}

missing_metrics <- setdiff(metric_columns, names(alpha))
if (length(missing_metrics) > 0) {
  stop("Alpha-diversity input is missing required metric column(s): ", paste(missing_metrics, collapse = ", "), call. = FALSE)
}

alpha$sampleID <- as.character(alpha$sampleID)
alpha$subjectID <- as.character(alpha$subjectID)

if (any(is.na(alpha$sampleID) | !nzchar(alpha$sampleID))) {
  stop("Alpha-diversity input contains missing sampleID values.", call. = FALSE)
}
if (anyDuplicated(alpha$sampleID) > 0) {
  duplicated_ids <- unique(alpha$sampleID[duplicated(alpha$sampleID)])
  stop("Alpha-diversity input contains duplicated sampleID values: ", paste(head(duplicated_ids, 10), collapse = ", "), call. = FALSE)
}
if (any(is.na(alpha$subjectID) | !nzchar(alpha$subjectID))) {
  stop("Alpha-diversity input contains missing subjectID values.", call. = FALSE)
}

for (metric in metric_columns) {
  alpha[[metric]] <- suppressWarnings(as.numeric(alpha[[metric]]))
}

if ("days_from_baseline" %in% names(alpha)) {
  alpha$days_from_baseline <- suppressWarnings(as.numeric(alpha$days_from_baseline))
}
if ("altitude" %in% names(alpha)) {
  alpha$altitude <- suppressWarnings(as.numeric(alpha$altitude))
}
if ("age" %in% names(alpha)) {
  alpha$age <- suppressWarnings(as.numeric(alpha$age))
}
if ("gender" %in% names(alpha)) {
  alpha$gender <- as.factor(alpha$gender)
}
if ("Timepoint" %in% names(alpha)) {
  alpha$Timepoint <- as.character(alpha$Timepoint)
}

natural_timepoint_order <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  numeric_tp <- x[grepl("^T[0-9]+$", x)]
  other_tp <- setdiff(x, numeric_tp)
  c(numeric_tp[order(as.integer(sub("^T", "", numeric_tp)))], sort(other_tp))
}

has_days <- "days_from_baseline" %in% names(alpha) && sum(!is.na(alpha$days_from_baseline)) >= 3
has_timepoint <- "Timepoint" %in% names(alpha) && sum(!is.na(alpha$Timepoint) & nzchar(alpha$Timepoint)) >= 3

if (has_days) {
  time_variable <- "days_from_baseline"
  alpha$alpha_lmm_time <- alpha$days_from_baseline
  alpha$alpha_lmm_time_source <- "days_from_baseline"
} else if (has_timepoint) {
  timepoint_levels <- natural_timepoint_order(alpha$Timepoint)
  alpha$Timepoint <- factor(alpha$Timepoint, levels = timepoint_levels, ordered = TRUE)
  alpha$timepoint_numeric <- as.numeric(alpha$Timepoint)
  time_variable <- "timepoint_numeric"
  alpha$alpha_lmm_time <- alpha$timepoint_numeric
  alpha$alpha_lmm_time_source <- "ordered_Timepoint"
} else {
  stop("No usable longitudinal time variable found. Need days_from_baseline or Timepoint.", call. = FALSE)
}

alpha$subjectID <- factor(alpha$subjectID)
if ("Timepoint" %in% names(alpha)) {
  tp_levels <- natural_timepoint_order(as.character(alpha$Timepoint))
  alpha$Timepoint <- factor(as.character(alpha$Timepoint), levels = tp_levels, ordered = TRUE)
}

diagnostic_lookup <- function(key, default = NA_character_) {
  if (!all(c("item", "value") %in% names(diagnostic))) {
    return(default)
  }
  hit <- diagnostic$value[diagnostic$item == key]
  if (length(hit) == 0) default else as.character(hit[[1]])
}

abundance_scale <- diagnostic_lookup("abundance_scale_decision", default = unique(alpha$abundance_scale_decision)[1])
chao1_ace_decision <- if ("chao1_ace_decision" %in% names(alpha)) unique(alpha$chao1_ace_decision)[1] else NA_character_

metric_availability <- lapply(metric_columns, function(metric) {
  values <- alpha[[metric]]
  non_missing <- sum(!is.na(values))
  finite_values <- values[is.finite(values)]
  negative_values <- sum(finite_values < 0)
  n_subjects <- dplyr::n_distinct(alpha$subjectID[!is.na(values)])
  repeated_subjects <- alpha %>%
    filter(!is.na(.data[[metric]]), !is.na(alpha_lmm_time)) %>%
    count(subjectID, name = "n") %>%
    filter(n >= 2) %>%
    nrow()

  reason <- "available"
  available <- TRUE
  if (non_missing < 10) {
    available <- FALSE
    reason <- "fewer_than_10_non_missing_samples"
  } else if (n_subjects < 3) {
    available <- FALSE
    reason <- "fewer_than_3_subjects"
  } else if (repeated_subjects < 2) {
    available <- FALSE
    reason <- "fewer_than_2_subjects_with_repeated_measurements"
  } else if (stats::var(values, na.rm = TRUE) <= 0) {
    available <- FALSE
    reason <- "zero_variance"
  } else if (negative_values > 0 && metric %in% c("Observed_richness", "Chao1", "ACE")) {
    available <- FALSE
    reason <- "negative_richness_values"
  }

  if (metric %in% c("Chao1", "ACE") && non_missing == 0 && identical(as.character(abundance_scale), "relative_abundance_like")) {
    reason <- "not_calculated_relative_or_non_count_input"
  }

  data.frame(
    metric = metric,
    available_for_lmm = available,
    reason = reason,
    non_missing_samples = non_missing,
    n_subjects = n_subjects,
    subjects_with_repeated_measurements = repeated_subjects,
    negative_values = negative_values,
    variance = if (length(finite_values) > 1) stats::var(finite_values) else NA_real_,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

impossible_rows <- lapply(c("Observed_richness", "Chao1", "ACE"), function(metric) {
  if (!metric %in% names(alpha)) {
    return(0L)
  }
  sum(alpha[[metric]] < 0, na.rm = TRUE)
})
names(impossible_rows) <- c("Observed_richness", "Chao1", "ACE")

input_diagnostics <- data.frame(
  item = c(
    "analysis_input",
    "alpha_diversity_reused",
    "time_variable",
    "time_source",
    "n_samples",
    "n_subjects",
    "n_timepoints",
    "abundance_scale_decision",
    "chao1_ace_decision",
    "missing_subject_ids",
    "duplicated_sample_ids",
    "metrics_available_for_lmm",
    paste0("negative_", names(impossible_rows), "_values")
  ),
  value = c(
    input_source,
    as.character(alpha_reused),
    time_variable,
    unique(alpha$alpha_lmm_time_source)[1],
    as.character(nrow(alpha)),
    as.character(dplyr::n_distinct(alpha$subjectID)),
    if ("Timepoint" %in% names(alpha)) as.character(dplyr::n_distinct(alpha$Timepoint, na.rm = TRUE)) else "NA",
    as.character(abundance_scale),
    as.character(chao1_ace_decision),
    as.character(sum(is.na(alpha$subjectID) | !nzchar(as.character(alpha$subjectID)))),
    as.character(anyDuplicated(alpha$sampleID)),
    paste(metric_availability$metric[metric_availability$available_for_lmm], collapse = "; "),
    as.character(unlist(impossible_rows))
  ),
  stringsAsFactors = FALSE
)

prepared_path <- file.path(tables_dir, "alpha_lmm_input.csv")
availability_path <- file.path(tables_dir, "alpha_lmm_metric_availability.csv")
diagnostics_path <- file.path(tables_dir, "alpha_lmm_input_diagnostics.csv")

utils::write.csv(alpha, prepared_path, row.names = FALSE)
utils::write.csv(metric_availability, availability_path, row.names = FALSE)
utils::write.csv(input_diagnostics, diagnostics_path, row.names = FALSE)

cat("Prepared input:", prepared_path, "\n")
cat("Metric availability:", availability_path, "\n")
cat("Input diagnostics:", diagnostics_path, "\n")
cat("Alpha-diversity LMM input preparation finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
