# ============================================
# validate_formula.R - Model Formula Validation
# ============================================
# Checks that config$model$formula is compatible with the processed data before
# fitting models.

library(yaml)
library(here)

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

extract_random_terms <- function(formula_text) {
  matches <- gregexpr("\\([^()]*\\|[^()]*\\)", formula_text, perl = TRUE)
  terms <- regmatches(formula_text, matches)[[1]]
  if (length(terms) == 1 && identical(terms, character(0))) character(0) else terms
}

parse_random_term <- function(term) {
  inner <- gsub("^\\(|\\)$", "", term)
  pieces <- strsplit(inner, "\\|", fixed = FALSE)[[1]]
  if (length(pieces) != 2) {
    return(list(term = term, slopes = character(0), group = NA_character_))
  }

  lhs <- trimws(pieces[1])
  group <- trimws(pieces[2])
  slope_parts <- trimws(strsplit(lhs, "\\+")[[1]])
  slopes <- setdiff(slope_parts, c("1", "0", "-1", ""))
  list(term = term, slopes = slopes, group = group)
}

formula_without_random_terms <- function(formula_text, random_terms) {
  fixed_text <- formula_text
  for (term in random_terms) {
    fixed_text <- gsub(paste0("\\s*\\+\\s*", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", term)), "", fixed_text)
    fixed_text <- gsub(gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", term), "1", fixed_text)
  }
  fixed_text
}

config <- yaml::read_yaml(here("config", "config.yml"))
processed_path <- resolve_path(config$data$processed %||% "data/processed/cleaned_data.rds")
report_path <- here("data", "processed", "formula_validation_report.txt")

dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

sink(report_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Formula Validation Report\n")
cat("=========================\n\n")
cat("Processed file:", processed_path, "\n")
cat("Formula report:", report_path, "\n\n")

if (!file.exists(processed_path)) {
  stop("Processed data not found: ", processed_path, call. = FALSE)
}

df <- readRDS(processed_path)

formula_text <- config$model$formula %||% ""
if (!nzchar(formula_text)) {
  stop("config$model$formula is empty.", call. = FALSE)
}

model_formula <- as.formula(formula_text)
formula_variables <- all.vars(model_formula)
missing_formula_variables <- setdiff(formula_variables, names(df))

random_terms <- extract_random_terms(formula_text)
parsed_random_terms <- lapply(random_terms, parse_random_term)
random_groups <- unique(stats::na.omit(vapply(parsed_random_terms, `[[`, character(1), "group")))
random_slopes <- unique(unlist(lapply(parsed_random_terms, `[[`, "slopes"), use.names = FALSE))

expected_group <- config$model$subject %||% "Subject"
fixed_formula_text <- formula_without_random_terms(formula_text, random_terms)
fixed_variables <- all.vars(as.formula(fixed_formula_text))

cat("Configured Formula\n")
cat("------------------\n")
cat("Formula:", formula_text, "\n")
cat("Family:", config$model$family %||% "<missing>", "\n")
cat("Optimizer:", config$model$optimizer %||% "<missing>", "\n")
cat("Formula variables:", paste(formula_variables, collapse = ", "), "\n")
cat("Fixed-effect variables:", paste(fixed_variables, collapse = ", "), "\n")
cat("Random-effect terms:", if (length(random_terms) > 0) paste(random_terms, collapse = ", ") else "<none>", "\n\n")

cat("Column Checks\n")
cat("-------------\n")
all_formula_columns_exist <- length(missing_formula_variables) == 0
cat(check_item(all_formula_columns_exist, "Every formula variable exists in processed data"), "\n")
if (!all_formula_columns_exist) {
  cat("Missing formula variables:", paste(missing_formula_variables, collapse = ", "), "\n")
}

random_group_exists <- expected_group %in% names(df)
cat(check_item(random_group_exists, paste("Random-effect grouping variable exists:", expected_group)), "\n")

if (length(random_groups) > 0) {
  cat("Random groups found in formula:", paste(random_groups, collapse = ", "), "\n")
}

uses_group <- "Group" %in% formula_variables
if (uses_group) {
  cat("[ ] Formula does not use Group - formula contains Group, but no true treatment/condition Group has been identified.\n")
} else {
  cat("[x] Formula does not use Group before a true treatment/condition variable is identified.\n")
}

cat("\nRandom Slope Checks\n")
cat("-------------------\n")
if (length(random_slopes) == 0) {
  cat("No random slope detected. Current random-effect structure is random intercept only.\n")
  random_slope_ok <- TRUE
} else {
  random_slope_ok <- TRUE
  for (slope in random_slopes) {
    exists <- slope %in% names(df)
    numeric <- exists && is.numeric(df[[slope]])
    cat(check_item(exists, paste("Random slope exists:", slope)), "\n")
    cat(check_item(numeric, paste("Random slope is numeric:", slope)), "\n")
    random_slope_ok <- random_slope_ok && exists && numeric
  }
}

cat("\nRepeated Observation Checks\n")
cat("---------------------------\n")
if (random_group_exists) {
  rows_per_subject <- sort(table(df[[expected_group]]), decreasing = TRUE)
  sampleids_per_subject <- stats::aggregate(sampleID ~ Subject, df, function(x) length(unique(x)))
  timepoints_per_subject <- stats::aggregate(Timepoint ~ Subject, df, function(x) length(unique(x)))
  days_per_subject <- stats::aggregate(days_from_baseline_subject ~ Subject, df, function(x) length(unique(x)))

  cat("Rows per Subject:\n")
  print(rows_per_subject)
  cat("\nUnique sampleID per Subject:\n")
  print(sampleids_per_subject)
  cat("\nUnique Timepoint per Subject:\n")
  print(timepoints_per_subject)
  cat("\nUnique days_from_baseline_subject per Subject:\n")
  print(days_per_subject)

  subjects_one_sample <- sampleids_per_subject$Subject[sampleids_per_subject$sampleID <= 1]
  enough_repeats_overall <- any(sampleids_per_subject$sampleID > 1)
  every_subject_repeated <- all(sampleids_per_subject$sampleID > 1)

  if (length(subjects_one_sample) > 0) {
    cat("\nWarning: subjects with only one sample/timepoint:", paste(subjects_one_sample, collapse = ", "), "\n")
  }
} else {
  rows_per_subject <- integer(0)
  sampleids_per_subject <- data.frame()
  timepoints_per_subject <- data.frame()
  days_per_subject <- data.frame()
  subjects_one_sample <- character(0)
  enough_repeats_overall <- FALSE
  every_subject_repeated <- FALSE
}

cat("\nComplexity Assessment\n")
cat("---------------------\n")
has_interaction <- grepl(":", formula_text, fixed = TRUE) || grepl("\\*", formula_text)
has_random_slope <- length(random_slopes) > 0
current_formula_is_simple <- !uses_group && !has_interaction && !has_random_slope &&
  identical(formula_text, "Y ~ days_from_baseline_subject + (1 | Subject)")

if (current_formula_is_simple) {
  cat("The configured formula is the recommended simple starting formula.\n")
} else {
  cat("The configured formula is more complex than the recommended starting formula or differs from it.\n")
}

if (uses_group) {
  cat("Warning: formulas containing Group should not be used yet because no true treatment/condition Group exists.\n")
}
if (has_random_slope) {
  cat("Warning: random-slope formulas should be tested only after checking repeated observations and model convergence.\n")
}
if (length(subjects_one_sample) > 0) {
  cat("Warning: random-effect estimation should account for subjects with only one sample/timepoint, currently:", paste(subjects_one_sample, collapse = ", "), "\n")
}

cat("\nRecommendation\n")
cat("--------------\n")
cat("Safest starting formula:\n")
cat('formula: "Y ~ days_from_baseline_subject + (1 | Subject)"\n')
cat('family: "gaussian"\n')
cat('optimizer: "bobyqa"\n\n')
cat("A later candidate after convergence and repeat-structure checks:\n")
cat("Y ~ days_from_baseline_subject + (1 + days_from_baseline_subject | Subject)\n")

formula_not_too_complex <- current_formula_is_simple

cat("\nChecklist Summary\n")
cat("-----------------\n")
cat(check_item(all_formula_columns_exist, "Every variable in config$model$formula exists in processed data"), "\n")
cat(check_item(random_group_exists, "Random-effect grouping variable exists"), "\n")
cat(check_item(random_slope_ok, "Random slope variable is numeric if modeled"), "\n")
cat(check_item(enough_repeats_overall, "There are enough repeated observations overall per random-effect group"), "\n")
cat(check_item(formula_not_too_complex, "Formula is not too complex for the starting analysis"), "\n")
cat(check_item(!uses_group, "Formula does not include Group before a true Group variable is identified"), "\n")

cat("\nFormula validation report saved to data/processed/formula_validation_report.txt\n")
