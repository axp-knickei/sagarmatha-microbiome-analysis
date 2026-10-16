# Fit sample-level longitudinal LMMs for alpha-diversity metrics.

required_packages <- c("dplyr", "lme4")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
})

output_dir <- Sys.getenv("BLACKICE_ALPHA_LMM_OUTPUT_DIR", unset = "results/alpha_diversity_lmm")
tables_dir <- file.path(output_dir, "tables")
logs_dir <- file.path(output_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(logs_dir, "alpha_lmm_fit.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Alpha-diversity LMM fitting started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

input_path <- file.path(tables_dir, "alpha_lmm_input.csv")
availability_path <- file.path(tables_dir, "alpha_lmm_metric_availability.csv")
if (!file.exists(input_path)) {
  stop("Missing prepared input: ", input_path, ". Run 01_prepare_alpha_lmm_input.R first.", call. = FALSE)
}
if (!file.exists(availability_path)) {
  stop("Missing metric availability table: ", availability_path, ". Run 01_prepare_alpha_lmm_input.R first.", call. = FALSE)
}

alpha <- utils::read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
availability <- utils::read.csv(availability_path, check.names = FALSE, stringsAsFactors = FALSE)

metric_columns <- c("Shannon", "Simpson", "Observed_richness", "Chao1", "ACE")
metric_labels <- c(
  Shannon = "Shannon diversity",
  Simpson = "Simpson diversity",
  Observed_richness = "Observed richness",
  Chao1 = "Chao1 richness",
  ACE = "ACE richness"
)

alpha$subjectID <- factor(alpha$subjectID)
if ("gender" %in% names(alpha)) alpha$gender <- factor(alpha$gender)
if ("Timepoint" %in% names(alpha)) alpha$Timepoint <- factor(alpha$Timepoint, ordered = TRUE)
for (metric in metric_columns) {
  if (metric %in% names(alpha)) alpha[[metric]] <- suppressWarnings(as.numeric(alpha[[metric]]))
}
numeric_columns <- intersect(c("days_from_baseline", "alpha_lmm_time", "timepoint_numeric", "altitude", "age"), names(alpha))
for (column in numeric_columns) {
  alpha[[column]] <- suppressWarnings(as.numeric(alpha[[column]]))
}

time_variable <- if ("days_from_baseline" %in% names(alpha) && sum(!is.na(alpha$days_from_baseline)) >= 3) {
  "days_from_baseline"
} else {
  "alpha_lmm_time"
}

fit_lmer_model <- function(formula, data) {
  warnings <- character()
  fit <- withCallingHandlers(
    {
      if (requireNamespace("lmerTest", quietly = TRUE)) {
        lmerTest::lmer(formula, data = data, REML = TRUE, control = lme4::lmerControl(optimizer = "bobyqa"))
      } else {
        lme4::lmer(formula, data = data, REML = TRUE, control = lme4::lmerControl(optimizer = "bobyqa"))
      }
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = unique(warnings))
}

convergence_messages <- function(fit) {
  messages <- fit@optinfo$conv$lme4$messages
  if (is.null(messages)) "" else paste(messages, collapse = "; ")
}

extract_fixed_effects <- function(fit, metric, model_type, formula_text, n_samples, n_subjects, warnings) {
  coef_table <- as.data.frame(summary(fit)$coefficients)
  coef_table$term <- rownames(coef_table)
  names(coef_table) <- sub("^Std\\. Error$", "std_error", names(coef_table))
  names(coef_table) <- sub("^t value$", "t_statistic", names(coef_table))
  names(coef_table) <- sub("^Pr\\(>\\|t\\|\\)$", "p_value", names(coef_table))
  names(coef_table) <- sub("^df$", "degrees_of_freedom", names(coef_table))
  names(coef_table) <- sub("^Estimate$", "estimate", names(coef_table))

  if (!"degrees_of_freedom" %in% names(coef_table)) coef_table$degrees_of_freedom <- NA_real_
  if (!"p_value" %in% names(coef_table)) {
    coef_table$p_value <- 2 * stats::pnorm(-abs(coef_table$t_statistic))
    p_method <- "normal_approximation_lme4_fallback"
  } else {
    p_method <- if (requireNamespace("lmerTest", quietly = TRUE)) "lmerTest_Satterthwaite" else "model_summary"
  }

  ci <- tryCatch(
    as.data.frame(suppressMessages(confint(fit, method = "Wald"))),
    error = function(e) data.frame(`2.5 %` = NA_real_, `97.5 %` = NA_real_, check.names = FALSE)
  )
  ci$term <- rownames(ci)
  names(ci) <- sub("^2\\.5 %$", "conf_low", names(ci))
  names(ci) <- sub("^97\\.5 %$", "conf_high", names(ci))

  coef_table %>%
    left_join(ci %>% select(term, conf_low, conf_high), by = "term") %>%
    transmute(
      metric = metric,
      metric_label = unname(metric_labels[metric]),
      model_type = model_type,
      term = term,
      estimate = estimate,
      standard_error = std_error,
      degrees_of_freedom = degrees_of_freedom,
      t_statistic = t_statistic,
      p_value = p_value,
      confidence_interval_low = conf_low,
      confidence_interval_high = conf_high,
      p_value_method = p_method,
      model_formula = formula_text,
      n_samples = n_samples,
      n_subjects = n_subjects,
      convergence_warnings = paste(warnings, collapse = "; ")
    )
}

extract_random_effects <- function(fit, metric, model_type) {
  as.data.frame(lme4::VarCorr(fit)) %>%
    transmute(
      metric = metric,
      model_type = model_type,
      group = grp,
      effect_1 = var1,
      effect_2 = var2,
      variance = vcov,
      standard_deviation = sdcor
    )
}

extract_fit_stats <- function(fit, metric, model_type, formula_text, n_samples, n_subjects, warnings, notes = "") {
  data.frame(
    metric = metric,
    metric_label = unname(metric_labels[metric]),
    model_type = model_type,
    model_formula = formula_text,
    n_samples = n_samples,
    n_subjects = n_subjects,
    AIC = stats::AIC(fit),
    BIC = stats::BIC(fit),
    logLik = as.numeric(stats::logLik(fit)),
    REML_criterion = tryCatch(lme4::REMLcrit(fit), error = function(e) NA_real_),
    singular_fit = lme4::isSingular(fit, tol = 1e-4),
    convergence_warnings = paste(c(warnings, convergence_messages(fit)), collapse = "; "),
    notes = notes,
    stringsAsFactors = FALSE
  )
}

make_predictions <- function(fit, metric, model_type, source_data, time_variable, n_grid = 100) {
  time_range <- range(source_data[[time_variable]], na.rm = TRUE)
  newdata <- data.frame(time_value = seq(time_range[1], time_range[2], length.out = n_grid))
  names(newdata) <- time_variable

  if ("subjectID" %in% all.vars(formula(fit))) {
    newdata$subjectID <- source_data$subjectID[which(!is.na(source_data$subjectID))[1]]
  }
  if ("altitude" %in% all.vars(formula(fit))) newdata$altitude <- stats::median(source_data$altitude, na.rm = TRUE)
  if ("age" %in% all.vars(formula(fit))) newdata$age <- stats::median(source_data$age, na.rm = TRUE)
  if ("gender" %in% all.vars(formula(fit))) {
    gender_mode <- names(sort(table(source_data$gender), decreasing = TRUE))[1]
    newdata$gender <- factor(gender_mode, levels = levels(source_data$gender))
  }

  fixed_term_labels <- attr(stats::terms(fit), "term.labels")
  fixed_formula <- if (length(fixed_term_labels) > 0) {
    stats::reformulate(fixed_term_labels)
  } else {
    stats::as.formula("~ 1")
  }
  design <- stats::model.matrix(fixed_formula, newdata)
  beta <- lme4::fixef(fit)
  vcov_beta <- as.matrix(stats::vcov(fit))
  common_cols <- intersect(colnames(design), names(beta))
  design <- design[, common_cols, drop = FALSE]
  beta <- beta[common_cols]
  vcov_beta <- vcov_beta[common_cols, common_cols, drop = FALSE]
  fit_values <- as.numeric(design %*% beta)
  se_values <- sqrt(pmax(0, diag(design %*% vcov_beta %*% t(design))))

  data.frame(
    metric = metric,
    metric_label = unname(metric_labels[metric]),
    model_type = model_type,
    time_variable = time_variable,
    time_value = newdata[[time_variable]],
    estimate = fit_values,
    standard_error = se_values,
    confidence_interval_low = fit_values - 1.96 * se_values,
    confidence_interval_high = fit_values + 1.96 * se_values,
    stringsAsFactors = FALSE
  )
}

adjusted_skip_reason <- function(data, time_variable) {
  required <- c("altitude", "age", "gender")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) return(paste0("missing_covariates: ", paste(missing, collapse = ", ")))
  complete <- data %>% filter(!is.na(.data[[time_variable]]), !is.na(altitude), !is.na(age), !is.na(gender))
  if (nrow(complete) < 30) return("fewer_than_30_complete_cases")
  if (dplyr::n_distinct(complete$subjectID) < 5) return("fewer_than_5_subjects_with_complete_covariates")
  if (dplyr::n_distinct(complete$gender) < 2) return("fewer_than_2_gender_levels")
  cor_time_alt <- suppressWarnings(stats::cor(complete[[time_variable]], complete$altitude, use = "complete.obs"))
  if (is.finite(cor_time_alt) && abs(cor_time_alt) > 0.85) {
    return(paste0("days_from_baseline_altitude_collinearity_r=", round(cor_time_alt, 3)))
  }
  "attempt"
}

fixed_effects <- list()
random_effects <- list()
fit_statistics <- list()
predictions <- list()
model_summary <- list()

for (metric in metric_columns) {
  availability_row <- availability %>% filter(metric == !!metric)
  metric_available <- nrow(availability_row) == 1 && isTRUE(as.logical(availability_row$available_for_lmm))
  metric_reason <- if (nrow(availability_row) == 1) availability_row$reason else "not_found_in_availability_table"

  if (!metric_available) {
    model_summary[[metric]] <- data.frame(
      metric = metric,
      metric_label = unname(metric_labels[metric]),
      status = "not_modeled",
      reason = metric_reason,
      time_variable = time_variable,
      main_model_formula = NA_character_,
      n_samples = sum(!is.na(alpha[[metric]])),
      n_subjects = dplyr::n_distinct(alpha$subjectID[!is.na(alpha[[metric]])]),
      adjusted_model_status = "not_attempted",
      adjusted_model_reason = metric_reason,
      stringsAsFactors = FALSE
    )
    next
  }

  metric_data <- alpha %>%
    filter(!is.na(.data[[metric]]), !is.na(.data[[time_variable]]), !is.na(subjectID)) %>%
    mutate(metric_value = .data[[metric]], subjectID = factor(subjectID))

  main_formula_text <- paste0("metric_value ~ ", time_variable, " + (1 | subjectID)")
  main_formula <- stats::as.formula(main_formula_text)
  main_result <- tryCatch(
    fit_lmer_model(main_formula, metric_data),
    error = function(e) list(fit = NULL, warnings = conditionMessage(e))
  )

  if (is.null(main_result$fit)) {
    model_summary[[metric]] <- data.frame(
      metric = metric,
      metric_label = unname(metric_labels[metric]),
      status = "model_failed",
      reason = paste(main_result$warnings, collapse = "; "),
      time_variable = time_variable,
      main_model_formula = main_formula_text,
      n_samples = nrow(metric_data),
      n_subjects = dplyr::n_distinct(metric_data$subjectID),
      adjusted_model_status = "not_attempted",
      adjusted_model_reason = "main_model_failed",
      stringsAsFactors = FALSE
    )
    next
  }

  fixed_effects[[paste(metric, "main", sep = "_")]] <- extract_fixed_effects(
    main_result$fit, metric, "main", main_formula_text, nrow(metric_data), dplyr::n_distinct(metric_data$subjectID), main_result$warnings
  )
  random_effects[[paste(metric, "main", sep = "_")]] <- extract_random_effects(main_result$fit, metric, "main")
  fit_statistics[[paste(metric, "main", sep = "_")]] <- extract_fit_stats(
    main_result$fit, metric, "main", main_formula_text, nrow(metric_data), dplyr::n_distinct(metric_data$subjectID), main_result$warnings
  )
  predictions[[paste(metric, "main", sep = "_")]] <- make_predictions(main_result$fit, metric, "main", metric_data, time_variable)

  adjusted_reason <- if (time_variable == "days_from_baseline") adjusted_skip_reason(metric_data, time_variable) else "adjusted_model_requires_days_from_baseline"
  adjusted_status <- "not_attempted"
  if (identical(adjusted_reason, "attempt")) {
    adjusted_data <- metric_data %>% filter(!is.na(altitude), !is.na(age), !is.na(gender))
    adjusted_formula_text <- paste0("metric_value ~ ", time_variable, " + altitude + age + gender + (1 | subjectID)")
    adjusted_formula <- stats::as.formula(adjusted_formula_text)
    adjusted_result <- tryCatch(
      fit_lmer_model(adjusted_formula, adjusted_data),
      error = function(e) list(fit = NULL, warnings = conditionMessage(e))
    )
    if (!is.null(adjusted_result$fit) && !lme4::isSingular(adjusted_result$fit, tol = 1e-4)) {
      adjusted_status <- "modeled"
      adjusted_reason <- "available_covariates_and_non_singular_fit"
      fixed_effects[[paste(metric, "adjusted", sep = "_")]] <- extract_fixed_effects(
        adjusted_result$fit, metric, "adjusted", adjusted_formula_text, nrow(adjusted_data), dplyr::n_distinct(adjusted_data$subjectID), adjusted_result$warnings
      )
      random_effects[[paste(metric, "adjusted", sep = "_")]] <- extract_random_effects(adjusted_result$fit, metric, "adjusted")
      fit_statistics[[paste(metric, "adjusted", sep = "_")]] <- extract_fit_stats(
        adjusted_result$fit, metric, "adjusted", adjusted_formula_text, nrow(adjusted_data), dplyr::n_distinct(adjusted_data$subjectID), adjusted_result$warnings
      )
    } else {
      adjusted_status <- "skipped"
      adjusted_reason <- if (is.null(adjusted_result$fit)) {
        paste("adjusted_model_failed:", paste(adjusted_result$warnings, collapse = "; "))
      } else {
        "adjusted_model_singular"
      }
    }
  } else {
    adjusted_status <- "skipped"
  }

  model_summary[[metric]] <- data.frame(
    metric = metric,
    metric_label = unname(metric_labels[metric]),
    status = if (lme4::isSingular(main_result$fit, tol = 1e-4)) "modeled_singular" else "modeled",
    reason = if (lme4::isSingular(main_result$fit, tol = 1e-4)) "main_model_singular_random_intercept" else "main_model_fit",
    time_variable = time_variable,
    main_model_formula = main_formula_text,
    n_samples = nrow(metric_data),
    n_subjects = dplyr::n_distinct(metric_data$subjectID),
    adjusted_model_status = adjusted_status,
    adjusted_model_reason = adjusted_reason,
    stringsAsFactors = FALSE
  )
}

empty_fixed <- data.frame(
  metric = character(), metric_label = character(), model_type = character(), term = character(),
  estimate = numeric(), standard_error = numeric(), degrees_of_freedom = numeric(), t_statistic = numeric(),
  p_value = numeric(), confidence_interval_low = numeric(), confidence_interval_high = numeric(),
  p_value_method = character(), model_formula = character(), n_samples = integer(), n_subjects = integer(),
  convergence_warnings = character(), stringsAsFactors = FALSE
)
empty_random <- data.frame(
  metric = character(), model_type = character(), group = character(), effect_1 = character(),
  effect_2 = character(), variance = numeric(), standard_deviation = numeric(), stringsAsFactors = FALSE
)
empty_fit <- data.frame(
  metric = character(), metric_label = character(), model_type = character(), model_formula = character(),
  n_samples = integer(), n_subjects = integer(), AIC = numeric(), BIC = numeric(), logLik = numeric(),
  REML_criterion = numeric(), singular_fit = logical(), convergence_warnings = character(), notes = character(),
  stringsAsFactors = FALSE
)
empty_prediction <- data.frame(
  metric = character(), metric_label = character(), model_type = character(), time_variable = character(),
  time_value = numeric(), estimate = numeric(), standard_error = numeric(),
  confidence_interval_low = numeric(), confidence_interval_high = numeric(), stringsAsFactors = FALSE
)

model_summary_df <- bind_rows(model_summary)
fixed_effects_df <- if (length(fixed_effects) > 0) bind_rows(fixed_effects) else empty_fixed
random_effects_df <- if (length(random_effects) > 0) bind_rows(random_effects) else empty_random
fit_statistics_df <- if (length(fit_statistics) > 0) bind_rows(fit_statistics) else empty_fit
predictions_df <- if (length(predictions) > 0) bind_rows(predictions) else empty_prediction

utils::write.csv(model_summary_df, file.path(tables_dir, "alpha_lmm_model_summary.csv"), row.names = FALSE)
utils::write.csv(fixed_effects_df, file.path(tables_dir, "alpha_lmm_fixed_effects.csv"), row.names = FALSE)
utils::write.csv(random_effects_df, file.path(tables_dir, "alpha_lmm_random_effects.csv"), row.names = FALSE)
utils::write.csv(fit_statistics_df, file.path(tables_dir, "alpha_lmm_model_fit_statistics.csv"), row.names = FALSE)
utils::write.csv(predictions_df, file.path(tables_dir, "alpha_lmm_predictions.csv"), row.names = FALSE)

cat("Modeled metrics:", paste(model_summary_df$metric[model_summary_df$status %in% c("modeled", "modeled_singular")], collapse = ", "), "\n")
cat("Alpha-diversity LMM fitting finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
