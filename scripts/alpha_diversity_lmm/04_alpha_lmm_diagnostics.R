# Generate alpha-diversity LMM diagnostics and companion report.

required_packages <- c("dplyr", "ggplot2", "lme4", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lme4)
})

output_dir <- Sys.getenv("BLACKICE_ALPHA_LMM_OUTPUT_DIR", unset = "results/alpha_diversity_lmm")
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
diagnostics_dir <- file.path(output_dir, "diagnostics")
logs_dir <- file.path(output_dir, "logs")
report_dir <- file.path("docs", "reports")

dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(logs_dir, "alpha_lmm_diagnostics.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Alpha-diversity LMM diagnostics started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

required_inputs <- file.path(tables_dir, c(
  "alpha_lmm_input.csv",
  "alpha_lmm_model_summary.csv",
  "alpha_lmm_fixed_effects.csv",
  "alpha_lmm_model_fit_statistics.csv",
  "alpha_lmm_predictions.csv",
  "alpha_lmm_metric_availability.csv",
  "alpha_lmm_input_diagnostics.csv"
))
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required diagnostics input(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)
}

alpha <- utils::read.csv(file.path(tables_dir, "alpha_lmm_input.csv"), check.names = FALSE, stringsAsFactors = FALSE)
model_summary <- utils::read.csv(file.path(tables_dir, "alpha_lmm_model_summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
fixed_effects <- utils::read.csv(file.path(tables_dir, "alpha_lmm_fixed_effects.csv"), check.names = FALSE, stringsAsFactors = FALSE)
fit_statistics <- utils::read.csv(file.path(tables_dir, "alpha_lmm_model_fit_statistics.csv"), check.names = FALSE, stringsAsFactors = FALSE)
availability <- utils::read.csv(file.path(tables_dir, "alpha_lmm_metric_availability.csv"), check.names = FALSE, stringsAsFactors = FALSE)
input_diagnostics <- utils::read.csv(file.path(tables_dir, "alpha_lmm_input_diagnostics.csv"), check.names = FALSE, stringsAsFactors = FALSE)

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
for (metric in metric_columns) alpha[[metric]] <- suppressWarnings(as.numeric(alpha[[metric]]))
for (column in intersect(c("days_from_baseline", "alpha_lmm_time", "altitude", "age"), names(alpha))) {
  alpha[[column]] <- suppressWarnings(as.numeric(alpha[[column]]))
}

time_variable <- unique(model_summary$time_variable[model_summary$status %in% c("modeled", "modeled_singular")])[1]
if (is.na(time_variable) || !nzchar(time_variable)) time_variable <- "alpha_lmm_time"

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

theme_diag <- function() {
  theme_classic(base_size = 9, base_family = "sans") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

save_diag_plot <- function(plot, stem, width = 5.2, height = 4.0) {
  ggsave(file.path(diagnostics_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(diagnostics_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 450)
  ggsave(file.path(diagnostics_dir, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite::svglite)
}

diagnostic_rows <- list()

for (metric in metric_columns) {
  summary_row <- model_summary %>% filter(metric == !!metric)
  if (nrow(summary_row) != 1 || !summary_row$status %in% c("modeled", "modeled_singular")) {
    diagnostic_rows[[metric]] <- data.frame(
      metric = metric,
      metric_label = unname(metric_labels[metric]),
      model_type = "main",
      diagnostic_status = "not_modeled",
      n_samples = sum(!is.na(alpha[[metric]])),
      n_subjects = dplyr::n_distinct(alpha$subjectID[!is.na(alpha[[metric]])]),
      singular_fit = NA,
      residual_mean = NA_real_,
      residual_sd = NA_real_,
      shapiro_wilk_p_value = NA_real_,
      random_intercept_variance = NA_real_,
      convergence_warnings = NA_character_,
      notes = if (nrow(summary_row) == 1) summary_row$reason else "missing_model_summary",
      stringsAsFactors = FALSE
    )
    next
  }

  metric_data <- alpha %>%
    filter(!is.na(.data[[metric]]), !is.na(.data[[time_variable]]), !is.na(subjectID)) %>%
    mutate(metric_value = .data[[metric]], subjectID = factor(subjectID))

  formula_text <- paste0("metric_value ~ ", time_variable, " + (1 | subjectID)")
  fit_result <- tryCatch(
    fit_lmer_model(stats::as.formula(formula_text), metric_data),
    error = function(e) list(fit = NULL, warnings = conditionMessage(e))
  )

  if (is.null(fit_result$fit)) {
    diagnostic_rows[[metric]] <- data.frame(
      metric = metric,
      metric_label = unname(metric_labels[metric]),
      model_type = "main",
      diagnostic_status = "diagnostic_refit_failed",
      n_samples = nrow(metric_data),
      n_subjects = dplyr::n_distinct(metric_data$subjectID),
      singular_fit = NA,
      residual_mean = NA_real_,
      residual_sd = NA_real_,
      shapiro_wilk_p_value = NA_real_,
      random_intercept_variance = NA_real_,
      convergence_warnings = paste(fit_result$warnings, collapse = "; "),
      notes = "Could not refit model for diagnostics.",
      stringsAsFactors = FALSE
    )
    next
  }

  fit <- fit_result$fit
  residual_data <- data.frame(
    fitted = stats::fitted(fit),
    residual = stats::resid(fit),
    sampleID = metric_data$sampleID,
    subjectID = metric_data$subjectID,
    stringsAsFactors = FALSE
  )

  resid_plot <- ggplot(residual_data, aes(x = fitted, y = residual)) +
    geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
    geom_point(color = "#0072B2", alpha = 0.75, size = 1.7) +
    labs(
      x = "Fitted value",
      y = "Residual",
      title = paste0(unname(metric_labels[metric]), ": residuals versus fitted values")
    ) +
    theme_diag()
  save_diag_plot(resid_plot, paste0("alpha_lmm_residuals_vs_fitted_", tolower(metric)))

  qq_data <- data.frame(
    theoretical = stats::qqnorm(residual_data$residual, plot.it = FALSE)$x,
    sample = stats::qqnorm(residual_data$residual, plot.it = FALSE)$y
  )
  qq_plot <- ggplot(qq_data, aes(x = theoretical, y = sample)) +
    geom_abline(slope = 1, intercept = 0, color = "gray55", linewidth = 0.3) +
    geom_point(color = "#0072B2", alpha = 0.75, size = 1.7) +
    labs(
      x = "Theoretical quantiles",
      y = "Sample residual quantiles",
      title = paste0(unname(metric_labels[metric]), ": residual Q-Q plot")
    ) +
    theme_diag()
  save_diag_plot(qq_plot, paste0("alpha_lmm_qq_", tolower(metric)))

  random_intercepts <- data.frame(
    subjectID = rownames(lme4::ranef(fit)$subjectID),
    random_intercept = lme4::ranef(fit)$subjectID[[1]],
    stringsAsFactors = FALSE
  )
  random_plot <- ggplot(random_intercepts, aes(x = reorder(subjectID, random_intercept), y = random_intercept)) +
    geom_hline(yintercept = 0, color = "gray55", linewidth = 0.3) +
    geom_point(color = "#004B75", size = 1.8) +
    coord_flip() +
    labs(
      x = "Subject",
      y = "Random intercept",
      title = paste0(unname(metric_labels[metric]), ": subject random intercepts")
    ) +
    theme_diag()
  save_diag_plot(random_plot, paste0("alpha_lmm_random_intercepts_", tolower(metric)), width = 5.2, height = 3.6)

  varcorr <- as.data.frame(lme4::VarCorr(fit))
  random_var <- varcorr$vcov[varcorr$grp == "subjectID" & is.na(varcorr$var2)][1]
  shapiro_p <- if (length(residual_data$residual) >= 3 && length(residual_data$residual) <= 5000) {
    stats::shapiro.test(residual_data$residual)$p.value
  } else {
    NA_real_
  }

  diagnostic_rows[[metric]] <- data.frame(
    metric = metric,
    metric_label = unname(metric_labels[metric]),
    model_type = "main",
    diagnostic_status = "diagnostics_generated",
    n_samples = nrow(metric_data),
    n_subjects = dplyr::n_distinct(metric_data$subjectID),
    singular_fit = lme4::isSingular(fit, tol = 1e-4),
    residual_mean = mean(residual_data$residual, na.rm = TRUE),
    residual_sd = stats::sd(residual_data$residual, na.rm = TRUE),
    shapiro_wilk_p_value = shapiro_p,
    random_intercept_variance = random_var,
    convergence_warnings = paste(fit_result$warnings, collapse = "; "),
    notes = "Inspect residual, Q-Q, and random-intercept plots before biological interpretation.",
    stringsAsFactors = FALSE
  )
}

diagnostics_summary <- bind_rows(diagnostic_rows)
diagnostics_path <- file.path(tables_dir, "alpha_lmm_diagnostics_summary.csv")
utils::write.csv(diagnostics_summary, diagnostics_path, row.names = FALSE)

get_diag <- function(key, default = "NA") {
  hit <- input_diagnostics$value[input_diagnostics$item == key]
  if (length(hit) == 0) default else as.character(hit[[1]])
}

main_time_effects <- fixed_effects %>%
  filter(model_type == "main", term == get_diag("time_variable")) %>%
  transmute(
    metric = metric_label,
    estimate = signif(estimate, 4),
    standard_error = signif(standard_error, 4),
    p_value = signif(p_value, 4),
    n_samples = n_samples,
    n_subjects = n_subjects
  )

format_table_lines <- function(df) {
  if (nrow(df) == 0) {
    return("- No modeled time effects were available.")
  }
  apply(df, 1, function(row) {
    paste0(
      "- ", row[["metric"]],
      ": estimate ", row[["estimate"]],
      ", SE ", row[["standard_error"]],
      ", p ", row[["p_value"]],
      " (n = ", row[["n_samples"]], " samples; ", row[["n_subjects"]], " subjects)"
    )
  })
}

table_files <- sort(basename(list.files(tables_dir, pattern = "^alpha_lmm_.*\\.csv$", full.names = FALSE)))
figure_files <- sort(basename(list.files(figures_dir, pattern = "^alpha_lmm_.*\\.(pdf|png|svg)$", full.names = FALSE)))
diagnostic_files <- sort(basename(list.files(diagnostics_dir, pattern = "^alpha_lmm_.*\\.(pdf|png|svg)$", full.names = FALSE)))

report_lines <- c(
  "# Alpha-Diversity LMM Analysis Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Inputs",
  "",
  paste0("- Primary input table: `", get_diag("analysis_input"), "`"),
  paste0("- Alpha-diversity reused: ", get_diag("alpha_diversity_reused")),
  paste0("- Time variable: `", get_diag("time_variable"), "` (source: ", get_diag("time_source"), ")"),
  paste0("- Samples: ", get_diag("n_samples")),
  paste0("- Subjects: ", get_diag("n_subjects")),
  paste0("- Timepoints: ", get_diag("n_timepoints")),
  paste0("- Abundance scale decision from alpha-diversity workflow: ", get_diag("abundance_scale_decision")),
  "",
  "## Metrics Modeled",
  "",
  paste0("- Available for LMM: ", get_diag("metrics_available_for_lmm")),
  "- Chao1 and ACE are retained in the availability table even when not modeled.",
  "",
  "## Model Formulas",
  "",
  paste0("- Main model: `metric_value ~ ", get_diag("time_variable"), " + (1 | subjectID)`"),
  "- Optional adjusted model: `metric_value ~ days_from_baseline + altitude + age + gender + (1 | subjectID)` when metadata completeness and collinearity checks allow.",
  "",
  "## Main Fixed-Effect Results",
  "",
  format_table_lines(main_time_effects),
  "",
  "## Adjusted Model Status",
  "",
  paste0(
    "- ", model_summary$metric_label, ": ",
    model_summary$adjusted_model_status, " (", model_summary$adjusted_model_reason, ")"
  ),
  "",
  "## Generated Tables",
  "",
  paste0("- `results/alpha_diversity_lmm/tables/", table_files, "`"),
  "",
  "## Generated Figures",
  "",
  paste0("- `results/alpha_diversity_lmm/figures/", figure_files, "`"),
  "",
  "## Diagnostic Figures",
  "",
  paste0("- `results/alpha_diversity_lmm/diagnostics/", diagnostic_files, "`"),
  "",
  "## Diagnostics Summary",
  "",
  paste0(
    "- ", diagnostics_summary$metric_label, ": ",
    diagnostics_summary$diagnostic_status,
    "; singular fit = ", diagnostics_summary$singular_fit,
    "; residual SD = ", signif(diagnostics_summary$residual_sd, 4),
    "; Shapiro-Wilk p = ", signif(diagnostics_summary$shapiro_wilk_p_value, 4)
  ),
  "",
  "## Scientific Interpretation Cautions",
  "",
  "- Alpha-diversity LMMs are sample-level models, not per-MAG models.",
  "- Shannon and Simpson are continuous diversity indices and are usually suitable for LMM if residual diagnostics are acceptable.",
  "- Observed richness, Chao1, and ACE are richness-like or count-derived metrics; LMMs for these metrics should be treated as exploratory unless diagnostics support the approximation.",
  "- Chao1 and ACE should ideally be estimated from raw count-like abundance data. The alpha-diversity workflow classified the gut MAG abundance table as relative-abundance-like, so Chao1 and ACE were not modeled when unavailable.",
  "- A significant time effect does not identify which MAGs changed; it only indicates a change in within-sample diversity.",
  "- If altitude and days from baseline are strongly correlated, adjusted model terms should not be interpreted as independent environmental effects.",
  "",
  "## Local Execution",
  "",
  "```bash",
  "Rscript scripts/alpha_diversity_lmm/01_prepare_alpha_lmm_input.R",
  "Rscript scripts/alpha_diversity_lmm/02_fit_alpha_diversity_lmm.R",
  "Rscript scripts/alpha_diversity_lmm/03_plot_alpha_diversity_lmm_trends.R",
  "Rscript scripts/alpha_diversity_lmm/04_alpha_lmm_diagnostics.R",
  "```",
  "",
  "## HPC Submission",
  "",
  "```bash",
  "qsub hpc/job_alpha_diversity_lmm.pbs",
  "```"
)

report_path <- file.path(report_dir, "alpha_diversity_lmm_report.md")
writeLines(report_lines, report_path)

session_path <- file.path(logs_dir, "alpha_lmm_session_info.txt")
utils::capture.output(utils::sessionInfo(), file = session_path)

cat("Diagnostics summary:", diagnostics_path, "\n")
cat("Report:", report_path, "\n")
cat("Session info:", session_path, "\n")
cat("Alpha-diversity LMM diagnostics finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
