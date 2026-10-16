#!/usr/bin/env Rscript
# Purpose: Quantify within-subject distance-to-baseline and adjacent-window
# changes, then fit the configured longitudinal distance models.
# Inputs: Prepared metadata and distance matrices from stages 02 and 03.
# Outputs: Longitudinal summaries, diagnostics, predictions, and fitted-model
# objects under the configured statistics directory.
# Dependencies: Shared beta-diversity helpers and mixed-model packages declared
# in environment.yml.
# Usage: Rscript analysis/beta_diversity/06_run_longitudinal_analysis.R --config config/beta_diversity.yml
# Workflow: Stage 06; supplies longitudinal evidence to figures and reporting.

source("analysis/beta_diversity/utils/io_helpers.R")
source("analysis/beta_diversity/utils/provenance_helpers.R")
args <- parse_beta_args(); config <- load_beta_config(args$config); ensure_beta_dirs(config)
require_beta_packages("lme4")

with_stage_log(config, "06_run_longitudinal_analysis", {
  data <- readRDS(beta_output_dir(config, "prepared_data", "prepared_data.rds")); metadata <- data$metadata
  files <- list.files(beta_output_dir(config, "distances"), pattern = "^distance_.*\\.rds$", full.names = TRUE)
  files <- files[!grepl("_transformed\\.rds$", files)]
  baseline_all <- list(); adjacent_all <- list(); models <- list(); model_rows <- list(); diagnostic_rows <- list(); prediction_rows <- list()
  for (file in files) {
    metric <- sub("\\.rds$", "", sub("^distance_", "", basename(file))); d <- as.matrix(readRDS(file))
    md <- metadata[match(rownames(d), metadata$sample_id), , drop = FALSE]
    baseline <- do.call(rbind, lapply(split(md, md$subject_id), function(subject) {
      subject <- subject[order(subject$sample_date, subject$sample_id), , drop = FALSE]
      ref <- subject$sample_id[1]
      data.frame(metric = metric, sample_id = subject$sample_id, subject_id = subject$subject_id,
                 timepoint = subject$timepoint, days_from_subject_baseline = subject$days_from_subject_baseline,
                 baseline_sample_id = ref, distance_to_baseline = d[subject$sample_id, ref])
    }))
    baseline_all[[metric]] <- baseline
    model_data <- baseline[baseline$days_from_subject_baseline > 0, , drop = FALSE]
    if (nrow(model_data) >= 8 && length(unique(model_data$subject_id)) >= 3) {
      fit <- lme4::lmer(distance_to_baseline ~ days_from_subject_baseline + (1 | subject_id), data = model_data, REML = TRUE)
      models[[metric]] <- fit; estimate <- summary(fit)$coefficients["days_from_subject_baseline", ]
      ci <- suppressMessages(confint(fit, parm = "days_from_subject_baseline", method = "Wald"))
      model_rows[[metric]] <- data.frame(metric = metric, model = "distance_to_baseline ~ days_from_subject_baseline + (1|subject_id)",
        term = "days_from_subject_baseline", estimate = estimate["Estimate"], std_error = estimate["Std. Error"],
        conf_low = ci[1], conf_high = ci[2], sample_count = nrow(model_data), subject_count = length(unique(model_data$subject_id)),
        singular = lme4::isSingular(fit), interpretation_caution = "Exploratory trend; residual and influence diagnostics required before inference.")
      residual <- stats::residuals(fit); fitted <- stats::fitted(fit)
      diagnostic_rows[[metric]] <- data.frame(metric = metric, residual_mean = mean(residual), residual_sd = stats::sd(residual),
        shapiro_w = if (length(residual) >= 3 && length(residual) <= 5000) stats::shapiro.test(residual)$statistic else NA_real_,
        shapiro_p = if (length(residual) >= 3 && length(residual) <= 5000) stats::shapiro.test(residual)$p.value else NA_real_,
        abs_residual_fitted_correlation = suppressWarnings(stats::cor(abs(residual), fitted)), singular = lme4::isSingular(fit))
      days <- seq(0, max(model_data$days_from_subject_baseline), length.out = 100)
      fixed <- lme4::fixef(fit); vc <- as.matrix(stats::vcov(fit)); design <- cbind(1, days)
      estimate_line <- as.vector(design %*% fixed); se_line <- sqrt(rowSums((design %*% vc) * design))
      prediction_rows[[metric]] <- data.frame(metric = metric, days_from_subject_baseline = days,
        estimate = estimate_line, conf_low = estimate_line - 1.96 * se_line, conf_high = estimate_line + 1.96 * se_line)
    }
    adjacent <- do.call(rbind, lapply(split(md, md$subject_id), function(subject) {
      subject <- subject[order(subject$sample_date, subject$sample_id), , drop = FALSE]
      if (nrow(subject) < 2) return(NULL)
      data.frame(metric = metric, subject_id = subject$subject_id[-1], from_sample = head(subject$sample_id, -1),
                 to_sample = tail(subject$sample_id, -1), from_timepoint = head(subject$timepoint, -1),
                 to_timepoint = tail(subject$timepoint, -1), transition_days = diff(subject$sample_date),
                 adjacent_distance = mapply(function(a, b) d[a, b], head(subject$sample_id, -1), tail(subject$sample_id, -1)))
    }))
    adjacent_all[[metric]] <- adjacent
  }
  write_tsv(do.call(rbind, baseline_all), beta_output_dir(config, "statistics", "distance_to_baseline.tsv"))
  write_tsv(do.call(rbind, adjacent_all), beta_output_dir(config, "statistics", "adjacent_window_distances.tsv"))
  write_tsv(if (length(model_rows)) do.call(rbind, model_rows) else data.frame(status = "SKIPPED_INSUFFICIENT_DATA"),
            beta_output_dir(config, "statistics", "longitudinal_lmm_results.tsv"))
  write_tsv(if (length(diagnostic_rows)) do.call(rbind, diagnostic_rows) else data.frame(status = "SKIPPED_INSUFFICIENT_DATA"),
            beta_output_dir(config, "statistics", "longitudinal_lmm_diagnostics.tsv"))
  write_tsv(if (length(prediction_rows)) do.call(rbind, prediction_rows) else data.frame(status = "SKIPPED_INSUFFICIENT_DATA"),
            beta_output_dir(config, "statistics", "longitudinal_lmm_estimated_trends.tsv"))
  saveRDS(models, beta_output_dir(config, "statistics", "longitudinal_lmm_models.rds"))
})
