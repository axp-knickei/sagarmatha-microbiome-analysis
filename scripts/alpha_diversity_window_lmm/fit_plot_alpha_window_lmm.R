# Fit and plot alpha-diversity LMMs across predefined biological windows.
#
# This is a sample-level alpha-diversity workflow, not a per-MAG abundance
# workflow. Sample alpha-diversity values are averaged within each subject and
# predefined W1-W6 biological window, then modeled as:
#   alpha_value ~ window_index + (1 | subjectID)

required_packages <- c("dplyr", "tidyr", "ggplot2", "lme4", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
})

repo_root <- getwd()
alpha_path <- Sys.getenv(
  "BLACKICE_ALPHA_METADATA_PATH",
  unset = "results/alpha_diversity/tables/alpha_diversity_with_metadata.csv"
)
windows_path <- Sys.getenv(
  "BLACKICE_WINDOW_DEFINITION_PATH",
  unset = "results/metadata_expedition_windows/proposed_biological_windows.csv"
)
output_dir <- Sys.getenv("BLACKICE_ALPHA_WINDOW_LMM_OUTPUT_DIR", unset = "results/alpha_diversity_window_lmm")

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
logs_dir <- file.path(output_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("docs", "reports"), recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(logs_dir, "alpha_window_lmm.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Alpha-diversity window LMM workflow started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Working directory:", repo_root, "\n")

required_files <- c(alpha_path, windows_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
}

alpha <- utils::read.csv(alpha_path, check.names = FALSE, stringsAsFactors = FALSE)
windows <- utils::read.csv(windows_path, check.names = FALSE, stringsAsFactors = FALSE)

required_alpha_columns <- c("sampleID", "subjectID", "Timepoint")
required_window_columns <- c("window_label", "timepoints", "forward_timepoint_path")
missing_alpha_columns <- setdiff(required_alpha_columns, names(alpha))
missing_window_columns <- setdiff(required_window_columns, names(windows))
if (length(missing_alpha_columns) > 0) {
  stop("Alpha-diversity table is missing required columns: ", paste(missing_alpha_columns, collapse = ", "), call. = FALSE)
}
if (length(missing_window_columns) > 0) {
  stop("Window table is missing required columns: ", paste(missing_window_columns, collapse = ", "), call. = FALSE)
}

metric_columns <- c("Shannon", "Simpson", "Observed_richness", "Chao1", "ACE")
metric_labels <- c(
  Shannon = "Shannon diversity",
  Simpson = "Simpson diversity",
  Observed_richness = "Observed richness",
  Chao1 = "Chao1 richness",
  ACE = "ACE richness"
)
available_metric_columns <- intersect(metric_columns, names(alpha))
if (length(available_metric_columns) == 0) {
  stop("No expected alpha-diversity metric columns were found.", call. = FALSE)
}

for (metric in available_metric_columns) {
  alpha[[metric]] <- suppressWarnings(as.numeric(alpha[[metric]]))
}

alpha <- alpha %>%
  mutate(
    sampleID = as.character(sampleID),
    subjectID = as.character(subjectID),
    Timepoint = as.character(Timepoint)
  )

if (anyDuplicated(alpha$sampleID) > 0) {
  stop("Alpha-diversity table contains duplicated sampleID values.", call. = FALSE)
}
if (any(is.na(alpha$subjectID) | !nzchar(alpha$subjectID))) {
  stop("Alpha-diversity table contains missing subjectID values.", call. = FALSE)
}

split_timepoints <- function(x) {
  trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
}

short_label <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("to", "to", x, fixed = TRUE)
  x
}

windows <- windows %>%
  mutate(
    window_index = paste0("W", seq_len(n())),
    window_display = paste0(window_index, ": ", window_label),
    window_axis_label = paste0(window_index, "\n", short_label(window_label)),
    timepoint_list = strsplit(gsub("\\s+", "", timepoints), ","),
    expected_timepoints = vapply(timepoint_list, length, integer(1))
  )

window_membership <- lapply(seq_len(nrow(windows)), function(i) {
  data.frame(
    window_index = windows$window_index[[i]],
    window_label = windows$window_label[[i]],
    window_display = windows$window_display[[i]],
    window_axis_label = windows$window_axis_label[[i]],
    window_order = i,
    timepoint = windows$timepoint_list[[i]],
    expected_timepoints = windows$expected_timepoints[[i]],
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

window_definitions <- windows %>%
  select(
    window_index,
    window_label,
    forward_timepoint_path,
    timepoints,
    date_span,
    sampling_altitude_span,
    n_complete_subjects_all_three_timepoints,
    lmm_sample_subject_support,
    biological_interpretation,
    figure_informed_altitude_exposure_interpretation,
    caveats
  )
utils::write.csv(window_definitions, file.path(tables_dir, "alpha_window_definitions.csv"), row.names = FALSE)

alpha_long <- alpha %>%
  select(sampleID, subjectID, Timepoint, all_of(available_metric_columns)) %>%
  pivot_longer(cols = all_of(available_metric_columns), names_to = "metric", values_to = "alpha_value") %>%
  mutate(metric_label = unname(metric_labels[metric])) %>%
  inner_join(window_membership, by = c("Timepoint" = "timepoint"), relationship = "many-to-many")

subject_window_alpha <- alpha_long %>%
  group_by(subjectID, metric, metric_label, window_index, window_label, window_display, window_axis_label, window_order, expected_timepoints) %>%
  summarise(
    alpha_value = mean(alpha_value, na.rm = TRUE),
    n_samples_in_window = sum(!is.na(alpha_value)),
    complete_window = n_samples_in_window >= dplyr::first(expected_timepoints),
    .groups = "drop"
  ) %>%
  filter(is.finite(alpha_value)) %>%
  mutate(
    window_index = factor(window_index, levels = windows$window_index, ordered = TRUE),
    window_axis_label = factor(window_axis_label, levels = windows$window_axis_label, ordered = TRUE)
  )

utils::write.csv(subject_window_alpha, file.path(tables_dir, "alpha_diversity_by_subject_window.csv"), row.names = FALSE)

metric_availability <- subject_window_alpha %>%
  group_by(metric, metric_label) %>%
  summarise(
    n_subject_window_values = n(),
    n_subjects = n_distinct(subjectID),
    n_windows_observed = n_distinct(window_index),
    available_for_lmm = n_subject_window_values >= 10 && n_subjects >= 3 && n_windows_observed >= 2 && stats::var(alpha_value, na.rm = TRUE) > 0,
    reason = case_when(
      n_subject_window_values < 10 ~ "fewer_than_10_subject_window_values",
      n_subjects < 3 ~ "fewer_than_3_subjects",
      n_windows_observed < 2 ~ "fewer_than_2_windows",
      stats::var(alpha_value, na.rm = TRUE) <= 0 ~ "zero_variance",
      TRUE ~ "available"
    ),
    .groups = "drop"
  )

missing_metrics <- setdiff(metric_columns, unique(metric_availability$metric))
if (length(missing_metrics) > 0) {
  metric_availability <- bind_rows(
    metric_availability,
    data.frame(
      metric = missing_metrics,
      metric_label = unname(metric_labels[missing_metrics]),
      n_subject_window_values = 0L,
      n_subjects = 0L,
      n_windows_observed = 0L,
      available_for_lmm = FALSE,
      reason = "metric_missing_or_all_NA",
      stringsAsFactors = FALSE
    )
  )
}
utils::write.csv(metric_availability, file.path(tables_dir, "alpha_window_metric_availability.csv"), row.names = FALSE)

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

window_pairs <- utils::combn(windows$window_index, 2, simplify = FALSE)
window_pair_table <- lapply(window_pairs, function(pair) {
  before <- windows %>% filter(window_index == pair[[1]])
  after <- windows %>% filter(window_index == pair[[2]])
  data.frame(
    comparison = paste(pair, collapse = "_vs_"),
    comparison_label = paste(pair, collapse = " vs "),
    before_window = pair[[1]],
    after_window = pair[[2]],
    before_label = before$window_label,
    after_label = after$window_label,
    before_axis_label = before$window_axis_label,
    after_axis_label = after$window_axis_label,
    adjacent_pair = as.integer(sub("^W", "", pair[[2]])) - as.integer(sub("^W", "", pair[[1]])) == 1,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()
utils::write.csv(window_pair_table, file.path(tables_dir, "alpha_window_pair_definitions.csv"), row.names = FALSE)

fixed_effects <- list()
estimated_means <- list()
pairwise_lmm <- list()
model_summary <- list()

for (metric in metric_columns) {
  available_row <- metric_availability %>% filter(metric == !!metric)
  available <- nrow(available_row) == 1 && isTRUE(available_row$available_for_lmm)
  if (!available) {
    model_summary[[metric]] <- data.frame(
      metric = metric,
      metric_label = unname(metric_labels[metric]),
      status = "not_modeled",
      reason = if (nrow(available_row) == 1) available_row$reason else "missing_availability_row",
      n_subject_window_values = if (nrow(available_row) == 1) available_row$n_subject_window_values else 0L,
      n_subjects = if (nrow(available_row) == 1) available_row$n_subjects else 0L,
      n_windows_observed = if (nrow(available_row) == 1) available_row$n_windows_observed else 0L,
      model_formula = NA_character_,
      singular_fit = NA,
      convergence_warnings = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }

  metric_data <- subject_window_alpha %>%
    filter(metric == !!metric) %>%
    mutate(
      subjectID = factor(subjectID),
      window_index = factor(as.character(window_index), levels = windows$window_index, ordered = FALSE),
      alpha_value = as.numeric(alpha_value)
    )

  formula_text <- "alpha_value ~ window_index + (1 | subjectID)"
  fit_result <- tryCatch(
    fit_lmer_model(stats::as.formula(formula_text), metric_data),
    error = function(e) list(fit = NULL, warnings = conditionMessage(e))
  )

  if (is.null(fit_result$fit)) {
    model_summary[[metric]] <- data.frame(
      metric = metric,
      metric_label = unname(metric_labels[metric]),
      status = "model_failed",
      reason = paste(fit_result$warnings, collapse = "; "),
      n_subject_window_values = nrow(metric_data),
      n_subjects = n_distinct(metric_data$subjectID),
      n_windows_observed = n_distinct(metric_data$window_index),
      model_formula = formula_text,
      singular_fit = NA,
      convergence_warnings = paste(fit_result$warnings, collapse = "; "),
      stringsAsFactors = FALSE
    )
    next
  }

  fit <- fit_result$fit
  coef_table <- as.data.frame(summary(fit)$coefficients)
  coef_table$term <- rownames(coef_table)
  names(coef_table) <- sub("^Estimate$", "estimate", names(coef_table))
  names(coef_table) <- sub("^Std\\. Error$", "standard_error", names(coef_table))
  names(coef_table) <- sub("^df$", "degrees_of_freedom", names(coef_table))
  names(coef_table) <- sub("^t value$", "t_statistic", names(coef_table))
  names(coef_table) <- sub("^Pr\\(>\\|t\\|\\)$", "p_value", names(coef_table))
  if (!"degrees_of_freedom" %in% names(coef_table)) coef_table$degrees_of_freedom <- NA_real_
  if (!"p_value" %in% names(coef_table)) {
    coef_table$p_value <- 2 * stats::pnorm(-abs(coef_table$t_statistic))
    p_method <- "normal_approximation_lme4_fallback"
  } else {
    p_method <- "lmerTest_Satterthwaite"
  }

  fixed_effects[[metric]] <- coef_table %>%
    transmute(
      metric = metric,
      metric_label = unname(metric_labels[metric]),
      term = term,
      estimate = estimate,
      standard_error = standard_error,
      degrees_of_freedom = degrees_of_freedom,
      t_statistic = t_statistic,
      p_value = p_value,
      p_value_method = p_method,
      model_formula = formula_text,
      n_subject_window_values = nrow(metric_data),
      n_subjects = n_distinct(metric_data$subjectID)
    )

  newdata <- data.frame(window_index = factor(windows$window_index, levels = windows$window_index, ordered = FALSE))
  design <- stats::model.matrix(stats::delete.response(stats::terms(fit, fixed.only = TRUE)), newdata)
  beta <- lme4::fixef(fit)
  vcov_beta <- as.matrix(stats::vcov(fit))
  fit_values <- as.numeric(design %*% beta)
  se_values <- sqrt(pmax(0, diag(design %*% vcov_beta %*% t(design))))

  estimated_means[[metric]] <- data.frame(
    metric = metric,
    metric_label = unname(metric_labels[metric]),
    window_index = windows$window_index,
    window_label = windows$window_label,
    window_axis_label = windows$window_axis_label,
    estimate = fit_values,
    standard_error = se_values,
    confidence_interval_low = fit_values - 1.96 * se_values,
    confidence_interval_high = fit_values + 1.96 * se_values,
    stringsAsFactors = FALSE
  )

  pairwise_lmm[[metric]] <- lapply(seq_len(nrow(window_pair_table)), function(i) {
    pair <- window_pair_table[i, ]
    before_i <- match(pair$before_window, windows$window_index)
    after_i <- match(pair$after_window, windows$window_index)
    contrast <- design[after_i, , drop = FALSE] - design[before_i, , drop = FALSE]
    estimate <- as.numeric(contrast %*% beta)
    se <- sqrt(as.numeric(contrast %*% vcov_beta %*% t(contrast)))
    statistic <- estimate / se
    p_value <- 2 * stats::pnorm(-abs(statistic))
    data.frame(
      metric = metric,
      metric_label = unname(metric_labels[metric]),
      comparison = pair$comparison,
      comparison_label = pair$comparison_label,
      before_window = pair$before_window,
      after_window = pair$after_window,
      before_label = pair$before_label,
      after_label = pair$after_label,
      adjacent_pair = pair$adjacent_pair,
      estimate_after_minus_before = estimate,
      standard_error = se,
      statistic = statistic,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    mutate(p_adj_bh_within_metric = stats::p.adjust(p_value, method = "BH"))

  model_summary[[metric]] <- data.frame(
    metric = metric,
    metric_label = unname(metric_labels[metric]),
    status = if (lme4::isSingular(fit, tol = 1e-4)) "modeled_singular" else "modeled",
    reason = if (lme4::isSingular(fit, tol = 1e-4)) "singular_random_intercept" else "model_fit",
    n_subject_window_values = nrow(metric_data),
    n_subjects = n_distinct(metric_data$subjectID),
    n_windows_observed = n_distinct(metric_data$window_index),
    model_formula = formula_text,
    singular_fit = lme4::isSingular(fit, tol = 1e-4),
    convergence_warnings = paste(fit_result$warnings, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

fixed_effects_df <- if (length(fixed_effects) > 0) bind_rows(fixed_effects) else data.frame()
estimated_means_df <- if (length(estimated_means) > 0) bind_rows(estimated_means) else data.frame()
pairwise_lmm_df <- if (length(pairwise_lmm) > 0) bind_rows(pairwise_lmm) else data.frame()
model_summary_df <- bind_rows(model_summary)

utils::write.csv(model_summary_df, file.path(tables_dir, "alpha_window_lmm_model_summary.csv"), row.names = FALSE)
utils::write.csv(fixed_effects_df, file.path(tables_dir, "alpha_window_lmm_fixed_effects.csv"), row.names = FALSE)
utils::write.csv(estimated_means_df, file.path(tables_dir, "alpha_window_lmm_estimated_means.csv"), row.names = FALSE)
utils::write.csv(pairwise_lmm_df, file.path(tables_dir, "alpha_window_lmm_pairwise_contrasts.csv"), row.names = FALSE)

pair_observed <- lapply(seq_len(nrow(window_pair_table)), function(i) {
  pair <- window_pair_table[i, ]
  before <- subject_window_alpha %>%
    filter(as.character(window_index) == pair$before_window) %>%
    rename(before_value = alpha_value, before_complete_window = complete_window)
  after <- subject_window_alpha %>%
    filter(as.character(window_index) == pair$after_window) %>%
    rename(after_value = alpha_value, after_complete_window = complete_window)

  before %>%
    select(subjectID, metric, metric_label, before_value, before_complete_window) %>%
    inner_join(
      after %>% select(subjectID, metric, metric_label, after_value, after_complete_window),
      by = c("subjectID", "metric", "metric_label")
    ) %>%
    mutate(
      comparison = pair$comparison,
      comparison_label = pair$comparison_label,
      before_window = pair$before_window,
      after_window = pair$after_window,
      before_axis_label = pair$before_axis_label,
      after_axis_label = pair$after_axis_label,
      adjacent_pair = pair$adjacent_pair,
      delta = after_value - before_value
    )
}) %>%
  bind_rows()

utils::write.csv(pair_observed, file.path(tables_dir, "alpha_window_pair_observed_values.csv"), row.names = FALSE)

paired_tests <- pair_observed %>%
  group_by(metric, metric_label, comparison, comparison_label, before_window, after_window, adjacent_pair) %>%
  summarise(
    n_paired_subjects = n_distinct(subjectID),
    mean_before = mean(before_value, na.rm = TRUE),
    mean_after = mean(after_value, na.rm = TRUE),
    mean_delta = mean(delta, na.rm = TRUE),
    median_delta = stats::median(delta, na.rm = TRUE),
    paired_t_p_value = if (n() >= 3 && stats::sd(delta, na.rm = TRUE) > 0) stats::t.test(after_value, before_value, paired = TRUE)$p.value else NA_real_,
    wilcoxon_p_value = if (n() >= 3 && any(delta != 0, na.rm = TRUE)) {
      suppressWarnings(stats::wilcox.test(after_value, before_value, paired = TRUE, exact = FALSE)$p.value)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  group_by(metric) %>%
  mutate(
    paired_t_p_adj_bh_within_metric = stats::p.adjust(paired_t_p_value, method = "BH"),
    wilcoxon_p_adj_bh_within_metric = stats::p.adjust(wilcoxon_p_value, method = "BH")
  ) %>%
  ungroup()

utils::write.csv(paired_tests, file.path(tables_dir, "alpha_window_paired_tests.csv"), row.names = FALSE)

theme_journal <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black", size = rel(0.9)),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.spacing = unit(0.85, "lines"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "gray25", hjust = 0),
      plot.caption = element_text(color = "gray35", hjust = 0, size = rel(0.82)),
      plot.margin = margin(8, 8, 8, 8)
    )
}

save_plot <- function(plot, stem, width, height) {
  ggplot2::ggsave(file.path(figures_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(file.path(figures_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 600)
  ggplot2::ggsave(file.path(figures_dir, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite::svglite)
}

available_modeled_metrics <- model_summary_df$metric[model_summary_df$status %in% c("modeled", "modeled_singular")]
plot_subject_window <- subject_window_alpha %>%
  filter(metric %in% available_modeled_metrics) %>%
  mutate(
    metric_label = factor(metric_label, levels = unname(metric_labels[metric_columns])),
    window_axis_label = factor(as.character(window_axis_label), levels = windows$window_axis_label, ordered = TRUE)
  )
plot_estimated_means <- estimated_means_df %>%
  mutate(
    metric_label = factor(metric_label, levels = unname(metric_labels[metric_columns])),
    window_axis_label = factor(window_axis_label, levels = windows$window_axis_label, ordered = TRUE)
  )

if (nrow(plot_subject_window) > 0) {
  window_lmm_plot <- ggplot() +
    geom_line(
      data = plot_subject_window,
      aes(x = window_axis_label, y = alpha_value, group = subjectID),
      color = "#8C8C8C",
      linewidth = 0.25,
      alpha = 0.45
    ) +
    geom_point(
      data = plot_subject_window,
      aes(x = window_axis_label, y = alpha_value),
      color = "#0072B2",
      size = 1.15,
      alpha = 0.62
    ) +
    geom_errorbar(
      data = plot_estimated_means,
      aes(x = window_axis_label, ymin = confidence_interval_low, ymax = confidence_interval_high),
      width = 0.12,
      color = "#004B75",
      linewidth = 0.45
    ) +
    geom_point(
      data = plot_estimated_means,
      aes(x = window_axis_label, y = estimate),
      color = "#004B75",
      size = 2.1
    ) +
    geom_line(
      data = plot_estimated_means,
      aes(x = window_axis_label, y = estimate, group = metric),
      color = "#004B75",
      linewidth = 0.65
    ) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 1) +
    labs(
      x = "Biological window label",
      y = "Alpha diversity",
      title = "Alpha diversity across predefined biological windows",
      subtitle = "Points and thin lines show subject-level window means; bold points and intervals show LMM estimated means with 95% CI.",
      caption = "Window values are means across available samples in each W1-W6 biological window. Chao1 and ACE are omitted when unavailable."
    ) +
    theme_journal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_plot(window_lmm_plot, "alpha_window_lmm_estimated_means_by_window", 8.6, 9.0)
}

pair_plot_data <- pair_observed %>%
  filter(metric %in% available_modeled_metrics) %>%
  select(subjectID, metric, metric_label, comparison, comparison_label, before_window, after_window, before_axis_label, after_axis_label, before_value, after_value) %>%
  pivot_longer(cols = c(before_value, after_value), names_to = "position", values_to = "alpha_value") %>%
  mutate(
    pair_subject_id = paste(subjectID, comparison, metric, sep = "__"),
    window_axis_label = ifelse(position == "before_value", before_axis_label, after_axis_label),
    window_tick_label = factor(
      ifelse(position == "before_value", before_window, after_window),
      levels = windows$window_index,
      ordered = TRUE
    ),
    before_after = factor(ifelse(position == "before_value", "Before", "After"), levels = c("Before", "After")),
    metric_label = factor(metric_label, levels = unname(metric_labels[metric_columns]))
  )

make_pair_plot <- function(data, title, subtitle) {
  ggplot(data, aes(x = window_tick_label, y = alpha_value, group = pair_subject_id)) +
    geom_line(color = "#8C8C8C", alpha = 0.42, linewidth = 0.24) +
    geom_point(color = "#0072B2", alpha = 0.68, size = 1.05) +
    stat_summary(aes(group = 1), fun = mean, geom = "line", color = "#D55E00", linewidth = 0.8) +
    stat_summary(aes(group = 1), fun = mean, geom = "point", color = "#D55E00", size = 1.9) +
    facet_grid(metric_label ~ comparison_label, scales = "free", space = "free_x") +
    labs(
      x = "Biological window label",
      y = "Alpha diversity",
      title = title,
      subtitle = subtitle,
      caption = "Thin paired lines connect the same subject between two biological-window means. Orange summaries show observed paired means; LMM pairwise estimates are in the contrast table."
    ) +
    theme_journal(base_size = 8)
}

if (nrow(pair_plot_data) > 0) {
  all_pair_plot <- make_pair_plot(
    pair_plot_data,
    "Paired alpha-diversity comparisons across all biological-window pairs",
    "All W1-W6 pairwise comparisons: W1 vs W2 through W5 vs W6."
  )
  save_plot(all_pair_plot, "alpha_window_lmm_paired_all_window_pairs", 14.5, 8.6)

  adjacent_pair_plot <- make_pair_plot(
    pair_plot_data %>% filter(comparison %in% window_pair_table$comparison[window_pair_table$adjacent_pair]),
    "Paired alpha-diversity comparisons across adjacent biological windows",
    "Adjacent comparisons only: W1 vs W2, W2 vs W3, W3 vs W4, W4 vs W5, W5 vs W6."
  )
  save_plot(adjacent_pair_plot, "alpha_window_lmm_paired_adjacent_window_pairs", 10.5, 8.0)

  for (metric in available_modeled_metrics) {
    metric_plot <- pair_plot_data %>%
      filter(metric == !!metric) %>%
      ggplot(aes(x = window_tick_label, y = alpha_value, group = pair_subject_id)) +
      geom_line(color = "#8C8C8C", alpha = 0.42, linewidth = 0.24) +
      geom_point(color = "#0072B2", alpha = 0.68, size = 1.05) +
      stat_summary(aes(group = 1), fun = mean, geom = "line", color = "#D55E00", linewidth = 0.8) +
      stat_summary(aes(group = 1), fun = mean, geom = "point", color = "#D55E00", size = 1.9) +
      facet_wrap(~ comparison_label, scales = "free", ncol = 5) +
      labs(
        x = "Biological window label",
        y = "Alpha diversity",
        title = paste0(unname(metric_labels[metric]), " paired comparisons across biological-window pairs"),
        subtitle = "All W1-W6 pairwise comparisons.",
        caption = "Thin paired lines connect subject-level window means. Orange summaries show observed paired means."
      ) +
      theme_journal(base_size = 8)
    save_plot(metric_plot, paste0("alpha_window_lmm_paired_", tolower(metric), "_all_window_pairs"), 13.0, 6.8)
  }
}

table_files <- sort(basename(list.files(tables_dir, pattern = "\\.csv$")))
figure_files <- sort(basename(list.files(figures_dir, pattern = "\\.(pdf|png|svg)$")))
modeled_metrics <- model_summary_df %>% filter(status %in% c("modeled", "modeled_singular"))
not_modeled_metrics <- model_summary_df %>% filter(!status %in% c("modeled", "modeled_singular"))

report_lines <- c(
  "# Alpha-Diversity Biological-Window LMM Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Inputs",
  "",
  paste0("- Alpha-diversity table: `", alpha_path, "`"),
  paste0("- Window definition table: `", windows_path, "`"),
  "",
  "## Window Labels",
  "",
  paste0("- ", windows$window_index, ": `", windows$window_label, "` (", windows$forward_timepoint_path, ")"),
  "",
  "## Model",
  "",
  "- Window-level subject means were modeled with `alpha_value ~ window_index + (1 | subjectID)`.",
  "- Pairwise LMM contrasts compare estimated window means for all W1-W6 combinations.",
  "",
  "## Metrics",
  "",
  paste0("- Modeled: ", paste(modeled_metrics$metric_label, collapse = "; ")),
  if (nrow(not_modeled_metrics) > 0) paste0("- Not modeled: ", paste(paste0(not_modeled_metrics$metric_label, " (", not_modeled_metrics$reason, ")"), collapse = "; ")) else "- Not modeled: none",
  "",
  "## Generated Tables",
  "",
  paste0("- `results/alpha_diversity_window_lmm/tables/", table_files, "`"),
  "",
  "## Generated Figures",
  "",
  paste0("- `results/alpha_diversity_window_lmm/figures/", figure_files, "`"),
  "",
  "## HPC Command",
  "",
  "```bash",
  "qsub hpc/job_alpha_diversity_window_lmm.pbs",
  "```",
  "",
  "## Cautions",
  "",
  "- This is alpha-diversity analysis by biological window, not per-MAG abundance modeling.",
  "- Subject-window alpha values are means across available samples in each window; W1 and W6 have limited complete-subject support.",
  "- Chao1 and ACE remain unavailable when upstream alpha-diversity results did not calculate them from relative-abundance-like input.",
  "- All pairwise window contrasts should be interpreted as a defined comparison family, not as independent discovery claims."
)

writeLines(report_lines, file.path("docs", "reports", "alpha_diversity_window_lmm_report.md"))
utils::capture.output(utils::sessionInfo(), file = file.path(logs_dir, "session_info.txt"))

cat("Modeled metrics:", paste(modeled_metrics$metric, collapse = ", "), "\n")
cat("Output directory:", output_dir, "\n")
cat("Alpha-diversity window LMM workflow finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
