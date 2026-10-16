# Plot observed subject trajectories plus model-estimated alpha-diversity LMM trends.

required_packages <- c("dplyr", "tidyr", "ggplot2", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

output_dir <- Sys.getenv("BLACKICE_ALPHA_LMM_OUTPUT_DIR", unset = "results/alpha_diversity_lmm")
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
logs_dir <- file.path(output_dir, "logs")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(logs_dir, "alpha_lmm_plot.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Alpha-diversity LMM plotting started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

input_path <- file.path(tables_dir, "alpha_lmm_input.csv")
summary_path <- file.path(tables_dir, "alpha_lmm_model_summary.csv")
predictions_path <- file.path(tables_dir, "alpha_lmm_predictions.csv")
for (path in c(input_path, summary_path, predictions_path)) {
  if (!file.exists(path)) stop("Missing required plotting input: ", path, call. = FALSE)
}

alpha <- utils::read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
model_summary <- utils::read.csv(summary_path, check.names = FALSE, stringsAsFactors = FALSE)
predictions <- utils::read.csv(predictions_path, check.names = FALSE, stringsAsFactors = FALSE)

metric_columns <- c("Shannon", "Simpson", "Observed_richness", "Chao1", "ACE")
metric_labels <- c(
  Shannon = "Shannon diversity",
  Simpson = "Simpson diversity",
  Observed_richness = "Observed richness",
  Chao1 = "Chao1 richness",
  ACE = "ACE richness"
)
filename_stems <- c(
  Shannon = "alpha_lmm_observed_model_trend_shannon",
  Simpson = "alpha_lmm_observed_model_trend_simpson",
  Observed_richness = "alpha_lmm_observed_model_trend_observed_richness",
  Chao1 = "alpha_lmm_observed_model_trend_chao1",
  ACE = "alpha_lmm_observed_model_trend_ace"
)

for (metric in metric_columns) {
  alpha[[metric]] <- suppressWarnings(as.numeric(alpha[[metric]]))
}
alpha$subjectID <- factor(alpha$subjectID)
if ("Timepoint" %in% names(alpha)) alpha$Timepoint <- factor(alpha$Timepoint, ordered = TRUE)
if ("days_from_baseline" %in% names(alpha)) alpha$days_from_baseline <- suppressWarnings(as.numeric(alpha$days_from_baseline))
if ("alpha_lmm_time" %in% names(alpha)) alpha$alpha_lmm_time <- suppressWarnings(as.numeric(alpha$alpha_lmm_time))

time_variable <- unique(predictions$time_variable[predictions$model_type == "main"])[1]
if (is.na(time_variable) || !nzchar(time_variable)) time_variable <- "alpha_lmm_time"
x_label <- if (identical(time_variable, "days_from_baseline")) "Days from baseline" else "Ordered timepoint"

theme_journal <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      legend.position = "none",
      panel.spacing = unit(1.1, "lines"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "gray25", hjust = 0),
      plot.caption = element_text(color = "gray35", hjust = 0, size = rel(0.82)),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_plot <- function(plot, filename, width, height) {
  ggplot2::ggsave(file.path(figures_dir, paste0(filename, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(file.path(figures_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 600)
  ggplot2::ggsave(file.path(figures_dir, paste0(filename, ".svg")), plot, width = width, height = height, device = svglite::svglite)
}

placeholder_plot <- function(metric, reason) {
  ggplot() +
    annotate("text", x = 0, y = 0.08, label = unname(metric_labels[metric]), fontface = "bold", size = 4.2) +
    annotate("text", x = 0, y = -0.08, label = paste("Not modeled:", reason), size = 3.2, color = "gray30") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(title = paste0(unname(metric_labels[metric]), " LMM trend unavailable")) +
    theme_void(base_family = "sans") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.margin = margin(15, 15, 15, 15))
}

plot_metric <- function(metric) {
  status_row <- model_summary %>% filter(metric == !!metric)
  status <- if (nrow(status_row) == 1) status_row$status else "not_modeled"
  reason <- if (nrow(status_row) == 1) status_row$reason else "missing_model_summary"
  if (!status %in% c("modeled", "modeled_singular")) {
    return(placeholder_plot(metric, reason))
  }

  observed <- alpha %>%
    filter(!is.na(.data[[metric]]), !is.na(.data[[time_variable]])) %>%
    transmute(subjectID = subjectID, time_value = .data[[time_variable]], metric_value = .data[[metric]])
  predicted <- predictions %>%
    filter(metric == !!metric, model_type == "main")

  ggplot() +
    geom_line(
      data = observed,
      aes(x = time_value, y = metric_value, group = subjectID),
      color = "#8C8C8C",
      linewidth = 0.32,
      alpha = 0.52
    ) +
    geom_point(
      data = observed,
      aes(x = time_value, y = metric_value),
      color = "#0072B2",
      size = 1.7,
      alpha = 0.68
    ) +
    geom_ribbon(
      data = predicted,
      aes(x = time_value, ymin = confidence_interval_low, ymax = confidence_interval_high),
      fill = "#0072B2",
      alpha = 0.16
    ) +
    geom_line(
      data = predicted,
      aes(x = time_value, y = estimate),
      color = "#004B75",
      linewidth = 0.9
    ) +
    labs(
      x = x_label,
      y = unname(metric_labels[metric]),
      title = paste0(unname(metric_labels[metric]), " over expedition time"),
      subtitle = "Observed subject trajectories with population-level LMM trend",
      caption = "Thin grey lines are subject trajectories; blue points are observed samples; bold line and ribbon are fixed-effect trend and 95% CI."
    ) +
    theme_journal()
}

for (metric in metric_columns) {
  plot <- plot_metric(metric)
  save_plot(plot, filename_stems[[metric]], 6.8, 4.4)
}

available_metrics <- model_summary %>%
  filter(status %in% c("modeled", "modeled_singular")) %>%
  pull(metric)

if (length(available_metrics) > 0) {
  observed_long <- alpha %>%
    select(sampleID, subjectID, all_of(time_variable), all_of(available_metrics)) %>%
    pivot_longer(cols = all_of(available_metrics), names_to = "metric", values_to = "metric_value") %>%
    filter(!is.na(metric_value), !is.na(.data[[time_variable]])) %>%
    mutate(
      metric_label = factor(metric_labels[metric], levels = metric_labels[metric_columns]),
      time_value = .data[[time_variable]]
    )

  predicted_long <- predictions %>%
    filter(model_type == "main", metric %in% available_metrics) %>%
    mutate(metric_label = factor(metric_labels[metric], levels = metric_labels[metric_columns]))

  combined <- ggplot() +
    geom_line(
      data = observed_long,
      aes(x = time_value, y = metric_value, group = subjectID),
      color = "#8C8C8C",
      linewidth = 0.25,
      alpha = 0.46
    ) +
    geom_point(
      data = observed_long,
      aes(x = time_value, y = metric_value),
      color = "#0072B2",
      size = 1.15,
      alpha = 0.62
    ) +
    geom_ribbon(
      data = predicted_long,
      aes(x = time_value, ymin = confidence_interval_low, ymax = confidence_interval_high),
      fill = "#0072B2",
      alpha = 0.15
    ) +
    geom_line(
      data = predicted_long,
      aes(x = time_value, y = estimate),
      color = "#004B75",
      linewidth = 0.75
    ) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
    labs(
      x = x_label,
      y = "Alpha-diversity metric value",
      title = "Alpha-diversity LMM trends across expedition time",
      subtitle = "Observed sample trajectories plus fixed-effect population trends",
      caption = "Chao1 and ACE panels are omitted when not calculated from relative-abundance-like input."
    ) +
    theme_journal(base_size = 8.5)
} else {
  combined <- placeholder_plot("Shannon", "no metrics were available for LMM modeling")
}

save_plot(combined, "alpha_lmm_observed_model_trend_all_metrics", 8.0, 6.2)

cat("Alpha-diversity LMM plotting finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
