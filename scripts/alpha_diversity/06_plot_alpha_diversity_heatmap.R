# Create publication-style subject x timepoint alpha-diversity heatmaps.

required_packages <- c("dplyr", "tidyr", "ggplot2", "svglite", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

repo_root <- getwd()
output_dir <- Sys.getenv("BLACKICE_ALPHA_OUTPUT_DIR", unset = "results/alpha_diversity")
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
heatmap_tables_dir <- file.path(tables_dir, "heatmaps")
heatmap_figures_dir <- file.path(figures_dir, "heatmaps")
logs_dir <- file.path(output_dir, "logs")
report_dir <- file.path("docs", "reports")

dir.create(heatmap_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(heatmap_figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(logs_dir, "alpha_diversity_heatmap.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Alpha-diversity heatmap workflow started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Working directory:", repo_root, "\n")

alpha_path <- Sys.getenv(
  "BLACKICE_ALPHA_HEATMAP_INPUT",
  unset = file.path(tables_dir, "alpha_diversity_with_metadata.csv")
)
metadata_path <- Sys.getenv("BLACKICE_METADATA_PATH", unset = file.path("data", "raw", "sample_info.csv"))
diagnostic_path <- file.path(tables_dir, "alpha_diversity_diagnostic_summary.csv")

if (!file.exists(alpha_path)) {
  fallback_paths <- c(
    file.path(tables_dir, "alpha_diversity_vegan.csv"),
    file.path(tables_dir, "alpha_diversity_scikitbio.csv")
  )
  fallback_paths <- fallback_paths[file.exists(fallback_paths)]
  if (length(fallback_paths) == 0) {
    stop("No alpha-diversity input table found. Expected: ", alpha_path, call. = FALSE)
  }
  alpha_path <- fallback_paths[[1]]
}

alpha <- utils::read.csv(alpha_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata_used <- FALSE
metadata <- NULL

required_core <- c("sampleID")
missing_core <- setdiff(required_core, names(alpha))
if (length(missing_core) > 0) {
  stop("Alpha-diversity table is missing required column(s): ", paste(missing_core, collapse = ", "), call. = FALSE)
}

if (!all(c("subjectID", "Timepoint") %in% names(alpha))) {
  if (!file.exists(metadata_path)) {
    stop(
      "Alpha-diversity table lacks subjectID/Timepoint and metadata file was not found: ",
      metadata_path,
      call. = FALSE
    )
  }
  metadata <- utils::read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"sampleID" %in% names(metadata)) {
    stop("Metadata table must contain sampleID for joining.", call. = FALSE)
  }
  if (!"subjectID" %in% names(metadata) && "group" %in% names(metadata)) {
    metadata$subjectID <- metadata$group
  }
  required_metadata <- c("sampleID", "subjectID", "Timepoint")
  missing_metadata <- setdiff(required_metadata, names(metadata))
  if (length(missing_metadata) > 0) {
    stop("Metadata table is missing required column(s): ", paste(missing_metadata, collapse = ", "), call. = FALSE)
  }
  alpha <- alpha %>%
    left_join(metadata, by = "sampleID", suffix = c("", "_metadata"))
  metadata_used <- TRUE
} else if (file.exists(metadata_path)) {
  metadata <- utils::read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
}

required_columns <- c("sampleID", "subjectID", "Timepoint")
missing_columns <- setdiff(required_columns, names(alpha))
if (length(missing_columns) > 0) {
  stop("Alpha-diversity plotting table is missing required column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

metric_columns <- c("Shannon", "Simpson", "Observed_richness", "Chao1", "ACE")
metric_labels <- c(
  Shannon = "Shannon diversity",
  Simpson = "Simpson diversity",
  Observed_richness = "Observed richness",
  Chao1 = "Chao1 richness",
  ACE = "ACE richness"
)
metric_file_labels <- c(
  Shannon = "shannon",
  Simpson = "simpson",
  Observed_richness = "observed_richness",
  Chao1 = "chao1",
  ACE = "ace"
)
missing_metric_columns <- setdiff(metric_columns, names(alpha))
if (length(missing_metric_columns) == length(metric_columns)) {
  stop("No expected alpha-diversity metric columns were found.", call. = FALSE)
}

for (metric in intersect(metric_columns, names(alpha))) {
  alpha[[metric]] <- suppressWarnings(as.numeric(alpha[[metric]]))
}

alpha <- alpha %>%
  mutate(
    sampleID = as.character(sampleID),
    subjectID = as.character(subjectID),
    Timepoint = as.character(Timepoint),
    days_from_baseline = if ("days_from_baseline" %in% names(.)) suppressWarnings(as.numeric(days_from_baseline)) else NA_real_,
    date_for_order = if ("date_converted" %in% names(.)) as.Date(date_converted) else as.Date(NA)
  )

if (!"date_converted" %in% names(alpha) && "date" %in% names(alpha)) {
  parsed_date <- suppressWarnings(as.Date(as.character(alpha$date), format = "%Y%m%d"))
  alpha$date_for_order <- parsed_date
}

if (any(is.na(alpha$sampleID) | !nzchar(alpha$sampleID))) {
  stop("Alpha-diversity table contains missing sampleID values.", call. = FALSE)
}
if (anyDuplicated(alpha$sampleID) > 0) {
  stop("Alpha-diversity table contains duplicated sampleID values.", call. = FALSE)
}
if (any(is.na(alpha$subjectID) | !nzchar(alpha$subjectID))) {
  stop("Alpha-diversity table contains missing subjectID values.", call. = FALSE)
}
if (any(is.na(alpha$Timepoint) | !nzchar(alpha$Timepoint))) {
  stop("Alpha-diversity table contains missing Timepoint values.", call. = FALSE)
}

natural_subject_order <- function(x) {
  x <- unique(as.character(x))
  prefix <- sub("[0-9]+$", "", x)
  suffix <- suppressWarnings(as.integer(sub("^.*?([0-9]+)$", "\\1", x)))
  x[order(prefix, suffix, x, na.last = TRUE)]
}

infer_timepoint_order <- function(df) {
  summary <- df %>%
    group_by(Timepoint) %>%
    summarise(
      median_days = suppressWarnings(stats::median(days_from_baseline, na.rm = TRUE)),
      median_date = suppressWarnings(stats::median(as.numeric(date_for_order), na.rm = TRUE)),
      numeric_suffix = ifelse(grepl("^T[0-9]+$", first(Timepoint)), as.integer(sub("^T", "", first(Timepoint))), NA_integer_),
      .groups = "drop"
    ) %>%
    mutate(
      median_days = ifelse(is.infinite(median_days) | is.nan(median_days), NA_real_, median_days),
      median_date = ifelse(is.infinite(median_date) | is.nan(median_date), NA_real_, median_date)
    )
  if (all(is.na(summary$median_days)) && all(is.na(summary$median_date)) && all(is.na(summary$numeric_suffix))) {
    stop("Could not infer biological timepoint order from Timepoint, days_from_baseline, or date.", call. = FALSE)
  }
  if (!all(is.na(summary$median_date))) {
    return(summary %>%
      arrange(ifelse(is.na(median_date), Inf, median_date), ifelse(is.na(numeric_suffix), Inf, numeric_suffix), Timepoint) %>%
      pull(Timepoint))
  }
  if (!all(is.na(summary$median_days))) {
    return(summary %>%
      arrange(ifelse(is.na(median_days), Inf, median_days), ifelse(is.na(numeric_suffix), Inf, numeric_suffix), Timepoint) %>%
      pull(Timepoint))
  }
  summary %>%
    arrange(ifelse(is.na(numeric_suffix), Inf, numeric_suffix), Timepoint) %>%
    pull(Timepoint)
}

subject_levels <- natural_subject_order(alpha$subjectID)
timepoint_levels <- infer_timepoint_order(alpha)

metric_availability <- lapply(metric_columns, function(metric) {
  if (!metric %in% names(alpha)) {
    return(data.frame(metric = metric, metric_label = metric_labels[[metric]], available = FALSE, reason = "missing_column", non_missing_values = 0L))
  }
  n_non_missing <- sum(!is.na(alpha[[metric]]))
  data.frame(
    metric = metric,
    metric_label = metric_labels[[metric]],
    available = n_non_missing > 0,
    reason = if (n_non_missing > 0) "available" else "all_values_missing",
    non_missing_values = n_non_missing,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()
utils::write.csv(metric_availability, file.path(heatmap_tables_dir, "alpha_diversity_heatmap_metric_availability.csv"), row.names = FALSE)

available_metrics <- metric_availability$metric[metric_availability$available]
if (length(available_metrics) == 0) {
  stop("No alpha-diversity metrics contain non-missing values for heatmap plotting.", call. = FALSE)
}

grid <- expand.grid(
  subjectID = subject_levels,
  Timepoint = timepoint_levels,
  metric = available_metrics,
  stringsAsFactors = FALSE
)

alpha_long <- alpha %>%
  select(sampleID, subjectID, Timepoint, any_of(c("days_from_baseline", "date", "date_converted", "altitude", "age", "gender")), all_of(available_metrics)) %>%
  pivot_longer(cols = all_of(available_metrics), names_to = "metric", values_to = "value") %>%
  mutate(metric_label = unname(metric_labels[metric]))

plot_long <- grid %>%
  left_join(alpha_long, by = c("subjectID", "Timepoint", "metric")) %>%
  mutate(
    metric = factor(metric, levels = available_metrics),
    metric_label = unname(metric_labels[as.character(metric)]),
    subjectID = factor(subjectID, levels = rev(subject_levels)),
    Timepoint = factor(Timepoint, levels = timepoint_levels, ordered = TRUE),
    missing_status = ifelse(is.na(sampleID) | is.na(value), "missing", "observed")
  ) %>%
  group_by(subjectID, metric) %>%
  mutate(
    subject_non_missing_values = sum(!is.na(value)),
    subject_sd = stats::sd(value, na.rm = TRUE),
    value_zscore = ifelse(subject_non_missing_values >= 2 & is.finite(subject_sd) & subject_sd > 0, (value - mean(value, na.rm = TRUE)) / subject_sd, NA_real_),
    zscore_status = case_when(
      missing_status == "missing" ~ "missing",
      subject_non_missing_values < 2 ~ "too_few_values",
      !is.finite(subject_sd) | subject_sd <= 0 ~ "zero_or_invalid_sd",
      TRUE ~ "calculated"
    )
  ) %>%
  ungroup()

utils::write.csv(plot_long, file.path(heatmap_tables_dir, "alpha_diversity_heatmap_long_table.csv"), row.names = FALSE)

matrix_paths <- character()
for (metric in available_metrics) {
  metric_name <- metric_file_labels[[metric]]
  metric_data <- plot_long %>% filter(as.character(metric) == !!metric)

  raw_matrix <- metric_data %>%
    select(subjectID, Timepoint, value) %>%
    mutate(subjectID = as.character(subjectID), Timepoint = as.character(Timepoint)) %>%
    pivot_wider(names_from = Timepoint, values_from = value) %>%
    arrange(factor(subjectID, levels = rev(subject_levels)))
  z_matrix <- metric_data %>%
    select(subjectID, Timepoint, value_zscore) %>%
    mutate(subjectID = as.character(subjectID), Timepoint = as.character(Timepoint)) %>%
    pivot_wider(names_from = Timepoint, values_from = value_zscore) %>%
    arrange(factor(subjectID, levels = rev(subject_levels)))

  raw_path <- file.path(heatmap_tables_dir, paste0("alpha_diversity_heatmap_matrix_", metric_name, "_raw.csv"))
  z_path <- file.path(heatmap_tables_dir, paste0("alpha_diversity_heatmap_matrix_", metric_name, "_zscore.csv"))
  utils::write.csv(raw_matrix, raw_path, row.names = FALSE)
  utils::write.csv(z_matrix, z_path, row.names = FALSE)
  matrix_paths <- c(matrix_paths, raw_path, z_path)
}

zscore_diagnostics <- plot_long %>%
  group_by(metric, metric_label, subjectID) %>%
  summarise(
    non_missing_values = first(subject_non_missing_values),
    subject_sd = first(subject_sd),
    zscore_status = paste(sort(unique(zscore_status)), collapse = "; "),
    .groups = "drop"
  )
utils::write.csv(zscore_diagnostics, file.path(heatmap_tables_dir, "alpha_diversity_heatmap_zscore_diagnostics.csv"), row.names = FALSE)

theme_heatmap <- function(base_size = 8) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.text = element_text(face = "bold", color = "black"),
      legend.position = "right",
      legend.title = element_text(size = rel(0.9)),
      legend.text = element_text(size = rel(0.82)),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "gray25", hjust = 0),
      plot.caption = element_text(color = "gray35", hjust = 0, size = rel(0.82)),
      plot.margin = margin(8, 8, 8, 8)
    )
}

save_plot <- function(plot, stem, width, height) {
  ggplot2::ggsave(file.path(heatmap_figures_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(file.path(heatmap_figures_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 600)
  ggplot2::ggsave(file.path(heatmap_figures_dir, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite::svglite)
}

for (metric in available_metrics) {
  metric_name <- metric_file_labels[[metric]]
  metric_label <- metric_labels[[metric]]
  metric_data <- plot_long %>% filter(as.character(metric) == !!metric)

  raw_plot <- ggplot(metric_data, aes(x = Timepoint, y = subjectID, fill = value)) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_viridis_c(
      option = "D",
      direction = 1,
      na.value = "#D9D9D9",
      name = metric_label
    ) +
    labs(
      x = "Timepoint",
      y = "Subject",
      title = paste0(metric_label, " across subjects and timepoints"),
      subtitle = "Raw alpha-diversity values; grey tiles mark missing subject-timepoint combinations.",
      caption = "Rows are subjects ordered naturally; columns are ordered by expedition timing."
    ) +
    theme_heatmap()
  save_plot(raw_plot, paste0("alpha_diversity_heatmap_", metric_name, "_raw"), 7.4, 4.8)

  z_plot <- ggplot(metric_data, aes(x = Timepoint, y = subjectID, fill = value_zscore)) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      na.value = "#D9D9D9",
      name = "Within-subject z"
    ) +
    labs(
      x = "Timepoint",
      y = "Subject",
      title = paste0(metric_label, " within-subject z-score heatmap"),
      subtitle = "Blue indicates lower-than-subject-average diversity; red indicates higher-than-subject-average diversity.",
      caption = "Grey tiles mark missing values or z-scores that could not be calculated safely."
    ) +
    theme_heatmap()
  save_plot(z_plot, paste0("alpha_diversity_heatmap_", metric_name, "_zscore"), 7.4, 4.8)
}

combined_data <- plot_long %>%
  mutate(metric_label = factor(metric_label, levels = unname(metric_labels[available_metrics])))

combined_raw <- ggplot(combined_data, aes(x = Timepoint, y = subjectID, fill = value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_viridis_c(option = "D", direction = 1, na.value = "#D9D9D9", name = "Raw value") +
  facet_wrap(~ metric_label, ncol = 1, scales = "free_x") +
  labs(
    x = "Timepoint",
    y = "Subject",
    title = "Alpha-diversity heatmap matrix across metrics",
    subtitle = "Raw values; interpret color intensity within each metric, not across metrics with different units.",
    caption = "Grey tiles mark missing subject-timepoint combinations or missing metric values."
  ) +
  theme_heatmap(base_size = 7.5)
save_plot(combined_raw, "alpha_diversity_heatmap_multimetric_raw", 7.8, 9.2)

combined_z <- ggplot(combined_data, aes(x = Timepoint, y = subjectID, fill = value_zscore)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "#D9D9D9", name = "Within-subject z") +
  facet_wrap(~ metric_label, ncol = 1, scales = "free_x") +
  labs(
    x = "Timepoint",
    y = "Subject",
    title = "Within-subject alpha-diversity z-score heatmap matrix",
    subtitle = "Blue indicates lower-than-subject-average diversity; red indicates higher-than-subject-average diversity.",
    caption = "Grey tiles mark missing values or z-scores that could not be calculated safely."
  ) +
  theme_heatmap(base_size = 7.5)
save_plot(combined_z, "alpha_diversity_heatmap_multimetric_zscore", 7.8, 9.2)

diagnostic <- if (file.exists(diagnostic_path)) {
  utils::read.csv(diagnostic_path, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  data.frame(item = character(), value = character(), stringsAsFactors = FALSE)
}
diag_value <- function(key, default = "NA") {
  if (!all(c("item", "value") %in% names(diagnostic))) return(default)
  hit <- diagnostic$value[diagnostic$item == key]
  if (length(hit) == 0) default else as.character(hit[[1]])
}

n_expected <- length(subject_levels) * length(timepoint_levels)
observed_grid <- plot_long %>%
  distinct(subjectID, Timepoint, sampleID) %>%
  filter(!is.na(sampleID))
n_observed <- nrow(observed_grid)
n_missing <- n_expected - n_observed

inspection_summary <- data.frame(
  item = c(
    "alpha_diversity_table",
    "metadata_file",
    "metadata_joined_in_script",
    "metadata_columns_detected",
    "subjects",
    "timepoints",
    "timepoint_order",
    "observed_subject_timepoint_samples",
    "missing_subject_timepoint_combinations",
    "metrics_plotted",
    "metrics_not_plotted",
    "abundance_scale_decision"
  ),
  value = c(
    alpha_path,
    if (file.exists(metadata_path)) metadata_path else "not_used_missing",
    as.character(metadata_used),
    paste(intersect(c("sampleID", "subjectID", "group", "Timepoint", "days_from_baseline", "date", "date_converted", "altitude", "gender", "age"), names(alpha)), collapse = "; "),
    as.character(length(subject_levels)),
    as.character(length(timepoint_levels)),
    paste(timepoint_levels, collapse = "; "),
    as.character(n_observed),
    as.character(n_missing),
    paste(unname(metric_labels[available_metrics]), collapse = "; "),
    paste(unname(metric_labels[setdiff(metric_columns, available_metrics)]), collapse = "; "),
    diag_value("abundance_scale_decision", default = unique(alpha$abundance_scale_decision)[1])
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(inspection_summary, file.path(heatmap_tables_dir, "alpha_diversity_heatmap_diagnostic_summary.csv"), row.names = FALSE)

figure_files <- sort(basename(list.files(heatmap_figures_dir, pattern = "\\.(pdf|png|svg)$")))
matrix_files <- sort(basename(list.files(heatmap_tables_dir, pattern = "^alpha_diversity_heatmap_matrix_.*\\.csv$")))

report_lines <- c(
  "# Alpha-Diversity Heatmap Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Inputs",
  "",
  paste0("- Alpha-diversity table: `", alpha_path, "`"),
  paste0("- Metadata file: `", if (file.exists(metadata_path)) metadata_path else "not used", "`"),
  paste0("- Metadata joined in heatmap script: ", metadata_used),
  paste0("- Detected metadata columns: ", inspection_summary$value[inspection_summary$item == "metadata_columns_detected"]),
  "",
  "## Matrix Dimensions",
  "",
  paste0("- Subjects: ", length(subject_levels)),
  paste0("- Timepoints: ", length(timepoint_levels), " (", paste(timepoint_levels, collapse = ", "), ")"),
  paste0("- Observed subject-timepoint samples: ", n_observed),
  paste0("- Missing subject-timepoint combinations: ", n_missing),
  "",
  "## Metrics",
  "",
  paste0("- Plotted: ", paste(unname(metric_labels[available_metrics]), collapse = "; ")),
  paste0("- Not plotted: ", paste(unname(metric_labels[setdiff(metric_columns, available_metrics)]), collapse = "; ")),
  "",
  "## Generated Matrix Tables",
  "",
  paste0("- `results/alpha_diversity/tables/heatmaps/", matrix_files, "`"),
  "",
  "## Generated Figures",
  "",
  paste0("- `results/alpha_diversity/figures/heatmaps/", figure_files, "`"),
  "",
  "## Interpretation",
  "",
  "- Raw heatmaps show the original alpha-diversity metric values for each observed subject-timepoint sample.",
  "- Row-scaled z-score heatmaps show whether each sample is higher or lower than that subject's own longitudinal mean for that metric.",
  "- Grey tiles represent missing subject-timepoint combinations, missing metric values, or z-scores that could not be calculated safely.",
  "",
  "## Cautions",
  "",
  "- These heatmaps are descriptive visualizations, not statistical tests.",
  "- Formal inference should use longitudinal models with subject-level random effects.",
  "- Chao1 and ACE should be interpreted cautiously when the abundance table is relative-abundance-like rather than raw count-like.",
  "- Extra labels `Ta` and `Tb` are retained and ordered by observed expedition timing because they are present in the dataset.",
  "",
  "## Local Run",
  "",
  "```bash",
  "Rscript scripts/alpha_diversity/06_plot_alpha_diversity_heatmap.R",
  "```",
  "",
  "## HPC Run",
  "",
  "```bash",
  "qsub hpc/job_alpha_diversity_heatmap.pbs",
  "```"
)
writeLines(report_lines, file.path(report_dir, "alpha_diversity_heatmap_report.md"))
utils::capture.output(utils::sessionInfo(), file = file.path(logs_dir, "alpha_diversity_heatmap_session_info.txt"))

cat("Metrics plotted:", paste(available_metrics, collapse = ", "), "\n")
cat("Heatmap figures:", heatmap_figures_dir, "\n")
cat("Heatmap tables:", heatmap_tables_dir, "\n")
cat("Alpha-diversity heatmap workflow finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
