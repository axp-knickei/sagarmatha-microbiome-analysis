# Create publication-oriented alpha-diversity figures and companion report.

required_packages <- c("dplyr", "tidyr", "ggplot2", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

library(dplyr)
library(tidyr)
library(ggplot2)

output_dir <- Sys.getenv("BLACKICE_ALPHA_OUTPUT_DIR", unset = "results/alpha_diversity")
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
logs_dir <- file.path(output_dir, "logs")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("docs/reports", recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(tables_dir, "alpha_diversity_with_metadata.csv")
diagnostic_path <- file.path(tables_dir, "alpha_diversity_diagnostic_summary.csv")
if (!file.exists(input_path)) stop("Missing alpha-diversity metadata table. Run 04_compare_alpha_diversity_tools.R first.", call. = FALSE)

alpha <- utils::read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
diagnostic <- if (file.exists(diagnostic_path)) utils::read.csv(diagnostic_path, stringsAsFactors = FALSE) else data.frame(item = character(), value = character())

metrics <- c("Shannon", "Simpson", "Observed_richness", "Chao1", "ACE")
metric_labels <- c(
  Shannon = "Shannon diversity",
  Simpson = "Simpson diversity",
  Observed_richness = "Observed richness",
  Chao1 = "Chao1 richness",
  ACE = "ACE richness"
)

timepoint_levels <- if ("Timepoint" %in% names(alpha)) {
  observed_timepoints <- unique(as.character(alpha$Timepoint[!is.na(alpha$Timepoint) & nzchar(as.character(alpha$Timepoint))]))
  numeric_timepoints <- observed_timepoints[grepl("^T[0-9]+$", observed_timepoints)]
  other_timepoints <- setdiff(observed_timepoints, numeric_timepoints)
  c(numeric_timepoints[order(as.integer(sub("^T", "", numeric_timepoints)))], sort(other_timepoints))
} else {
  "All samples"
}

alpha_long <- alpha %>%
  pivot_longer(cols = all_of(metrics), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = metrics, labels = metric_labels[metrics]),
    Timepoint = if ("Timepoint" %in% names(.)) factor(Timepoint, levels = timepoint_levels, ordered = TRUE) else factor("All samples"),
    subjectID = if ("subjectID" %in% names(.)) as.character(subjectID) else sampleID
  )

group_counts <- alpha_long %>% filter(!is.na(value)) %>% count(Timepoint, metric, name = "n")
min_group_n <- if (nrow(group_counts) > 0) min(group_counts$n) else 0
violin_caution <- min_group_n < 5

theme_nature <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      legend.position = "right",
      legend.title = element_blank(),
      panel.spacing = unit(1.1, "lines"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "gray25"),
      plot.margin = margin(7, 7, 7, 7)
    )
}

palette_subject <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442",
  "#56B4E9", "#E69F00", "#000000", "#999999", "#882255",
  "#44AA99", "#AA4499", "#332288"
)

save_plot <- function(plot, filename, width, height) {
  ggsave(file.path(figures_dir, paste0(filename, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(figures_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 450)
  ggsave(file.path(figures_dir, paste0(filename, ".svg")), plot, width = width, height = height, device = svglite::svglite)
}

primary_metric <- "Shannon diversity"

boxplot_jitter <- alpha_long %>%
  filter(metric == primary_metric, !is.na(value)) %>%
  ggplot(aes(x = Timepoint, y = value)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35, fill = "white", color = "black") +
  geom_point(aes(color = subjectID), position = position_jitter(width = 0.12, height = 0), size = 1.9, alpha = 0.9) +
  scale_color_manual(values = rep(palette_subject, length.out = length(unique(alpha_long$subjectID)))) +
  labs(x = "Timepoint", y = "Shannon diversity", title = "Gut MAG alpha diversity across expedition timepoints") +
  theme_nature()
save_plot(boxplot_jitter, "alpha_diversity_boxplot_jitter_by_timepoint", 7.2, 4.6)

violin_plot <- alpha_long %>%
  filter(metric == primary_metric, !is.na(value)) %>%
  ggplot(aes(x = Timepoint, y = value)) +
  geom_violin(fill = "#DCEEF8", color = "#4D4D4D", linewidth = 0.25, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.3, fill = "white") +
  geom_point(aes(color = subjectID), position = position_jitter(width = 0.1, height = 0), size = 1.6, alpha = 0.85) +
  scale_color_manual(values = rep(palette_subject, length.out = length(unique(alpha_long$subjectID)))) +
  labs(
    x = "Timepoint",
    y = "Shannon diversity",
    title = "Distribution of gut MAG alpha diversity by timepoint",
    subtitle = if (violin_caution) "Caution: at least one timepoint has fewer than five samples." else NULL
  ) +
  theme_nature()
save_plot(violin_plot, "alpha_diversity_violin_boxplot_jitter_by_timepoint", 7.2, 4.6)

x_var <- if ("days_from_baseline" %in% names(alpha) && any(!is.na(alpha$days_from_baseline))) "days_from_baseline" else "Timepoint"
spaghetti_data <- alpha_long %>% filter(metric == primary_metric, !is.na(value))
spaghetti <- ggplot(spaghetti_data, aes(x = .data[[x_var]], y = value, group = subjectID, color = subjectID)) +
  geom_line(linewidth = 0.45, alpha = 0.75) +
  geom_point(size = 1.9, alpha = 0.95) +
  scale_color_manual(values = rep(palette_subject, length.out = length(unique(alpha_long$subjectID)))) +
  labs(
    x = if (x_var == "days_from_baseline") "Days from subject baseline" else "Timepoint",
    y = "Shannon diversity",
    title = "Longitudinal gut MAG alpha diversity trajectories"
  ) +
  theme_nature()
save_plot(spaghetti, "alpha_diversity_longitudinal_spaghetti_by_subject", 7.4, 4.8)

facet_panel <- alpha_long %>%
  ggplot(aes(x = Timepoint, y = value)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.25, fill = "white", na.rm = TRUE) +
  geom_point(aes(color = subjectID), position = position_jitter(width = 0.11, height = 0), size = 1.1, alpha = 0.85, na.rm = TRUE) +
  geom_text(
    data = data.frame(
      Timepoint = factor(timepoint_levels[ceiling(length(timepoint_levels) / 2)], levels = timepoint_levels, ordered = TRUE),
      value = 0,
      metric = factor(c("Chao1 richness", "ACE richness"), levels = metric_labels[metrics]),
      label = "Not calculated\nrelative-abundance input"
    ),
    aes(x = Timepoint, y = value, label = label),
    inherit.aes = FALSE,
    size = 3,
    color = "gray30",
    lineheight = 0.95
  ) +
  scale_color_manual(values = rep(palette_subject, length.out = length(unique(alpha_long$subjectID)))) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2, drop = FALSE) +
  labs(x = "Timepoint", y = "Alpha-diversity metric", title = "Gut MAG alpha-diversity metrics") +
  theme_nature(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(facet_panel, "alpha_diversity_faceted_multimetric_panel", 8.1, 7.2)

if ("days_from_baseline" %in% names(alpha) && sum(!is.na(alpha$days_from_baseline)) >= 5) {
  trend <- alpha_long %>%
    filter(metric == primary_metric, !is.na(value), !is.na(days_from_baseline)) %>%
    ggplot(aes(x = days_from_baseline, y = value)) +
    geom_point(aes(color = subjectID), size = 1.9, alpha = 0.9) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, color = "black", fill = "gray80", linewidth = 0.55) +
    scale_color_manual(values = rep(palette_subject, length.out = length(unique(alpha_long$subjectID)))) +
    labs(x = "Days from subject baseline", y = "Shannon diversity", title = "Estimated alpha-diversity trend over expedition time") +
    theme_nature()
  save_plot(trend, "alpha_diversity_model_trend_days_from_baseline", 7.2, 4.6)
}

heatmap <- alpha_long %>%
  filter(metric == primary_metric, !is.na(value)) %>%
  ggplot(aes(x = Timepoint, y = subjectID, fill = value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_viridis_c(option = "C", name = "Shannon") +
  labs(x = "Timepoint", y = "Subject", title = "Subject-by-timepoint alpha-diversity heatmap") +
  theme_nature()
save_plot(heatmap, "alpha_diversity_subject_timepoint_heatmap", 6.8, 4.8)

get_diag <- function(key, default = "NA") {
  hit <- diagnostic$value[diagnostic$item == key]
  if (length(hit) == 0) default else hit[[1]]
}

figure_files <- sort(basename(list.files(figures_dir, pattern = "^alpha_diversity_.*\\.(pdf|png|svg)$", full.names = TRUE)))
report_lines <- c(
  "# Alpha-Diversity Analysis Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Inputs",
  "",
  paste0("- Abundance table: `", get_diag("abundance_file"), "`"),
  paste0("- Metadata table: `", get_diag("metadata_file"), "`"),
  paste0("- Taxonomy table: `", get_diag("taxonomy_file"), "`"),
  paste0("- Detected orientation: ", get_diag("abundance_orientation")),
  paste0("- Samples: ", get_diag("prepared_samples")),
  paste0("- MAG/features: ", get_diag("prepared_MAG_features")),
  paste0("- Samples matched with metadata: ", get_diag("matched_samples")),
  paste0("- Detected abundance scale: ", get_diag("abundance_scale_decision")),
  "",
  "## Metrics and Tools",
  "",
  "- Metrics: Shannon diversity, Simpson diversity, observed richness, Chao1 richness, ACE richness.",
  "- Tools: R `vegan` and Python `scikit-bio`.",
  "- Primary plotting table: `results/alpha_diversity/tables/alpha_diversity_with_metadata.csv`.",
  "",
  "## Scientific Caution",
  "",
  "The gut MAG abundance table was detected as relative-abundance-like continuous data rather than raw counts. Chao1 and ACE are count-based richness estimators, so the primary workflow records those columns as `NA` instead of treating relative abundance as ideal count data.",
  "",
  "## Generated Tables",
  "",
  "- `results/alpha_diversity/tables/alpha_diversity_vegan.csv`",
  "- `results/alpha_diversity/tables/alpha_diversity_scikitbio.csv`",
  "- `results/alpha_diversity/tables/alpha_diversity_tool_comparison.csv`",
  "- `results/alpha_diversity/tables/alpha_diversity_with_metadata.csv`",
  "- `results/alpha_diversity/tables/alpha_diversity_diagnostic_summary.csv`",
  "",
  "## Generated Figures",
  "",
  paste0("- `results/alpha_diversity/figures/", figure_files, "`"),
  "",
  "## Figure Notes",
  "",
  "- Boxplot plus jittered points compares Shannon diversity across timepoints while showing each sample.",
  "- Violin plus boxplot plus jittered points shows the same grouping with a distribution layer; interpret cautiously where group sample size is small.",
  "- Longitudinal spaghetti plot connects samples from the same subject across expedition time.",
  "- Faceted multi-metric panel shows all requested metrics; Chao1 and ACE panels are annotated as not calculated when the input is not count-like.",
  "- Model trend plot shows a smoothed Shannon-diversity trend over days from baseline when that variable is available.",
  "- Subject-by-timepoint heatmap shows individual trajectories and missing timepoints.",
  "",
  "## Running the Workflow",
  "",
  "Local:",
  "",
  "```bash",
  "make alpha-diversity",
  "```",
  "",
  "or run each step directly:",
  "",
  "```bash",
  "Rscript scripts/alpha_diversity/01_prepare_alpha_diversity_inputs.R",
  "Rscript scripts/alpha_diversity/02_calculate_alpha_diversity_vegan.R",
  "python scripts/alpha_diversity/03_calculate_alpha_diversity_scikitbio.py",
  "Rscript scripts/alpha_diversity/04_compare_alpha_diversity_tools.R",
  "Rscript scripts/alpha_diversity/05_plot_alpha_diversity_publication.R",
  "```",
  "",
  "HPC:",
  "",
  "```bash",
  "qsub hpc/job_alpha_diversity.pbs",
  "```",
  "",
  "With an external HPC raw-data directory:",
  "",
  "```bash",
  "qsub -v BLACKICE_RAW_DIR=/path/to/raw-data hpc/job_alpha_diversity.pbs",
  "```",
  "",
  "The raw-data directory must contain the configured abundance, metadata, and taxonomy files. Defaults are `gut.csv`, `sample_info.csv`, and `GDTB.csv`."
)

writeLines(report_lines, "docs/reports/alpha_diversity_analysis_report.md")
capture.output(utils::sessionInfo(), file = file.path(logs_dir, "alpha_diversity_plot_session_info.txt"))
cat("Wrote alpha-diversity figures and report.\n")
