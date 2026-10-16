#!/usr/bin/env Rscript

# Purpose: Summarize and visualize the complete matched-timepoint companion
# results without selecting or hiding timepoint comparisons.
# Inputs: Tables written by run_matched_timepoint_companion.R.
# Outputs: coverage/evidence figures and a concise generated Markdown report.
# Dependencies: ggplot2, yaml.
# Invocation: run from the repository root after the companion analysis.
# Relationship: reporting only; it does not rerun NetMoss2 or paired inference.

suppressPackageStartupMessages({
  library(ggplot2)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x)) || identical(x, "")) y else x
}

usage <- function() {
  cat(
    "Usage: Rscript scripts/netmoss/summarize_matched_timepoint_companion.R [options]\n",
    "  --config PATH       Default: config/netmoss_matched_timepoint_companion.yml\n",
    "  --output-dir PATH   Override configured output directory\n",
    "  --force             Replace existing named report outputs\n",
    "  --help\n",
    sep = ""
  )
}

parse_args <- function(arguments) {
  parsed <- list(
    config = "config/netmoss_matched_timepoint_companion.yml",
    output_dir = NULL,
    force = FALSE
  )
  index <- 1L
  while (index <= length(arguments)) {
    argument <- arguments[[index]]
    if (argument == "--help") {
      usage()
      quit(status = 0)
    } else if (argument == "--force") {
      parsed$force <- TRUE
      index <- index + 1L
    } else if (argument %in% c("--config", "--output-dir")) {
      if (index == length(arguments)) stop("Missing value after ", argument, call. = FALSE)
      key <- if (argument == "--config") "config" else "output_dir"
      parsed[[key]] <- arguments[[index + 1L]]
      index <- index + 2L
    } else {
      stop("Unknown argument: ", argument, call. = FALSE)
    }
  }
  parsed
}

assert_input <- function(path) {
  if (!file.exists(path)) stop("Required companion output not found: ", path, call. = FALSE)
}

assert_output <- function(path, force) {
  if (file.exists(path) && !force) {
    stop("Refusing to overwrite ", path, ". Use --force.", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

arguments <- parse_args(commandArgs(trailingOnly = TRUE))
if (!file.exists(arguments$config)) stop("Configuration not found: ", arguments$config, call. = FALSE)
config <- yaml::read_yaml(arguments$config)
output_dir <- arguments$output_dir %||% config$outputs$output_dir

pair_path <- file.path(output_dir, "tables", "timepoint_pair_support.csv")
genus_path <- file.path(output_dir, "tables", "all_genus_pair_results.csv")
participant_path <- file.path(output_dir, "tables", "all_participant_genus_changes.csv")
manifest_path <- file.path(output_dir, "reproducibility", "analysis_manifest.yml")
invisible(lapply(c(pair_path, genus_path, participant_path, manifest_path), assert_input))

pair_support <- utils::read.csv(pair_path, stringsAsFactors = FALSE)
genus_results <- utils::read.csv(genus_path, stringsAsFactors = FALSE)
participant_results <- utils::read.csv(participant_path, stringsAsFactors = FALSE)
manifest <- yaml::read_yaml(manifest_path)
timepoint_order <- unlist(manifest$timepoints %||% config$design$timepoint_order)

figure_dir <- file.path(output_dir, "figures")
report_path <- file.path(output_dir, "report", "matched_timepoint_companion_report.md")
figure_paths <- file.path(figure_dir, c(
  "matched_participant_support.pdf",
  "matched_participant_support.svg",
  "top_ranked_evidence_map.pdf",
  "top_ranked_evidence_map.svg",
  "top_ranked_participant_changes.pdf",
  "top_ranked_participant_changes.svg"
))
invisible(lapply(c(report_path, figure_paths), assert_output, force = arguments$force))

pair_support$timepoint_a <- factor(pair_support$timepoint_a, levels = timepoint_order)
pair_support$timepoint_b <- factor(pair_support$timepoint_b, levels = rev(timepoint_order))
coverage_plot <- ggplot(pair_support, aes(timepoint_a, timepoint_b)) +
  geom_tile(aes(fill = matched_participants), color = "white", linewidth = 0.4) +
  geom_text(aes(label = matched_participants), size = 3.2) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  labs(
    title = "Matched participants available for every timepoint pair",
    x = "Earlier configured timepoint",
    y = "Later configured timepoint",
    fill = "Participants"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(figure_paths[[1]], coverage_plot, width = 8.5, height = 6.5)
ggsave(figure_paths[[2]], coverage_plot, width = 8.5, height = 6.5)

rankable <- genus_results[is.finite(genus_results$global_bh_q_value), , drop = FALSE]
rankable <- rankable[order(
  rankable$global_bh_q_value,
  -abs(rankable$mean_clr_change),
  rankable$pair_index,
  rankable$genus
), , drop = FALSE]
top_ranked <- head(rankable, 30L)
if (!nrow(top_ranked)) stop("No inferential genus-pair results are available to summarize.", call. = FALSE)
top_ranked$evidence_label <- paste(top_ranked$genus, top_ranked$pair_id, sep = " | ")
top_ranked$evidence_label <- factor(top_ranked$evidence_label, levels = rev(top_ranked$evidence_label))
top_ranked$minus_log10_q <- -log10(pmax(top_ranked$global_bh_q_value, .Machine$double.xmin))

evidence_plot <- ggplot(top_ranked, aes(mean_clr_change, evidence_label)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.5) +
  geom_point(aes(size = minus_log10_q, color = direction_concordance), alpha = 0.85) +
  scale_color_viridis_c(option = "D") +
  labs(
    title = "Top 30 fully retained genus-timepoint results",
    subtitle = "Ranking is diagnostic; candidate status still requires every predeclared safeguard",
    x = "Mean paired CLR change",
    y = NULL,
    size = "-log10(global q)",
    color = "Direction\nconcordance"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_line(color = "grey92"))

ggsave(figure_paths[[3]], evidence_plot, width = 11, height = 9)
ggsave(figure_paths[[4]], evidence_plot, width = 11, height = 9)

top_result <- top_ranked[1, , drop = FALSE]
top_participants <- participant_results[
  participant_results$pair_id == top_result$pair_id &
    participant_results$genus == top_result$genus,
  ,
  drop = FALSE
]
top_long <- rbind(
  data.frame(
    subject_id = top_participants$subject_id,
    timepoint = top_participants$timepoint_a,
    clr = top_participants$clr_a
  ),
  data.frame(
    subject_id = top_participants$subject_id,
    timepoint = top_participants$timepoint_b,
    clr = top_participants$clr_b
  )
)
top_long$timepoint <- factor(
  top_long$timepoint,
  levels = c(top_result$timepoint_a, top_result$timepoint_b)
)
participant_plot <- ggplot(top_long, aes(timepoint, clr, group = subject_id)) +
  geom_line(alpha = 0.55, color = "grey45") +
  geom_point(aes(color = subject_id), size = 2.2) +
  labs(
    title = paste0("Highest-ranked diagnostic result: ", top_result$genus),
    subtitle = paste0(top_result$pair_id, "; status = ", top_result$interpretation_status),
    x = NULL,
    y = "CLR abundance",
    color = "Participant"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

ggsave(figure_paths[[5]], participant_plot, width = 8.5, height = 5.5)
ggsave(figure_paths[[6]], participant_plot, width = 8.5, height = 5.5)

candidate_count <- sum(genus_results$predeclared_candidate, na.rm = TRUE)
descriptive_pairs <- sum(!pair_support$inferential_pair)
report_lines <- c(
  "# Matched-timepoint NetMoss companion report",
  "",
  "## Status",
  "",
  "This is an abundance-based paired companion analysis. It is not NetMoss, does not estimate participant-level networks, and does not validate microbial biomarkers.",
  "",
  "## Complete analysis accounting",
  "",
  paste0("- Timepoints: ", length(timepoint_order)),
  paste0("- Enumerated timepoint pairs: ", nrow(pair_support)),
  paste0("- Pairs below the inferential participant minimum: ", descriptive_pairs),
  paste0("- Tested genus–pair hypotheses: ", manifest$tested_genus_pair_hypotheses),
  paste0("- Predeclared candidates passing every safeguard: ", candidate_count),
  paste0("- Timepoint pairs meeting the configured NetMoss profile minimum: ",
         sum(pair_support$netmoss_network_feasible)),
  "",
  "All comparisons remain in the output tables. High scores were not used to choose which comparisons to retain.",
  "",
  "## Interpretation boundary",
  "",
  "The top-ranked diagnostic figure is included to inspect behavior, not to declare a discovery. Biological interpretation requires the global BH result, bootstrap interval, direction concordance, leave-one-participant-out sign stability, adequate detection, and independent validation.",
  "",
  "## Figures",
  "",
  "- `figures/matched_participant_support.pdf`",
  "- `figures/top_ranked_evidence_map.pdf`",
  "- `figures/top_ranked_participant_changes.pdf`"
)
writeLines(report_lines, report_path)

cat("Matched-timepoint companion summary completed.\n")
cat("Report: ", normalizePath(report_path, mustWork = TRUE), "\n", sep = "")
