#!/usr/bin/env Rscript

# Purpose: Create reader-first figures and a Markdown report from the NetMoss
# altitude stability tables. This script does not refit networks or change any
# stability classification.

suppressPackageStartupMessages({
  library(ggplot2)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x)) || identical(x, "")) y else x
}

parse_args <- function(arguments) {
  parsed <- list(
    config = "config/netmoss_altitude_stability.yml",
    output_dir = NULL,
    force = FALSE
  )
  index <- 1L
  while (index <= length(arguments)) {
    argument <- arguments[[index]]
    if (argument == "--force") {
      parsed$force <- TRUE
      index <- index + 1L
    } else if (argument == "--help") {
      cat("Usage: Rscript scripts/netmoss/summarize_altitude_stability.R [--config PATH] [--output-dir PATH] [--force]\n")
      quit(status = 0)
    } else if (argument %in% c("--config", "--output-dir")) {
      if (index == length(arguments)) stop("Missing value after ", argument, call. = FALSE)
      parsed[[if (argument == "--config") "config" else "output_dir"]] <- arguments[[index + 1L]]
      index <- index + 2L
    } else stop("Unknown argument: ", argument, call. = FALSE)
  }
  parsed
}

assert_input <- function(path) {
  if (!file.exists(path)) stop("Required stability output not found: ", path, call. = FALSE)
}

assert_output <- function(path, force) {
  if (file.exists(path) && !force) stop("Refusing to overwrite ", path, ". Use --force.", call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

save_pair <- function(plot, stem, directory, width, height) {
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(directory, paste0(stem, ".svg")), plot, width = width, height = height)
}

arguments <- parse_args(commandArgs(trailingOnly = TRUE))
if (!file.exists(arguments$config)) stop("Configuration not found: ", arguments$config, call. = FALSE)
config <- yaml::read_yaml(arguments$config)
output_dir <- arguments$output_dir %||% config$outputs$output_dir
paths <- file.path(output_dir, "tables", c(
  "scenario_summary.csv", "replicate_summary.csv", "taxon_rank_stability.csv",
  "participant_influence.csv", "edge_stability.csv", "cross_rank_concordance.csv"
))
invisible(lapply(paths, assert_input))
scenario <- utils::read.csv(paths[[1]], stringsAsFactors = FALSE)
replicates <- utils::read.csv(paths[[2]], stringsAsFactors = FALSE)
taxa <- utils::read.csv(paths[[3]], stringsAsFactors = FALSE)
influence <- utils::read.csv(paths[[4]], stringsAsFactors = FALSE)
edges <- utils::read.csv(paths[[5]], stringsAsFactors = FALSE)
cross_rank <- utils::read.csv(paths[[6]], stringsAsFactors = FALSE)

figure_dir <- file.path(output_dir, "figures")
report_path <- file.path(output_dir, "report", "netmoss_altitude_stability_report.md")
figure_stems <- c("scenario_network_stability", "primary_taxon_rank_stability", "leave_one_participant_influence")
figure_paths <- unlist(lapply(figure_stems, function(stem) file.path(figure_dir, paste0(stem, c(".pdf", ".svg")))))
invisible(lapply(c(report_path, figure_paths), assert_output, force = arguments$force))

completed <- scenario[scenario$status == "completed_internal_stability", , drop = FALSE]
if (!nrow(completed)) stop("No completed stability scenarios are available to summarize.", call. = FALSE)
scenario_long <- rbind(
  data.frame(
    scenario_id = completed$scenario_id,
    taxonomic_rank = completed$taxonomic_rank,
    group = config$design$group_a_label,
    median_similarity = completed$median_bootstrap_network_similarity_a,
    stable = completed$operational_network_stability
  ),
  data.frame(
    scenario_id = completed$scenario_id,
    taxonomic_rank = completed$taxonomic_rank,
    group = config$design$group_b_label,
    median_similarity = completed$median_bootstrap_network_similarity_b,
    stable = completed$operational_network_stability
  )
)
scenario_long <- scenario_long[is.finite(scenario_long$median_similarity), , drop = FALSE]
if (nrow(scenario_long)) {
  scenario_plot <- ggplot(scenario_long, aes(median_similarity, scenario_id, color = group, shape = stable)) +
    geom_vline(xintercept = as.numeric(config$stability$network_similarity_threshold), linetype = 2, color = "grey45") +
    geom_point(size = 3) +
    facet_wrap(~taxonomic_rank, scales = "free_y") +
    scale_x_continuous(limits = c(-1, 1)) +
    labs(
      title = "Bootstrap similarity of each reference association network",
      subtitle = "Dashed line is the predeclared operational threshold; passing is internal reproducibility only",
      x = "Median Spearman correlation of network edges",
      y = NULL,
      color = "Network",
      shape = "Scenario passes"
    ) +
    theme_minimal(base_size = 10)
} else {
  scenario_plot <- ggplot() + annotate("text", x = 0, y = 0, label = "No successful bootstrap networks") +
    xlim(-1, 1) + ylim(-1, 1) + theme_void() + labs(title = "Bootstrap network stability unavailable")
}
save_pair(scenario_plot, figure_stems[[1]], figure_dir, 10, 6.5)

primary_taxa <- taxa[taxa$scenario_id %in% c("genus_matched_p80", "family_matched_p50"), , drop = FALSE]
primary_taxa <- primary_taxa[is.finite(primary_taxa$bootstrap_top10_frequency), , drop = FALSE]
primary_taxa <- primary_taxa[order(primary_taxa$taxonomic_rank, primary_taxa$reference_rank), , drop = FALSE]
if (nrow(primary_taxa)) {
  primary_taxa <- do.call(rbind, lapply(split(primary_taxa, primary_taxa$taxonomic_rank), head, 15L))
  primary_taxa$label <- factor(primary_taxa$taxon, levels = rev(unique(primary_taxa$taxon)))
  taxon_plot <- ggplot(primary_taxa, aes(bootstrap_top10_frequency, label, color = operationally_stable_taxon)) +
    geom_vline(xintercept = as.numeric(config$stability$selection_frequency_threshold), linetype = 2, color = "grey45") +
    geom_point(size = 2.6) +
    facet_wrap(~taxonomic_rank, scales = "free_y") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      title = "Top-ranked taxa are retained only when resampling is consistent",
      x = "Bootstrap frequency in NetMoss top 10",
      y = NULL,
      color = "Operationally stable"
    ) +
    theme_minimal(base_size = 10)
} else {
  taxon_plot <- ggplot() + annotate("text", x = 0, y = 0, label = "No successful primary bootstrap rankings") +
    xlim(-1, 1) + ylim(-1, 1) + theme_void() + labs(title = "Primary NetMoss rank stability unavailable")
}
save_pair(taxon_plot, figure_stems[[2]], figure_dir, 10, 7)

primary_influence <- influence[influence$scenario_id %in% c("genus_matched_p80", "family_matched_p50"), , drop = FALSE]
primary_influence <- primary_influence[is.finite(primary_influence$mean_network_similarity), , drop = FALSE]
primary_influence$participant_label <- factor(
  primary_influence$omitted_participant,
  levels = rev(sort(unique(primary_influence$omitted_participant)))
)
if (nrow(primary_influence)) {
  influence_plot <- ggplot(primary_influence, aes(mean_network_similarity, participant_label, color = scenario_id)) +
    geom_point(size = 2.5) +
    labs(
      title = "Effect of omitting each matched participant",
      subtitle = "Lower similarity identifies participants with greater network influence",
      x = "Mean similarity of the two leave-one-out networks",
      y = "Omitted participant",
      color = "Primary scenario"
    ) +
    theme_minimal(base_size = 10)
} else {
  influence_plot <- ggplot() + annotate("text", x = 0, y = 0, label = "No successful primary leave-one-out networks") +
    xlim(-1, 1) + ylim(-1, 1) + theme_void() + labs(title = "Participant influence unavailable")
}
save_pair(influence_plot, figure_stems[[3]], figure_dir, 9, 6)

primary <- scenario[scenario$scenario_id %in% c("genus_matched_p80", "family_matched_p50"), , drop = FALSE]
not_estimable <- scenario[!scenario$estimable, , drop = FALSE]
failed <- replicates[!replicates$success, , drop = FALSE]
stable_edge_count <- sum(edges$operationally_stable_edge, na.rm = TRUE)
stable_taxon_count <- sum(taxa$operationally_stable_taxon, na.rm = TRUE)
cross_rank_count <- sum(cross_rank$cross_rank_status == "stable_at_genus_and_parent_family", na.rm = TRUE)
primary_lines <- unlist(lapply(seq_len(nrow(primary)), function(index) {
  row <- primary[index, ]
  paste0(
    "- `", row$scenario_id, "`: ", if (isTRUE(row$operational_network_stability)) "passes" else "does not pass",
    "; bootstrap success ", sprintf("%.1f%%", 100 * row$bootstrap_success_fraction),
    "; median network similarities ", sprintf("%.3f", row$median_bootstrap_network_similarity_a),
    " and ", sprintf("%.3f", row$median_bootstrap_network_similarity_b), "."
  )
}))
report_lines <- c(
  "# NetMoss altitude network-stability report",
  "",
  "## Answer-first status",
  "",
  primary_lines,
  "",
  "Passing these thresholds means internal reproducibility under this dataset and resampling design. It is not independent validation and does not establish a causal altitude effect.",
  "",
  "## Complete accounting",
  "",
  paste0("- Configured scenarios: ", nrow(scenario)),
  paste0("- Scenarios below the unchanged NetMoss minimum module size: ", nrow(not_estimable)),
  paste0("- Failed resampling fits: ", nrow(failed)),
  paste0("- Operationally stable reference edges across all scenarios: ", stable_edge_count),
  paste0("- Operationally stable reference top taxa across all scenarios: ", stable_taxon_count),
  paste0("- Genus results stable together with their parent family: ", cross_rank_count),
  "",
  "Participants without a low/initial profile are excluded only from the matched primary cohort and remain in the all-participant sensitivity analysis.",
  "",
  "## Interpretation boundary",
  "",
  "The low/initial state combines baseline, initial reach, and recovery. Differences may reflect expedition stage, calendar time, diet, exertion, illness, medication, or recovery as well as altitude. Unstable results are unsuitable for biological claims; stable results remain exploratory high-altitude expedition-period associations.",
  "",
  "## Principal files",
  "",
  "- `tables/scenario_summary.csv`",
  "- `tables/edge_stability.csv`",
  "- `tables/taxon_rank_stability.csv`",
  "- `tables/participant_influence.csv`",
  "- `tables/cross_rank_concordance.csv`"
)
writeLines(report_lines, report_path)

cat("NetMoss altitude stability summary completed.\n")
cat("Report: ", normalizePath(report_path, mustWork = TRUE), "\n", sep = "")
