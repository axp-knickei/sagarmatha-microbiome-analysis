# ============================================
# summarize_hpc_candidate_signals.R
# ============================================
# Creates exploratory candidate-signal summaries from the current HPC per-MAG
# LMM results. Raw p < 0.05 MAGs are treated as candidates only, not as
# FDR-significant discoveries.

parse_args <- function(args) {
  parsed <- list(
    input = "data/processed/cleaned_data.rds",
    results_csv = "",
    output_dir = "results/hpc_per_mag_lmertest/candidates",
    candidate_p = 0.05,
    top_n = 20L
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      cat(
        "Usage: Rscript scripts/summarize_hpc_candidate_signals.R [options]\n\n",
        "Options:\n",
        "  --input PATH          Processed RDS input. Default: data/processed/cleaned_data.rds\n",
        "  --results-csv PATH    HPC result CSV. Default: combined CSV if present, else single-run CSV.\n",
        "  --output-dir PATH     Output directory. Default: results/hpc_per_mag_lmertest/candidates\n",
        "  --candidate-p VALUE   Raw p-value threshold for candidates. Default: 0.05\n",
        "  --top-n N             Number of top candidates for trajectory faceting. Default: 20\n",
        "  --help                Show this help.\n",
        sep = ""
      )
      quit(status = 0)
    } else if (arg %in% c("--input", "--results-csv", "--output-dir", "--candidate-p", "--top-n")) {
      if (i == length(args)) {
        stop("Missing value for argument: ", arg, call. = FALSE)
      }
      value <- args[[i + 1]]
      if (arg == "--input") {
        parsed$input <- value
      } else if (arg == "--results-csv") {
        parsed$results_csv <- value
      } else if (arg == "--output-dir") {
        parsed$output_dir <- sub("/+$", "", value)
      } else if (arg == "--candidate-p") {
        parsed$candidate_p <- as.numeric(value)
      } else if (arg == "--top-n") {
        parsed$top_n <- as.integer(value)
      }
      i <- i + 2
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }

  if (!is.finite(parsed$candidate_p) || parsed$candidate_p <= 0 || parsed$candidate_p >= 1) {
    stop("--candidate-p must be a number between 0 and 1.", call. = FALSE)
  }
  if (is.na(parsed$top_n) || parsed$top_n < 1) {
    stop("--top-n must be a positive integer.", call. = FALSE)
  }

  parsed
}

default_results_csv <- function() {
  candidates <- c(
    "results/hpc_per_mag_lmertest/tables/per_mag_lmertest_hpc_results_combined.csv",
    "results/hpc_per_mag_lmertest/tables/per_mag_lmertest_hpc_results.csv"
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(
      "No default HPC result CSV found. Provide --results-csv or run the HPC workflow first.",
      call. = FALSE
    )
  }
  existing[[1]]
}

safe_filename <- function(x) {
  safe <- gsub("[^A-Za-z0-9._-]+", "_", x)
  safe <- gsub("_+", "_", safe)
  safe <- gsub("^_|_$", "", safe)
  ifelse(nchar(safe) == 0, "unnamed_MAG", safe)
}

select_first_existing <- function(names_to_check, available_names) {
  found <- intersect(names_to_check, available_names)
  if (length(found) == 0) character(0) else found[[1]]
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!nzchar(args$results_csv)) {
  args$results_csv <- default_results_csv()
}

suppressPackageStartupMessages({
  library(ggplot2)
})

if (!file.exists(args$input)) {
  stop("Input file not found: ", args$input, call. = FALSE)
}
if (!file.exists(args$results_csv)) {
  stop("HPC results CSV not found: ", args$results_csv, call. = FALSE)
}

tables_dir <- file.path(args$output_dir, "tables")
figures_dir <- file.path(args$output_dir, "figures")
logs_dir <- file.path(args$output_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

results <- utils::read.csv(args$results_csv, stringsAsFactors = FALSE, check.names = FALSE)
df <- readRDS(args$input)

required_result_cols <- c("MAG", "status", "primary_analysis", "time_estimate", "time_p_value")
missing_result_cols <- setdiff(required_result_cols, names(results))
if (length(missing_result_cols) > 0) {
  stop("Result CSV is missing required columns: ", paste(missing_result_cols, collapse = ", "), call. = FALSE)
}

required_data_cols <- c("MAG", "Y", "Subject", "Timepoint")
missing_data_cols <- setdiff(required_data_cols, names(df))
if (length(missing_data_cols) > 0) {
  stop("Processed data is missing required columns: ", paste(missing_data_cols, collapse = ", "), call. = FALSE)
}

q_col <- select_first_existing(c("time_p_adj_bh_primary", "time_p_adj_bh_all"), names(results))
if (length(q_col) == 0) {
  results$time_p_adj_bh_primary <- NA_real_
  q_col <- "time_p_adj_bh_primary"
}

results$abs_time_estimate <- abs(results$time_estimate)
results$interpretation_class <- ifelse(
  results$primary_analysis %in% TRUE & results[[q_col]] < 0.05,
  "FDR_significant",
  ifelse(
    results$primary_analysis %in% TRUE & results$time_p_value < args$candidate_p,
    "nominal_candidate",
    "not_candidate"
  )
)

taxonomy_cols <- intersect(
  c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "taxonomy", "Taxonomy"),
  names(df)
)
if (length(taxonomy_cols) > 0) {
  taxonomy_map <- unique(df[, c("MAG", taxonomy_cols), drop = FALSE])
  taxonomy_map <- taxonomy_map[!duplicated(taxonomy_map$MAG), , drop = FALSE]
  results <- merge(results, taxonomy_map, by = "MAG", all.x = TRUE, sort = FALSE)
}

candidate_results <- results[
  results$primary_analysis %in% TRUE &
    !is.na(results$time_p_value) &
    results$time_p_value < args$candidate_p,
  ,
  drop = FALSE
]
candidate_results <- candidate_results[order(candidate_results$time_p_value, -candidate_results$abs_time_estimate), , drop = FALSE]

effect_ranked <- results[
  results$primary_analysis %in% TRUE & is.finite(results$time_estimate),
  ,
  drop = FALSE
]
effect_ranked <- effect_ranked[order(-effect_ranked$abs_time_estimate, effect_ranked$time_p_value), , drop = FALSE]

utils::write.csv(
  candidate_results,
  file.path(tables_dir, "nominal_candidate_mags_raw_p_lt_0_05.csv"),
  row.names = FALSE
)
utils::write.csv(
  effect_ranked,
  file.path(tables_dir, "effect_size_ranked_primary_mags.csv"),
  row.names = FALSE
)

if (length(taxonomy_cols) > 0 && nrow(candidate_results) > 0) {
  taxonomy_summary <- do.call(
    rbind,
    lapply(taxonomy_cols, function(col) {
      counts <- sort(table(candidate_results[[col]], useNA = "ifany"), decreasing = TRUE)
      data.frame(
        taxonomy_level = col,
        taxonomy_value = names(counts),
        candidate_mag_count = as.integer(counts),
        stringsAsFactors = FALSE
      )
    })
  )
  utils::write.csv(
    taxonomy_summary,
    file.path(tables_dir, "candidate_taxonomy_summary.csv"),
    row.names = FALSE
  )
}

candidate_mags <- candidate_results$MAG
plot_df <- df[df$MAG %in% candidate_mags, , drop = FALSE]
if (nrow(plot_df) > 0) {
  positive_y <- df$Y[df$Y > 0]
  pseudocount <- if (length(positive_y) > 0) min(positive_y, na.rm = TRUE) / 2 else NA_real_
  plot_df$response_log10_pseudocount <- log10(plot_df$Y + pseudocount)

  time_summary <- aggregate(
    response_log10_pseudocount ~ MAG + Timepoint,
    data = plot_df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(time_summary)[names(time_summary) == "response_log10_pseudocount"] <- "mean_response"

  heatmap_path <- file.path(figures_dir, "candidate_heatmap_mean_log10_response.png")
  heatmap_plot <- ggplot(time_summary, aes(x = Timepoint, y = MAG, fill = mean_response)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(option = "C", na.value = "grey90") +
    labs(
      title = "Exploratory candidate MAG heatmap",
      subtitle = "Mean log10(Y + pseudocount) by timepoint; candidates are raw p < 0.05 only",
      x = "Timepoint",
      y = "MAG",
      fill = "Mean response"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank()
  )
  ggsave(heatmap_path, heatmap_plot, width = 9, height = max(5, min(14, 0.18 * length(candidate_mags))), dpi = 300)
  ggsave(
    sub("\\.png$", ".svg", heatmap_path),
    heatmap_plot,
    width = 9,
    height = max(5, min(14, 0.18 * length(candidate_mags))),
    device = svglite::svglite
  )

  top_mags <- utils::head(candidate_mags, args$top_n)
  trajectory_df <- plot_df[plot_df$MAG %in% top_mags, , drop = FALSE]
  trajectory_summary <- aggregate(
    response_log10_pseudocount ~ MAG + Timepoint,
    data = trajectory_df,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  )
  expanded <- data.frame(
    MAG = trajectory_summary$MAG,
    Timepoint = trajectory_summary$Timepoint,
    mean_response = trajectory_summary$response_log10_pseudocount[, "mean"],
    se_response = trajectory_summary$response_log10_pseudocount[, "se"],
    stringsAsFactors = FALSE
  )

  trajectory_path <- file.path(figures_dir, "top_candidate_trajectories.png")
  trajectory_plot <- ggplot(expanded, aes(x = Timepoint, y = mean_response, group = MAG)) +
    geom_line(color = "#2b6cb0", linewidth = 0.5) +
    geom_point(color = "#1a365d", size = 1.4) +
    geom_errorbar(aes(ymin = mean_response - se_response, ymax = mean_response + se_response), width = 0.15, linewidth = 0.25) +
    facet_wrap(~MAG, scales = "free_y") +
    labs(
      title = "Exploratory trajectories for top nominal candidate MAGs",
      subtitle = "Facets ranked by raw p-value; error bars show standard error",
      x = "Timepoint",
      y = "Mean log10(Y + pseudocount)"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(trajectory_path, trajectory_plot, width = 11, height = max(5, ceiling(length(top_mags) / 4) * 2.2), dpi = 300)
  ggsave(
    sub("\\.png$", ".svg", trajectory_path),
    trajectory_plot,
    width = 11,
    height = max(5, ceiling(length(top_mags) / 4) * 2.2),
    device = svglite::svglite
  )

  direction_table <- candidate_results[, intersect(
    c("MAG", "time_estimate", "time_p_value", q_col, "zero_fraction", "nonzero_sample_count", "subjects_with_nonzero", taxonomy_cols),
    names(candidate_results)
  ), drop = FALSE]
  direction_table$effect_direction <- ifelse(direction_table$time_estimate > 0, "increasing_over_time", "decreasing_over_time")
  utils::write.csv(
    direction_table,
    file.path(tables_dir, "candidate_direction_summary.csv"),
    row.names = FALSE
  )
}

summary_path <- file.path(logs_dir, "candidate_signal_summary.txt")
summary_connection <- file(summary_path, open = "wt")
sink(summary_connection, split = TRUE)

cat("HPC Candidate Signal Summary\n")
cat("============================\n\n")
cat("Processed input:", args$input, "\n")
cat("HPC results CSV:", args$results_csv, "\n")
cat("Output directory:", args$output_dir, "\n")
cat("Candidate rule: primary_analysis == TRUE and raw time_p_value <", args$candidate_p, "\n\n")

cat("Counts:\n")
cat("Total MAG result rows:", nrow(results), "\n")
cat("Primary-analysis MAGs:", sum(results$primary_analysis %in% TRUE, na.rm = TRUE), "\n")
cat("FDR-significant MAGs:", sum(results$primary_analysis %in% TRUE & results[[q_col]] < 0.05, na.rm = TRUE), "\n")
cat("Nominal candidate MAGs:", nrow(candidate_results), "\n\n")

cat("Interpretation:\n")
cat("Nominal candidates are exploratory signals only. Use these outputs for heatmaps, trajectories, taxonomy summaries, effect-size ranking, and follow-up stage-based tests. Do not report raw p < 0.05 MAGs as significant discoveries unless they pass the configured FDR criterion.\n\n")

cat("Outputs:\n")
cat("- ", file.path(tables_dir, "nominal_candidate_mags_raw_p_lt_0_05.csv"), "\n", sep = "")
cat("- ", file.path(tables_dir, "effect_size_ranked_primary_mags.csv"), "\n", sep = "")
if (length(taxonomy_cols) > 0 && nrow(candidate_results) > 0) {
  cat("- ", file.path(tables_dir, "candidate_taxonomy_summary.csv"), "\n", sep = "")
}
if (nrow(plot_df) > 0) {
  cat("- ", file.path(tables_dir, "candidate_direction_summary.csv"), "\n", sep = "")
  cat("- ", file.path(figures_dir, "candidate_heatmap_mean_log10_response.png"), "\n", sep = "")
  cat("- ", file.path(figures_dir, "top_candidate_trajectories.png"), "\n", sep = "")
}

sink()
close(summary_connection)

cat("Saved candidate summaries to ", args$output_dir, "\n", sep = "")
