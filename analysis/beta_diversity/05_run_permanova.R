#!/usr/bin/env Rscript
# Purpose: Run configured PERMANOVA, dispersion, and pairwise tests with the
# workflow's restricted-permutation rules.
# Inputs: Prepared metadata and distance matrices from stages 02 and 03.
# Outputs: Statistical result and dispersion TSV files under the configured
# results directory.
# Dependencies: Shared beta-diversity permutation/validation helpers and vegan.
# Usage: Rscript analysis/beta_diversity/05_run_permanova.R --config config/beta_diversity.yml
# Workflow: Stage 05; results are summarized by figures and the run report.

source("analysis/beta_diversity/utils/io_helpers.R")
source("analysis/beta_diversity/utils/provenance_helpers.R")
source("analysis/beta_diversity/utils/permutation_helpers.R")
args <- parse_beta_args(); config <- load_beta_config(args$config); ensure_beta_dirs(config)
require_beta_packages(c("vegan", "permute"))

with_stage_log(config, "05_run_permanova", {
  data <- readRDS(beta_output_dir(config, "prepared_data", "prepared_data.rds")); metadata <- data$metadata
  nperm <- config$analysis$permutation_count
  if (grepl("test", basename(config$.config_path), ignore.case = TRUE)) nperm <- config$analysis$smoke_test_permutations
  files <- list.files(beta_output_dir(config, "distances"), pattern = "^distance_.*\\.rds$", full.names = TRUE)
  files <- files[!grepl("_transformed\\.rds$", files)]
  all_results <- list(); dispersion_results <- list(); pairwise_results <- list()
  for (file in files) {
    metric <- sub("\\.rds$", "", sub("^distance_", "", basename(file))); distance <- readRDS(file)
    metadata <- metadata[match(attr(distance, "Labels"), metadata$sample_id), , drop = FALSE]
    control <- restricted_permutations(metadata, "subject_id", nperm)
    set.seed(config$analysis$random_seed)
    fit <- vegan::adonis2(distance ~ subject_id + timepoint, data = metadata, permutations = control, by = "margin")
    result <- tidy_adonis(fit, metric, "distance ~ subject_id + timepoint", "subject fixed blocking term; within-subject restricted permutations", nperm, nrow(metadata))
    result$hypothesis_valid <- result$hypothesis != "subject_id"
    result$interpretation_caution <- ifelse(result$hypothesis == "subject_id",
      "NOT TESTED: within-subject permutations cannot test a between-subject blocking term.",
      ifelse(result$hypothesis == "timepoint", "Primary within-subject timepoint hypothesis.", "Model bookkeeping row."))
    result$p_value[result$hypothesis == "subject_id"] <- NA_real_
    all_results[[metric]] <- result
    group <- droplevels(metadata$timepoint)
    set.seed(config$analysis$random_seed)
    bd <- vegan::betadisper(distance, group, bias.adjust = TRUE)
    bd_test <- vegan::permutest(bd, permutations = control)
    tab <- as.data.frame(bd_test$tab); tab$hypothesis <- rownames(tab)
    dispersion_results[[metric]] <- data.frame(metric = metric, hypothesis = tab$hypothesis,
      df = tab$Df, sum_squares = tab$`Sum Sq`, statistic = tab$F,
      p_value = tab$`Pr(>F)`, sample_count = nrow(metadata), permutations = nperm)
    write_tsv(data.frame(sample_id = names(bd$distances), timepoint = group,
                         distance_to_centroid = bd$distances),
              beta_output_dir(config, "statistics", paste0("dispersion_distances_", metric, ".tsv")))
    global_p <- result$p_value[result$hypothesis == "timepoint"]
    if (length(global_p) && is.finite(global_p) && global_p < 0.05) {
      pairs <- combn(levels(group), 2, simplify = FALSE)
      for (pair in pairs) {
        keep <- group %in% pair; submeta <- droplevels(metadata[keep, , drop = FALSE]); subdist <- as.dist(as.matrix(distance)[keep, keep])
        if (length(unique(submeta$subject_id)) < 2 || any(table(submeta$timepoint) < 2)) next
        subcontrol <- restricted_permutations(submeta, "subject_id", nperm)
        subfit <- vegan::adonis2(subdist ~ subject_id + timepoint, data = submeta, permutations = subcontrol, by = "margin")
        subres <- tidy_adonis(subfit, metric, paste(pair, collapse = " vs "), "subject fixed blocking term; within-subject restricted permutations", nperm, nrow(submeta))
        subres$contrast <- paste(pair, collapse = " vs "); pairwise_results[[length(pairwise_results) + 1]] <- subres[subres$hypothesis == "timepoint", ]
      }
    }
  }
  permanova <- do.call(rbind, all_results); dispersion <- do.call(rbind, dispersion_results)
  primary <- permanova$metric %in% unlist(config$analysis$primary_metrics) & permanova$hypothesis == "timepoint"
  permanova$p_adjusted_bh <- NA_real_; permanova$p_adjusted_bh[primary] <- p.adjust(permanova$p_value[primary], "BH")
  dispersion$p_adjusted_bh <- p.adjust(dispersion$p_value, "BH")
  write_tsv(permanova, beta_output_dir(config, "statistics", "permanova_results.tsv"))
  write_tsv(dispersion, beta_output_dir(config, "statistics", "dispersion_results.tsv"))
  primary_tests <- permanova[permanova$hypothesis == "timepoint", , drop = FALSE]
  dispersion_model <- dispersion[dispersion$hypothesis == "Groups", c("metric", "statistic", "p_value", "p_adjusted_bh"), drop = FALSE]
  names(dispersion_model)[-1] <- c("dispersion_statistic", "dispersion_p_value", "dispersion_p_adjusted_bh")
  complete <- merge(primary_tests, dispersion_model, by = "metric", all.x = TRUE, sort = FALSE)
  complete$interpretation_caution <- "PERMANOVA may reflect centroid location, dispersion, or both; inspect the paired dispersion result."
  write_tsv(complete, beta_output_dir(config, "statistics", "complete_statistical_results.tsv"))
  if (length(pairwise_results)) {
    pairwise <- do.call(rbind, pairwise_results); pairwise$p_adjusted_bh <- ave(pairwise$p_value, pairwise$metric, FUN = function(x) p.adjust(x, "BH"))
  } else pairwise <- data.frame(status = "SKIPPED_NO_SIGNIFICANT_GLOBAL_TEST")
  write_tsv(pairwise, beta_output_dir(config, "statistics", "pairwise_permanova_results.tsv"))
  beta_log(config, "05_run_permanova", "COMPLETED_PERMANOVA", paste("metrics=", length(all_results), "permutations=", nperm))
})
