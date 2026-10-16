#!/usr/bin/env Rscript
source("analysis/functional_community_dynamics/R/common.R")

args <- parse_dynamics_args()
cfg <- read_dynamics_config(args$config)
assert_dynamics_packages()
paths <- dynamics_paths(cfg)
marker <- file.path(paths$qc, "02_run_multivariate_analysis.complete")
guard_dynamics_stage(marker, args, cfg)
prepared_file <- file.path(paths$intermediate, "prepared_profiles.rds")
if (!file.exists(prepared_file)) stop("Stage 01 prepared profiles are missing.", call. = FALSE)
prepared <- readRDS(prepared_file)
metadata <- copy(prepared$metadata)
layer_metadata <- copy(prepared$layer_metadata)
metrics <- configured_metrics(cfg)
nperm <- as.integer(Sys.getenv("FUNCDYN_PERMUTATIONS", unset = cfg$analysis$permutation_count))
set.seed(as.integer(cfg$analysis$random_seed))

global_rows <- list()
distance_objects <- list()
pcoa_rows <- list()

for (layer in names(prepared$matrices)) {
  matrix <- prepared$matrices[[layer]]
  current_layer <- layer
  role <- layer_metadata[layer == current_layer, role]
  for (metric in metrics$metric) {
    current_metric <- metric
    key <- paste(layer, metric, sep = "__")
    status <- "fitted"
    result <- tryCatch({
      distance <- profile_distance(matrix, metric)
      if (any(!is.finite(as.numeric(distance)))) stop("Distance contains non-finite values.")
      attr(distance, "Labels") <- rownames(matrix)
      distance_objects[[key]] <- distance
      pcoa_rows[[key]] <- make_pcoa(distance, metadata, layer, metric)
      model_data <- data.frame(subject_id = factor(metadata$subject_id), timepoint = factor(metadata$timepoint))
      permutations <- restricted_permutations(metadata, nperm)
      fit <- vegan::adonis2(distance ~ subject_id + timepoint, data = model_data, permutations = permutations, by = "terms")
      timepoint_row <- tidy_adonis(fit, layer, metric, nperm, nrow(metadata))[hypothesis == "timepoint"]
      dispersion <- vegan::betadisper(distance, group = model_data$timepoint, type = "median")
      dispersion_test <- vegan::permutest(dispersion, permutations = permutations)
      dtab <- dispersion_test$tab
      data.table(
        layer = layer, role = role, metric = metric,
        metric_role = metrics[metric == current_metric, metric_role],
        samples = nrow(metadata), subjects = uniqueN(metadata$subject_id), timepoints = uniqueN(metadata$timepoint),
        permutations = nperm, permanova_df = timepoint_row$df, permanova_r_squared = timepoint_row$r_squared,
        permanova_f = timepoint_row$statistic, permanova_p_value = timepoint_row$p_value,
        dispersion_f = dtab[1L, "F"], dispersion_p_value = dtab[1L, "Pr(>F)"],
        distance_min = min(as.numeric(distance)), distance_median = median(as.numeric(distance)),
        distance_max = max(as.numeric(distance)), status = status
      )
    }, error = function(e) {
      data.table(
        layer = layer, role = role, metric = metric,
        metric_role = metrics[metric == current_metric, metric_role], samples = nrow(metadata),
        subjects = uniqueN(metadata$subject_id), timepoints = uniqueN(metadata$timepoint), permutations = nperm,
        permanova_df = NA_real_, permanova_r_squared = NA_real_, permanova_f = NA_real_, permanova_p_value = NA_real_,
        dispersion_f = NA_real_, dispersion_p_value = NA_real_, distance_min = NA_real_,
        distance_median = NA_real_, distance_max = NA_real_, status = paste0("error: ", conditionMessage(e))
      )
    })
    global_rows[[key]] <- result
  }
}

global <- rbindlist(global_rows, fill = TRUE)
global[, permanova_p_adj_bh_functional_family := NA_real_]
global[layer != "taxonomic_MAG" & role != "descriptive", permanova_p_adj_bh_functional_family := safe_bh(permanova_p_value)]
global[, dispersion_p_adj_bh := safe_bh(dispersion_p_value)]
global[, `:=`(
  supported_global = layer != "taxonomic_MAG" & role != "descriptive" & status == "fitted" & permanova_p_adj_bh_functional_family < as.numeric(cfg$analysis$fdr_threshold),
  dispersion_supported = status == "fitted" & dispersion_p_adj_bh < as.numeric(cfg$analysis$fdr_threshold)
)]
global[, inference_note := fifelse(
  dispersion_supported,
  "PERMANOVA location result requires caution because dispersion also differs.",
  "No BH-supported dispersion difference detected for this layer-metric combination."
)]
write_dynamics_table(global, file.path(paths$multivariate, "Table_27_global_PERMANOVA_and_dispersion_results.csv"))

pair_rows <- list()
functional_supported <- global[supported_global == TRUE, .(layer, metric)]
timepoints <- timepoint_levels(cfg)
pair_grid <- as.data.table(t(combn(timepoints, 2L)))
setnames(pair_grid, c("before_timepoint", "after_timepoint"))

for (i in seq_len(nrow(functional_supported))) {
  layer <- functional_supported$layer[[i]]
  metric <- functional_supported$metric[[i]]
  full_distance <- distance_objects[[paste(layer, metric, sep = "__")]]
  for (j in seq_len(nrow(pair_grid))) {
    before <- pair_grid$before_timepoint[[j]]
    after <- pair_grid$after_timepoint[[j]]
    paired_subjects <- metadata[, .(has_before = any(timepoint == before), has_after = any(timepoint == after)), by = subject_id][has_before & has_after, subject_id]
    subset_meta <- metadata[subject_id %in% paired_subjects & as.character(timepoint) %in% c(before, after)]
    if (length(paired_subjects) < 2L || nrow(subset_meta) < 4L) {
      pair_rows[[length(pair_rows) + 1L]] <- data.table(layer = layer, metric = metric, before_timepoint = before, after_timepoint = after, paired_subjects = length(paired_subjects), samples = nrow(subset_meta), r_squared = NA_real_, statistic = NA_real_, p_value = NA_real_, status = "insufficient_paired_subjects")
      next
    }
    subset_distance <- as.dist(as.matrix(full_distance)[subset_meta$sample_id, subset_meta$sample_id, drop = FALSE])
    model_data <- data.frame(subject_id = factor(subset_meta$subject_id), timepoint = factor(subset_meta$timepoint, levels = c(before, after)))
    fitted <- tryCatch({
      fit <- vegan::adonis2(subset_distance ~ subject_id + timepoint, data = model_data, permutations = restricted_permutations(subset_meta, nperm), by = "terms")
      row <- tidy_adonis(fit, layer, metric, nperm, nrow(subset_meta))[hypothesis == "timepoint"]
      data.table(layer = layer, metric = metric, before_timepoint = before, after_timepoint = after, paired_subjects = length(paired_subjects), samples = nrow(subset_meta), r_squared = row$r_squared, statistic = row$statistic, p_value = row$p_value, status = "fitted")
    }, error = function(e) data.table(layer = layer, metric = metric, before_timepoint = before, after_timepoint = after, paired_subjects = length(paired_subjects), samples = nrow(subset_meta), r_squared = NA_real_, statistic = NA_real_, p_value = NA_real_, status = paste0("error: ", conditionMessage(e))))
    pair_rows[[length(pair_rows) + 1L]] <- fitted
  }
}

pairwise <- if (length(pair_rows)) rbindlist(pair_rows, fill = TRUE) else data.table(
  layer = character(), metric = character(), before_timepoint = character(), after_timepoint = character(),
  paired_subjects = integer(), samples = integer(), r_squared = numeric(), statistic = numeric(), p_value = numeric(), status = character()
)
if (nrow(pairwise)) pairwise[, p_adj_bh_across_45_contrasts := safe_bh(p_value), by = .(layer, metric)] else pairwise[, p_adj_bh_across_45_contrasts := numeric()]
pairwise[, significant_bh := status == "fitted" & p_adj_bh_across_45_contrasts < as.numeric(cfg$analysis$fdr_threshold)]
write_dynamics_table(pairwise, file.path(paths$multivariate, "Table_28_supported_functional_timepoint_contrasts.csv"))
saveRDS(distance_objects, file.path(paths$intermediate, "distance_objects.rds"))
write_dynamics_table(rbindlist(pcoa_rows, fill = TRUE), file.path(paths$reproducibility, "Figure_16_plot_data.csv"))
write_dynamics_session(cfg, "02_run_multivariate_analysis")
writeLines(c(
  "Stage 02 complete", paste("Global layer-metric tests:", nrow(global)),
  paste("BH-supported functional global tests:", sum(global$supported_global)),
  paste("Pairwise tests fitted:", sum(pairwise$status == "fitted")),
  "PERMANOVA uses participant-restricted permutations; dispersion is reported alongside every global test."
), marker)
