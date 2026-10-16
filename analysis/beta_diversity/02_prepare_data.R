#!/usr/bin/env Rscript
# Purpose: Filter validated beta-diversity inputs and construct aligned,
# canonical community, metadata, and taxonomy objects.
# Inputs: prepared_data/validated_inputs.rds plus filtering settings from the
# selected beta-diversity YAML configuration.
# Outputs: Filtering/QC tables, canonical TSV files, and
# prepared_data/prepared_data.rds under the configured results directory.
# Dependencies: Shared beta-diversity helpers and environment.yml packages.
# Usage: Rscript analysis/beta_diversity/02_prepare_data.R --config config/beta_diversity.yml
# Workflow: Stage 02; consumed by distance and downstream analysis stages.

source("analysis/beta_diversity/utils/io_helpers.R")
source("analysis/beta_diversity/utils/provenance_helpers.R")
args <- parse_beta_args(); config <- load_beta_config(args$config); ensure_beta_dirs(config)

with_stage_log(config, "02_prepare_data", {
  input <- beta_output_dir(config, "prepared_data", "validated_inputs.rds")
  if (!file.exists(input)) stop("Run 01_validate_inputs.R first.", call. = FALSE)
  data <- readRDS(input); matrix <- data$community_matrix
  min_library <- config$data$minimum_sample_library_size %||% 0
  min_prevalence <- config$data$minimum_mag_prevalence %||% 0
  min_total <- config$data$minimum_total_abundance %||% 0
  sample_reason <- ifelse(rowSums(matrix) == 0, "all_zero", ifelse(rowSums(matrix) < min_library, "below_minimum_library_size", NA))
  feature_reason <- ifelse(colSums(matrix) == 0, "all_zero", ifelse(colMeans(matrix > 0) < min_prevalence, "below_prevalence", ifelse(colSums(matrix) < min_total, "below_total_abundance", NA)))
  excluded_samples <- data.frame(sample_id = rownames(matrix)[!is.na(sample_reason)], reason = na.omit(sample_reason))
  excluded_features <- data.frame(feature_id = colnames(matrix)[!is.na(feature_reason)], reason = na.omit(feature_reason))
  keep_samples <- is.na(sample_reason); keep_features <- is.na(feature_reason)
  filtered <- matrix[keep_samples, keep_features, drop = FALSE]
  if (nrow(filtered) < 3 || ncol(filtered) < 2) stop("Filtering leaves insufficient data for beta diversity.", call. = FALSE)
  metadata <- data$metadata[match(rownames(filtered), data$metadata$sample_id), , drop = FALSE]
  taxonomy <- data$taxonomy[data$taxonomy[[config$columns$taxonomy_feature_id]] %in% colnames(filtered), , drop = FALSE]
  write_tsv(excluded_samples, beta_output_dir(config, "qc", "excluded_samples.tsv"))
  write_tsv(excluded_features, beta_output_dir(config, "qc", "excluded_features.tsv"))
  write_matrix_tsv(filtered, beta_output_dir(config, "prepared_data", "community_matrix_filtered.tsv"))
  write_tsv(metadata, beta_output_dir(config, "prepared_data", "metadata_canonical.tsv"))
  write_tsv(taxonomy, beta_output_dir(config, "prepared_data", "taxonomy_filtered.tsv"))
  sensitivity <- data.frame(
    prevalence_threshold = unlist(config$data$sensitivity_prevalence),
    retained_features = vapply(unlist(config$data$sensitivity_prevalence), function(threshold) sum(colMeans(matrix > 0) >= threshold & colSums(matrix) >= min_total), integer(1))
  )
  write_tsv(sensitivity, beta_output_dir(config, "qc", "prevalence_sensitivity_summary.tsv"))
  summary <- c("# Filtering summary", "", paste("Generated:", timestamp_beta()), "",
               sprintf("- Before: %d samples x %d MAGs.", nrow(matrix), ncol(matrix)),
               sprintf("- After: %d samples x %d MAGs.", nrow(filtered), ncol(filtered)),
               sprintf("- Minimum library size: %g.", min_library), sprintf("- Minimum prevalence: %g.", min_prevalence),
               sprintf("- Minimum total abundance: %g.", min_total),
               sprintf("- Excluded samples: %d; excluded MAGs: %d.", nrow(excluded_samples), nrow(excluded_features)),
               paste0("- Prevalence sensitivity retained features: ", paste(paste0(sensitivity$prevalence_threshold, "=", sensitivity$retained_features), collapse = ", "), "."),
               "- Exclusions are listed in the adjacent TSV files; no row was silently removed.")
  writeLines(summary, beta_output_dir(config, "qc", "filtering_summary.md"))
  output <- beta_output_dir(config, "prepared_data", "prepared_data.rds")
  assert_writable_output(output, args$force)
  saveRDS(list(community_matrix = filtered, metadata = metadata, taxonomy = taxonomy,
               exclusions = list(samples = excluded_samples, features = excluded_features)), output)
  beta_log(config, "02_prepare_data", "COMPLETED_FILTERING", sprintf("samples=%d features=%d", nrow(filtered), ncol(filtered)))
})
