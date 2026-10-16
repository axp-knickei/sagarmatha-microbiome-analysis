#!/usr/bin/env Rscript
# Purpose: Validate configured beta-diversity inputs, identifiers, and data
# alignment before filtering or calculating ecological distances.
# Inputs: Raw abundance, metadata, and taxonomy files named by the selected
# beta-diversity YAML configuration.
# Outputs: QC tables and prepared_data/validated_inputs.rds under the configured
# results directory.
# Dependencies: Shared helpers in analysis/beta_diversity/utils/ and packages
# declared in environment.yml.
# Usage: Rscript analysis/beta_diversity/01_validate_inputs.R --config config/beta_diversity.yml
# Workflow: Stage 01; consumed by 02_prepare_data.R.

source("analysis/beta_diversity/utils/io_helpers.R")
source("analysis/beta_diversity/utils/provenance_helpers.R")
source("analysis/beta_diversity/utils/validation_helpers.R")
args <- parse_beta_args(); config <- load_beta_config(args$config); ensure_beta_dirs(config)

with_stage_log(config, "01_validate_inputs", {
  abundance <- read_delimited_table(config$inputs$abundance_file)
  metadata <- read_delimited_table(config$inputs$metadata_file)
  taxonomy <- read_delimited_table(config$inputs$taxonomy_file)
  checked <- validate_community(abundance, metadata, taxonomy, config)
  metadata <- canonicalize_metadata(metadata, config)
  metadata <- metadata[match(checked$sample_ids, metadata$sample_id), , drop = FALSE]
  matrix <- checked$matrix
  sample_qc <- data.frame(sample_id = rownames(matrix), library_size = rowSums(matrix),
                          detected_features = rowSums(matrix > 0), zero_fraction = rowMeans(matrix == 0),
                          all_zero = rowSums(matrix) == 0)
  feature_qc <- data.frame(feature_id = colnames(matrix), total_abundance = colSums(matrix),
                           prevalence_n = colSums(matrix > 0), prevalence_fraction = colMeans(matrix > 0),
                           zero_fraction = colMeans(matrix == 0), all_zero = colSums(matrix) == 0)
  id_report <- data.frame(id_type = character(), identifier = character())
  if (length(checked$metadata_only)) {
    id_report <- rbind(id_report, data.frame(id_type = rep("metadata_without_abundance", length(checked$metadata_only)), identifier = checked$metadata_only))
  }
  if (length(checked$taxonomy_missing)) {
    id_report <- rbind(id_report, data.frame(id_type = rep("abundance_without_taxonomy", length(checked$taxonomy_missing)), identifier = checked$taxonomy_missing))
  }
  write_tsv(sample_qc, beta_output_dir(config, "qc", "sample_qc.tsv"))
  write_tsv(feature_qc, beta_output_dir(config, "qc", "feature_qc.tsv"))
  write_tsv(id_report, beta_output_dir(config, "qc", "id_matching_report.tsv"))
  output <- beta_output_dir(config, "prepared_data", "validated_inputs.rds")
  assert_writable_output(output, args$force)
  saveRDS(list(community_matrix = matrix, metadata = metadata, taxonomy = taxonomy), output)
  record_run_provenance(config, args)
  beta_log(config, "01_validate_inputs", "COMPLETED_VALIDATION",
           sprintf("samples=%d features=%d metadata_only=%d zeros=%.6f", nrow(matrix), ncol(matrix), length(checked$metadata_only), mean(matrix == 0)))
})
