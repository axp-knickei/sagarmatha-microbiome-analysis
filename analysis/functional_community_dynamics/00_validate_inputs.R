#!/usr/bin/env Rscript
source("analysis/functional_community_dynamics/R/common.R")

args <- parse_dynamics_args()
cfg <- read_dynamics_config(args$config)
assert_dynamics_packages()
paths <- dynamics_paths(cfg)
marker <- file.path(paths$qc, "00_validate_inputs.complete")
guard_dynamics_stage(marker, args, cfg)
inputs <- resolve_dynamics_inputs(cfg)

required_paths <- c(inputs$abundance_file, inputs$metadata_file, inputs$mag_functional_results_dir)
if (any(!file.exists(required_paths))) {
  stop("Missing required input(s): ", paste(required_paths[!file.exists(required_paths)], collapse = ", "), call. = FALSE)
}

metadata_all <- prepare_dynamics_metadata(inputs$metadata_file, cfg)
sample_col <- cfg$columns$sample_id
abundance <- fread(inputs$abundance_file, select = sample_col)
if (!sample_col %in% names(abundance)) stop("Abundance input lacks sample column: ", sample_col, call. = FALSE)
metadata <- metadata_all[sample_id %in% abundance[[sample_col]]]

layers <- configured_layers(cfg)
audit <- rbindlist(lapply(seq_len(nrow(layers)), function(i) {
  file <- file.path(inputs$mag_functional_results_dir, layers$file[[i]])
  if (!file.exists(file)) stop("Functional-potential input not found: ", file, call. = FALSE)
  preview <- read_dynamics_table(file, nrows = 1000L)
  absent <- setdiff(c("sampleID", "term", "potential", "Subject", "Timepoint"), names(preview))
  if (length(absent)) stop("Functional input is missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  data.table(
    check = paste0(layers$layer[[i]], "_input"), value = normalizePath(file), status = "pass",
    checksum_md5 = checksum_file(file)
  )
}))

expected <- as.integer(cfg$analysis$expected_sample_count)
if (isTRUE(cfg$analysis$enforce_expected_sample_count) && nrow(metadata) != expected) {
  stop("Expected ", expected, " included samples but found ", nrow(metadata), ".", call. = FALSE)
}

summary <- rbind(
  data.table(check = "included_samples", value = nrow(metadata), status = "pass", checksum_md5 = ""),
  data.table(check = "metadata_T1_T10_rows_before_abundance_intersection", value = nrow(metadata_all), status = "informational", checksum_md5 = ""),
  data.table(check = "metadata_rows_excluded_without_abundance", value = nrow(metadata_all) - nrow(metadata), status = "informational", checksum_md5 = ""),
  data.table(check = "included_subjects", value = uniqueN(metadata$subject_id), status = "pass", checksum_md5 = ""),
  data.table(check = "included_timepoints", value = paste(sort(unique(as.character(metadata$timepoint))), collapse = ","), status = "pass", checksum_md5 = ""),
  data.table(check = "excluded_timepoints", value = paste(unlist(cfg$analysis$excluded_timepoints), collapse = ","), status = "decision", checksum_md5 = ""),
  data.table(check = "abundance_input", value = normalizePath(inputs$abundance_file), status = "pass", checksum_md5 = checksum_file(inputs$abundance_file)),
  data.table(check = "metadata_input", value = normalizePath(inputs$metadata_file), status = "pass", checksum_md5 = checksum_file(inputs$metadata_file)),
  audit,
  fill = TRUE
)
write_dynamics_table(summary, file.path(paths$qc, "Table_26_preflight_input_validation.csv"))
write_dynamics_session(cfg, "00_validate_inputs")
writeLines(c(
  "Stage 00 complete", paste("Included samples:", nrow(metadata)),
  paste("Included subjects:", uniqueN(metadata$subject_id)),
  "Ta and Tb are excluded. Functional values remain predicted genomic potential."
), marker)
