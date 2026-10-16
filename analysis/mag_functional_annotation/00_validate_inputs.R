#!/usr/bin/env Rscript

# Purpose: Validate configuration, file structure, canonical joins, timepoint
# decisions, and the unnormalized gut abundance scale before full catalogue work.
# Inputs: Full MMseqs membership and eggNOG annotation TSVs plus gut.csv,
# sample_info.csv, GDTB.csv, and config/mag_functional_annotation.yml.
# Outputs: Numbered QC tables and a completion marker under the configured output.
# Dependencies: R, data.table, yaml, and analysis/mag_functional_annotation/R/common.R.
# Usage: Rscript analysis/mag_functional_annotation/00_validate_inputs.R --config config/mag_functional_annotation.yml
# Workflow: Mandatory gate before any biological analysis; reads only small TSV previews.

source("analysis/mag_functional_annotation/R/common.R")
args <- parse_stage_args()
cfg <- read_magfunc_config(args$config)
assert_packages()
inputs <- resolve_input_paths(cfg)
paths <- output_paths(cfg)
guard_stage(file.path(paths$qc, "00_validate_inputs.complete"), args, cfg)

missing_files <- names(inputs)[!vapply(inputs, file.exists, logical(1))]
if (length(missing_files) > 0) {
  detail <- paste(paste0(missing_files, "=", unlist(inputs[missing_files])), collapse = "; ")
  stop("Required input files are missing: ", detail, call. = FALSE)
}

preview_target <- 1000L
cluster_lines <- count_file_lines(inputs$cluster_file)
annotation_lines <- count_file_lines(inputs$annotation_file)
annotation_file_header <- annotation_header(inputs$annotation_file)

cluster_preview <- read_cluster(inputs$cluster_file, nrows = preview_target)
annotation_columns <- unique(c(
  cfg$columns$annotation_query,
  "seed_ortholog", "evalue", "score", "Description", "Preferred_name",
  vapply(ontology_specs(cfg), `[[`, character(1), "column")
))
annotation_preview <- read_annotation(inputs$annotation_file, cfg, annotation_columns, nrows = preview_target)

expected_cluster_preview_rows <- min(preview_target, cluster_lines)
if (!isTRUE(all.equal(nrow(cluster_preview), expected_cluster_preview_rows))) {
  stop("Cluster preview ended early: expected ", expected_cluster_preview_rows,
       " rows but parsed ", nrow(cluster_preview), ".", call. = FALSE)
}
if (annotation_lines - annotation_file_header$line >= preview_target && nrow(annotation_preview) != preview_target) {
  stop("Annotation preview ended early: expected ", preview_target,
       " data rows but parsed ", nrow(annotation_preview), ".", call. = FALSE)
}

gut <- data.table::fread(inputs$abundance_file)
metadata <- data.table::fread(inputs$metadata_file)
taxonomy <- data.table::fread(inputs$taxonomy_file)

sample_col <- cfg$columns$sample_id
tax_mag_col <- cfg$columns$taxonomy_mag_id
required_gut <- sample_col
required_metadata <- unlist(cfg$columns[c("sample_id", "subject_id", "timepoint", "date", "altitude")], use.names = FALSE)
if (!required_gut %in% names(gut)) stop("gut.csv is missing sample identifier column: ", sample_col, call. = FALSE)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0) stop("Metadata is missing columns: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
if (!tax_mag_col %in% names(taxonomy)) stop("Taxonomy is missing MAG key column: ", tax_mag_col, call. = FALSE)

gut_mags <- setdiff(names(gut), sample_col)
taxonomy_mags <- as.character(taxonomy[[tax_mag_col]])
if (anyDuplicated(gut[[sample_col]])) stop("gut.csv contains duplicate sample IDs.", call. = FALSE)
if (anyDuplicated(gut_mags)) stop("gut.csv contains duplicate MAG feature names.", call. = FALSE)
if (anyDuplicated(taxonomy_mags)) stop("Taxonomy contains duplicate Genome keys.", call. = FALSE)
missing_gut_taxonomy <- setdiff(gut_mags, taxonomy_mags)
if (isTRUE(cfg$identifiers$require_all_gut_mags_in_taxonomy) && length(missing_gut_taxonomy) > 0) {
  stop(length(missing_gut_taxonomy), " gut MAGs are absent from taxonomy.", call. = FALSE)
}

metadata_timepoints <- sort(unique(as.character(metadata[[cfg$columns$timepoint]])))
excluded <- unlist(cfg$analysis$excluded_timepoints)
included <- unlist(cfg$analysis$ordered_timepoints)
unexpected <- setdiff(metadata_timepoints, c(excluded, included))
if (length(unexpected) > 0) stop("Unexpected metadata timepoints: ", paste(unexpected, collapse = ", "), call. = FALSE)

abundance <- as.matrix(gut[, ..gut_mags])
storage.mode(abundance) <- "double"
if (any(!is.finite(abundance)) || any(abundance < 0)) stop("gut.csv MAG values must be finite and non-negative.", call. = FALSE)
row_sums <- rowSums(abundance)

reported_cluster <- as.numeric(cfg$quality_control$reported_cluster_lines)
reported_annotation <- as.numeric(cfg$quality_control$reported_annotation_lines)
if (isTRUE(cfg$quality_control$enforce_reported_line_counts)) {
  if (!isTRUE(all.equal(cluster_lines, reported_cluster))) stop("Cluster line count differs from the configured reported count.", call. = FALSE)
  if (!isTRUE(all.equal(annotation_lines, reported_annotation))) stop("Annotation line count differs from the configured reported count.", call. = FALSE)
}

summary <- data.table::data.table(
  check = c(
    "cluster_file_lines", "reported_cluster_lines", "annotation_file_lines", "reported_annotation_lines",
    "cluster_preview_rows", "cluster_preview_columns", "annotation_header_fields",
    "annotation_preview_rows", "annotation_preview_columns", "gut_samples", "gut_mags",
    "taxonomy_mags", "gut_mags_missing_taxonomy", "metadata_rows", "metadata_subjects",
    "excluded_timepoints_present", "included_timepoints_present", "gut_row_sum_min",
    "gut_row_sum_median", "gut_row_sum_max", "gut_row_normalization",
    "protein_clustering_scope", "mmseqs_min_sequence_identity",
    "mmseqs_coverage_threshold", "mmseqs_coverage_mode", "mmseqs_parameter_status"
  ),
  value = as.character(c(
    cluster_lines, reported_cluster, annotation_lines, reported_annotation,
    nrow(cluster_preview), ncol(cluster_preview), length(annotation_file_header$names),
    nrow(annotation_preview), ncol(annotation_preview), nrow(gut), length(gut_mags),
    length(taxonomy_mags), length(missing_gut_taxonomy), nrow(metadata),
    data.table::uniqueN(metadata[[cfg$columns$subject_id]]),
    paste(intersect(excluded, metadata_timepoints), collapse = ","),
    paste(intersect(included, metadata_timepoints), collapse = ","),
    min(row_sums), stats::median(row_sums), max(row_sums), cfg$analysis$abundance_row_normalization,
    cfg$clustering_provenance$scope,
    cfg$clustering_provenance$mmseqs_min_sequence_identity %||% "not_provided",
    cfg$clustering_provenance$mmseqs_coverage_threshold %||% "not_provided",
    cfg$clustering_provenance$mmseqs_coverage_mode %||% "not_provided",
    cfg$clustering_provenance$status
  )),
  status = c(
    "informational", "reference", "informational", "reference", "pass", "pass", "pass",
    "pass", "pass", "pass", "pass",
    "pass", if (length(missing_gut_taxonomy) == 0) "pass" else "fail", "pass", "pass",
    "decision", "decision", "informational", "informational", "informational", "decision",
    "collaborator_report", "unresolved", "unresolved", "unresolved", "unresolved"
  )
)

table_path <- file.path(paths$qc, "Table_01_input_validation_summary.csv")
write_table(summary, table_path)
write_table(data.table(sampleID = gut[[sample_col]], abundance_sum_as_provided = row_sums), file.path(paths$qc, "gut_abundance_row_sums_as_provided.csv"))
writeLines(c(
  "MAG functional annotation input validation passed.",
  paste("Configuration:", cfg$.config_path),
  paste("Ta/Tb exclusion:", paste(excluded, collapse = ", ")),
  "gut.csv values will be used as provided; no row normalization is permitted.",
  "Full-file uniqueness and join checks are performed in stage 01."
), file.path(paths$qc, "00_validate_inputs.complete"))

register_results(cfg, data.table::data.table(
  result_number = "Table 01",
  category = "quality_control",
  question = "Are the configured inputs structurally ready for the full workflow?",
  artifact_type = "table",
  path = relative_to_root(table_path),
  description = "Input paths, reported line counts, canonical gut-taxonomy coverage, timepoint decision, and unnormalized abundance-scale checks.",
  stage = "00_validate_inputs"
))
write_session_info(cfg, "00_validate_inputs")
message("Stage 00 completed: ", table_path)
