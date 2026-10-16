#!/usr/bin/env Rscript

# Purpose: Read the complete MMseqs and eggNOG tables, audit uniqueness and
# joins, parse MAG keys, and build compact reusable catalogue intermediates.
# Inputs: Stage-00-validated files configured in mag_functional_annotation.yml.
# Outputs: QC tables, cluster summaries, representative annotations, term
# tables, and representative-to-MAG member counts.
# Dependencies: R, data.table, yaml, and common.R.
# Usage: Rscript analysis/mag_functional_annotation/01_prepare_catalogue.R --config config/mag_functional_annotation.yml
# Workflow: Full-data HPC preparation; biological stages consume its RDS files.

source("analysis/mag_functional_annotation/R/common.R")
args <- parse_stage_args()
cfg <- read_magfunc_config(args$config)
assert_packages()
inputs <- resolve_input_paths(cfg)
paths <- output_paths(cfg)
if (!file.exists(file.path(paths$qc, "00_validate_inputs.complete"))) stop("Stage 00 completion marker is missing.", call. = FALSE)
guard_stage(file.path(paths$qc, "01_prepare_catalogue.complete"), args, cfg)

message("Reading complete MMseqs cluster table...")
cluster_raw <- read_cluster(inputs$cluster_file)
raw_cluster_rows <- nrow(cluster_raw)
cluster <- unique(cluster_raw, by = c("representative_id", "member_id"))
duplicate_pairs <- raw_cluster_rows - nrow(cluster)
rm(cluster_raw); invisible(gc())

member_assignments <- cluster[, .(n_representatives = data.table::uniqueN(representative_id)), by = member_id]
multi_members <- member_assignments[n_representatives > 1L]
if (isTRUE(cfg$quality_control$fail_on_multi_assigned_members) && nrow(multi_members) > 0) {
  write_table(head(multi_members, cfg$quality_control$diagnostic_example_limit), file.path(paths$qc, "multi_assigned_member_examples.csv"))
  stop(nrow(multi_members), " member proteins map to more than one representative.", call. = FALSE)
}
rm(member_assignments); invisible(gc())

pattern <- cfg$identifiers$mag_pattern
cluster[, representative_mag := parse_mag_id(representative_id, pattern)]
cluster[, member_mag := parse_mag_id(member_id, pattern)]
unparsed_representatives <- unique(cluster[is.na(representative_mag), representative_id])
unparsed_members <- unique(cluster[is.na(member_mag), member_id])
limit <- as.integer(cfg$quality_control$diagnostic_example_limit)
if (length(unparsed_representatives) > 0) write_table(data.table(representative_id = head(unparsed_representatives, limit)), file.path(paths$qc, "unparsed_representative_examples.csv"))
if (length(unparsed_members) > 0) write_table(data.table(member_id = head(unparsed_members, limit)), file.path(paths$qc, "unparsed_member_examples.csv"))
if (isTRUE(cfg$identifiers$require_all_member_ids_parse) && length(unparsed_members) > 0) {
  stop(length(unparsed_members), " unique member protein IDs could not be parsed into MAG IDs.", call. = FALSE)
}

gut <- data.table::fread(inputs$abundance_file, nrows = 0L)
gut_mags <- setdiff(names(gut), cfg$columns$sample_id)
taxonomy <- data.table::fread(inputs$taxonomy_file)
taxonomy_mags <- as.character(taxonomy[[cfg$columns$taxonomy_mag_id]])
all_parsed_mags <- sort(unique(stats::na.omit(c(cluster$representative_mag, cluster$member_mag))))
parsed_not_taxonomy <- setdiff(all_parsed_mags, taxonomy_mags)
if (length(parsed_not_taxonomy) > 0) write_table(data.table(MAG = head(parsed_not_taxonomy, limit)), file.path(paths$qc, "parsed_mags_missing_taxonomy_examples.csv"))
if (isTRUE(cfg$identifiers$require_all_parsed_mags_in_taxonomy) && length(parsed_not_taxonomy) > 0) {
  stop(length(parsed_not_taxonomy), " parsed MAG IDs are absent from GDTB.csv.", call. = FALSE)
}

cluster_mag <- cluster[!is.na(member_mag), .(member_gene_count = .N), by = .(representative_id, member_mag)]
cluster_size <- cluster[, .(
  cluster_size = .N,
  representative_mag = representative_mag[[1]],
  self_membership = any(representative_id == member_id)
), by = representative_id]
cluster_mag_summary <- cluster_mag[, .(
  n_member_mags = data.table::uniqueN(member_mag),
  gut_member_gene_count = sum(member_gene_count[member_mag %in% gut_mags]),
  n_gut_member_mags = data.table::uniqueN(member_mag[member_mag %in% gut_mags])
), by = representative_id]
cluster_summary <- merge(cluster_size, cluster_mag_summary, by = "representative_id", all.x = TRUE)
cluster_summary[is.na(gut_member_gene_count), gut_member_gene_count := 0L]
cluster_summary[is.na(n_gut_member_mags), n_gut_member_mags := 0L]

message("Reading complete eggNOG annotation table...")
specs <- ontology_specs(cfg)
annotation_columns <- unique(c(
  cfg$columns$annotation_query,
  "seed_ortholog", "evalue", "score", "Description", "Preferred_name",
  vapply(specs, `[[`, character(1), "column")
))
annotation_raw <- read_annotation(inputs$annotation_file, cfg, annotation_columns)
raw_annotation_rows <- nrow(annotation_raw)
annotation <- unique(annotation_raw)
exact_annotation_duplicates <- raw_annotation_rows - nrow(annotation)
rm(annotation_raw); invisible(gc())

query_col <- cfg$columns$annotation_query
query_counts <- annotation[, .N, by = .(query = get(query_col))]
conflicting_queries <- query_counts[N > 1L]
if (nrow(conflicting_queries) > 0) {
  write_table(head(conflicting_queries, limit), file.path(paths$qc, "conflicting_annotation_query_examples.csv"))
  if (isTRUE(cfg$quality_control$fail_on_conflicting_annotation_queries)) {
    stop(nrow(conflicting_queries), " annotation queries retain conflicting duplicate rows after exact deduplication.", call. = FALSE)
  }
  data.table::setorderv(annotation, query_col)
  annotation <- annotation[, .SD[1L], by = query_col]
}

coverage_fields <- unique(c("seed_ortholog", "Description", "Preferred_name", vapply(specs, `[[`, character(1), "column")))
functional_fields <- setdiff(coverage_fields, "seed_ortholog")
present_matrix <- lapply(functional_fields, function(field) !is.na(annotation[[field]]) & nzchar(trimws(annotation[[field]])))
annotation[, annotation_field_count := Reduce(`+`, present_matrix)]
annotation[, has_any_annotation := annotation_field_count > 0L]
unknown_pattern <- cfg$analysis$unknown_description_pattern
annotation[, description_unknown := !is.na(Description) & grepl(unknown_pattern, Description, ignore.case = TRUE, perl = TRUE)]
annotation[, preferred_name_missing := is.na(Preferred_name) | !nzchar(trimws(Preferred_name))]
annotation_summary <- annotation[, .(
  query = get(query_col),
  has_annotation_row = TRUE,
  has_any_annotation,
  annotation_field_count,
  Description,
  Preferred_name,
  description_unknown,
  preferred_name_missing,
  seed_ortholog,
  evalue = suppressWarnings(as.numeric(evalue)),
  score = suppressWarnings(as.numeric(score))
)]

field_coverage <- data.table::rbindlist(lapply(coverage_fields, function(field) {
  present <- !is.na(annotation[[field]]) & nzchar(trimws(annotation[[field]]))
  data.table(field = field, annotated_queries = sum(present), total_queries = length(present), coverage = mean(present))
}))
write_table(field_coverage, file.path(paths$qc, "annotation_field_coverage.csv"))

term_dir <- file.path(paths$intermediate, "term_tables")
dir.create(term_dir, recursive = TRUE, showWarnings = FALSE)
term_counts <- data.table::rbindlist(lapply(specs, function(spec) {
  message("Normalizing annotation terms: ", spec$id)
  terms <- split_annotation_terms(annotation, query_col, spec)
  save_object(terms, file.path(term_dir, paste0(spec$id, ".rds")))
  data.table(ontology = spec$id, representative_term_pairs = nrow(terms), distinct_terms = data.table::uniqueN(terms$term), annotated_representatives = data.table::uniqueN(terms$query))
}))
write_table(term_counts, file.path(paths$qc, "annotation_term_table_counts.csv"))
rm(annotation); invisible(gc())

cluster_annotation <- merge(
  cluster_summary,
  annotation_summary,
  by.x = "representative_id",
  by.y = "query",
  all.x = TRUE
)
logical_cols <- c("has_annotation_row", "has_any_annotation", "description_unknown", "preferred_name_missing")
for (field in logical_cols) cluster_annotation[is.na(get(field)), (field) := FALSE]
cluster_annotation[, unknown_or_unannotated := !has_annotation_row | !has_any_annotation | description_unknown]

write_table(
  cluster_annotation[, .(
    representative_id, cluster_size, has_annotation_row,
    has_any_functional_annotation = has_any_annotation,
    functional_annotation_field_count = annotation_field_count,
    description_unknown, unknown_or_unannotated
  )],
  file.path(paths$qc, "representative_annotation_status.csv.gz")
)
write_table(
  cluster_annotation[unknown_or_unannotated == TRUE, .(
    representative_id, cluster_size, has_annotation_row,
    has_any_functional_annotation = has_any_annotation,
    functional_annotation_field_count = annotation_field_count,
    Description, Preferred_name, description_unknown
  )],
  file.path(paths$qc, "representatives_unknown_or_unannotated.csv.gz")
)

cluster_representatives <- unique(cluster_summary$representative_id)
annotation_queries <- unique(annotation_summary$query)
reps_without_annotation <- setdiff(cluster_representatives, annotation_queries)
annotations_without_cluster <- setdiff(annotation_queries, cluster_representatives)
if (length(reps_without_annotation) > 0) write_table(data.table(representative_id = head(reps_without_annotation, limit)), file.path(paths$qc, "representatives_without_annotation_examples.csv"))
if (length(annotations_without_cluster) > 0) write_table(data.table(query = head(annotations_without_cluster, limit)), file.path(paths$qc, "annotations_without_cluster_examples.csv"))

cross_mag <- cluster_summary[n_member_mags > 1L][order(-n_member_mags, -cluster_size)]
write_table(head(cross_mag, limit), file.path(paths$qc, "cross_mag_cluster_examples.csv"))

qc_summary <- data.table::data.table(
  metric = c(
    "raw_cluster_rows", "unique_representative_member_pairs", "exact_duplicate_cluster_pairs",
    "unique_representatives", "unique_member_proteins", "multi_assigned_members",
    "self_membership_representatives", "unparsed_representative_ids", "unparsed_member_ids",
    "unique_parsed_mags", "parsed_mags_matching_taxonomy", "parsed_mags_matching_gut",
    "gut_mags_represented_by_members", "cross_mag_representative_clusters",
    "raw_annotation_rows", "unique_annotation_rows", "exact_duplicate_annotation_rows",
    "unique_annotation_queries", "conflicting_annotation_queries", "representatives_with_annotation_row",
    "representatives_without_annotation_row", "representatives_with_functional_annotation",
    "representatives_unknown_or_unannotated", "annotation_queries_without_cluster"
  ),
  value = c(
    raw_cluster_rows, nrow(cluster), duplicate_pairs,
    length(cluster_representatives), data.table::uniqueN(cluster$member_id), nrow(multi_members),
    sum(cluster_summary$self_membership), length(unparsed_representatives), length(unparsed_members),
    length(all_parsed_mags), sum(all_parsed_mags %in% taxonomy_mags), sum(all_parsed_mags %in% gut_mags),
    data.table::uniqueN(cluster$member_mag[cluster$member_mag %in% gut_mags]), nrow(cross_mag),
    raw_annotation_rows, nrow(annotation_summary), exact_annotation_duplicates,
    data.table::uniqueN(annotation_summary$query), nrow(conflicting_queries),
    sum(cluster_representatives %in% annotation_queries), length(reps_without_annotation),
    sum(cluster_annotation$has_any_annotation), sum(cluster_annotation$unknown_or_unannotated),
    length(annotations_without_cluster)
  )
)
qc_summary[, percent_of_relevant_denominator := NA_real_]
qc_summary[metric == "parsed_mags_matching_taxonomy", percent_of_relevant_denominator := value / length(all_parsed_mags) * 100]
qc_summary[metric == "parsed_mags_matching_gut", percent_of_relevant_denominator := value / length(all_parsed_mags) * 100]
qc_summary[metric == "representatives_with_annotation_row", percent_of_relevant_denominator := value / length(cluster_representatives) * 100]
qc_summary[metric == "representatives_without_annotation_row", percent_of_relevant_denominator := value / length(cluster_representatives) * 100]
qc_summary[metric == "representatives_with_functional_annotation", percent_of_relevant_denominator := value / length(cluster_representatives) * 100]
qc_summary[metric == "representatives_unknown_or_unannotated", percent_of_relevant_denominator := value / length(cluster_representatives) * 100]

qc_path <- file.path(paths$qc, "Table_02_full_catalogue_qc_summary.csv")
write_table(qc_summary, qc_path)
write_table(cluster_summary[, .(
  clusters = .N,
  minimum = min(cluster_size),
  q25 = stats::quantile(cluster_size, 0.25),
  median = stats::median(cluster_size),
  mean = mean(cluster_size),
  q75 = stats::quantile(cluster_size, 0.75),
  p95 = stats::quantile(cluster_size, 0.95),
  p99 = stats::quantile(cluster_size, 0.99),
  maximum = max(cluster_size),
  singleton_fraction = mean(cluster_size == 1L)
)], file.path(paths$qc, "cluster_size_summary.csv"))

save_object(cluster_summary, file.path(paths$intermediate, "cluster_summary.rds"))
save_object(cluster_mag, file.path(paths$intermediate, "cluster_mag_member_counts.rds"))
save_object(annotation_summary, file.path(paths$intermediate, "annotation_summary.rds"))
save_object(cluster_annotation, file.path(paths$intermediate, "cluster_annotation_summary.rds"))
save_object(gut_mags, file.path(paths$intermediate, "gut_mag_ids.rds"))
writeLines("Stage 01 complete", file.path(paths$qc, "01_prepare_catalogue.complete"))

register_results(cfg, data.table::data.table(
  result_number = "Table 02",
  category = "quality_control",
  question = "Do cluster membership, annotation queries, MAG parsing, and joins satisfy the analysis contract?",
  artifact_type = "table",
  path = relative_to_root(qc_path),
  description = "Complete-file counts, uniqueness, parsing, taxonomy/gut coverage, cross-MAG clustering, and representative annotation coverage.",
  stage = "01_prepare_catalogue"
))
write_session_info(cfg, "01_prepare_catalogue")
message("Stage 01 completed: ", qc_path)
