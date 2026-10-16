#!/usr/bin/env Rscript

# Purpose: Audit the protein-cluster and MAG-level evidence behind the two
# selected functional candidates using reviewed, offline curation metadata.
# Inputs: Stage-01 catalogue intermediates, Stage-08 carrier membership,
# Stage-09 completion marker, taxonomy, and tracked candidate curation.
# Outputs: Numbered Tables 23-25, a completion marker, session information,
# and updated follow-up result indexes.
# Dependencies: R, data.table, yaml, common.R, and followup.R.
# Usage: Rscript analysis/mag_functional_annotation/10_audit_candidate_annotations.R --config config/mag_functional_annotation.yml
# Workflow: Protein-evidence audit only; it does not run a new statistical
# screen or establish genomic neighborhood, pathway completeness, or activity.

source("analysis/mag_functional_annotation/R/common.R")
source("analysis/mag_functional_annotation/R/followup.R")
args <- parse_stage_args()
cfg <- read_magfunc_config(args$config)
assert_packages()
paths <- output_paths(cfg)
inputs <- resolve_input_paths(cfg)

if (!file.exists(file.path(paths$qc, "09_analyze_shared_component_drivers.complete"))) {
  stop("Stage 09 completion marker is missing.", call. = FALSE)
}
marker <- file.path(paths$qc, "10_audit_candidate_annotations.complete")
guard_stage(marker, args, cfg)

candidates <- configured_followup_candidates(cfg)
curation <- configured_candidate_curation(cfg)
if (nrow(candidates) != 2L) {
  stop("The candidate annotation audit currently requires exactly two configured candidates.", call. = FALSE)
}
if (!setequal(candidates$candidate_id, curation$candidate_id)) {
  stop("Candidate identities and curation metadata do not match.", call. = FALSE)
}
candidate_ids <- candidates$candidate_id

selected_representatives <- data.table::rbindlist(lapply(seq_len(nrow(candidates)), function(i) {
  candidate <- candidates[i]
  terms <- read_object(file.path(
    paths$intermediate, "term_tables", paste0(candidate$ontology[[1L]], ".rds")
  ))
  terms <- data.table::as.data.table(terms)[term == candidate$term[[1L]]]
  if (nrow(terms) == 0L) {
    stop("No representative protein clusters were found for ", candidate$candidate_id[[1L]], ".", call. = FALSE)
  }
  unique(terms[, .(
    candidate_id = candidate$candidate_id[[1L]],
    ontology = candidate$ontology[[1L]],
    term = candidate$term[[1L]],
    representative_id = query
  )])
}))
if (anyDuplicated(selected_representatives, by = c("candidate_id", "representative_id"))) {
  stop("Candidate-representative relations must be unique.", call. = FALSE)
}

cluster_summary <- data.table::as.data.table(read_object(file.path(
  paths$intermediate, "cluster_summary.rds"
)))
annotation_summary <- data.table::as.data.table(read_object(file.path(
  paths$intermediate, "annotation_summary.rds"
)))
cluster_mag <- data.table::as.data.table(read_object(file.path(
  paths$intermediate, "cluster_mag_member_counts.rds"
)))
gut_mags <- as.character(read_object(file.path(paths$intermediate, "gut_mag_ids.rds")))

representative_labels <- selected_representatives[, .(
  candidate_terms_on_representative = paste(sort(unique(candidate_id)), collapse = ";"),
  candidate_terms_on_representative_count = data.table::uniqueN(candidate_id)
), by = representative_id]
same_representative_ids <- representative_labels[
  candidate_terms_on_representative_count == nrow(candidates), representative_id
]

table23_data <- merge(
  selected_representatives,
  cluster_summary,
  by = "representative_id",
  all.x = TRUE
)
table23_data <- merge(
  table23_data,
  annotation_summary,
  by.x = "representative_id",
  by.y = "query",
  all.x = TRUE
)
table23_data <- merge(table23_data, representative_labels, by = "representative_id", all.x = TRUE)
if (anyNA(table23_data$cluster_size) || anyNA(table23_data$representative_mag)) {
  stop("One or more configured candidate representatives are missing from Stage-01 cluster summaries.", call. = FALSE)
}
if (anyNA(table23_data$has_annotation_row) || any(!table23_data$has_annotation_row)) {
  stop("Every configured candidate representative must retain its annotation row.", call. = FALSE)
}
table23_data[, `:=`(
  total_member_gene_count = as.integer(cluster_size),
  carries_both_candidate_terms_same_representative =
    candidate_terms_on_representative_count == nrow(candidates)
)]
table23_data <- table23_data[, .(
  candidate_id, ontology, term, representative_id, representative_mag,
  Description, Preferred_name, seed_ortholog, score, evalue,
  cluster_size, total_member_gene_count, n_member_mags,
  gut_member_gene_count, n_gut_member_mags,
  candidate_terms_on_representative,
  carries_both_candidate_terms_same_representative
)]
data.table::setorder(table23_data, candidate_id, representative_id)

gut_cluster_mag <- cluster_mag[member_mag %in% gut_mags]
carrier_evidence <- merge(
  gut_cluster_mag,
  selected_representatives,
  by = "representative_id",
  allow.cartesian = TRUE
)
carrier_summary <- carrier_evidence[, .(
  representative_clusters = data.table::uniqueN(representative_id),
  member_gene_count = sum(member_gene_count)
), by = .(candidate_id, MAG = member_mag)]
if (!setequal(unique(carrier_summary$candidate_id), candidate_ids)) {
  stop("Candidate carrier evidence is incomplete.", call. = FALSE)
}

table15_path <- file.path(
  paths$carrier_decomposition, "Table_15_candidate_carrier_set_membership.csv"
)
table15 <- read_table(table15_path)
required_table15 <- c("MAG", "candidate_ids", "candidates_carried")
missing_table15 <- setdiff(required_table15, names(table15))
if (length(missing_table15) > 0L) {
  stop("Table 15 is missing columns: ", paste(missing_table15, collapse = ", "), call. = FALSE)
}

union_mags <- sort(unique(carrier_summary$MAG))
table24_data <- data.table::data.table(
  MAG = union_mags,
  candidate_1_id = candidate_ids[[1L]],
  candidate_2_id = candidate_ids[[2L]]
)
for (candidate_index in seq_along(candidate_ids)) {
  candidate_key <- candidate_ids[[candidate_index]]
  candidate_carriers <- carrier_summary[candidate_id == candidate_key, .(
    MAG,
    carries_candidate = TRUE,
    representative_clusters,
    member_gene_count
  )]
  data.table::setnames(
    candidate_carriers,
    c("carries_candidate", "representative_clusters", "member_gene_count"),
    paste0(
      "candidate_", candidate_index,
      c("_carried", "_representative_clusters", "_member_gene_count")
    )
  )
  table24_data <- merge(table24_data, candidate_carriers, by = "MAG", all.x = TRUE)
  carried_column <- paste0("candidate_", candidate_index, "_carried")
  representative_column <- paste0("candidate_", candidate_index, "_representative_clusters")
  member_column <- paste0("candidate_", candidate_index, "_member_gene_count")
  table24_data[is.na(get(carried_column)), (carried_column) := FALSE]
  table24_data[is.na(get(representative_column)), (representative_column) := 0L]
  table24_data[is.na(get(member_column)), (member_column) := 0L]
}
table24_data[, shared_carrier := candidate_1_carried & candidate_2_carried]
same_representative_mags <- unique(gut_cluster_mag[
  representative_id %in% same_representative_ids, member_mag
])
table24_data[, both_terms_on_same_representative_cluster := MAG %in% same_representative_mags]

taxonomy <- data.table::fread(inputs$taxonomy_file)
taxonomy_key <- cfg$columns$taxonomy_mag_id
if (anyDuplicated(taxonomy[[taxonomy_key]])) {
  stop("Taxonomy MAG identifiers must be unique for Stage 10.", call. = FALSE)
}
if (taxonomy_key != "MAG") data.table::setnames(taxonomy, taxonomy_key, "MAG")
if ("ordeer" %in% names(taxonomy) && !"order" %in% names(taxonomy)) taxonomy[, order := ordeer]
taxonomy_columns <- intersect(
  c("MAG", "domain", "phylum", "class", "order", "family", "genus", "species"),
  names(taxonomy)
)
table24_data <- merge(table24_data, taxonomy[, ..taxonomy_columns], by = "MAG", all.x = TRUE)
if (anyNA(table24_data$phylum) || anyNA(table24_data$genus)) {
  stop("Every candidate-union MAG must retain its taxonomy join.", call. = FALSE)
}
data.table::setorder(table24_data, -shared_carrier, MAG)

table15_candidate_lists <- strsplit(as.character(table15$candidate_ids), ";", fixed = TRUE)
for (candidate_index in seq_along(candidate_ids)) {
  candidate_key <- candidate_ids[[candidate_index]]
  expected_mags <- sort(table15$MAG[vapply(
    table15_candidate_lists, function(ids) candidate_key %in% ids, logical(1)
  )])
  observed_mags <- sort(table24_data[
    get(paste0("candidate_", candidate_index, "_carried")), MAG
  ])
  if (!identical(expected_mags, observed_mags)) {
    stop("Stage-10 carrier MAGs do not reconcile with Table 15 for ", candidate_key, ".", call. = FALSE)
  }
}
if (!setequal(table15$MAG, table24_data$MAG) ||
    sum(table15$candidates_carried == 2L) != sum(table24_data$shared_carrier)) {
  stop("Stage-10 union/shared carrier totals do not reconcile with Table 15.", call. = FALSE)
}

summarize_labels <- function(values, maximum_labels = 10L) {
  values <- trimws(as.character(values))
  values <- values[!is.na(values) & nzchar(values) & values != "-"]
  if (length(values) == 0L) return("")
  counts <- sort(table(values), decreasing = TRUE)
  keep <- seq_len(min(length(counts), maximum_labels))
  paste0(names(counts)[keep], " (", as.integer(counts[keep]), ")", collapse = "; ")
}

table25_data <- data.table::rbindlist(lapply(seq_len(nrow(curation)), function(i) {
  curated <- curation[i]
  candidate_rows <- table23_data[candidate_id == curated$candidate_id[[1L]]]
  carried_column <- paste0("candidate_", i, "_carried")
  carrier_count <- sum(table24_data[[carried_column]])
  same_representative_count <- sum(candidate_rows$carries_both_candidate_terms_same_representative)
  data.table::data.table(
    candidate_id = curated$candidate_id[[1L]],
    ontology = curated$ontology[[1L]],
    term = curated$term[[1L]],
    accepted_name = curated$accepted_name[[1L]],
    gene_symbol = curated$gene_symbol[[1L]],
    entity_type = curated$entity_type[[1L]],
    reaction = curated$reaction[[1L]],
    authoritative_source_names = curated$authoritative_source_names[[1L]],
    authoritative_source_urls = curated$authoritative_source_urls[[1L]],
    review_date = curated$review_date[[1L]],
    representative_protein_clusters = nrow(candidate_rows),
    representatives_with_description = sum(
      !is.na(candidate_rows$Description) & nzchar(trimws(candidate_rows$Description))
    ),
    unique_project_description_labels = data.table::uniqueN(
      candidate_rows[!is.na(Description) & nzchar(trimws(Description)), Description]
    ),
    top_project_description_labels = summarize_labels(candidate_rows$Description),
    top_project_preferred_name_labels = summarize_labels(candidate_rows$Preferred_name),
    preferred_name_matching_curated_symbol = sum(
      tolower(candidate_rows$Preferred_name) == tolower(curated$gene_symbol[[1L]]),
      na.rm = TRUE
    ),
    carrier_MAGs = carrier_count,
    shared_carrier_MAGs = sum(table24_data$shared_carrier),
    candidate_specific_carrier_MAGs = carrier_count - sum(table24_data$shared_carrier),
    carrier_union_MAGs = nrow(table24_data),
    representative_clusters_carrying_both_terms = length(same_representative_ids),
    carrier_MAGs_with_both_terms_on_same_representative = length(same_representative_mags),
    protein_overlap_status = if (length(same_representative_ids) > 0L) {
      "At least one representative cluster carries both configured terms"
    } else {
      "No representative cluster carries both configured terms"
    },
    relationship_statement = curated$relationship_statement[[1L]],
    pathway_claim_readiness = if (length(same_representative_ids) > 0L) {
      "Protein-level co-annotation is present, but locus and pathway completeness remain unresolved"
    } else {
      "MAG-level co-carriage only; not ready for a complete Eut or Pdu pathway claim"
    },
    candidate_claim_boundary = curated$candidate_claim_boundary[[1L]],
    relationship_claim_boundary = curated$relationship_claim_boundary[[1L]]
  )
}))
data.table::setorder(table25_data, candidate_id)

table23 <- file.path(paths$candidate_curation, "Table_23_candidate_representative_protein_clusters.csv")
table24 <- file.path(paths$candidate_curation, "Table_24_candidate_carrier_MAG_concordance.csv")
table25 <- file.path(paths$candidate_curation, "Table_25_candidate_curation_summary.csv")
write_table(table23_data, table23)
write_table(table24_data, table24)
write_table(table25_data, table25)

table_records <- data.table::data.table(
  result_number = sprintf("Table %02d", 23:25),
  category = "candidate_annotation_audit",
  question = c(
    "Which representative protein clusters carry each selected candidate annotation?",
    "How do candidate carriage and same-representative evidence reconcile across the 73-MAG union?",
    "What biological identity and claim boundary are supported for each selected candidate?"
  ),
  artifact_type = "table",
  path = vapply(c(table23, table24, table25), relative_to_root, character(1)),
  description = c(
    "Candidate-representative eggNOG labels, alignment scores, cluster sizes, and total/gut MAG membership.",
    "Exact candidate-union MAG carriage, representative counts, taxonomy, shared-carrier status, and same-representative overlap.",
    "Reviewed candidate identities, project annotation-label distributions, carrier reconciliation, and pathway-claim readiness."
  ),
  stage = "10_audit_candidate_annotations"
)
register_results(cfg, table_records)

registry <- read_table(file.path(paths$intermediate, "results_registry.csv"))
expected_numbers <- c(sprintf("Table %02d", 8:25), sprintf("Figure %02d", 7:15))
followup_registry <- registry[result_number %in% expected_numbers]
missing_numbers <- setdiff(expected_numbers, unique(followup_registry$result_number))
if (length(missing_numbers) > 0L) {
  stop("Follow-up registry is missing: ", paste(missing_numbers, collapse = ", "), call. = FALSE)
}
data.table::setorder(followup_registry, result_number, path)
followup_csv <- file.path(paths$root, "FOLLOWUP_RESULTS_INDEX.csv")
followup_md <- file.path(paths$root, "FOLLOWUP_RESULTS_INDEX.md")
write_table(followup_registry, followup_csv)
escape_md <- function(x) gsub("\\|", "\\\\|", as.character(x))
lines <- c(
  "# MAG Functional Annotation Follow-up Results Index", "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
  "## Analysis contract", "",
  "- The two candidates and shared component were selected from the same data; follow-up p-values are not independent replication.",
  "- `gut.csv` values remain unnormalized, and Ta/Tb remain excluded.",
  "- Stage 10 audits protein-cluster and MAG-carriage evidence; it does not run a new statistical screen.",
  "- Shared MAG carriage does not establish that both terms occur on the same representative protein or genomic locus.",
  "- Results describe predicted genomic potential, not expression, activity, pathway completion, metabolite production, adaptation, or causal altitude effects.",
  "- MMseqs identity, coverage, mode, version, and command remain an external provenance gap.", "",
  "Number | Category | Type | Path | Purpose", "--- | --- | --- | --- | ---"
)
for (row_index in seq_len(nrow(followup_registry))) {
  lines <- c(lines, paste(
    escape_md(followup_registry$result_number[[row_index]]),
    escape_md(followup_registry$category[[row_index]]),
    escape_md(followup_registry$artifact_type[[row_index]]),
    paste0("`", escape_md(followup_registry$path[[row_index]]), "`"),
    escape_md(followup_registry$description[[row_index]]),
    sep = " | "
  ))
}
lines <- c(lines, "", "## Interpretation order", "",
  "1. Read Table 25 for reviewed identities and explicit biological claim boundaries.",
  "2. Use Tables 23-24 to audit representative annotations, carrier reconciliation, and same-representative evidence.",
  "3. Read Table 11 before interpreting the robustness of either original selected candidate.",
  "4. Read Table 17 to establish whether temporal evidence resides in shared or candidate-specific carriers.",
  "5. Use Tables 18-22 and Figures 13-15 for shared-component contrasts, participant influence, and numerical drivers.",
  "6. Treat all follow-up findings as selection-aware same-data explanation, not independent replication or pathway completion."
)
writeLines(lines, followup_md)

writeLines(c(
  "Stage 10 complete",
  paste0(candidate_ids[[1L]], " representative protein clusters: ", sum(table23_data$candidate_id == candidate_ids[[1L]])),
  paste0(candidate_ids[[2L]], " representative protein clusters: ", sum(table23_data$candidate_id == candidate_ids[[2L]])),
  paste0(candidate_ids[[1L]], " carrier MAGs: ", sum(table24_data$candidate_1_carried)),
  paste0(candidate_ids[[2L]], " carrier MAGs: ", sum(table24_data$candidate_2_carried)),
  paste0("Shared carrier MAGs: ", sum(table24_data$shared_carrier)),
  paste0("Carrier union MAGs: ", nrow(table24_data)),
  paste0("Representative clusters carrying both terms: ", length(same_representative_ids)),
  paste0("Carrier MAGs with both terms on the same representative: ", length(same_representative_mags)),
  "Relationship claim: co-carriage compatible with bacterial-microcompartment metabolism.",
  "No locus, operon, complete-pathway, expression, activity, metabolite-production, adaptation, or causal claim is supported."
), marker)
write_session_info(cfg, "10_audit_candidate_annotations")
message("Stage 10 completed: ", table25)
