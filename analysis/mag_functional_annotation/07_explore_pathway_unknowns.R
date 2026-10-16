#!/usr/bin/env Rscript

# Purpose: Extend selected-candidate interpretation with functional co-carriage
# and audit whether unknown singleton proteins can be distinguished from likely
# technical fragmentation with the currently available inputs.
# Inputs: Stage-01/02 intermediates, Stage-06 outputs, taxonomy, and optional
# protein-level or MAG-quality diagnostic tables.
# Outputs: Numbered Tables 12-14, Figure 10, and FOLLOWUP_RESULTS_INDEX files.
# Dependencies: R, data.table, ggplot2, yaml, common.R, and followup.R.
# Usage: Rscript analysis/mag_functional_annotation/07_explore_pathway_unknowns.R --config config/mag_functional_annotation.yml
# Workflow: Explicitly exploratory; co-carriage is not pathway completeness and
# unknown status is not evidence of biological novelty.

source("analysis/mag_functional_annotation/R/common.R")
source("analysis/mag_functional_annotation/R/followup.R")
args <- parse_stage_args()
cfg <- read_magfunc_config(args$config)
assert_packages()
paths <- output_paths(cfg)
inputs <- resolve_input_paths(cfg)
if (!file.exists(file.path(paths$qc, "06_validate_selected_candidates.complete"))) {
  stop("Stage 06 completion marker is missing.", call. = FALSE)
}
marker <- file.path(paths$qc, "07_explore_pathway_unknowns.complete")
guard_stage(marker, args, cfg)

candidates <- configured_followup_candidates(cfg)
gut_mags <- read_object(file.path(paths$intermediate, "gut_mag_ids.rds"))
gut_mags <- sort(unique(as.character(gut_mags)))
candidate_carrier_table <- read_table(file.path(paths$targeted_validation, "Table_08_candidate_carrier_MAG_taxonomy.csv"))
co_ontologies <- unlist(cfg$followup$co_carriage_ontologies)
minimum_co_carriage_mags <- as.integer(cfg$followup$minimum_co_carriage_mags)

message("Testing candidate-associated functional co-carriage across gut MAGs...")
co_carriage_rows <- list()
for (candidate_index in seq_len(nrow(candidates))) {
  candidate <- candidates[candidate_index]
  candidate_id <- candidate$candidate_id[[1L]]
  candidate_carriers <- unique(candidate_carrier_table[candidate_id == candidate$candidate_id[[1L]], MAG])
  candidate_carriers <- intersect(candidate_carriers, gut_mags)
  noncandidate_mags <- setdiff(gut_mags, candidate_carriers)
  if (length(candidate_carriers) == 0L || length(noncandidate_mags) == 0L) {
    stop("Candidate carrier/noncarrier comparison is unavailable for ", candidate_id, ".", call. = FALSE)
  }
  for (target_ontology in co_ontologies) {
    mag_function_path <- file.path(paths$intermediate, "mag_function_tables", paste0(target_ontology, ".rds"))
    if (!file.exists(mag_function_path)) {
      warning("Skipping unavailable co-carriage ontology: ", target_ontology, call. = FALSE)
      next
    }
    presence <- unique(data.table::as.data.table(read_object(mag_function_path))[MAG %in% gut_mags, .(MAG, term)])
    term_counts <- presence[, .(all_term_carrier_mags = data.table::uniqueN(MAG)), by = term]
    eligible_terms <- term_counts[all_term_carrier_mags >= minimum_co_carriage_mags, term]
    presence <- presence[term %in% eligible_terms]
    if (nrow(presence) == 0L) next
    results <- presence[, {
      term_carriers <- unique(MAG)
      a <- length(intersect(candidate_carriers, term_carriers))
      b <- length(candidate_carriers) - a
      c <- length(intersect(noncandidate_mags, term_carriers))
      d <- length(noncandidate_mags) - c
      test <- tryCatch(stats::fisher.test(matrix(c(a, b, c, d), nrow = 2L)), error = function(e) NULL)
      list(
        candidate_carrier_mags_with_term = a,
        candidate_carrier_mags_without_term = b,
        noncandidate_mags_with_term = c,
        noncandidate_mags_without_term = d,
        candidate_carrier_fraction = a / length(candidate_carriers),
        noncandidate_carrier_fraction = c / length(noncandidate_mags),
        odds_ratio = if (is.null(test)) NA_real_ else unname(test$estimate),
        lower_95 = if (is.null(test)) NA_real_ else unname(test$conf.int[[1L]]),
        upper_95 = if (is.null(test)) NA_real_ else unname(test$conf.int[[2L]]),
        p_value = if (is.null(test)) NA_real_ else test$p.value
      )
    }, by = term]
    results[, `:=`(
      candidate_id = candidate_id,
      candidate_ontology = candidate$ontology[[1L]],
      candidate_term = candidate$term[[1L]],
      comparison_ontology = target_ontology,
      candidate_carrier_mags = length(candidate_carriers),
      noncandidate_mags = length(noncandidate_mags),
      same_as_candidate_term = target_ontology == candidate$ontology[[1L]] & term == candidate$term[[1L]],
      analysis_scope = "exploratory MAG-level co-carriage; not pathway completeness or activity"
    )]
    results[, p_adj_bh := safe_bh(p_value)]
    results[, significant_bh := is.finite(p_adj_bh) & p_adj_bh <= as.numeric(cfg$analysis$fdr_threshold)]
    co_carriage_rows[[paste(candidate_id, target_ontology, sep = "_")]] <- results
  }
}
co_carriage <- data.table::rbindlist(co_carriage_rows, fill = TRUE)
if (nrow(co_carriage) > 0L) {
  data.table::setorder(co_carriage, candidate_id, comparison_ontology, p_adj_bh, p_value, term, na.last = TRUE)
}

message("Auditing unknown singleton proteins with available catalogue metadata...")
cluster_annotation <- data.table::as.data.table(read_object(file.path(paths$intermediate, "cluster_annotation_summary.rds")))
singletons <- cluster_annotation[cluster_size == 1L]
provenance <- parse_mag_construction_provenance(singletons$representative_mag)
singletons <- cbind(singletons, provenance[, -"MAG"])
taxonomy <- data.table::fread(inputs$taxonomy_file)
tax_id <- cfg$columns$taxonomy_mag_id
if ("ordeer" %in% names(taxonomy) && !"order" %in% names(taxonomy)) taxonomy[, order := ordeer]
if (anyDuplicated(taxonomy[[tax_id]])) stop("Taxonomy MAG identifiers must be unique for singleton joins.", call. = FALSE)
singletons <- merge(singletons, taxonomy, by.x = "representative_mag", by.y = tax_id, all.x = TRUE)
singletons[, gut_selected_MAG := representative_mag %in% gut_mags]

wilson_interval <- function(successes, total, z = 1.96) {
  if (!is.finite(total) || total <= 0) return(c(NA_real_, NA_real_))
  proportion <- successes / total
  denominator <- 1 + z^2 / total
  center <- (proportion + z^2 / (2 * total)) / denominator
  half <- z * sqrt((proportion * (1 - proportion) + z^2 / (4 * total)) / total) / denominator
  c(max(0, center - half), min(1, center + half))
}

summarize_unknown_dimension <- function(data, scope, dimension, column = NULL) {
  if (is.null(column)) {
    out <- data[, .(
      singleton_clusters = .N,
      unknown_or_unannotated_singletons = sum(unknown_or_unannotated),
      unknown_or_unannotated_fraction = mean(unknown_or_unannotated),
      annotation_row_fraction = mean(has_annotation_row),
      any_annotation_fraction = mean(has_any_annotation),
      description_unknown_fraction = mean(description_unknown),
      median_eggnog_score_when_present = stats::median(score, na.rm = TRUE)
    )]
    out[, level := "all"]
  } else {
    work <- data[!is.na(get(column)) & nzchar(as.character(get(column)))]
    out <- work[, .(
      singleton_clusters = .N,
      unknown_or_unannotated_singletons = sum(unknown_or_unannotated),
      unknown_or_unannotated_fraction = mean(unknown_or_unannotated),
      annotation_row_fraction = mean(has_annotation_row),
      any_annotation_fraction = mean(has_any_annotation),
      description_unknown_fraction = mean(description_unknown),
      median_eggnog_score_when_present = stats::median(score, na.rm = TRUE)
    ), by = .(level = as.character(get(column)))]
  }
  if (nrow(out) == 0L) {
    out[, `:=`(
      lower_95_wilson = numeric(),
      upper_95_wilson = numeric(),
      scope = character(),
      dimension = character(),
      interpretation = character()
    )]
    return(out)
  }
  intervals <- t(mapply(wilson_interval, out$unknown_or_unannotated_singletons, out$singleton_clusters))
  out[, `:=`(
    lower_95_wilson = intervals[, 1L],
    upper_95_wilson = intervals[, 2L],
    scope = scope,
    dimension = dimension,
    interpretation = "descriptive catalogue fraction; protein clusters are not biological replicates"
  )]
  out
}

scopes <- list(global_catalogue = singletons, gut_selected_MAGs = singletons[gut_selected_MAG == TRUE])
dimensions <- list(
  overall = NULL,
  construction_project = "construction_project",
  construction_specimen = "construction_specimen",
  construction_tier = "construction_tier",
  phylum = "phylum",
  genus = "genus"
)
unknown_rows <- list()
for (scope_name in names(scopes)) {
  for (dimension_name in names(dimensions)) {
    unknown_rows[[paste(scope_name, dimension_name, sep = "_")]] <- summarize_unknown_dimension(
      scopes[[scope_name]], scope_name, dimension_name, dimensions[[dimension_name]]
    )
  }
}
unknown_summary <- data.table::rbindlist(unknown_rows, fill = TRUE)
unknown_summary <- unknown_summary[
  dimension == "overall" |
    singleton_clusters >= as.integer(cfg$followup$minimum_unknown_group_clusters)
]
data.table::setcolorder(unknown_summary, c(
  "scope", "dimension", "level", "singleton_clusters",
  "unknown_or_unannotated_singletons", "unknown_or_unannotated_fraction",
  "lower_95_wilson", "upper_95_wilson"
))

protein_metrics_path <- resolve_optional_followup_input(
  cfg, "protein_metrics_file", cfg$followup$optional_input_environment_overrides$protein_metrics_file
)
mag_quality_path <- resolve_optional_followup_input(
  cfg, "mag_quality_file", cfg$followup$optional_input_environment_overrides$mag_quality_file
)
readiness <- data.table::data.table(
  evidence_component = c(
    "protein_length_amino_acids", "partial_gene_or_contig_edge_status",
    "MAG_completeness_and_contamination", "MMseqs_identity_and_coverage",
    "authoritative_pathway_membership"
  ),
  required_for = c(
    "distinguish short predictions from candidate novelty",
    "test gene-calling or contig-edge fragmentation",
    "test whether unknown fractions track MAG quality",
    "interpret cluster redundancy and cross-MAG sharing",
    "test biochemical pathway completeness rather than co-carriage"
  ),
  source_path = c(
    protein_metrics_path, protein_metrics_path, mag_quality_path,
    NA_character_, NA_character_
  ),
  status = c(
    ifelse(!is.na(protein_metrics_path) && file.exists(protein_metrics_path), "available", "not_provided"),
    ifelse(!is.na(protein_metrics_path) && file.exists(protein_metrics_path), "available", "not_provided"),
    ifelse(!is.na(mag_quality_path) && file.exists(mag_quality_path), "available", "not_provided"),
    ifelse(is.null(cfg$clustering_provenance$mmseqs_min_sequence_identity), "not_provided", "available"),
    "not_provided"
  ),
  consequence = c(
    "Novelty-versus-short-protein interpretation is unresolved without this input.",
    "Fragmentation interpretation is unresolved without this input.",
    "MAG-quality confounding is unresolved without this input.",
    "Sequence-cluster boundaries cannot receive firm biological interpretation.",
    "Current results describe co-carriage only; they do not establish complete pathways."
  )
)

fragmentation_model <- data.table::data.table(
  term = character(), estimate = numeric(), standard_error = numeric(),
  z_value = numeric(), p_value = numeric(), status = character()
)
if (!is.na(protein_metrics_path) && file.exists(protein_metrics_path)) {
  protein_metrics <- data.table::fread(protein_metrics_path)
  required <- c("representative_id", "protein_length_aa", "partial_gene", "contig_edge_distance_bp")
  absent <- setdiff(required, names(protein_metrics))
  if (length(absent) > 0L) {
    readiness[evidence_component %in% c("protein_length_amino_acids", "partial_gene_or_contig_edge_status"),
      `:=`(status = "invalid_schema", consequence = paste("Missing columns:", paste(absent, collapse = ", ")))]
  } else {
    diagnostic <- merge(
      singletons[, .(representative_id, unknown_or_unannotated)],
      protein_metrics[, ..required],
      by = "representative_id",
      all = FALSE
    )
    diagnostic[, partial_gene := as.logical(partial_gene)]
    diagnostic <- diagnostic[
      is.finite(protein_length_aa) & protein_length_aa > 0 &
        is.finite(contig_edge_distance_bp) & contig_edge_distance_bp >= 0 &
        !is.na(partial_gene)
    ]
    if (nrow(diagnostic) > 0L && data.table::uniqueN(diagnostic$unknown_or_unannotated) == 2L) {
      model <- stats::glm(
        unknown_or_unannotated ~ scale(log10(protein_length_aa)) + partial_gene +
          scale(log10(contig_edge_distance_bp + 1)),
        family = stats::binomial(), data = diagnostic
      )
      coefficients <- summary(model)$coefficients
      fragmentation_model <- data.table::data.table(
        term = rownames(coefficients),
        estimate = coefficients[, "Estimate"],
        standard_error = coefficients[, "Std. Error"],
        z_value = coefficients[, "z value"],
        p_value = coefficients[, "Pr(>|z|)"],
        status = "fitted catalogue-level diagnostic; clusters are not biological replicates"
      )
    } else {
      fragmentation_model <- data.table::data.table(
        term = NA_character_, estimate = NA_real_, standard_error = NA_real_,
        z_value = NA_real_, p_value = NA_real_, status = "insufficient matched protein metrics"
      )
    }
  }
}

mag_quality_summary <- data.table::data.table(
  quality_group = character(), singleton_clusters = integer(),
  unknown_or_unannotated_fraction = numeric(), status = character()
)
if (!is.na(mag_quality_path) && file.exists(mag_quality_path)) {
  mag_quality <- data.table::fread(mag_quality_path)
  required <- c("MAG", "completeness", "contamination")
  absent <- setdiff(required, names(mag_quality))
  if (length(absent) > 0L) {
    readiness[evidence_component == "MAG_completeness_and_contamination",
      `:=`(status = "invalid_schema", consequence = paste("Missing columns:", paste(absent, collapse = ", ")))]
  } else {
    quality_diagnostic <- merge(
      singletons[, .(MAG = representative_mag, unknown_or_unannotated)],
      mag_quality[, ..required], by = "MAG", all = FALSE
    )
    quality_diagnostic <- quality_diagnostic[is.finite(completeness) & is.finite(contamination)]
    if (nrow(quality_diagnostic) > 0L) {
      quality_diagnostic[, quality_score := completeness - 5 * contamination]
      quality_diagnostic[, quality_group := cut(
        quality_score,
        breaks = unique(stats::quantile(quality_score, probs = seq(0, 1, 0.25), na.rm = TRUE)),
        include.lowest = TRUE
      )]
      mag_quality_summary <- quality_diagnostic[!is.na(quality_group), .(
        singleton_clusters = .N,
        unknown_or_unannotated_fraction = mean(unknown_or_unannotated),
        status = "descriptive quality-score quartile; clusters are not biological replicates"
      ), by = quality_group]
    }
  }
}

table12 <- file.path(paths$exploratory_followup, "Table_12_candidate_function_co_carriage.csv.gz")
table13 <- file.path(paths$exploratory_followup, "Table_13_unknown_singleton_catalogue_audit.csv")
table14 <- file.path(paths$exploratory_followup, "Table_14_novelty_fragmentation_evidence_readiness.csv")
write_table(co_carriage, table12)
write_table(unknown_summary, table13)
write_table(readiness, table14)
write_table(fragmentation_model, file.path(paths$exploratory_followup, "optional_protein_fragmentation_model.csv"))
write_table(mag_quality_summary, file.path(paths$exploratory_followup, "optional_MAG_quality_unknown_summary.csv"))

figure_data <- unknown_summary[
  scope == "gut_selected_MAGs" &
    dimension %in% c("construction_specimen", "construction_tier", "phylum")
]
largest_phyla <- figure_data[dimension == "phylum"][order(-singleton_clusters), head(level, 10L)]
figure_data <- figure_data[dimension != "phylum" | level %in% largest_phyla]
if (nrow(figure_data) == 0L) stop("No eligible gut singleton groups are available for Figure 10.", call. = FALSE)
figure_data[, level_label := paste0(level, " (n=", scales::comma(singleton_clusters), ")")]
figure_data[, level_label := factor(level_label, levels = rev(unique(level_label[order(dimension, unknown_or_unannotated_fraction)])))]
dimension_labels <- c(
  construction_specimen = "Construction-source specimen",
  construction_tier = "Retained bin label",
  phylum = "GTDB phylum"
)
figure_data[, dimension_label := unname(dimension_labels[dimension])]
p10 <- ggplot2::ggplot(figure_data, ggplot2::aes(unknown_or_unannotated_fraction, level_label)) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = lower_95_wilson, xmax = upper_95_wilson),
    orientation = "y", width = 0, linewidth = 0.45, colour = magfunc_palette[["muted"]]
  ) +
  ggplot2::geom_point(shape = 21, fill = "white", colour = magfunc_palette[["blue"]], size = 2) +
  ggplot2::facet_wrap(~dimension_label, scales = "free_y", ncol = 1) +
  ggplot2::scale_x_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  ggplot2::labs(
    title = "Unknown or unannotated singleton clusters among gut-selected MAGs",
    subtitle = "Descriptive fractions with Wilson intervals; clusters are catalogue entries, not biological replicates",
    x = "Unknown or unannotated fraction", y = NULL
  ) + theme_nature()
files10 <- save_publication_plot(p10, "Figure_10_unknown_singleton_catalogue_audit", cfg, width = 7.4, height = 8.0)
write_table(figure_data, file.path(paths$reproducibility, "Figure_10_plot_data.csv"))

chart_map <- data.table::data.table(
  figure = sprintf("Figure %02d", 7:10),
  analytical_question = c(
    "Are selected candidate trajectories distributed across participants and timepoints?",
    "Which genera contribute to each selected candidate potential over time?",
    "Does omitting one participant materially alter each selected candidate result?",
    "Which gut-MAG groups contain larger descriptive fractions of unknown singleton clusters?"
  ),
  chart_family = c(
    "longitudinal trend", "small-multiple composition trend",
    "influence diagnostic", "uncertainty comparison"
  ),
  visual_encoding = c(
    "participant lines plus descriptive timepoint means",
    "top-three genera plus Other on a common 0-100% share scale",
    "raw omnibus p-value after each leave-one-participant-out refit on a log scale",
    "group fraction with 95% Wilson interval"
  ),
  denominator_or_unit = c(
    "sum of unnormalized gut.csv values for carrier MAGs",
    "within-sample candidate potential",
    "same-data selected candidate model; diagnostic only",
    "singleton protein clusters within each displayed gut-MAG group"
  ),
  palette_policy = c(
    "grey participants; orange descriptive mean",
    "single blue series; position and facet labels carry identity",
    "single blue series; point shape carries model diagnostic",
    "single blue estimate with muted interval"
  )
)
write_table(chart_map, file.path(paths$reproducibility, "followup_publication_chart_map.csv"))

table_records <- data.table::data.table(
  result_number = sprintf("Table %02d", 12:14),
  category = "pathway_unknown_exploration",
  question = c(
    "Which functions are enriched among MAGs carrying the selected candidates?",
    "How are unknown singleton clusters distributed across gut-MAG construction and taxonomic groups?",
    "Can current evidence distinguish protein novelty from assembly or gene-prediction fragmentation?"
  ),
  artifact_type = "table",
  path = vapply(c(table12, table13, table14), relative_to_root, character(1)),
  description = c(
    "Exploratory MAG-level functional co-carriage enrichment; not pathway completeness or activity.",
    "Descriptive unknown/unannotated singleton fractions with explicit catalogue denominators.",
    "Availability and schema readiness for protein length, partial-gene, contig-edge, MAG-quality, MMseqs, and pathway evidence."
  ),
  stage = "07_explore_pathway_unknowns"
)
figure_record <- data.table::data.table(
  result_number = "Figure 10",
  category = "pathway_unknown_exploration",
  question = "Which gut-MAG groups have higher descriptive fractions of unknown singleton clusters?",
  artifact_type = "figure",
  path = vapply(files10, relative_to_root, character(1)),
  description = "Descriptive unknown/unannotated singleton fractions; intervals do not make clusters biological replicates.",
  stage = "07_explore_pathway_unknowns"
)
register_results(cfg, data.table::rbindlist(list(table_records, figure_record), fill = TRUE))

registry <- data.table::fread(file.path(paths$intermediate, "results_registry.csv"))
followup_registry <- registry[
  result_number %in% c(sprintf("Table %02d", 8:14), sprintf("Figure %02d", 7:10))
]
data.table::setorder(followup_registry, result_number, path)
followup_csv <- file.path(paths$root, "FOLLOWUP_RESULTS_INDEX.csv")
followup_md <- file.path(paths$root, "FOLLOWUP_RESULTS_INDEX.md")
write_table(followup_registry, followup_csv)
escape_md <- function(x) gsub("\\|", "\\\\|", as.character(x))
lines <- c(
  "# MAG Functional Annotation Follow-up Results Index", "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
  "## Analysis contract", "",
  "- The two candidates were selected from the same Stage-03 screen; follow-up p-values are not independent replication.",
  "- The unnormalized gut.csv analysis remains primary; row-sum normalization is sensitivity-only.",
  "- Co-carriage does not establish pathway completeness, activity, or metabolite production.",
  "- Unknown singleton status does not establish biological novelty.", "",
  "Number | Category | Type | Path | Purpose", "--- | --- | --- | --- | ---"
)
for (i in seq_len(nrow(followup_registry))) {
  lines <- c(lines, paste(
    escape_md(followup_registry$result_number[[i]]),
    escape_md(followup_registry$category[[i]]),
    escape_md(followup_registry$artifact_type[[i]]),
    paste0("`", escape_md(followup_registry$path[[i]]), "`"),
    escape_md(followup_registry$description[[i]]),
    sep = " | "
  ))
}
lines <- c(lines, "", "## Interpretation order", "",
  "1. Read Table 11 before interpreting the robustness of either selected signal.",
  "2. Read Table 09 and Figure 09 to identify participant influence.",
  "3. Use Tables 08 and 10 with Figures 07-08 for carrier and timepoint interpretation.",
  "4. Treat Table 12 as exploratory co-carriage, not pathway completeness.",
  "5. Read Table 14 before making any novelty-versus-fragmentation statement from Table 13 or Figure 10."
)
writeLines(lines, followup_md)

writeLines(c(
  "Stage 07 complete",
  "Functional co-carriage is exploratory and does not establish a complete pathway.",
  "Unknown singleton status alone cannot distinguish novelty from assembly or gene-prediction fragmentation.",
  "Consult Table 14 before making novelty claims."
), marker)
write_session_info(cfg, "07_explore_pathway_unknowns")
message("Stage 07 completed: ", followup_md)
