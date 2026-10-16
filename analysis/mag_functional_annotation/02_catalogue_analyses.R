#!/usr/bin/env Rscript

# Purpose: Run catalogue-level functional convergence, cluster-size/unknown,
# and gut-MAG taxonomic distribution analyses from stage-01 intermediates.
# Inputs: Cluster, annotation, term, MAG-membership, gut, and taxonomy files.
# Outputs: Numbered result tables and MAG-by-function matrices for longitudinal work.
# Dependencies: R, data.table, yaml, splines, and common.R.
# Usage: Rscript analysis/mag_functional_annotation/02_catalogue_analyses.R --config config/mag_functional_annotation.yml
# Workflow: Implements priority questions 1-3 without treating member counts as abundance.

source("analysis/mag_functional_annotation/R/common.R")
args <- parse_stage_args()
cfg <- read_magfunc_config(args$config)
assert_packages()
inputs <- resolve_input_paths(cfg)
paths <- output_paths(cfg)
if (!file.exists(file.path(paths$qc, "01_prepare_catalogue.complete"))) stop("Stage 01 completion marker is missing.", call. = FALSE)
guard_stage(file.path(paths$qc, "02_catalogue_analyses.complete"), args, cfg)

cluster_annotation <- read_object(file.path(paths$intermediate, "cluster_annotation_summary.rds"))
cluster_mag <- read_object(file.path(paths$intermediate, "cluster_mag_member_counts.rds"))
gut_mags <- read_object(file.path(paths$intermediate, "gut_mag_ids.rds"))
specs <- ontology_specs(cfg)

message("Calculating functional convergence and member-gene representation...")
convergence_tables <- lapply(specs, function(spec) {
  terms <- read_object(file.path(paths$intermediate, "term_tables", paste0(spec$id, ".rds")))
  joined <- merge(
    terms,
    cluster_annotation[, .(query = representative_id, cluster_size, gut_member_gene_count, n_member_mags, n_gut_member_mags)],
    by = "query",
    all.x = FALSE,
    all.y = FALSE
  )
  result <- joined[, .(
    cluster_richness = data.table::uniqueN(query),
    member_gene_representation = sum(cluster_size),
    mean_members_per_cluster = as.numeric(mean(cluster_size)),
    median_members_per_cluster = as.numeric(stats::median(cluster_size)),
    gut_cluster_richness = data.table::uniqueN(query[gut_member_gene_count > 0]),
    gut_member_gene_representation = sum(gut_member_gene_count),
    median_member_mags_per_cluster = as.numeric(stats::median(n_member_mags)),
    maximum_member_mags_per_cluster = as.numeric(max(n_member_mags))
  ), by = term]
  result[, ontology := spec$id]
  data.table::setcolorder(result, c("ontology", "term", setdiff(names(result), c("ontology", "term"))))
  data.table::setorder(result, -cluster_richness, -member_gene_representation, term)
  write_table(result, file.path(paths$convergence, paste0("functional_convergence_", spec$id, ".csv.gz")))
  result
})
convergence <- data.table::rbindlist(convergence_tables, fill = TRUE)
convergence_path <- file.path(paths$convergence, "Table_03_functional_convergence_all_ontologies.csv.gz")
write_table(convergence, convergence_path)

message("Testing annotation and unknown-protein patterns across cluster size...")
unknown_data <- copy(cluster_annotation)
unknown_data[, log10_cluster_size := log10(cluster_size)]
unknown_data[, size_group := cut(
  cluster_size,
  breaks = c(-Inf, 1, 5, 20, 100, Inf),
  labels = c("1", "2-5", "6-20", "21-100", ">100"),
  right = TRUE
)]
unknown_summary <- unknown_data[, .(
  clusters = .N,
  annotation_row_fraction = mean(has_annotation_row),
  any_annotation_fraction = mean(has_any_annotation),
  unknown_or_unannotated_fraction = mean(unknown_or_unannotated),
  description_unknown_fraction = mean(description_unknown),
  median_annotation_fields = as.numeric(stats::median(annotation_field_count, na.rm = TRUE))
), by = size_group]

fit_binomial_spline <- function(outcome) {
  grouped <- unknown_data[, .(
    successes = sum(get(outcome), na.rm = TRUE),
    trials = sum(!is.na(get(outcome)))
  ), by = log10_cluster_size]
  grouped[, failures := trials - successes]
  distinct_sizes <- nrow(grouped)
  rhs <- if (distinct_sizes >= 3L) {
    paste0("splines::ns(log10_cluster_size, df = ", min(4L, distinct_sizes - 1L), ")")
  } else {
    "log10_cluster_size"
  }
  formula <- stats::as.formula(paste0("cbind(successes, failures) ~ ", rhs))
  tryCatch(stats::glm(formula, family = stats::binomial(), data = grouped), error = function(e) e)
}
models <- list(
  annotation_probability = fit_binomial_spline("has_annotation_row"),
  unknown_probability = fit_binomial_spline("unknown_or_unannotated")
)
model_status <- data.table::rbindlist(lapply(names(models), function(name) {
  model <- models[[name]]
  data.table::data.table(
    model = name,
    status = if (inherits(model, "error")) conditionMessage(model) else "fitted"
  )
}))
write_table(model_status, file.path(paths$unknown, "unknown_annotation_model_status.csv"))
failed_models <- model_status[status != "fitted"]
if (nrow(failed_models) > 0L) {
  stop(
    "Cluster-size annotation models failed: ",
    paste0(failed_models$model, " [", failed_models$status, "]", collapse = "; "),
    ". See unknown_annotation_model_status.csv.",
    call. = FALSE
  )
}
grid <- data.table(log10_cluster_size = seq(min(unknown_data$log10_cluster_size), max(unknown_data$log10_cluster_size), length.out = 200L))
prediction_tables <- lapply(names(models), function(name) {
  model <- models[[name]]
  pred <- stats::predict(model, newdata = grid, type = "link", se.fit = TRUE)
  data.table(
    model = name,
    log10_cluster_size = grid$log10_cluster_size,
    cluster_size = 10^grid$log10_cluster_size,
    probability = stats::plogis(pred$fit),
    lower = stats::plogis(pred$fit - 1.96 * pred$se.fit),
    upper = stats::plogis(pred$fit + 1.96 * pred$se.fit),
    status = "fitted"
  )
})
unknown_predictions <- data.table::rbindlist(prediction_tables, fill = TRUE)

singleton_model <- tryCatch(
  stats::glm(unknown_or_unannotated ~ I(cluster_size == 1L), family = stats::binomial(), data = unknown_data),
  error = function(e) e
)
if (inherits(singleton_model, "error") || nrow(summary(singleton_model)$coefficients) < 2L) {
  failure_status <- if (inherits(singleton_model, "error")) conditionMessage(singleton_model) else "singleton and larger clusters were not both represented"
  singleton_result <- data.table(contrast = "singleton_vs_larger", log_odds = NA_real_, odds_ratio = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_, p_value = NA_real_, status = failure_status)
} else {
  coefficient <- summary(singleton_model)$coefficients[2L, ]
  singleton_result <- data.table(
    contrast = "singleton_vs_larger",
    log_odds = unname(coefficient[[1]]),
    odds_ratio = exp(unname(coefficient[[1]])),
    lower_95 = exp(unname(coefficient[[1]] - 1.96 * coefficient[[2]])),
    upper_95 = exp(unname(coefficient[[1]] + 1.96 * coefficient[[2]])),
    p_value = unname(coefficient[[4]]),
    status = "fitted"
  )
}
unknown_path <- file.path(paths$unknown, "Table_04_unknown_annotation_by_cluster_size.csv")
write_table(unknown_summary, unknown_path)
write_table(unknown_predictions, file.path(paths$unknown, "unknown_annotation_probability_predictions.csv"))
write_table(singleton_result, file.path(paths$unknown, "singleton_unknown_enrichment_model.csv"))

message("Building gut MAG-by-function and taxonomic distribution tables...")
taxonomy <- data.table::fread(inputs$taxonomy_file)
tax_id <- cfg$columns$taxonomy_mag_id
if ("ordeer" %in% names(taxonomy) && !"order" %in% names(taxonomy)) taxonomy[, order := ordeer]
taxonomy_gut <- taxonomy[get(tax_id) %in% gut_mags]
gut_cluster_mag <- cluster_mag[member_mag %in% gut_mags]
mag_function_dir <- file.path(paths$intermediate, "mag_function_tables")
dir.create(mag_function_dir, recursive = TRUE, showWarnings = FALSE)

taxonomy_summaries <- lapply(specs, function(spec) {
  message("Aggregating MAG functions and taxonomy: ", spec$id)
  terms <- read_object(file.path(paths$intermediate, "term_tables", paste0(spec$id, ".rds")))
  member_terms <- merge(
    gut_cluster_mag,
    terms,
    by.x = "representative_id",
    by.y = "query",
    all = FALSE,
    allow.cartesian = TRUE
  )
  mag_function <- member_terms[, .(
    member_gene_count = sum(member_gene_count),
    cluster_count = data.table::uniqueN(representative_id)
  ), by = .(MAG = member_mag, term)]
  mag_function[, ontology := spec$id]
  data.table::setcolorder(mag_function, c("ontology", "MAG", "term", "cluster_count", "member_gene_count"))
  save_object(mag_function, file.path(mag_function_dir, paste0(spec$id, ".rds")))
  write_table(mag_function, file.path(paths$taxonomy, paste0("mag_function_", spec$id, ".csv.gz")))

  annotated <- merge(mag_function, taxonomy_gut, by.x = "MAG", by.y = tax_id, all.x = TRUE)
  summary <- annotated[, .(
    carrier_mags = data.table::uniqueN(MAG),
    mag_cluster_occurrences = sum(cluster_count),
    member_gene_sequences = sum(member_gene_count)
  ), by = .(ontology, term, phylum, genus)]
  write_table(summary, file.path(paths$taxonomy, paste0("taxonomy_function_distribution_", spec$id, ".csv.gz")))
  summary
})
taxonomy_summary <- data.table::rbindlist(taxonomy_summaries, fill = TRUE)
taxonomy_path <- file.path(paths$taxonomy, "Table_05_taxonomic_function_distribution_all_ontologies.csv.gz")
write_table(taxonomy_summary, taxonomy_path)
writeLines("Stage 02 complete", file.path(paths$qc, "02_catalogue_analyses.complete"))

register_results(cfg, data.table::data.table(
  result_number = c("Table 03", "Table 04", "Table 05"),
  category = c("functional_convergence", "unknown_by_cluster_size", "taxonomic_distribution"),
  question = c(
    "How many distinct sequence clusters and member genes support each predicted function?",
    "Are small clusters more often unannotated or poorly characterized?",
    "Which gut MAG taxa carry each predicted function?"
  ),
  artifact_type = "table",
  path = vapply(c(convergence_path, unknown_path, taxonomy_path), relative_to_root, character(1)),
  description = c(
    "Cluster richness and member-gene representation; member counts are catalogue membership, not abundance.",
    "Annotation and unknown fractions across prespecified cluster-size groups, with separate model outputs.",
    "MAG carrier counts and member-gene sequence counts by function, phylum, and genus."
  ),
  stage = "02_catalogue_analyses"
))
write_session_info(cfg, "02_catalogue_analyses")
message("Stage 02 completed.")
