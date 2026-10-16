# Public replacement: neutral input-derived reporting; original preserved privately.
#!/usr/bin/env Rscript

# Purpose: Build a publication-quality four-panel summary of MAG eligibility,
# native Random Forest importance, and rank stability across subject-disjoint
# leave-one-subject-out (LOSO) models.
# Inputs: The five feature CSV files from one completed subject-validation run,
# its resolved configuration and performance outputs, and data/raw/GDTB.csv.
# Outputs: Standalone and combined figures, exact plot-ready tables, an exact
# taxonomy-join audit, caption, methods note, validation manifest, and session
# information below the selected run's figures/feature_stability directory.
# Dependencies: R packages declared in environment.yml; no model is refitted.
# Usage: Rscript scripts/plot_random_forest_feature_stability.R \
#   --feature-results-dir results/random_forest_subject_validation/development_20260724_150733/features \
#   --taxonomy-file data/raw/GDTB.csv \
#   --output-dir results/random_forest_subject_validation/development_20260724_150733/figures/feature_stability
# Relationship: Publication post-processing for a completed subject-aware
# Random Forest validation run. Native importance remains exploratory.

required_packages <- c(
  "dplyr", "ggplot2", "grid", "gridExtra", "ragg",
  "scales", "svglite", "viridisLite", "yaml"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    ". Use the repository environment-update workflow on a compute node.",
    call. = FALSE
  )
}

set.seed(20260716)

usage <- function() {
  cat(
    "Usage: Rscript scripts/plot_random_forest_feature_stability.R [options]\n\n",
    "Options:\n",
    "  --feature-results-dir PATH  Directory containing the five feature CSV files.\n",
    "  --taxonomy-file PATH        Exact MAG taxonomy CSV; default: data/raw/GDTB.csv.\n",
    "  --output-dir PATH           Figure and audit output directory.\n",
    "  --n-candidates INTEGER      Candidate MAG count; default: 20.\n",
    "  --output-prefix TEXT        Optional filename prefix.\n",
    "  --help                      Show this message.\n",
    sep = ""
  )
}

parse_args <- function(args) {
  defaults <- list(
    feature_results_dir = "results/random_forest_subject_validation/development_20260724_150733/features",
    taxonomy_file = "data/raw/GDTB.csv",
    output_dir = "results/random_forest_subject_validation/development_20260724_150733/figures/feature_stability",
    n_candidates = 20L,
    output_prefix = ""
  )
  mapping <- c(
    "--feature-results-dir" = "feature_results_dir",
    "--taxonomy-file" = "taxonomy_file",
    "--output-dir" = "output_dir",
    "--n-candidates" = "n_candidates",
    "--output-prefix" = "output_prefix"
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--help")) {
      usage()
      quit(status = 0L)
    }
    if (!(arg %in% names(mapping))) stop("Unknown option: ", arg, call. = FALSE)
    if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
    defaults[[mapping[[arg]]]] <- args[[i + 1L]]
    i <- i + 2L
  }
  defaults$n_candidates <- suppressWarnings(as.integer(defaults$n_candidates))
  if (!is.finite(defaults$n_candidates) || defaults$n_candidates < 1L) {
    stop("--n-candidates must be a positive integer.", call. = FALSE)
  }
  if (grepl("[/\\\\]", defaults$output_prefix)) {
    stop("--output-prefix must not contain a path separator.", call. = FALSE)
  }
  defaults
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
repo_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(repo_root, "environment.yml")) ||
    !file.exists(file.path(repo_root, "analysis", "random_forest_subject_validation", "run_random_forest_subject_validation.R"))) {
  stop("Run this script from the sagarmatha-microbiome-analysis repository root.", call. = FALSE)
}

normalize_existing <- function(path, label) {
  if (!file.exists(path) && !dir.exists(path)) stop(label, " does not exist: ", path, call. = FALSE)
  normalizePath(path, mustWork = TRUE)
}

feature_dir <- normalize_existing(args$feature_results_dir, "Feature-results directory")
taxonomy_path <- normalize_existing(args$taxonomy_file, "Taxonomy file")
output_dir <- normalizePath(args$output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

prefix <- if (nzchar(args$output_prefix)) paste0(args$output_prefix, "_") else ""
out_path <- function(name) file.path(output_dir, paste0(prefix, name))

input_names <- c(
  candidate = "candidate_MAG_summary.csv",
  eligibility = "feature_eligibility_by_fold.csv",
  importance = "feature_importance_by_fold.csv",
  stability = "feature_importance_stability.csv",
  selection = "feature_selection_frequency.csv"
)
input_paths <- file.path(feature_dir, input_names)
names(input_paths) <- names(input_names)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop("Missing feature input(s): ", paste(missing_inputs, collapse = "; "), call. = FALSE)
}

run_dir <- dirname(feature_dir)
resolved_config_path <- file.path(run_dir, "reports", "resolved_config.yml")
overall_metrics_path <- file.path(run_dir, "predictions", "loso_overall_metrics.csv")
per_class_metrics_path <- file.path(run_dir, "predictions", "loso_per_class_metrics.csv")
for (path in resolved_config_path) {
  if (!file.exists(path)) stop("Required run-level evidence is missing: ", path, call. = FALSE)
}

fmt_mtime <- function(path) format(file.info(path)$mtime, "%Y-%m-%d %H:%M:%S %Z", tz = "UTC")
cat("Random Forest MAG feature-stability plotting\n")
cat("Repository root: ", repo_root, "\n", sep = "")
cat("Feature-results directory: ", feature_dir, "\n", sep = "")
cat("Taxonomy file: ", taxonomy_path, "\n", sep = "")
cat("Output directory: ", output_dir, "\n", sep = "")
for (path in c(input_paths, taxonomy_path, resolved_config_path, overall_metrics_path, per_class_metrics_path)) {
  cat("Input: ", path, " | modified: ", fmt_mtime(path), "\n", sep = "")
}
cat("Package versions:\n")
for (pkg in required_packages) cat("  ", pkg, ": ", as.character(utils::packageVersion(pkg)), "\n", sep = "")

read_csv_checked <- function(path) {
  utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}
write_csv <- function(data, path) utils::write.csv(data, path, row.names = FALSE, na = "")
write_tsv <- function(data, path) utils::write.table(
  data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
require_columns <- function(data, required, source_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(source_name, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}
require_unique_key <- function(data, columns, source_name) {
  key <- do.call(paste, c(data[columns], sep = "\r"))
  duplicated_key <- duplicated(key) | duplicated(key, fromLast = TRUE)
  if (any(duplicated_key)) {
    examples <- unique(key[duplicated_key])[seq_len(min(5L, length(unique(key[duplicated_key]))))]
    stop(source_name, " contains duplicate key(s) for ", paste(columns, collapse = "+"),
         ": ", paste(examples, collapse = "; "), call. = FALSE)
  }
}
assert_same_set <- function(x, y, x_name, y_name) {
  if (!setequal(x, y)) {
    stop(x_name, " and ", y_name, " have different MAG identifier sets.", call. = FALSE)
  }
}
assert_close <- function(observed, expected, label, tolerance = 1e-10) {
  if (length(observed) != length(expected) || any(is.na(observed) != is.na(expected)) ||
      any(abs(observed - expected) > tolerance, na.rm = TRUE)) {
    stop(label, " does not agree with independently recomputed values.", call. = FALSE)
  }
}

candidate <- read_csv_checked(input_paths[["candidate"]])
eligibility <- read_csv_checked(input_paths[["eligibility"]])
importance <- read_csv_checked(input_paths[["importance"]])
stability <- read_csv_checked(input_paths[["stability"]])
selection <- read_csv_checked(input_paths[["selection"]])
taxonomy <- read_csv_checked(taxonomy_path)
# Performance summaries are context, not prerequisites for feature calculations.
read_optional_metrics <- function(path) {
  if (file.exists(path)) read_csv_checked(path) else data.frame()
}
overall_metrics <- read_optional_metrics(overall_metrics_path)
per_class_metrics <- read_optional_metrics(per_class_metrics_path)

require_columns(candidate, c(
  "feature", "importance_folds", "median_native_importance", "native_importance_iqr",
  "median_importance_rank", "importance_rank_iqr", "eligible_folds", "selected_folds",
  "selection_frequency_when_eligible"
), "candidate_MAG_summary.csv")
require_columns(eligibility, c(
  "MAG", "missing_fraction", "sample_prevalence", "subject_prevalence", "mean_abundance",
  "variance", "fold_id", "heldout_subject", "eligible", "selected"
), "feature_eligibility_by_fold.csv")
require_columns(importance, c(
  "feature", "fold_id", "heldout_subject", "MeanDecreaseAccuracy", "MeanDecreaseGini"
), "feature_importance_by_fold.csv")
require_columns(stability, c(
  "feature", "importance_folds", "median_native_importance", "native_importance_iqr",
  "median_importance_rank", "importance_rank_iqr"
), "feature_importance_stability.csv")
require_columns(selection, c(
  "feature", "eligible_folds", "selected_folds", "selection_frequency_when_eligible"
), "feature_selection_frequency.csv")
require_columns(taxonomy, "Genome", "GDTB.csv")
# Missing values/columns stay unavailable; malformed probabilities remain errors.
validate_probability <- function(x, label) {
  text <- trimws(as.character(x))
  missing <- is.na(x) | text == ""
  value <- if (is.numeric(x)) as.numeric(x) else suppressWarnings(as.numeric(text))
  if (any(!missing & (!is.finite(value) | value < 0 | value > 1))) {
    stop(label, " must contain probabilities in [0, 1] or missing values.", call. = FALSE)
  }
  value[missing] <- NA_real_
  value
}
if (nrow(overall_metrics) > 1L) stop("Overall performance must contain at most one pooled summary row.", call. = FALSE)
if (!nrow(overall_metrics)) overall_metrics <- data.frame(balanced_accuracy = NA_real_, macro_f1 = NA_real_)
for (metric in c("balanced_accuracy", "macro_f1")) {
  if (!metric %in% names(overall_metrics)) overall_metrics[[metric]] <- NA_real_
  overall_metrics[[metric]] <- validate_probability(overall_metrics[[metric]], metric)
}
if (!"class" %in% names(per_class_metrics) || !nrow(per_class_metrics)) {
  per_class_metrics <- data.frame(class = character(), recall = numeric())
} else {
  if (anyNA(per_class_metrics$class) || any(!nzchar(trimws(as.character(per_class_metrics$class)))) || anyDuplicated(per_class_metrics$class)) {
    stop("Performance class identifiers must be nonempty and unique.", call. = FALSE)
  }
  if (!"recall" %in% names(per_class_metrics)) per_class_metrics$recall <- NA_real_
  per_class_metrics$recall <- validate_probability(per_class_metrics$recall, "class recall")
}

id_trim_counts <- c(
  candidate = sum(candidate$feature != trimws(candidate$feature)),
  eligibility = sum(eligibility$MAG != trimws(eligibility$MAG)),
  importance = sum(importance$feature != trimws(importance$feature)),
  stability = sum(stability$feature != trimws(stability$feature)),
  selection = sum(selection$feature != trimws(selection$feature)),
  taxonomy = sum(taxonomy$Genome != trimws(taxonomy$Genome))
)
candidate$feature <- trimws(candidate$feature)
eligibility$MAG <- trimws(eligibility$MAG)
importance$feature <- trimws(importance$feature)
stability$feature <- trimws(stability$feature)
selection$feature <- trimws(selection$feature)
taxonomy$Genome <- trimws(taxonomy$Genome)

require_unique_key(candidate, "feature", "candidate_MAG_summary.csv")
require_unique_key(eligibility, c("MAG", "fold_id"), "feature_eligibility_by_fold.csv")
require_unique_key(importance, c("feature", "fold_id"), "feature_importance_by_fold.csv")
require_unique_key(stability, "feature", "feature_importance_stability.csv")
require_unique_key(selection, "feature", "feature_selection_frequency.csv")

feature_ids <- sort(unique(stability$feature))
fold_ids <- unique(eligibility$fold_id)
heldout_subjects <- unique(eligibility$heldout_subject)
subject_number <- suppressWarnings(as.integer(sub("^[^0-9]*", "", heldout_subjects)))
if (anyNA(subject_number)) stop("Held-out subject labels cannot be naturally ordered.", call. = FALSE)
heldout_subjects <- heldout_subjects[order(subject_number, heldout_subjects)]
fold_ids <- unique(eligibility$fold_id[match(heldout_subjects, eligibility$heldout_subject)])

if (length(feature_ids) != 441L) stop("Expected 441 unique MAGs; found ", length(feature_ids), ".", call. = FALSE)
if (length(heldout_subjects) != 13L) stop("Expected 13 held-out subjects; found ", length(heldout_subjects), ".", call. = FALSE)
if (nrow(eligibility) != length(feature_ids) * length(heldout_subjects)) {
  stop("Eligibility table is not the complete 441 x 13 MAG-fold grid.", call. = FALSE)
}
assert_same_set(feature_ids, eligibility$MAG, "Stability", "eligibility")
assert_same_set(feature_ids, importance$feature, "Stability", "importance")
assert_same_set(feature_ids, selection$feature, "Stability", "selection-frequency")
if (!all(candidate$feature %in% feature_ids)) stop("Candidate summary contains an unknown MAG.", call. = FALSE)

resolved_config <- yaml::read_yaml(resolved_config_path)
if (!is.null(resolved_config$preprocessing$top_n)) {
  stop("Expected preprocessing.top_n to be disabled (null); found: ", resolved_config$preprocessing$top_n, call. = FALSE)
}
if (any(eligibility$selected & !eligibility$eligible)) {
  stop("At least one MAG was selected when it was not eligible.", call. = FALSE)
}
if (!all(eligibility$selected == eligibility$eligible)) {
  stop("Not every eligible MAG was included even though top_n is disabled.", call. = FALSE)
}

computed_frequency <- stats::aggregate(
  cbind(eligible = as.integer(eligibility$eligible), selected = as.integer(eligibility$selected)),
  by = list(feature = eligibility$MAG), FUN = sum
)
computed_frequency$selection_frequency_when_eligible <- with(
  computed_frequency, ifelse(eligible == 0L, NA_real_, selected / eligible)
)
frequency_check <- merge(
  selection, computed_frequency, by = "feature", sort = FALSE,
  suffixes = c("_stored", "_computed")
)
assert_close(frequency_check$eligible_folds, frequency_check$eligible, "Eligible-fold counts")
assert_close(frequency_check$selected_folds, frequency_check$selected, "Selected-fold counts")
assert_close(
  frequency_check$selection_frequency_when_eligible_stored,
  frequency_check$selection_frequency_when_eligible_computed,
  "Selection frequencies"
)

importance$recomputed_importance_rank <- ave(
  importance$MeanDecreaseAccuracy, importance$fold_id,
  FUN = function(x) rank(-x, ties.method = "average")
)
importance_split <- split(importance, importance$feature)
recomputed_stability <- do.call(rbind, lapply(names(importance_split), function(feature_id) {
  x <- importance_split[[feature_id]]
  data.frame(
    feature = feature_id,
    importance_folds_recomputed = length(unique(x$fold_id)),
    importance_q25 = unname(stats::quantile(x$MeanDecreaseAccuracy, 0.25, type = 7)),
    median_importance_recomputed = stats::median(x$MeanDecreaseAccuracy),
    importance_q75 = unname(stats::quantile(x$MeanDecreaseAccuracy, 0.75, type = 7)),
    importance_iqr_recomputed = stats::IQR(x$MeanDecreaseAccuracy, type = 7),
    median_rank_recomputed = stats::median(x$recomputed_importance_rank),
    rank_iqr_recomputed = stats::IQR(x$recomputed_importance_rank, type = 7),
    stringsAsFactors = FALSE
  )
}))
rownames(recomputed_stability) <- NULL
stability_check <- merge(stability, recomputed_stability, by = "feature", sort = FALSE)
assert_close(stability_check$importance_folds, stability_check$importance_folds_recomputed, "Importance-fold counts")
assert_close(stability_check$median_native_importance, stability_check$median_importance_recomputed, "Median importance")
assert_close(stability_check$native_importance_iqr, stability_check$importance_iqr_recomputed, "Importance IQR")
assert_close(stability_check$median_importance_rank, stability_check$median_rank_recomputed, "Median importance rank")
assert_close(stability_check$importance_rank_iqr, stability_check$rank_iqr_recomputed, "Importance-rank IQR")

expected_candidate_order <- stability[order(
  stability$median_importance_rank, -stability$median_native_importance, stability$feature
), "feature"]
if (!identical(candidate$feature, utils::head(expected_candidate_order, nrow(candidate)))) {
  stop("candidate_MAG_summary.csv order does not match median rank then median importance.", call. = FALSE)
}
if (args$n_candidates > nrow(candidate)) {
  stop("Requested ", args$n_candidates, " candidates, but candidate_MAG_summary.csv contains only ",
       nrow(candidate), ".", call. = FALSE)
}

eligible_keys <- paste(eligibility$MAG[eligibility$eligible], eligibility$fold_id[eligibility$eligible], sep = "\r")
importance_keys <- paste(importance$feature, importance$fold_id, sep = "\r")
missing_importance_keys <- setdiff(eligible_keys, importance_keys)
unexpected_importance_keys <- setdiff(importance_keys, eligible_keys)
if (length(missing_importance_keys) > 0L) {
  stop("Eligible MAG-fold combination(s) lack importance: ",
       paste(utils::head(missing_importance_keys, 10L), collapse = "; "), call. = FALSE)
}
if (length(unexpected_importance_keys) > 0L) {
  stop("Importance exists for ineligible MAG-fold combination(s): ",
       paste(utils::head(unexpected_importance_keys, 10L), collapse = "; "), call. = FALSE)
}

duplicate_genomes <- unique(taxonomy$Genome[duplicated(taxonomy$Genome) | duplicated(taxonomy$Genome, fromLast = TRUE)])
nonidentical_duplicate_genomes <- character()
if (length(duplicate_genomes) > 0L) {
  nonidentical_duplicate_genomes <- duplicate_genomes[vapply(duplicate_genomes, function(genome) {
    rows <- taxonomy[taxonomy$Genome == genome, , drop = FALSE]
    nrow(unique(rows)) > 1L
  }, logical(1))]
}
if (length(nonidentical_duplicate_genomes) > 0L) {
  stop(
    "GDTB.csv contains duplicate Genome keys with non-identical taxonomy records: ",
    paste(nonidentical_duplicate_genomes, collapse = "; "), call. = FALSE
  )
}
taxonomy_rows_input <- nrow(taxonomy)
taxonomy <- taxonomy[!duplicated(taxonomy$Genome), , drop = FALSE]
taxonomy_rows_deduplicated <- nrow(taxonomy)
require_unique_key(taxonomy, "Genome", "GDTB.csv after exact-row deduplication")

rank_source <- c(
  domain = "domain", phylum = "phylum", class = "class",
  order = if ("order" %in% names(taxonomy)) "order" else if ("ordeer" %in% names(taxonomy)) "ordeer" else NA_character_,
  family = "family", genus = "genus", species = "species"
)
missing_rank_sources <- rank_source[is.na(rank_source) | !(rank_source %in% names(taxonomy))]
if (length(missing_rank_sources) > 0L) {
  stop("GDTB.csv lacks one or more required taxonomy ranks; observed columns: ",
       paste(names(taxonomy), collapse = ", "), call. = FALSE)
}

taxonomy_join <- data.frame(
  Genome = taxonomy$Genome,
  taxonomy_label = if ("label" %in% names(taxonomy)) taxonomy$label else NA_character_,
  domain = taxonomy[[rank_source[["domain"]]]],
  phylum = taxonomy[[rank_source[["phylum"]]]],
  class = taxonomy[[rank_source[["class"]]]],
  order = taxonomy[[rank_source[["order"]]]],
  family = taxonomy[[rank_source[["family"]]]],
  genus = taxonomy[[rank_source[["genus"]]]],
  species = taxonomy[[rank_source[["species"]]]],
  stringsAsFactors = FALSE
)
rank_columns <- c("domain", "phylum", "class", "order", "family", "genus", "species")
taxonomy_join$taxonomy_complete <- apply(taxonomy_join[rank_columns], 1L, function(x) {
  x <- trimws(as.character(x))
  paste(x[!is.na(x) & nzchar(x)], collapse = ";")
})

clean_rank <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^[a-z]__", "", x)
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}
for (rank_name in rank_columns) {
  taxonomy_join[[paste0(rank_name, "_display")]] <- clean_rank(taxonomy_join[[rank_name]])
}
taxonomy_join$display_taxon <- vapply(seq_len(nrow(taxonomy_join)), function(i) {
  choices <- taxonomy_join[i, paste0(c("species", "genus", "family", "order"), "_display"), drop = TRUE]
  chosen <- choices[!is.na(choices) & nzchar(choices)]
  if (length(chosen)) chosen[[1L]] else taxonomy_join$Genome[[i]]
}, character(1))

all_feature_ids <- sort(unique(c(
  candidate$feature, eligibility$MAG, importance$feature, stability$feature, selection$feature
)))
matched_index <- match(all_feature_ids, taxonomy_join$Genome)
unmatched_ids <- all_feature_ids[is.na(matched_index)]
matched_ids <- all_feature_ids[!is.na(matched_index)]
match_percentage <- 100 * length(matched_ids) / length(all_feature_ids)

join_taxonomy <- function(data, id_column, source_name) {
  before <- nrow(data)
  joined <- dplyr::left_join(data, taxonomy_join, by = stats::setNames("Genome", id_column))
  if (nrow(joined) != before) stop(source_name, " row count changed during taxonomy join.", call. = FALSE)
  joined
}

candidate_rows_before <- nrow(candidate)
eligibility_rows_before <- nrow(eligibility)
importance_rows_before <- nrow(importance)
stability_rows_before <- nrow(stability)
selection_rows_before <- nrow(selection)
if ("taxonomy" %in% names(candidate)) names(candidate)[names(candidate) == "taxonomy"] <- "source_taxonomy"
candidate_joined <- join_taxonomy(candidate, "feature", "candidate_MAG_summary.csv")
eligibility_joined <- join_taxonomy(eligibility, "MAG", "feature_eligibility_by_fold.csv")
importance_joined <- join_taxonomy(importance, "feature", "feature_importance_by_fold.csv")
stability_joined <- join_taxonomy(stability, "feature", "feature_importance_stability.csv")
selection_joined <- join_taxonomy(selection, "feature", "feature_selection_frequency.csv")

row_counts_unchanged <- identical(c(
  nrow(candidate_joined), nrow(eligibility_joined), nrow(importance_joined),
  nrow(stability_joined), nrow(selection_joined)
), c(
  candidate_rows_before, eligibility_rows_before, importance_rows_before,
  stability_rows_before, selection_rows_before
))
if (!row_counts_unchanged) stop("At least one feature-table row count changed after taxonomy joining.", call. = FALSE)

audit <- data.frame(
  check = c(
    "taxonomy_file", "taxonomy_file_mtime_utc", "taxonomy_rows_input",
    "taxonomy_rows_after_identical_deduplication", "duplicate_genome_values",
    "identical_duplicate_rows_removed", "nonidentical_duplicate_genomes",
    "total_unique_feature_MAGs", "matched_unique_MAGs", "unmatched_unique_MAGs",
    "taxonomy_match_percentage", "candidate_rows_input", "candidate_rows_output",
    "eligibility_rows_input", "eligibility_rows_output", "importance_rows_input",
    "importance_rows_output", "stability_rows_input", "stability_rows_output",
    "selection_rows_input", "selection_rows_output", "trimmed_identifier_values",
    "taxonomy_join_passed"
  ),
  value = c(
    taxonomy_path, fmt_mtime(taxonomy_path), taxonomy_rows_input,
    taxonomy_rows_deduplicated, length(duplicate_genomes),
    taxonomy_rows_input - taxonomy_rows_deduplicated, length(nonidentical_duplicate_genomes),
    length(all_feature_ids), length(matched_ids), length(unmatched_ids),
    sprintf("%.6f", match_percentage), candidate_rows_before, nrow(candidate_joined),
    eligibility_rows_before, nrow(eligibility_joined), importance_rows_before,
    nrow(importance_joined), stability_rows_before, nrow(stability_joined),
    selection_rows_before, nrow(selection_joined), sum(id_trim_counts),
    length(unmatched_ids) == 0L && row_counts_unchanged && length(nonidentical_duplicate_genomes) == 0L
  ),
  detail = c(
    "Repository-located file; project interface intentionally uses filename GDTB.csv.",
    "UTC modification time", "Rows before duplicate-key assessment",
    "Exact duplicate rows may be safely removed", "Number of duplicated Genome keys",
    "Removed only when every taxonomy field was identical", "Critical many-to-many risk",
    "Union across all five feature tables", "Exact full-identifier matches",
    if (length(unmatched_ids)) paste(unmatched_ids, collapse = ";") else "None",
    "Matched / total unique feature MAGs x 100", "Before exact join", "After exact join",
    "Before exact join", "After exact join", "Before exact join", "After exact join",
    "Before exact join", "After exact join", "Before exact join", "After exact join",
    paste(names(id_trim_counts), id_trim_counts, sep = "=", collapse = "; "),
    "Requires complete matching, stable row counts, and no non-identical duplicate taxonomy keys"
  ),
  stringsAsFactors = FALSE
)
write_csv(audit, out_path("taxonomy_join_audit.csv"))
write_tsv(audit, out_path("taxonomy_join_audit.tsv"))
write_csv(data.frame(MAG = unmatched_ids, stringsAsFactors = FALSE), out_path("taxonomy_unmatched_MAGs.csv"))

audit_report <- c(
  "Random Forest MAG taxonomy-join audit",
  "=====================================",
  paste0("Taxonomy file: ", taxonomy_path),
  paste0("Taxonomy file modified (UTC): ", fmt_mtime(taxonomy_path)),
  paste0("Taxonomy input rows: ", taxonomy_rows_input),
  paste0("Unique Genome keys after safe exact-row deduplication: ", taxonomy_rows_deduplicated),
  paste0("Duplicated Genome keys: ", length(duplicate_genomes)),
  paste0("Non-identical duplicate Genome keys: ", length(nonidentical_duplicate_genomes)),
  paste0("Unique feature MAGs before join: ", length(all_feature_ids)),
  paste0("Matched unique MAGs: ", length(matched_ids)),
  paste0("Unmatched unique MAGs: ", length(unmatched_ids)),
  paste0("Match percentage: ", sprintf("%.2f%%", match_percentage)),
  paste0("All feature-table row counts unchanged: ", row_counts_unchanged),
  paste0("Identifiers trimmed for accidental edge whitespace: ", sum(id_trim_counts)),
  "Join keys: feature = Genome for candidate, importance, stability and selection tables; MAG = Genome for eligibility.",
  "The source taxonomy uses the project-established filename GDTB.csv and raw column spelling ordeer; joined outputs expose this rank as order.",
  if (length(unmatched_ids)) paste0("Unmatched identifiers: ", paste(unmatched_ids, collapse = "; ")) else "Unmatched identifiers: none",
  paste0("AUDIT STATUS: ", if (length(unmatched_ids) == 0L && row_counts_unchanged) "PASS" else "FAIL")
)
writeLines(audit_report, out_path("taxonomy_join_report.txt"), useBytes = TRUE)
if (length(unmatched_ids) > 0L) {
  stop("Taxonomy join is incomplete; see taxonomy_unmatched_MAGs.csv and taxonomy_join_report.txt.", call. = FALSE)
}

make_unique_suffixes <- function(ids, minimum = 10L) {
  width <- min(minimum, max(nchar(ids)))
  suffix <- substring(ids, pmax(1L, nchar(ids) - width + 1L))
  while (anyDuplicated(suffix) && width < max(nchar(ids))) {
    width <- width + 1L
    suffix <- substring(ids, pmax(1L, nchar(ids) - width + 1L))
  }
  if (anyDuplicated(suffix)) suffix <- ids
  suffix
}

candidate_order <- utils::head(candidate_joined$feature, args$n_candidates)

rank_cutpoints <- unname(stats::quantile(
  stability_joined$importance_rank_iqr, probs = c(1 / 3, 2 / 3), type = 7, na.rm = TRUE
))
stability_joined$stability_category <- cut(
  stability_joined$importance_rank_iqr,
  breaks = c(-Inf, rank_cutpoints, Inf),
  labels = c("Relatively stable", "Intermediate", "Unstable"),
  include.lowest = TRUE
)
stability_joined <- dplyr::left_join(
  stability_joined, recomputed_stability[c(
    "feature", "importance_q25", "median_importance_recomputed", "importance_q75"
  )], by = "feature"
)
stability_joined <- dplyr::left_join(
  stability_joined,
  selection[c(
    "feature", "eligible_folds", "selected_folds", "selection_frequency_when_eligible"
  )],
  by = "feature"
)
stability_joined$short_MAG_suffix <- make_unique_suffixes(stability_joined$feature)
if (anyNA(stability_joined$taxonomy_label) || any(!nzchar(stability_joined$taxonomy_label))) {
  stop("One or more matched MAGs lack the project taxonomy label.", call. = FALSE)
}
if (anyDuplicated(stability_joined$taxonomy_label)) {
  stop("Project taxonomy labels are not unique across the 441 feature MAGs.", call. = FALSE)
}
stability_joined$display_label <- paste0(
  stability_joined$display_taxon, "\n[", stability_joined$taxonomy_label, "]"
)
stability_joined$is_candidate <- stability_joined$feature %in% candidate_order
stability_joined$is_top_three <- stability_joined$feature %in% utils::head(candidate_order, 3L)
candidate_plot <- stability_joined[match(candidate_order, stability_joined$feature), , drop = FALSE]
if (anyDuplicated(candidate_plot$display_label)) stop("Candidate display labels are not unique.", call. = FALSE)
candidate_label_map <- stats::setNames(candidate_plot$display_label, candidate_plot$feature)

subjects_table <- unique(eligibility_joined[c("fold_id", "heldout_subject")])
subjects_table$subject_number <- suppressWarnings(as.integer(sub("^[^0-9]*", "", subjects_table$heldout_subject)))
subjects_table <- subjects_table[order(subjects_table$subject_number, subjects_table$heldout_subject), , drop = FALSE]
if (anyDuplicated(subjects_table$heldout_subject) || anyDuplicated(subjects_table$fold_id)) {
  stop("Fold IDs and held-out subjects are not one-to-one.", call. = FALSE)
}

panel_a_grid <- expand.grid(
  feature = candidate_order,
  heldout_subject = subjects_table$heldout_subject,
  stringsAsFactors = FALSE
)
panel_a_eligibility <- eligibility_joined[c("MAG", "heldout_subject", "fold_id", "eligible", "selected")]
names(panel_a_eligibility)[names(panel_a_eligibility) == "MAG"] <- "feature"
panel_a_grid <- merge(
  panel_a_grid, panel_a_eligibility,
  by = c("feature", "heldout_subject"), all.x = TRUE, sort = FALSE
)
panel_a_grid <- merge(
  panel_a_grid,
  importance_joined[c("feature", "fold_id", "MeanDecreaseAccuracy", "MeanDecreaseGini")],
  by = c("feature", "fold_id"), all.x = TRUE, sort = FALSE
)
panel_a_grid$display_label <- unname(candidate_label_map[panel_a_grid$feature])
if (anyNA(panel_a_grid$eligible)) stop("Panel A lacks an eligibility record.", call. = FALSE)
if (any(panel_a_grid$eligible & is.na(panel_a_grid$MeanDecreaseAccuracy))) {
  bad <- panel_a_grid[panel_a_grid$eligible & is.na(panel_a_grid$MeanDecreaseAccuracy), ]
  stop("Panel A eligible MAG-fold combination lacks importance: ",
       paste(paste(bad$feature, bad$fold_id, sep = "/"), collapse = "; "), call. = FALSE)
}
panel_a_grid$importance_for_plot <- ifelse(
  panel_a_grid$eligible, panel_a_grid$MeanDecreaseAccuracy, NA_real_
)
panel_a_grid$heldout_label <- paste0("Held out: ", panel_a_grid$heldout_subject)
panel_a_grid$heldout_label <- factor(
  panel_a_grid$heldout_label,
  levels = paste0("Held out: ", subjects_table$heldout_subject)
)
panel_a_grid$display_label <- factor(
  panel_a_grid$display_label,
  levels = rev(candidate_plot$display_label)
)

candidate_plot$display_label <- factor(
  candidate_plot$display_label,
  levels = rev(candidate_plot$display_label)
)

eligible_counts <- table(selection_joined$eligible_folds)
eligible_groups <- c("13", "12", "11", "<=10")
panel_d_counts <- data.frame(
  eligible_fold_group = factor(eligible_groups, levels = eligible_groups),
  MAG_count = c(
    sum(selection_joined$eligible_folds == 13L),
    sum(selection_joined$eligible_folds == 12L),
    sum(selection_joined$eligible_folds == 11L),
    sum(selection_joined$eligible_folds <= 10L)
  )
)
affected_features <- selection_joined$feature[selection_joined$eligible_folds < length(heldout_subjects)]
if (!length(affected_features)) stop("No fold-variable MAG exists for Panel D.", call. = FALSE)
affected_labels <- stability_joined$display_label[match(affected_features, stability_joined$feature)]
panel_d_tiles <- eligibility_joined[eligibility_joined$MAG %in% affected_features, , drop = FALSE]
panel_d_tiles$display_label <- factor(
  stability_joined$display_label[match(panel_d_tiles$MAG, stability_joined$feature)],
  levels = rev(affected_labels)
)
panel_d_tiles$heldout_label <- factor(
  paste0("Held out: ", panel_d_tiles$heldout_subject),
  levels = paste0("Held out: ", subjects_table$heldout_subject)
)
panel_d_tiles$eligibility_status <- factor(
  ifelse(panel_d_tiles$eligible, "Eligible", "Ineligible"),
  levels = c("Eligible", "Ineligible")
)

thresholds <- resolved_config$preprocessing
ineligible_rows <- panel_d_tiles[!panel_d_tiles$eligible, , drop = FALSE]
ineligibility_reason <- vapply(seq_len(nrow(ineligible_rows)), function(i) {
  row <- ineligible_rows[i, ]
  reasons <- character()
  if (is.finite(row$missing_fraction) && row$missing_fraction > thresholds$max_missing_fraction) {
    reasons <- c(reasons, paste0("missing ", scales::percent(row$missing_fraction, accuracy = 0.1)))
  }
  if (is.finite(row$sample_prevalence) && row$sample_prevalence < thresholds$sample_prevalence) {
    reasons <- c(reasons, paste0("sample prevalence ", scales::percent(row$sample_prevalence, accuracy = 0.1),
                                " < ", scales::percent(thresholds$sample_prevalence, accuracy = 0.1)))
  }
  if (is.finite(row$subject_prevalence) && row$subject_prevalence < thresholds$subject_prevalence_min) {
    reasons <- c(reasons, paste0("subject prevalence ", row$subject_prevalence,
                                " < ", thresholds$subject_prevalence_min))
  }
  if (is.finite(row$mean_abundance) && row$mean_abundance < thresholds$min_mean_abundance) {
    reasons <- c(reasons, "mean abundance below configured minimum")
  }
  if (!is.finite(row$variance) || row$variance <= 0) reasons <- c(reasons, "non-positive variance")
  paste0(row$heldout_subject, ": ", if (length(reasons)) paste(reasons, collapse = ", ") else "reason not represented by audited filter columns")
}, character(1))
eligibility_annotation <- paste(ineligibility_reason, collapse = "; ")

balanced_accuracy <- overall_metrics$balanced_accuracy[[1L]]
macro_f1 <- overall_metrics$macro_f1[[1L]]
zero_recall_classes <- as.character(per_class_metrics$class[is.finite(per_class_metrics$recall) & per_class_metrics$recall == 0])
format_metric <- function(x, digits = 6L) {
  if (length(x) == 1L && is.finite(x)) sprintf(paste0("%.", digits, "f"), x) else "unavailable"
}
recall_available <- sum(is.finite(per_class_metrics$recall))
recall_missing <- nrow(per_class_metrics) - recall_available
recall_text <- if (!nrow(per_class_metrics)) {
  "Class recall unavailable (no identified class summaries supplied)."
} else if (!recall_available) {
  paste0("Class recall unavailable for all ", nrow(per_class_metrics), " supplied classes.")
} else {
  paste0(length(zero_recall_classes), "/", recall_available,
         " classes with available recall have zero recall",
         if (length(zero_recall_classes)) paste0(" (", paste(zero_recall_classes, collapse = ", "), ")") else "",
         if (recall_missing) paste0("; recall unavailable for ", recall_missing, " additional classes") else "", ".")
}
performance_context <- paste0("Pooled LOSO balanced accuracy = ", format_metric(balanced_accuracy),
                              "; macro F1 = ", format_metric(macro_f1), "; ", recall_text)

palette <- list(
  ink = "#252525", grey = "#D9D9D9", grid = "#E7E5E1", offwhite = "#FFF8E7",
  blue = "#2166AC", blue_dark = "#174A7E", orange = "#B35806",
  stable = "#0072B2", intermediate = "#E69F00", unstable = "#767676",
  paper = "#FCFBF8"
)
base_theme <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title = ggplot2::element_text(face = "bold", colour = palette$ink, size = base_size + 0.8),
      plot.subtitle = ggplot2::element_text(colour = "#555555", size = base_size - 0.5, lineheight = 1.05),
      plot.caption = ggplot2::element_text(colour = "#555555", size = base_size - 1.2, hjust = 0),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 3, colour = palette$ink),
      plot.tag.position = c(0.01, 0.99),
      axis.title = ggplot2::element_text(colour = palette$ink, size = base_size - 0.2),
      axis.text = ggplot2::element_text(colour = palette$ink, size = base_size - 1),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 0.4),
      legend.text = ggplot2::element_text(size = base_size - 0.8),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(8, 12, 8, 10)
    )
}
wrap_text <- function(x, width) paste(base::strwrap(x, width = width), collapse = "\n")

panel_a <- ggplot2::ggplot(
  panel_a_grid,
  ggplot2::aes(x = heldout_label, y = display_label, fill = importance_for_plot)
) +
  ggplot2::geom_tile(width = 0.94, height = 0.92) +
  ggplot2::scale_fill_gradient2(
    low = palette$blue, mid = palette$offwhite, high = palette$orange,
    midpoint = 0, na.value = palette$grey,
    name = "Mean decrease\nin accuracy",
    guide = ggplot2::guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(38, "mm"), barheight = grid::unit(3, "mm")
    )
  ) +
  ggplot2::labs(
    title = "Fold-specific native importance",
    subtitle = paste0(
      args$n_candidates, " highest-ranked exploratory MAGs\nGrey = ineligible",
      if (!any(!panel_a_grid$eligible)) " (none among displayed MAGs)" else ""
    ),
    x = NULL, y = NULL
  ) +
  base_theme(8.3) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 50, hjust = 1, vjust = 1, size = 6.6),
    axis.text.y = ggplot2::element_text(size = 6.3, lineheight = 0.95),
    axis.ticks = ggplot2::element_blank(),
    legend.margin = ggplot2::margin(t = -2),
    plot.margin = ggplot2::margin(10, 8, 8, 12)
  ) +
  ggplot2::coord_cartesian(clip = "off")

stability_colours <- c(
  "Relatively stable" = palette$stable,
  "Intermediate" = palette$intermediate,
  "Unstable" = palette$unstable
)
stability_shapes <- c("Relatively stable" = 21, "Intermediate" = 22, "Unstable" = 24)
panel_b <- ggplot2::ggplot(candidate_plot, ggplot2::aes(y = display_label)) +
  ggplot2::geom_vline(xintercept = 0, colour = "#777777", linetype = "22", linewidth = 0.4) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = importance_q25, xend = importance_q75,
      y = display_label, yend = display_label
    ),
    linewidth = 0.7, colour = "#5F5F5F"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = median_native_importance, fill = stability_category, shape = stability_category),
    size = 2.5, stroke = 0.45, colour = palette$ink
  ) +
  ggplot2::geom_point(
    data = candidate_plot[candidate_plot$is_top_three, , drop = FALSE],
    ggplot2::aes(x = median_native_importance, y = display_label),
    shape = 21, size = 4.0, stroke = 0.8, fill = NA, colour = palette$ink,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_fill_manual(values = stability_colours, drop = FALSE, name = "Rank-IQR tertile") +
  ggplot2::scale_shape_manual(values = stability_shapes, drop = FALSE, name = "Rank-IQR tertile") +
  ggplot2::labs(
    title = "Typical importance and fold variability",
    subtitle = paste0(
      "Points: medians; intervals: fold IQR.\nRank-IQR tertiles: stable <=",
      scales::number(rank_cutpoints[[1]], accuracy = 0.1), "; intermediate <=",
      scales::number(rank_cutpoints[[2]], accuracy = 0.1), "; unstable >",
      scales::number(rank_cutpoints[[2]], accuracy = 0.1),
      ".\nOutlined symbols identify the top three."
    ),
    x = "MeanDecreaseAccuracy across eligible folds", y = NULL
  ) +
  base_theme(8.3) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE, title.position = "top"),
    shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE, title.position = "top")
  ) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 6.3, lineheight = 0.95),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = ggplot2::margin(10, 12, 8, 8)
  ) +
  ggplot2::coord_cartesian(clip = "off")

noncandidate <- stability_joined[!stability_joined$is_candidate, , drop = FALSE]
candidate_landscape <- stability_joined[stability_joined$is_candidate, , drop = FALSE]
top_three <- stability_joined[stability_joined$is_top_three, , drop = FALSE]
top_three <- top_three[match(utils::head(candidate_order, 3L), top_three$feature), , drop = FALSE]
top_three$label_x <- 42
top_three$label_y <- c(92, 58, 26)
panel_c <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = noncandidate,
    ggplot2::aes(x = median_importance_rank, y = importance_rank_iqr),
    colour = "#AFAFAF", alpha = 0.36, size = 1.25
  ) +
  ggplot2::geom_point(
    data = candidate_landscape,
    ggplot2::aes(
      x = median_importance_rank, y = importance_rank_iqr,
      fill = median_native_importance, size = eligible_folds
    ),
    shape = 21, colour = palette$ink, stroke = 0.35
  ) +
  ggplot2::geom_segment(
    data = top_three,
    ggplot2::aes(
      x = median_importance_rank, y = importance_rank_iqr,
      xend = label_x - 1.5, yend = label_y
    ),
    linewidth = 0.3, colour = "#777777",
    show.legend = FALSE
  ) +
  ggplot2::geom_label(
    data = top_three,
    ggplot2::aes(x = label_x, y = label_y, label = display_label),
    size = 2.1, lineheight = 0.9, colour = palette$ink,
    fill = scales::alpha(palette$paper, 0.94),
    label.padding = grid::unit(0.12, "lines"), hjust = 0,
    show.legend = FALSE
  ) +
  ggplot2::annotate(
    "text", x = 110, y = 0,
    label = "More consistently high-ranking", hjust = 1, vjust = -0.25,
    size = 2.35, colour = palette$blue_dark
  ) +
  ggplot2::scale_fill_viridis_c(
    option = "C", direction = -1, name = "Median native\nimportance"
  ) +
  ggplot2::scale_size_continuous(
    range = c(2.0, 3.2), breaks = sort(unique(candidate_landscape$eligible_folds)),
    name = "Eligible folds"
  ) +
  ggplot2::labs(
    title = "Rank-stability landscape",
    subtitle = paste0("All ", length(feature_ids), " MAGs shown; coloured points are the ", args$n_candidates, " exploratory candidates."),
    x = "Median importance rank (better to the left)",
    y = "Importance-rank IQR (more stable below)"
  ) +
  base_theme(8.3) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(30, "mm"), barheight = grid::unit(3, "mm")
    ),
    size = ggplot2::guide_legend(title.position = "top")
  ) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(colour = palette$grid, linewidth = 0.25),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = ggplot2::margin(10, 18, 10, 12)
  ) +
  ggplot2::coord_cartesian(clip = "off")

panel_d_count_plot <- ggplot2::ggplot(
  panel_d_counts,
  ggplot2::aes(x = eligible_fold_group, y = MAG_count)
) +
  ggplot2::geom_col(fill = palette$blue, width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = MAG_count), vjust = -0.35,
    size = 2.7, colour = palette$ink
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.14)),
    breaks = scales::breaks_pretty(n = 3)
  ) +
  ggplot2::labs(
    title = "Eligibility across folds",
    subtitle = paste0(sum(selection_joined$eligible_folds == 13L), "/", length(feature_ids),
                      " MAGs eligible in all 13 folds"),
    x = "Eligible folds", y = "MAG count"
  ) +
  base_theme(8.0) +
  ggplot2::theme(
    legend.position = "none",
    plot.margin = ggplot2::margin(8, 10, 6, 8)
  )

panel_d_tile_plot <- ggplot2::ggplot(
  panel_d_tiles,
  ggplot2::aes(x = heldout_label, y = display_label, fill = eligibility_status)
) +
  ggplot2::geom_tile(width = 0.94, height = 0.7) +
  ggplot2::scale_fill_manual(
    values = c("Eligible" = palette$blue, "Ineligible" = palette$grey),
    drop = FALSE, name = "Filter status",
    guide = ggplot2::guide_legend(title.position = "top")
  ) +
  ggplot2::labs(
    title = "Fold-specific exclusions",
    subtitle = wrap_text(eligibility_annotation, 52L),
    x = NULL, y = NULL,
    caption = "Eligibility controls model entry after training-fold filtering; it is not feature importance."
  ) +
  base_theme(7.7) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 50, hjust = 1, vjust = 1, size = 6.1),
    axis.text.y = ggplot2::element_text(size = 6.1),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(5, 10, 8, 8)
  ) +
  ggplot2::coord_cartesian(clip = "off")

panel_d <- gridExtra::arrangeGrob(
  panel_d_count_plot, panel_d_tile_plot,
  ncol = 1, heights = c(0.48, 0.52)
)

panel_a_tagged <- panel_a + ggplot2::labs(tag = "A")
panel_b_tagged <- panel_b + ggplot2::labs(tag = "B")
panel_c_tagged <- panel_c + ggplot2::labs(tag = "C")
panel_d_tagged <- gridExtra::arrangeGrob(
  panel_d,
  top = grid::textGrob(
    "D", x = 0.01, hjust = 0,
    gp = grid::gpar(fontfamily = "sans", fontsize = 12, fontface = "bold", col = palette$ink)
  )
)

layout_matrix <- rbind(
  c(1, 2, 2, 2, 3, 3, 3, 6),
  c(1, 2, 2, 2, 3, 3, 3, 6),
  c(1, 4, 4, 4, 5, 5, 5, 6)
)
combined_body <- gridExtra::arrangeGrob(
  grobs = list(
    grid::nullGrob(), panel_a_tagged, panel_b_tagged,
    panel_c_tagged, panel_d_tagged, grid::nullGrob()
  ),
  layout_matrix = layout_matrix,
  widths = grid::unit.c(
    grid::unit(0.35, "in"), rep(grid::unit(1, "null"), 6L),
    grid::unit(0.08, "in")
  ),
  heights = c(1, 1, 1.15)
)
performance_note <- paste0(
  "Performance context: ", performance_context,
  " Native importance is an exploratory algorithm-specific ranking, not an effect size or significance test."
)
performance_note <- wrap_text(performance_note, 180L)
combined_figure <- gridExtra::arrangeGrob(
  combined_body,
  grid::textGrob(
    performance_note, x = 0.01, hjust = 0,
    gp = grid::gpar(fontfamily = "sans", fontsize = 7.2, col = "#4F4F4F")
  ),
  ncol = 1, heights = c(1, 0.075), padding = grid::unit(0.04, "in")
)

draw_grob <- function(grob) {
  grid::grid.newpage()
  grid::grid.draw(grob)
}
save_grob_bundle <- function(grob, stem, width, height, png_dpi = 400L, include_tiff = FALSE) {
  raster_override <- Sys.getenv("PUBLIC_RF_RASTER_DPI", unset = "")
  tiff_dpi <- 600L
  if (nzchar(raster_override)) {
    override <- suppressWarnings(as.integer(raster_override))
    if (!is.finite(override) || override <= 0L) stop("PUBLIC_RF_RASTER_DPI must be a positive integer.", call. = FALSE)
    png_dpi <- override
    tiff_dpi <- override
  }
  pdf_path <- paste0(stem, ".pdf")
  png_path <- paste0(stem, ".png")
  svg_path <- paste0(stem, ".svg")
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = "sans", onefile = TRUE)
  draw_grob(grob)
  grDevices::dev.off()
  ragg::agg_png(png_path, width = width, height = height, units = "in", res = png_dpi, background = "white")
  draw_grob(grob)
  grDevices::dev.off()
  svglite::svglite(svg_path, width = width, height = height, bg = "white")
  draw_grob(grob)
  grDevices::dev.off()
  cat("Exported ", basename(stem), ": ", width, " x ", height,
      " in; PNG ", png_dpi, " dpi; PDF/SVG vector\n", sep = "")
  if (include_tiff) {
    tiff_path <- paste0(stem, ".tiff")
    ragg::agg_tiff(
      tiff_path, width = width, height = height, units = "in", res = tiff_dpi,
      compression = "lzw", background = "white"
    )
    draw_grob(grob)
    grDevices::dev.off()
    cat("Exported ", basename(tiff_path), ": ", width, " x ", height,
        " in; ", tiff_dpi, " dpi; LZW compression\n", sep = "")
  }
}

save_grob_bundle(panel_a_tagged, out_path("panel_A_fold_importance_heatmap"), 8.2, 7.8)
save_grob_bundle(panel_b_tagged, out_path("panel_B_median_importance_iqr"), 8.2, 7.8)
save_grob_bundle(panel_c_tagged, out_path("panel_C_rank_stability_landscape"), 8.2, 5.8)
save_grob_bundle(panel_d_tagged, out_path("panel_D_feature_eligibility_qc"), 8.2, 5.8)
save_grob_bundle(combined_figure, out_path("figure_RF_MAG_feature_stability_ABCD"), 12.0, 11.2, include_tiff = TRUE)

candidate_output <- candidate_joined
candidate_output <- dplyr::left_join(
  candidate_output,
  stability_joined[c(
    "feature", "importance_q25", "importance_q75", "stability_category",
    "short_MAG_suffix", "display_label"
  )],
  by = "feature"
)
write_csv(candidate_output, out_path("candidate_MAG_summary_with_taxonomy.csv"))
write_csv(importance_joined, out_path("feature_importance_by_fold_with_taxonomy.csv"))
write_csv(eligibility_joined, out_path("feature_eligibility_by_fold_with_taxonomy.csv"))
write_csv(stability_joined, out_path("feature_importance_stability_with_taxonomy.csv"))
write_csv(selection_joined, out_path("feature_selection_frequency_with_taxonomy.csv"))
write_csv(panel_a_grid, out_path("panel_A_plot_data.csv"))
write_csv(candidate_plot, out_path("panel_B_plot_data.csv"))
write_csv(stability_joined, out_path("panel_C_plot_data.csv"))
write_csv(panel_d_counts, out_path("panel_D_eligibility_counts.csv"))
write_csv(panel_d_tiles, out_path("panel_D_eligibility_tiles.csv"))

caption <- paste0(
  "# Figure | Random Forest MAG eligibility, importance and stability\n\n",
  paste0("**Caption.** Fold-specific eligibility and native Random Forest importance of metagenome-assembled genomes across ", length(heldout_subjects), " leave-one-subject-out models. "),
  "**(A)** `MeanDecreaseAccuracy` for the ", args$n_candidates,
  paste0(" highest-ranked exploratory MAGs in each subject-held-out fold. Grey cells indicate fold-specific ineligibility; ", sum(!panel_a_grid$eligible), " displayed candidate-fold cells were ineligible. "),
  "**(B)** Median importance and interquartile range across eligible folds. Point colours and shapes divide the observed all-MAG rank-IQR distribution into tertiles (relatively stable ≤ ",
  scales::number(rank_cutpoints[[1]], accuracy = 0.1), ", intermediate > ",
  scales::number(rank_cutpoints[[1]], accuracy = 0.1), " to ≤ ",
  scales::number(rank_cutpoints[[2]], accuracy = 0.1), ", unstable > ",
  scales::number(rank_cutpoints[[2]], accuracy = 0.1), "); continuous rank IQR is retained in the source data. ",
  paste0("**(C)** Median importance rank versus rank interquartile range; lower values on both axes indicate more consistently high-ranking features. All ", length(feature_ids), " MAGs are shown, with the exploratory candidate set emphasized. "),
  "**(D)** Fold-specific eligibility summary. In total, ", sum(selection_joined$eligible_folds == 13L),
  "/", length(feature_ids), " MAGs were eligible in all 13 folds; ", eligibility_annotation, ". ",
  "Importance is an algorithm-specific predictive ranking measure and is not an effect size, abundance change, statistical-significance test or evidence of causation. ",
  "Pooled LOSO balanced accuracy = ", format_metric(balanced_accuracy, 3L),
  "; macro F1 = ", format_metric(macro_f1, 3L), "; ", recall_text,
  " These metrics are reported without a qualitative performance category. The displayed MAGs remain exploratory predictive candidates.\n"
)
writeLines(caption, out_path("figure_RF_MAG_feature_stability_caption.md"), useBytes = TRUE)

methods_note <- c(
  "# Methods note: Random Forest MAG feature stability",
  "",
  paste0("The supplied subject-validation outputs represent ", length(heldout_subjects), " leave-one-subject-out (LOSO) Random Forest models. Each model held out one participant and estimated eligibility using only the remaining participants. Missingness, sample prevalence, subject prevalence, mean abundance and non-zero finite variance were evaluated independently within each training fold."),
  "",
  paste0("All ", length(feature_ids), " MAGs were considered. Because `preprocessing.top_n` was null, every eligible MAG entered its fold-specific model; selection frequency therefore measures inclusion after filtering and is not interpreted as feature importance."),
  "",
  paste0("Taxonomy was joined by the complete unchanged identifier using `feature = GDTB.csv$Genome` (or `MAG = GDTB.csv$Genome` in the eligibility table). The exact join matched ", length(matched_ids), "/", length(all_feature_ids), " MAGs (", sprintf("%.2f%%", match_percentage), "). No taxonomy was inferred from MAG names. The project intentionally uses the filename `GDTB.csv`; its raw `ordeer` field was retained semantically as the order rank in joined outputs."),
  "",
  "`MeanDecreaseAccuracy` and `MeanDecreaseGini` are native Random Forest importance measures. They are algorithm-specific predictive rankings, can be affected by feature prevalence, variance and correlation, and are not biological effect estimates or significance tests.",
  "",
  "For every MAG, median importance, the 25th and 75th percentiles and IQR were calculated from eligible-fold `MeanDecreaseAccuracy` values. Importance ranks were recomputed within each fold, and their medians and IQRs were checked against the stored stability table. Rank-IQR display categories are descriptive tertiles of the observed 441-MAG distribution, not biological or inferential thresholds.",
  "",
  paste0(performance_context, " No qualitative discrimination category is assigned. Missing summary values are not inferred from other files. Native importance remains an exploratory predictive measure and does not establish statistical significance, biological effects or validated biomarkers.")
)
writeLines(methods_note, out_path("figure_RF_MAG_feature_stability_methods.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = out_path("figure_generation_session_info.txt"))

expected_figure_files <- c(
  paste0(prefix, "panel_A_fold_importance_heatmap.", c("pdf", "png", "svg")),
  paste0(prefix, "panel_B_median_importance_iqr.", c("pdf", "png", "svg")),
  paste0(prefix, "panel_C_rank_stability_landscape.", c("pdf", "png", "svg")),
  paste0(prefix, "panel_D_feature_eligibility_qc.", c("pdf", "png", "svg")),
  paste0(prefix, "figure_RF_MAG_feature_stability_ABCD.", c("pdf", "png", "svg", "tiff"))
)
expected_data_files <- paste0(prefix, c(
  "candidate_MAG_summary_with_taxonomy.csv",
  "feature_importance_by_fold_with_taxonomy.csv",
  "feature_eligibility_by_fold_with_taxonomy.csv",
  "feature_importance_stability_with_taxonomy.csv",
  "feature_selection_frequency_with_taxonomy.csv",
  "taxonomy_join_audit.csv", "taxonomy_join_audit.tsv", "taxonomy_join_report.txt",
  "taxonomy_unmatched_MAGs.csv", "panel_A_plot_data.csv", "panel_B_plot_data.csv",
  "panel_C_plot_data.csv", "panel_D_eligibility_counts.csv", "panel_D_eligibility_tiles.csv",
  "figure_RF_MAG_feature_stability_caption.md",
  "figure_RF_MAG_feature_stability_methods.md", "figure_generation_session_info.txt"
))
expected_paths <- file.path(output_dir, c(expected_figure_files, expected_data_files))

file_signature_ok <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "pdf") return(identical(readBin(path, "raw", n = 4L), charToRaw("%PDF")))
  if (ext == "png") {
    expected <- as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
    return(identical(readBin(path, "raw", n = 8L), expected))
  }
  if (ext %in% c("tif", "tiff")) {
    actual <- readBin(path, "raw", n = 4L)
    return(identical(actual, as.raw(c(0x49, 0x49, 0x2a, 0x00))) ||
           identical(actual, as.raw(c(0x4d, 0x4d, 0x00, 0x2a))))
  }
  if (ext == "svg") {
    header <- paste(readLines(path, n = 10L, warn = FALSE), collapse = " ")
    return(grepl("<svg", header, fixed = TRUE))
  }
  TRUE
}

exists_check <- file.exists(expected_paths)
nonempty_check <- exists_check & file.info(expected_paths)$size > 0
signature_check <- vapply(expected_paths, function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) return(FALSE)
  file_signature_ok(path)
}, logical(1))
validation <- data.frame(
  check = c(
    "expected_unique_MAGs", "expected_LOSO_folds", "complete_eligibility_grid",
    "expected_eligible_importance_rows", "selection_equals_eligibility_top_n_disabled",
    "candidate_order_verified", "eligible_candidate_fold_importance_complete",
    "negative_importance_retained", "taxonomy_join_complete", "taxonomy_join_row_counts_unchanged",
    "expected_outputs_exist", "expected_outputs_nonempty", "file_signatures_valid"
  ),
  passed = c(
    length(feature_ids) == 441L, length(heldout_subjects) == 13L,
    nrow(eligibility) == 5733L,
    nrow(importance) == sum(eligibility$eligible),
    all(eligibility$selected == eligibility$eligible) && is.null(resolved_config$preprocessing$top_n),
    identical(candidate$feature, utils::head(expected_candidate_order, nrow(candidate))),
    !any(panel_a_grid$eligible & is.na(panel_a_grid$MeanDecreaseAccuracy)),
    sum(importance_joined$MeanDecreaseAccuracy < 0, na.rm = TRUE) == sum(importance$MeanDecreaseAccuracy < 0, na.rm = TRUE),
    length(unmatched_ids) == 0L, row_counts_unchanged,
    all(exists_check), all(nonempty_check), all(signature_check)
  ),
  detail = c(
    paste(length(feature_ids), "unique MAGs"),
    paste(length(heldout_subjects), "held-out subjects"),
    paste(nrow(eligibility), "rows; expected", length(feature_ids) * length(heldout_subjects)),
    paste(nrow(importance), "importance rows for", sum(eligibility$eligible), "eligible combinations"),
    "Every eligible MAG entered its model; inclusion is not importance.",
    paste(nrow(candidate), "stored candidates ordered by median rank then median importance"),
    paste(sum(panel_a_grid$eligible), "eligible displayed combinations have importance"),
    paste(sum(importance$MeanDecreaseAccuracy < 0, na.rm = TRUE), "negative values retained"),
    paste(length(matched_ids), "of", length(all_feature_ids), "unique MAGs matched exactly"),
    "All five joined tables retain their input row counts.",
    paste(sum(exists_check), "of", length(expected_paths), "expected outputs exist"),
    paste(sum(nonempty_check), "of", length(expected_paths), "expected outputs are non-empty"),
    paste(sum(signature_check), "of", length(expected_paths), "expected outputs have valid format signatures")
  ),
  stringsAsFactors = FALSE
)
write_csv(validation, out_path("figure_generation_validation.csv"))
if (!all(validation$passed)) {
  failed <- validation$check[!validation$passed]
  stop("Figure-generation validation failed: ", paste(failed, collapse = ", "), call. = FALSE)
}

manifest <- data.frame(
  file = basename(c(expected_paths, out_path("figure_generation_validation.csv"))),
  bytes = file.info(c(expected_paths, out_path("figure_generation_validation.csv")))$size,
  stringsAsFactors = FALSE
)
write_csv(manifest, out_path("figure_output_manifest.csv"))

cat("\nVerified source counts:\n")
cat("  Unique MAGs: ", length(feature_ids), "\n", sep = "")
cat("  LOSO folds: ", length(heldout_subjects), "\n", sep = "")
cat("  Eligibility rows: ", nrow(eligibility), "\n", sep = "")
cat("  Eligible importance rows: ", nrow(importance), "\n", sep = "")
cat("  MAGs eligible in all 13 folds: ", sum(selection$eligible_folds == 13L), "\n", sep = "")
cat("  MAGs eligible in 11 folds: ", sum(selection$eligible_folds == 11L), "\n", sep = "")
cat("  Negative MeanDecreaseAccuracy values retained: ",
    sum(importance$MeanDecreaseAccuracy < 0, na.rm = TRUE), "\n", sep = "")
cat("  Taxonomy matches: ", length(matched_ids), "/", length(all_feature_ids),
    " (", sprintf("%.2f%%", match_percentage), ")\n", sep = "")
cat("  Balanced accuracy: ", format_metric(balanced_accuracy), "\n", sep = "")
cat("  Macro F1: ", format_metric(macro_f1), "\n", sep = "")
cat("  Zero-recall classes: ", paste(zero_recall_classes, collapse = ", "), "\n", sep = "")
cat("\nAll figure-generation validation checks passed.\n")
cat("Output manifest: ", out_path("figure_output_manifest.csv"), "\n", sep = "")
