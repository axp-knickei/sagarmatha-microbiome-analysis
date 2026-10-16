#!/usr/bin/env Rscript
source("analysis/functional_community_dynamics/R/common.R")

args <- parse_dynamics_args()
cfg <- read_dynamics_config(args$config)
assert_dynamics_packages()
paths <- dynamics_paths(cfg)
marker <- file.path(paths$qc, "01_prepare_profiles.complete")
guard_dynamics_stage(marker, args, cfg)
if (!file.exists(file.path(paths$qc, "00_validate_inputs.complete"))) stop("Stage 00 must complete first.", call. = FALSE)
inputs <- resolve_dynamics_inputs(cfg)
metadata_all <- prepare_dynamics_metadata(inputs$metadata_file, cfg)
abundance <- fread(inputs$abundance_file)
sample_col <- cfg$columns$sample_id
if (anyDuplicated(abundance[[sample_col]])) stop("Abundance sample identifiers are not unique.", call. = FALSE)
metadata <- metadata_all[sample_id %in% abundance[[sample_col]]]
setorder(metadata, subject_id, sample_date, sample_id)
sample_ids <- metadata$sample_id
expected <- as.integer(cfg$analysis$expected_sample_count)
if (isTRUE(cfg$analysis$enforce_expected_sample_count) && length(sample_ids) != expected) stop("Expected ", expected, " abundance-matched T1-T10 samples but found ", length(sample_ids), ".", call. = FALSE)
missing_abundance <- setdiff(sample_ids, abundance[[sample_col]])
if (length(missing_abundance)) stop("Included samples missing from abundance input: ", paste(missing_abundance, collapse = ", "), call. = FALSE)
abundance <- abundance[match(sample_ids, abundance[[sample_col]])]
mag_cols <- setdiff(names(abundance), sample_col)
taxonomic <- as.matrix(abundance[, ..mag_cols])
storage.mode(taxonomic) <- "double"
rownames(taxonomic) <- sample_ids
if (any(!is.finite(taxonomic) | taxonomic < 0)) stop("Taxonomic MAG matrix contains non-finite or negative values.", call. = FALSE)

layers <- configured_layers(cfg)
matrices <- list(taxonomic_MAG = taxonomic)
audit <- list()
audit[[1L]] <- data.table(
  layer = "taxonomic_MAG", role = "comparator", samples = nrow(taxonomic), subjects = uniqueN(metadata$subject_id),
  metadata_T1_T10_rows_before_abundance_intersection = nrow(metadata_all), metadata_rows_excluded_without_abundance = nrow(metadata_all) - nrow(metadata),
  timepoints = uniqueN(metadata$timepoint), raw_features = ncol(taxonomic), retained_features = ncol(taxonomic),
  filtered_features = 0L, zero_fraction = mean(taxonomic == 0), prevalence_threshold = NA_real_,
  row_sum_min = min(rowSums(taxonomic)), row_sum_median = median(rowSums(taxonomic)), row_sum_max = max(rowSums(taxonomic)),
  values_finite_nonnegative = TRUE, sample_order_identical = identical(rownames(taxonomic), sample_ids), metadata_identical = TRUE,
  source_checksum_md5 = checksum_file(inputs$abundance_file), retained_matrix_checksum_md5 = checksum_object(taxonomic),
  aitchison_zero_replacement_passed = all(is.finite(multiplicative_replace(taxonomic))) && all(multiplicative_replace(taxonomic) > 0),
  zero_replacement_fraction = 0.65,
  transformation_contract = "sample closure before Hellinger or Bray-Curtis; multiplicative zero replacement before CLR/Aitchison"
)

for (i in seq_len(nrow(layers))) {
  x <- read_dynamics_table(file.path(inputs$mag_functional_results_dir, layers$file[[i]]))
  x[, sampleID := as.character(sampleID)]
  x <- x[sampleID %in% sample_ids]
  if (!setequal(unique(x$sampleID), sample_ids)) stop("Functional layer ", layers$layer[[i]], " does not contain exactly the included samples.", call. = FALSE)
  if (any(!is.finite(x$potential) | x$potential < 0)) stop("Functional layer ", layers$layer[[i]], " contains non-finite or negative potential.", call. = FALSE)
  if (anyDuplicated(x[, .(sampleID, term)])) stop("Duplicate sample-term rows in layer ", layers$layer[[i]], ".", call. = FALSE)
  absent_metadata <- setdiff(c("Subject", "Timepoint"), names(x))
  if (length(absent_metadata)) stop("Functional layer ", layers$layer[[i]], " lacks metadata columns: ", paste(absent_metadata, collapse = ", "), call. = FALSE)
  observed_metadata <- unique(x[, .(sampleID, Subject = as.character(Subject), Timepoint = as.character(Timepoint))])
  expected_metadata <- metadata[, .(sampleID = sample_id, Subject = subject_id, Timepoint = as.character(timepoint))]
  setorder(observed_metadata, sampleID)
  setorder(expected_metadata, sampleID)
  metadata_identical <- identical(observed_metadata, expected_metadata)
  if (!metadata_identical) stop("Functional layer ", layers$layer[[i]], " has Subject/Timepoint metadata inconsistent with sample_info.csv.", call. = FALSE)
  wide <- dcast(x, sampleID ~ term, value.var = "potential", fill = 0)
  wide <- wide[match(sample_ids, sampleID)]
  raw_features <- ncol(wide) - 1L
  m <- as.matrix(wide[, -"sampleID"])
  storage.mode(m) <- "double"
  rownames(m) <- wide$sampleID
  prevalence <- colMeans(m > 0)
  keep <- prevalence >= as.numeric(cfg$analysis$minimum_feature_prevalence) & apply(m, 2L, stats::var) > 0
  m <- m[, keep, drop = FALSE]
  if (ncol(m) < as.integer(cfg$analysis$minimum_features_per_layer)) stop("Too few retained features in layer ", layers$layer[[i]], ".", call. = FALSE)
  if (any(rowSums(m) <= 0)) stop("Layer ", layers$layer[[i]], " has an all-zero retained profile.", call. = FALSE)
  matrices[[layers$layer[[i]]]] <- m
  audit[[length(audit) + 1L]] <- data.table(
    layer = layers$layer[[i]], role = layers$role[[i]], samples = nrow(m), subjects = uniqueN(metadata$subject_id),
    metadata_T1_T10_rows_before_abundance_intersection = nrow(metadata_all), metadata_rows_excluded_without_abundance = nrow(metadata_all) - nrow(metadata),
    timepoints = uniqueN(metadata$timepoint), raw_features = raw_features, retained_features = ncol(m),
    filtered_features = raw_features - ncol(m), zero_fraction = mean(m == 0),
    prevalence_threshold = as.numeric(cfg$analysis$minimum_feature_prevalence),
    row_sum_min = min(rowSums(m)), row_sum_median = median(rowSums(m)), row_sum_max = max(rowSums(m)),
    values_finite_nonnegative = TRUE, sample_order_identical = identical(rownames(m), sample_ids), metadata_identical = metadata_identical,
    source_checksum_md5 = checksum_file(file.path(inputs$mag_functional_results_dir, layers$file[[i]])), retained_matrix_checksum_md5 = checksum_object(m),
    aitchison_zero_replacement_passed = all(is.finite(multiplicative_replace(m))) && all(multiplicative_replace(m) > 0),
    zero_replacement_fraction = 0.65,
    transformation_contract = "predicted potential; sample closure before Hellinger or Bray-Curtis; multiplicative zero replacement before CLR/Aitchison"
  )
}

audit <- rbindlist(audit, fill = TRUE)
if (!all(audit$sample_order_identical)) stop("Matrix sample ordering validation failed.", call. = FALSE)
saveRDS(list(metadata = metadata, matrices = matrices, layer_metadata = rbind(data.table(layer = "taxonomic_MAG", role = "comparator"), layers[, .(layer, role)])), file.path(paths$intermediate, "prepared_profiles.rds"))
write_dynamics_table(audit, file.path(paths$prepared, "Table_26_functional_profile_input_and_filtering_audit.csv"))
write_dynamics_session(cfg, "01_prepare_profiles")
writeLines(c("Stage 01 complete", paste("Samples:", length(sample_ids)), paste("Matrices:", paste(names(matrices), collapse = ", ")), "All matrices use identical sample ordering."), marker)
