duplicates <- function(x) unique(x[duplicated(x)])

validate_community <- function(abundance, metadata, taxonomy, config) {
  sample_col <- config$columns$sample_id
  feature_id <- config$columns$taxonomy_feature_id
  if (!sample_col %in% names(abundance)) stop("FAILED_NO_SAMPLE_BY_MAG_MATRIX: sample ID column is absent: ", sample_col, call. = FALSE)
  model_signature <- c("estimate", "std.error", "p.value", "p_value", "coefficient", "time_p")
  if (sum(tolower(names(abundance)) %in% model_signature) >= 2) stop("FAILED_NO_SAMPLE_BY_MAG_MATRIX: model-summary columns detected; coefficients are not abundance profiles.", call. = FALSE)
  sample_ids <- as.character(abundance[[sample_col]])
  feature_cols <- setdiff(names(abundance), sample_col)
  if (length(duplicates(sample_ids))) stop("FAILED_DUPLICATE_SAMPLE_IDS: ", paste(head(duplicates(sample_ids)), collapse = ", "), call. = FALSE)
  if (length(duplicates(feature_cols))) stop("FAILED_DUPLICATE_MAG_IDS", call. = FALSE)
  numeric_ok <- vapply(abundance[feature_cols], is.numeric, logical(1))
  if (!all(numeric_ok)) stop("FAILED_NON_NUMERIC_FEATURES: ", paste(feature_cols[!numeric_ok], collapse = ", "), call. = FALSE)
  matrix <- as.matrix(abundance[feature_cols]); storage.mode(matrix) <- "double"
  rownames(matrix) <- sample_ids
  if (any(!is.finite(matrix))) stop("FAILED_NON_FINITE_ABUNDANCE", call. = FALSE)
  if (any(matrix < 0)) stop("FAILED_NEGATIVE_ABUNDANCE", call. = FALSE)
  metadata_ids <- as.character(metadata[[sample_col]])
  if (length(duplicates(metadata_ids))) stop("FAILED_DUPLICATE_SAMPLE_IDS in metadata", call. = FALSE)
  missing_metadata <- setdiff(sample_ids, metadata_ids)
  if (length(missing_metadata)) stop("FAILED_METADATA_ID_MISMATCH: abundance samples lack metadata: ", paste(head(missing_metadata), collapse = ", "), call. = FALSE)
  if (!feature_id %in% names(taxonomy)) stop("FAILED_TAXONOMY_FEATURE_ID", call. = FALSE)
  if (length(duplicates(as.character(taxonomy[[feature_id]])))) stop("FAILED_DUPLICATE_TAXONOMY_IDS", call. = FALSE)
  list(matrix = matrix, sample_ids = sample_ids, feature_ids = feature_cols,
       metadata_only = setdiff(metadata_ids, sample_ids), taxonomy_missing = setdiff(feature_cols, taxonomy[[feature_id]]))
}

validate_distance <- function(distance, expected_ids) {
  matrix <- as.matrix(distance)
  stopifnot(identical(rownames(matrix), expected_ids), identical(colnames(matrix), expected_ids))
  if (any(!is.finite(matrix))) stop("Distance contains NA, NaN, or infinite values.", call. = FALSE)
  if (max(abs(matrix - t(matrix))) > 1e-10) stop("Distance is not symmetric.", call. = FALSE)
  if (max(abs(diag(matrix))) > 1e-10) stop("Distance diagonal is not zero.", call. = FALSE)
  invisible(TRUE)
}
