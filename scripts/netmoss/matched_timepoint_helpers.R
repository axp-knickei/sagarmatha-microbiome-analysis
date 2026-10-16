# Shared functions for the matched-timepoint companion analysis.
#
# This file deliberately contains no NetMoss2 call. The companion workflow
# estimates paired abundance changes; it does not pretend that one participant
# or a small matched subset is a microbial association network.

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x)) || identical(x, "")) y else x
}

read_csv_checked <- function(path, label) {
  if (!file.exists(path)) stop(label, " file not found: ", path, call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

assert_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

assert_output <- function(path, force) {
  if (file.exists(path) && !force) {
    stop("Refusing to overwrite ", path, ". Use --force.", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

close_columns <- function(matrix) {
  totals <- colSums(matrix, na.rm = TRUE)
  if (any(!is.finite(totals) | totals <= 0)) {
    stop("One or more sample profiles have non-positive totals.", call. = FALSE)
  }
  sweep(matrix, 2, totals, "/")
}

prepare_genus_matrices <- function(abundance, metadata, taxonomy, config) {
  columns <- config$columns
  assert_columns(abundance, columns$sample_id, "Abundance")
  assert_columns(
    metadata,
    c(columns$sample_id, columns$subject_id, columns$timepoint),
    "Metadata"
  )
  assert_columns(
    taxonomy,
    c(columns$taxonomy_feature_id, columns$taxonomy_rank),
    "Taxonomy"
  )

  if (anyDuplicated(abundance[[columns$sample_id]])) {
    stop("Abundance contains duplicated sample IDs.", call. = FALSE)
  }
  if (anyDuplicated(metadata[[columns$sample_id]])) {
    stop("Metadata contains duplicated sample IDs.", call. = FALSE)
  }

  feature_ids <- setdiff(names(abundance), columns$sample_id)
  numeric_features <- vapply(abundance[feature_ids], is.numeric, logical(1))
  if (!all(numeric_features)) {
    stop(
      "Non-numeric abundance columns: ",
      paste(feature_ids[!numeric_features], collapse = ", "),
      call. = FALSE
    )
  }

  sample_ids <- as.character(abundance[[columns$sample_id]])
  metadata_ids <- as.character(metadata[[columns$sample_id]])
  missing_metadata <- setdiff(sample_ids, metadata_ids)
  if (length(missing_metadata)) {
    stop("Abundance samples lack metadata: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
  }

  taxonomy_index <- match(feature_ids, as.character(taxonomy[[columns$taxonomy_feature_id]]))
  if (anyNA(taxonomy_index)) {
    stop("Some abundance features lack an exact taxonomy match.", call. = FALSE)
  }
  genus_names <- trimws(as.character(taxonomy[[columns$taxonomy_rank]][taxonomy_index]))
  missing_genus <- is.na(genus_names) | !nzchar(genus_names)
  genus_names[missing_genus] <- feature_ids[missing_genus]

  sample_by_mag <- as.matrix(abundance[feature_ids])
  storage.mode(sample_by_mag) <- "double"
  rownames(sample_by_mag) <- sample_ids
  if (any(!is.finite(sample_by_mag)) || any(sample_by_mag < 0)) {
    stop("Abundance values must be finite and non-negative.", call. = FALSE)
  }

  genus_by_sample <- rowsum(t(sample_by_mag), group = genus_names, reorder = TRUE)
  relative_abundance <- close_columns(genus_by_sample)
  positive_values <- relative_abundance[relative_abundance > 0]
  if (!length(positive_values)) stop("No positive genus abundances were found.", call. = FALSE)
  zero_fraction <- as.numeric(
    config$transformation$zero_replacement_fraction_of_smallest_positive
  )
  zero_replacement <- min(positive_values) * zero_fraction
  replaced <- relative_abundance
  replaced[replaced == 0] <- zero_replacement
  replaced <- close_columns(replaced)
  clr <- log(replaced)
  clr <- sweep(clr, 2, colMeans(clr), "-")

  matched_metadata <- metadata[match(sample_ids, metadata_ids), , drop = FALSE]
  matched_metadata$sample_id_internal <- sample_ids
  matched_metadata$subject_id_internal <- as.character(
    matched_metadata[[columns$subject_id]]
  )
  matched_metadata$timepoint_internal <- as.character(
    matched_metadata[[columns$timepoint]]
  )

  duplicate_grain <- duplicated(
    matched_metadata[c("subject_id_internal", "timepoint_internal")]
  )
  if (any(duplicate_grain)) {
    stop("More than one abundance sample exists for a subject-timepoint pair.", call. = FALSE)
  }

  list(
    raw = genus_by_sample,
    relative = relative_abundance,
    clr = clr,
    metadata = matched_metadata,
    zero_replacement = zero_replacement,
    feature_count = length(feature_ids),
    genus_count = nrow(genus_by_sample)
  )
}

make_timepoint_pairs <- function(timepoints) {
  indices <- utils::combn(seq_along(timepoints), 2L)
  data.frame(
    pair_index = seq_len(ncol(indices)),
    timepoint_a = timepoints[indices[1L, ]],
    timepoint_b = timepoints[indices[2L, ]],
    pair_id = paste0(timepoints[indices[1L, ]], "_vs_", timepoints[indices[2L, ]]),
    stringsAsFactors = FALSE
  )
}

exact_sign_flip_pvalues <- function(changes) {
  participant_count <- ncol(changes)
  sign_grid <- as.matrix(
    expand.grid(rep(list(c(-1, 1)), participant_count), KEEP.OUT.ATTRS = FALSE)
  )
  null_means <- changes %*% t(sign_grid) / participant_count
  observed <- rowMeans(changes)
  rowMeans(abs(null_means) >= abs(observed) - sqrt(.Machine$double.eps))
}

bootstrap_mean_intervals <- function(changes, replicates, confidence_level) {
  participant_count <- ncol(changes)
  weights <- vapply(seq_len(replicates), function(unused) {
    tabulate(
      sample.int(participant_count, participant_count, replace = TRUE),
      nbins = participant_count
    ) / participant_count
  }, numeric(participant_count))
  bootstrap_means <- changes %*% weights
  alpha <- (1 - confidence_level) / 2
  intervals <- apply(
    bootstrap_means,
    1,
    stats::quantile,
    probs = c(alpha, 1 - alpha),
    names = FALSE,
    type = 8
  )
  if (is.null(dim(intervals))) intervals <- matrix(intervals, nrow = 2L)
  list(lower = intervals[1L, ], upper = intervals[2L, ])
}

leave_one_out_sign_stability <- function(changes) {
  participant_count <- ncol(changes)
  observed_sign <- sign(rowMeans(changes))
  if (participant_count <= 2L) return(rep(NA_real_, nrow(changes)))
  leave_one_out_means <- vapply(seq_len(participant_count), function(index) {
    rowMeans(changes[, -index, drop = FALSE])
  }, numeric(nrow(changes)))
  rowMeans(sign(leave_one_out_means) == observed_sign)
}

analyze_timepoint_pair <- function(pair, prepared, config, bootstrap_replicates) {
  metadata <- prepared$metadata
  metadata_a <- metadata[metadata$timepoint_internal == pair$timepoint_a, , drop = FALSE]
  metadata_b <- metadata[metadata$timepoint_internal == pair$timepoint_b, , drop = FALSE]
  subjects <- sort(intersect(metadata_a$subject_id_internal, metadata_b$subject_id_internal))
  participant_count <- length(subjects)

  sample_ids_for <- function(data, selected_subjects) {
    data$sample_id_internal[match(selected_subjects, data$subject_id_internal)]
  }
  ids_a <- sample_ids_for(metadata_a, subjects)
  ids_b <- sample_ids_for(metadata_b, subjects)
  raw_a <- prepared$raw[, ids_a, drop = FALSE]
  raw_b <- prepared$raw[, ids_b, drop = FALSE]
  relative_a <- prepared$relative[, ids_a, drop = FALSE]
  relative_b <- prepared$relative[, ids_b, drop = FALSE]
  clr_a <- prepared$clr[, ids_a, drop = FALSE]
  clr_b <- prepared$clr[, ids_b, drop = FALSE]
  clr_change <- clr_b - clr_a

  detected_participants <- rowSums(raw_a > 0 | raw_b > 0)
  minimum_pairs <- as.integer(config$design$minimum_pairs_for_inference)
  minimum_detected <- as.integer(config$design$minimum_detected_participants)
  inferential_pair <- participant_count >= minimum_pairs
  testable <- inferential_pair & detected_participants >= minimum_detected

  p_value <- rep(NA_real_, nrow(clr_change))
  if (inferential_pair && any(testable)) {
    p_value[testable] <- exact_sign_flip_pvalues(clr_change[testable, , drop = FALSE])
  }

  confidence_level <- as.numeric(config$inference$bootstrap_confidence_level)
  intervals <- bootstrap_mean_intervals(
    clr_change,
    replicates = bootstrap_replicates,
    confidence_level = confidence_level
  )
  relative_pseudocount <- prepared$zero_replacement
  log2_relative_change <- log2(
    (relative_b + relative_pseudocount) / (relative_a + relative_pseudocount)
  )
  concordance <- pmax(rowMeans(clr_change > 0), rowMeans(clr_change < 0))

  genus_results <- data.frame(
    pair_index = pair$pair_index,
    pair_id = pair$pair_id,
    timepoint_a = pair$timepoint_a,
    timepoint_b = pair$timepoint_b,
    matched_participants = participant_count,
    genus = rownames(clr_change),
    detected_participants = detected_participants,
    inferential_pair = inferential_pair,
    testable = testable,
    mean_clr_change = rowMeans(clr_change),
    median_clr_change = apply(clr_change, 1, stats::median),
    mean_log2_relative_change = rowMeans(log2_relative_change),
    direction_concordance = concordance,
    leave_one_out_sign_stability = leave_one_out_sign_stability(clr_change),
    bootstrap_ci_lower = intervals$lower,
    bootstrap_ci_upper = intervals$upper,
    exact_p_value = p_value,
    stringsAsFactors = FALSE
  )

  participant_results <- do.call(rbind, lapply(seq_along(subjects), function(index) {
    data.frame(
      pair_index = pair$pair_index,
      pair_id = pair$pair_id,
      timepoint_a = pair$timepoint_a,
      timepoint_b = pair$timepoint_b,
      subject_id = subjects[[index]],
      sample_id_a = ids_a[[index]],
      sample_id_b = ids_b[[index]],
      genus = rownames(clr_change),
      relative_abundance_a = relative_a[, index],
      relative_abundance_b = relative_b[, index],
      clr_a = clr_a[, index],
      clr_b = clr_b[, index],
      clr_change = clr_change[, index],
      direction = ifelse(
        clr_change[, index] > 0,
        "increase",
        ifelse(clr_change[, index] < 0, "decrease", "unchanged")
      ),
      stringsAsFactors = FALSE
    )
  }))

  pair_summary <- data.frame(
    pair_index = pair$pair_index,
    pair_id = pair$pair_id,
    timepoint_a = pair$timepoint_a,
    timepoint_b = pair$timepoint_b,
    matched_participants = participant_count,
    inferential_pair = inferential_pair,
    netmoss_network_feasible = participant_count >= as.integer(
      config$design$netmoss_network_minimum_profiles
    ),
    participant_ids = paste(subjects, collapse = ";"),
    stringsAsFactors = FALSE
  )

  list(
    pair_summary = pair_summary,
    genus_results = genus_results,
    participant_results = participant_results
  )
}
