# Shared preparation, resampling, and stability functions for the exploratory
# low/initial-versus-high-altitude NetMoss stability workflow.

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

validate_stability_config <- function(config) {
  group_a <- unlist(config$design$group_a_timepoints)
  group_b <- unlist(config$design$group_b_timepoints)
  excluded <- unlist(config$design$excluded_timepoints %||% character(0))
  if (length(intersect(group_a, group_b))) stop("Configured analysis groups overlap.", call. = FALSE)
  if (length(intersect(c(group_a, group_b), excluded))) {
    stop("An analyzed timepoint is also configured as excluded.", call. = FALSE)
  }
  if (!identical(config$taxonomic_analyses$prevalence_rule, "detected_in_both_groups")) {
    stop("Only the predeclared detected_in_both_groups prevalence rule is supported.", call. = FALSE)
  }
  thresholds <- c(
    unlist(config$taxonomic_analyses$genus$sensitivity_thresholds),
    unlist(config$taxonomic_analyses$family$sensitivity_thresholds),
    config$taxonomic_analyses$genus$primary_prevalence_threshold,
    config$taxonomic_analyses$family$primary_prevalence_threshold,
    config$stability$minimum_success_fraction,
    config$stability$selection_frequency_threshold,
    config$stability$network_similarity_threshold,
    config$stability$module_ari_threshold,
    config$stability$degenerate_module_fraction
  )
  if (any(!is.finite(thresholds)) || any(thresholds < 0 | thresholds > 1)) {
    stop("All prevalence and stability thresholds must be within [0, 1].", call. = FALSE)
  }
  if (as.integer(config$netmoss$minimum_module_size) < 2L) {
    stop("NetMoss minimum module size must be at least two.", call. = FALSE)
  }
  invisible(TRUE)
}

close_columns <- function(matrix) {
  totals <- colSums(matrix, na.rm = TRUE)
  if (any(!is.finite(totals) | totals <= 0)) {
    stop("One or more profiles have non-positive totals.", call. = FALSE)
  }
  sweep(matrix, 2, totals, "/")
}

aggregate_subject_profiles <- function(matrix, metadata, subjects, subject_column, method) {
  profiles <- vapply(subjects, function(subject) {
    ids <- metadata$sample_id_internal[as.character(metadata[[subject_column]]) == subject]
    values <- matrix[, ids, drop = FALSE]
    if (!ncol(values)) stop("No samples found for participant ", subject, call. = FALSE)
    if (method == "median") apply(values, 1, stats::median)
    else if (method == "mean") rowMeans(values)
    else stop("Unsupported participant aggregation: ", method, call. = FALSE)
  }, numeric(nrow(matrix)))
  if (is.null(dim(profiles))) profiles <- matrix(profiles, ncol = 1L)
  rownames(profiles) <- rownames(matrix)
  colnames(profiles) <- subjects
  profiles
}

prepare_rank_profiles <- function(abundance, metadata, taxonomy, config, rank) {
  columns <- config$columns
  assert_columns(abundance, columns$sample_id, "Abundance")
  assert_columns(metadata, c(columns$sample_id, columns$subject_id, columns$timepoint), "Metadata")
  assert_columns(taxonomy, c(columns$taxonomy_feature_id, rank), "Taxonomy")
  if (anyDuplicated(abundance[[columns$sample_id]])) stop("Duplicated abundance sample IDs.", call. = FALSE)
  if (anyDuplicated(metadata[[columns$sample_id]])) stop("Duplicated metadata sample IDs.", call. = FALSE)

  feature_ids <- setdiff(names(abundance), columns$sample_id)
  if (!all(vapply(abundance[feature_ids], is.numeric, logical(1)))) {
    stop("All abundance features must be numeric.", call. = FALSE)
  }
  sample_ids <- as.character(abundance[[columns$sample_id]])
  metadata_index <- match(sample_ids, as.character(metadata[[columns$sample_id]]))
  if (anyNA(metadata_index)) stop("One or more abundance samples lack metadata.", call. = FALSE)
  taxonomy_index <- match(feature_ids, as.character(taxonomy[[columns$taxonomy_feature_id]]))
  if (anyNA(taxonomy_index)) stop("One or more abundance features lack taxonomy.", call. = FALSE)

  labels <- trimws(as.character(taxonomy[[rank]][taxonomy_index]))
  missing_labels <- is.na(labels) | !nzchar(labels)
  labels[missing_labels] <- feature_ids[missing_labels]
  sample_matrix <- as.matrix(abundance[feature_ids])
  storage.mode(sample_matrix) <- "double"
  rownames(sample_matrix) <- sample_ids
  if (any(!is.finite(sample_matrix)) || any(sample_matrix < 0)) {
    stop("Abundance values must be finite and non-negative.", call. = FALSE)
  }
  rank_by_sample <- rowsum(t(sample_matrix), group = labels, reorder = TRUE)
  matched_metadata <- metadata[metadata_index, , drop = FALSE]
  matched_metadata$sample_id_internal <- sample_ids
  timepoint <- as.character(matched_metadata[[columns$timepoint]])
  excluded <- unlist(config$design$excluded_timepoints %||% character(0))
  keep <- !timepoint %in% excluded
  matched_metadata <- matched_metadata[keep, , drop = FALSE]
  rank_by_sample <- rank_by_sample[, matched_metadata$sample_id_internal, drop = FALSE]

  timepoint <- as.character(matched_metadata[[columns$timepoint]])
  subject <- as.character(matched_metadata[[columns$subject_id]])
  group_a_meta <- matched_metadata[timepoint %in% unlist(config$design$group_a_timepoints), , drop = FALSE]
  group_b_meta <- matched_metadata[timepoint %in% unlist(config$design$group_b_timepoints), , drop = FALSE]
  subjects_a <- sort(unique(as.character(group_a_meta[[columns$subject_id]])))
  subjects_b <- sort(unique(as.character(group_b_meta[[columns$subject_id]])))
  common <- sort(intersect(subjects_a, subjects_b))
  method <- config$design$participant_aggregation
  list(
    rank = rank,
    feature_count = nrow(rank_by_sample),
    matched_subjects = common,
    subjects_a = subjects_a,
    subjects_b = subjects_b,
    matched_a = aggregate_subject_profiles(rank_by_sample, group_a_meta, common, columns$subject_id, method),
    matched_b = aggregate_subject_profiles(rank_by_sample, group_b_meta, common, columns$subject_id, method),
    all_a = aggregate_subject_profiles(rank_by_sample, group_a_meta, subjects_a, columns$subject_id, method),
    all_b = aggregate_subject_profiles(rank_by_sample, group_b_meta, subjects_b, columns$subject_id, method),
    raw_sample_count_a = nrow(group_a_meta),
    raw_sample_count_b = nrow(group_b_meta)
  )
}

filter_rank_profiles <- function(group_a, group_b, threshold, detection = 0) {
  prevalence_a <- rowMeans(group_a > detection)
  prevalence_b <- rowMeans(group_b > detection)
  variable <- apply(cbind(group_a, group_b), 1, stats::var) > 0
  retained <- prevalence_a >= threshold & prevalence_b >= threshold & variable
  table <- data.frame(
    taxon = rownames(group_a),
    prevalence_group_a = prevalence_a,
    prevalence_group_b = prevalence_b,
    variable = variable,
    retained = retained,
    stringsAsFactors = FALSE
  )
  list(
    group_a = group_a[retained, , drop = FALSE],
    group_b = group_b[retained, , drop = FALSE],
    filtering = table
  )
}

build_network <- function(abundance, method) {
  if (method == "sparcc") {
    if (!requireNamespace("rsparcc", quietly = TRUE)) stop("rsparcc is required.", call. = FALSE)
    invisible(utils::capture.output(
      result <- suppressWarnings(suppressMessages(rsparcc::sparcc(t(abundance)))),
      type = "output"
    ))
    network <- result$CORR
  } else if (method == "spearman") {
    network <- stats::cor(t(abundance), method = "spearman", use = "pairwise.complete.obs")
  } else stop("Unsupported network method: ", method, call. = FALSE)
  network[!is.finite(network)] <- 0
  network <- (network + t(network)) / 2
  diag(network) <- 1
  rownames(network) <- rownames(abundance)
  colnames(network) <- rownames(abundance)
  network
}

normalize_netmoss_scores <- function(result, taxa) {
  scores <- as.data.frame(result[[1]], stringsAsFactors = FALSE)
  clean <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  cleaned <- clean(names(scores))
  score_index <- match(clean(c("NetMoss_Score", "NetMoss_score", "score")), cleaned, nomatch = 0L)
  score_index <- score_index[score_index > 0L]
  if (!nrow(scores) || !length(score_index)) {
    aligned <- stats::setNames(rep(0, length(taxa)), taxa)
    attr(aligned, "scored") <- stats::setNames(rep(FALSE, length(taxa)), taxa)
    return(aligned)
  }
  taxon_index <- match(clean(c("taxon_names", "taxon", "name")), cleaned, nomatch = 0L)
  taxon_index <- taxon_index[taxon_index > 0L]
  score_taxa <- if (length(taxon_index)) as.character(scores[[taxon_index[[1]]]]) else rownames(scores)
  values <- suppressWarnings(as.numeric(scores[[score_index[[1]]]]))
  aligned <- values[match(taxa, score_taxa)]
  scored <- is.finite(aligned)
  if (!any(scored)) stop("NetMoss did not return any aligned taxon scores.", call. = FALSE)
  aligned[!scored] <- 0
  aligned <- stats::setNames(aligned, taxa)
  attr(aligned, "scored") <- stats::setNames(scored, taxa)
  aligned
}

module_is_degenerate <- function(labels, fraction_threshold) {
  assigned <- labels[labels > 0]
  if (!length(assigned) || length(unique(assigned)) < 2L) return(TRUE)
  max(table(assigned)) / length(labels) >= fraction_threshold
}

adjusted_rand_index <- function(reference, observed) {
  if (length(reference) != length(observed) || length(reference) < 2L) return(NA_real_)
  contingency <- table(reference, observed)
  choose2 <- function(x) x * (x - 1) / 2
  sum_cells <- sum(choose2(contingency))
  sum_rows <- sum(choose2(rowSums(contingency)))
  sum_cols <- sum(choose2(colSums(contingency)))
  total <- choose2(length(reference))
  if (total == 0) return(NA_real_)
  expected <- sum_rows * sum_cols / total
  maximum <- (sum_rows + sum_cols) / 2
  if (maximum == expected) return(if (identical(reference, observed)) 1 else 0)
  (sum_cells - expected) / (maximum - expected)
}

run_netmoss_analysis <- function(group_a, group_b, network_method, config, seed) {
  set.seed(seed)
  group_a <- close_columns(group_a)
  group_b <- close_columns(group_b)
  network_a <- build_network(group_a, network_method)
  network_b <- build_network(group_b, network_method)
  taxa <- rownames(group_a)
  control <- data.frame(X = taxa, group_a, check.names = FALSE)
  case <- data.frame(X = taxa, group_b, check.names = FALSE)
  invisible(utils::capture.output(
    result <- suppressWarnings(suppressMessages(NetMoss2::NetMoss(
      case_dir = case,
      control_dir = control,
      net_case_dir = network_b,
      net_control_dir = network_a,
      scaled = isTRUE(config$netmoss$scaled),
      deepSplit = as.integer(config$netmoss$deep_split),
      minModuleSize = as.integer(config$netmoss$minimum_module_size),
      soft_pow = as.numeric(config$netmoss$soft_power)
    ))),
    type = "output"
  ))
  invisible(utils::capture.output(
    modules <- suppressWarnings(suppressMessages(NetMoss2:::divModule(
      case_union = result[[2]],
      control_union = result[[3]],
      deepSplit = as.integer(config$netmoss$deep_split),
      minModuleSize = as.integer(config$netmoss$minimum_module_size),
      soft_pow = as.numeric(config$netmoss$soft_power)
    ))),
    type = "output"
  ))
  module_b <- stats::setNames(as.integer(modules[[3]]), taxa)
  module_a <- stats::setNames(as.integer(modules[[4]]), taxa)
  normalized_scores <- normalize_netmoss_scores(result, taxa)
  list(
    network_a = network_a,
    network_b = network_b,
    scores = normalized_scores,
    scored = attr(normalized_scores, "scored"),
    module_a = module_a,
    module_b = module_b
  )
}

upper_triangle <- function(matrix) matrix[upper.tri(matrix)]

same_strong_edge <- function(reference, observed, threshold) {
  reference_values <- upper_triangle(reference)
  observed_values <- upper_triangle(observed)
  abs(observed_values) >= threshold & sign(observed_values) == sign(reference_values)
}

network_similarity <- function(reference, observed) {
  suppressWarnings(stats::cor(upper_triangle(reference), upper_triangle(observed), method = "spearman"))
}

make_resample_specifications <- function(group_a, group_b, bootstrap_replicates, seed) {
  subjects_a <- colnames(group_a)
  subjects_b <- colnames(group_b)
  common <- intersect(subjects_a, subjects_b)
  extra_a <- setdiff(subjects_a, common)
  extra_b <- setdiff(subjects_b, common)
  set.seed(seed)
  bootstrap <- lapply(seq_len(bootstrap_replicates), function(index) {
    common_draw <- sample(common, length(common), replace = TRUE)
    list(
      type = "bootstrap",
      id = sprintf("bootstrap_%04d", index),
      omitted = NA_character_,
      subjects_a = c(common_draw, if (length(extra_a)) sample(extra_a, length(extra_a), replace = TRUE)),
      subjects_b = c(common_draw, if (length(extra_b)) sample(extra_b, length(extra_b), replace = TRUE)),
      seed = seed + index
    )
  })
  union_subjects <- sort(union(subjects_a, subjects_b))
  leave_one_out <- lapply(seq_along(union_subjects), function(index) {
    omitted <- union_subjects[[index]]
    list(
      type = "leave_one_out",
      id = paste0("leave_out_", omitted),
      omitted = omitted,
      subjects_a = setdiff(subjects_a, omitted),
      subjects_b = setdiff(subjects_b, omitted),
      seed = seed + bootstrap_replicates + index
    )
  })
  c(bootstrap, leave_one_out)
}

select_resampled_profiles <- function(matrix, subjects, prefix) {
  selected <- matrix[, subjects, drop = FALSE]
  colnames(selected) <- sprintf("%s_draw_%03d", prefix, seq_along(subjects))
  selected
}

analyze_resample <- function(specification, group_a, group_b, reference, config, network_method) {
  tryCatch({
    observed <- run_netmoss_analysis(
      select_resampled_profiles(group_a, specification$subjects_a, "A"),
      select_resampled_profiles(group_b, specification$subjects_b, "B"),
      network_method,
      config,
      specification$seed
    )
    taxa <- names(reference$scores)
    ranks <- stats::setNames(rep(Inf, length(taxa)), taxa)
    observed_scored <- observed$scored[taxa] %in% TRUE
    ranks[observed_scored] <- rank(-observed$scores[taxa][observed_scored], ties.method = "min")
    top <- observed_scored & ranks <= as.integer(config$netmoss$top_taxa)
    degenerate_a <- module_is_degenerate(observed$module_a, config$stability$degenerate_module_fraction)
    degenerate_b <- module_is_degenerate(observed$module_b, config$stability$degenerate_module_fraction)
    list(
      success = TRUE,
      type = specification$type,
      id = specification$id,
      omitted = specification$omitted,
      error = NA_character_,
      network_similarity_a = network_similarity(reference$network_a, observed$network_a),
      network_similarity_b = network_similarity(reference$network_b, observed$network_b),
      score_rank_correlation = suppressWarnings(stats::cor(reference$scores, observed$scores[taxa], method = "spearman")),
      module_ari_a = if (degenerate_a) NA_real_ else adjusted_rand_index(reference$module_a, observed$module_a),
      module_ari_b = if (degenerate_b) NA_real_ else adjusted_rand_index(reference$module_b, observed$module_b),
      degenerate_module_a = degenerate_a,
      degenerate_module_b = degenerate_b,
      edge_same_a = same_strong_edge(reference$network_a, observed$network_a, config$network$edge_magnitude_threshold),
      edge_same_b = same_strong_edge(reference$network_b, observed$network_b, config$network$edge_magnitude_threshold),
      top_taxa = stats::setNames(top, taxa)
    )
  }, error = function(error) {
    list(
      success = FALSE,
      type = specification$type,
      id = specification$id,
      omitted = specification$omitted,
      error = conditionMessage(error)
    )
  })
}

edge_index <- function(taxa, condition, reference) {
  indices <- which(upper.tri(reference), arr.ind = TRUE)
  data.frame(
    condition = condition,
    taxon_1 = taxa[indices[, 1]],
    taxon_2 = taxa[indices[, 2]],
    reference_correlation = reference[indices],
    stringsAsFactors = FALSE
  )
}

operational_network_pass <- function(success_fraction, median_a, median_b, config) {
  is.finite(success_fraction) &&
    success_fraction >= as.numeric(config$stability$minimum_success_fraction) &&
    is.finite(median_a) && median_a >= as.numeric(config$stability$network_similarity_threshold) &&
    is.finite(median_b) && median_b >= as.numeric(config$stability$network_similarity_threshold)
}

operational_selection <- function(reference_selected, bootstrap_frequency, loo_count,
                                  requested_loo, network_pass, config) {
  loo_required <- ceiling(
    as.numeric(config$stability$selection_frequency_threshold) * requested_loo
  )
  network_pass & reference_selected &
    bootstrap_frequency >= as.numeric(config$stability$selection_frequency_threshold) &
    loo_count >= loo_required
}
