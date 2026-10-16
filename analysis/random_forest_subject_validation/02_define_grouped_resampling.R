#!/usr/bin/env Rscript
# Purpose: Define deterministic leave-one-subject-out and repeated grouped-fold
# split constructors that keep subjects disjoint across analysis and assessment.
# Inputs: Canonical sample metadata and resampling settings from validation config.
# Outputs: Functions when sourced; a split-validation message when run directly.
# Dependencies: R/common.R and dplyr helpers loaded there.
# Usage: Source from validation stages, or run directly to validate LOSO splits.
# Workflow: Shared resampling component used by baseline, LOSO, and grouped-CV stages.

source("analysis/random_forest_subject_validation/R/common.R")

make_loso_splits <- function(metadata, subject_col) {
  subjects <- sort(unique(as.character(metadata[[subject_col]])))
  if (length(subjects) < 2L) stop("LOSO requires at least two subjects.", call. = FALSE)
  lapply(subjects, function(s) list(id = paste0("LOSO_", s), holdout_subject = s,
    analysis = which(as.character(metadata[[subject_col]]) != s),
    assessment = which(as.character(metadata[[subject_col]]) == s)))
}

subject_class_matrix <- function(metadata, subject_col, outcome_col, levels) {
  table(factor(metadata[[subject_col]], levels = sort(unique(metadata[[subject_col]]))),
        factor(as.character(metadata[[outcome_col]]), levels = levels))
}

make_grouped_folds <- function(metadata, subject_col, outcome_col, levels, k, repeats, seed) {
  subjects <- sort(unique(as.character(metadata[[subject_col]])))
  counts <- subject_class_matrix(metadata, subject_col, outcome_col, levels)
  k <- min(k, length(subjects))
  if (k < 2L) stop("Grouped CV requires at least two folds.", call. = FALSE)
  if (repeats < 1L) stop("Grouped CV repeats must be positive.", call. = FALSE)
  output <- list(); z <- 1L
  for (r in seq_len(repeats)) {
    set.seed(seed + r)
    shuffled <- sample(subjects)
    order_subjects <- shuffled[order(rowSums(counts[shuffled, , drop = FALSE]), decreasing = TRUE)]
    fold_counts <- matrix(0, k, length(levels)); assignment <- setNames(integer(length(subjects)), subjects)
    for (s in order_subjects) {
      scores <- vapply(seq_len(k), function(f) sum((fold_counts[f, ] + counts[s, ])^2) + sum(fold_counts[f, ]), numeric(1))
      chosen <- which(scores == min(scores))[1]
      assignment[s] <- chosen; fold_counts[chosen, ] <- fold_counts[chosen, ] + counts[s, ]
    }
    if (any(assignment == 0L) || any(tabulate(assignment, nbins = k) == 0L)) {
      stop("Grouped fold construction produced an invalid or empty fold.", call. = FALSE)
    }
    for (f in seq_len(k)) {
      assess_subjects <- names(assignment)[assignment == f]
      output[[z]] <- list(id = sprintf("Repeat%02d_Fold%02d", r, f), repeat_id = r, fold = f,
        analysis = which(!metadata[[subject_col]] %in% assess_subjects), assessment = which(metadata[[subject_col]] %in% assess_subjects))
      z <- z + 1L
    }
  }
  output
}

if (sys.nframe() == 0L) {
  cfg <- read_validation_config(); d <- read_long_data(cfg); m <- sample_metadata(d, cfg)
  splits <- make_loso_splits(m, cfg$data$subject_id)
  stopifnot(length(splits) == n_distinct(m[[cfg$data$subject_id]]))
  cat("LOSO splits validated:", length(splits), "\n")
}
