suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(yaml)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

read_validation_config <- function(path = "config/random_forest_subject_validation.yml") {
  if (!file.exists(path)) stop("Configuration file not found: ", path, call. = FALSE)
  cfg <- yaml::read_yaml(path)
  mode_override <- Sys.getenv("RF_VALIDATION_MODE", unset = "")
  if (nzchar(mode_override)) cfg$analysis$mode <- mode_override
  cfg$.config_path <- path
  validate_validation_config(cfg)
  cfg
}

validate_validation_config <- function(cfg) {
  modes <- c("smoke", "development", "production")
  if (!cfg$analysis$mode %in% modes) stop("analysis.mode must be smoke, development, or production.", call. = FALSE)
  if (!is.numeric(cfg$analysis$seed) || length(cfg$analysis$seed) != 1L) stop("analysis.seed must be one number.", call. = FALSE)
  baseline_repeats <- unlist(cfg$analysis$baseline_repeats[c("smoke", "development", "production")])
  if (length(baseline_repeats) != 3L || any(!is.finite(baseline_repeats)) || any(baseline_repeats < 1)) {
    stop("analysis.baseline_repeats must define positive smoke, development, and production values.", call. = FALSE)
  }
  required_data <- c("input", "existing_run", "sample_id", "subject_id", "feature_id", "abundance", "outcome", "altitude_levels")
  missing <- required_data[vapply(required_data, function(x) is.null(cfg$data[[x]]), logical(1))]
  if (length(missing)) stop("Missing data configuration: ", paste(missing, collapse = ", "), call. = FALSE)
  levels <- as.character(unlist(cfg$data$altitude_levels))
  if (length(levels) < 2L || anyDuplicated(levels)) stop("data.altitude_levels must contain unique class labels.", call. = FALSE)
  numeric_levels <- suppressWarnings(as.numeric(levels))
  if (anyNA(numeric_levels) || is.unsorted(numeric_levels, strictly = TRUE)) {
    stop("data.altitude_levels must be unique numeric labels in increasing order.", call. = FALSE)
  }
  numeric_unit <- function(x, name, allow_zero = FALSE) {
    lower <- if (allow_zero) 0 else 0
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || (!allow_zero && x == 0)) {
      stop(name, " must be one finite ", if (allow_zero) "non-negative" else "positive", " number.", call. = FALSE)
    }
  }
  numeric_unit(cfg$preprocessing$pseudocount, "preprocessing.pseudocount")
  for (name in c("max_missing_fraction", "sample_prevalence")) {
    x <- cfg$preprocessing[[name]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1) stop("preprocessing.", name, " must be in [0, 1].", call. = FALSE)
  }
  numeric_unit(cfg$preprocessing$subject_prevalence_min, "preprocessing.subject_prevalence_min")
  if (!is.numeric(cfg$preprocessing$min_mean_abundance) || length(cfg$preprocessing$min_mean_abundance) != 1L ||
      !is.finite(cfg$preprocessing$min_mean_abundance) || cfg$preprocessing$min_mean_abundance < 0) {
    stop("preprocessing.min_mean_abundance must be one finite non-negative number.", call. = FALSE)
  }
  if (!cfg$preprocessing$transform %in% c("none", "log10_pseudocount")) stop("Unsupported preprocessing.transform.", call. = FALSE)
  if (!identical(cfg$model$implementation, "randomForest")) stop("Only model.implementation=randomForest is supported.", call. = FALSE)
  if (!identical(cfg$model$tuning, "fixed")) stop("Only fixed hyperparameters are implemented; nested tuning is not implemented.", call. = FALSE)
  if (!isTRUE(cfg$model$replace) && cfg$model$sample_fraction > 1) stop("sample_fraction cannot exceed 1 without replacement.", call. = FALSE)
  if (!is.numeric(cfg$model$sample_fraction) || cfg$model$sample_fraction <= 0 || cfg$model$sample_fraction > 1) stop("model.sample_fraction must be in (0, 1].", call. = FALSE)
  for (name in c("smoke_trees", "development_trees", "production_trees", "min_node_size")) numeric_unit(cfg$model[[name]], paste0("model.", name))
  numeric_unit(cfg$grouped_cv$folds, "grouped_cv.folds")
  numeric_unit(cfg$grouped_cv$repeats, "grouped_cv.repeats")
  for (mode in modes) {
    sequence <- unlist(cfg$tree_stability[[mode]])
    if (!length(sequence) || any(!is.finite(sequence)) || any(sequence < 1)) stop("tree_stability.", mode, " must contain positive tree counts.", call. = FALSE)
  }
  if (isTRUE(cfg$outcomes$numeric_regression) || isTRUE(cfg$outcomes$phase) || isTRUE(cfg$outcomes$baseline_change)) {
    stop("Only four-class classification plus ordinal metrics are implemented; disabled outcome flags must remain false.", call. = FALSE)
  }
  invisible(TRUE)
}

resolve_input <- function(cfg, allow_fallback = TRUE) {
  if (file.exists(cfg$data$input)) return(cfg$data$input)
  fallback <- cfg$data$smoke_input_fallback %||% ""
  if (allow_fallback && identical(cfg$analysis$mode, "smoke") && nzchar(fallback) && file.exists(fallback)) {
    warning("Configured input is absent; using documented smoke fallback: ", fallback, call. = FALSE)
    return(fallback)
  }
  stop("Processed input is unavailable: ", cfg$data$input, call. = FALSE)
}

make_run_directory <- function(cfg, run_id = NULL) {
  root <- cfg$results$output_root
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  run_id <- run_id %||% paste0(cfg$analysis$mode, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  candidate <- file.path(root, run_id)
  if (dir.exists(candidate)) stop("Refusing to reuse output directory: ", candidate, call. = FALSE)
  dir.create(candidate, recursive = TRUE)
  for (d in c("audit", "predictions", "comparisons", "features", "reports", "figures/png", "figures/pdf", "figures/svg")) {
    dir.create(file.path(candidate, d), recursive = TRUE, showWarnings = FALSE)
  }
  candidate
}

read_long_data <- function(cfg) {
  path <- resolve_input(cfg)
  x <- readRDS(path)
  required <- unlist(cfg$data[c("sample_id", "subject_id", "feature_id", "abundance", "outcome")], use.names = FALSE)
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Input is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  validate_long_data(x, cfg)
  attr(x, "input_path") <- path
  x
}

validate_long_data <- function(x, cfg) {
  sid <- cfg$data$sample_id; subject <- cfg$data$subject_id; fid <- cfg$data$feature_id
  value <- cfg$data$abundance; outcome <- cfg$data$outcome
  blank <- function(z) is.na(z) | !nzchar(trimws(as.character(z)))
  if (any(blank(x[[sid]]))) stop("Input contains missing or blank sample identifiers.", call. = FALSE)
  if (any(blank(x[[subject]]))) stop("Input contains missing or blank subject identifiers.", call. = FALSE)
  if (any(blank(x[[fid]]))) stop("Input contains missing or blank feature identifiers.", call. = FALSE)
  if (any(is.na(x[[outcome]]))) stop("Input contains missing outcomes.", call. = FALSE)
  observed_levels <- sort(unique(as.character(x[[outcome]])))
  expected_levels <- sort(as.character(unlist(cfg$data$altitude_levels)))
  if (!identical(observed_levels, expected_levels)) stop("Observed outcome classes do not match configured altitude_levels.", call. = FALSE)
  if (!is.numeric(x[[value]])) stop("Abundance column must be numeric.", call. = FALSE)
  if (any(!is.finite(x[[value]]), na.rm = TRUE)) stop("Abundance contains infinite values.", call. = FALSE)
  if (any(x[[value]] < 0, na.rm = TRUE)) stop("Abundance contains negative values.", call. = FALSE)
  key <- paste(x[[sid]], x[[fid]], sep = "\r")
  if (anyDuplicated(key)) stop("Input contains duplicate sampleID-MAG combinations.", call. = FALSE)
  sample_features <- table(x[[sid]])
  n_features <- length(unique(x[[fid]]))
  if (any(sample_features != n_features)) {
    stop("The sample-by-feature grid is incomplete; absent combinations are not assumed to be biological zeros.", call. = FALSE)
  }
  feature_sets <- tapply(x[[fid]], x[[sid]], function(z) paste(sort(z), collapse = "\r"))
  if (length(unique(feature_sets)) != 1L) stop("Samples do not contain identical feature identifiers.", call. = FALSE)
  metadata_cols <- setdiff(names(x), c(fid, value))
  inconsistent <- vapply(split(seq_len(nrow(x)), x[[sid]]), function(index) {
    any(vapply(x[index, metadata_cols, drop = FALSE], function(z) length(unique(z)) != 1L, logical(1)))
  }, logical(1))
  if (any(inconsistent)) stop("Sample metadata are inconsistent within sampleID: ", paste(names(inconsistent)[inconsistent], collapse = ", "), call. = FALSE)
  sample_subjects <- unique(x[c(sid, subject)])
  if (anyDuplicated(sample_subjects[[sid]])) stop("A sampleID maps to more than one subject.", call. = FALSE)
  invisible(TRUE)
}

sample_metadata <- function(x, cfg) {
  sid <- cfg$data$sample_id
  fid <- cfg$data$feature_id
  val <- cfg$data$abundance
  x |> distinct(.data[[sid]], .keep_all = TRUE) |> select(-all_of(c(fid, val)))
}

wide_abundance <- function(x, cfg) {
  sid <- cfg$data$sample_id; fid <- cfg$data$feature_id; val <- cfg$data$abundance
  x |>
    transmute(!!sid := as.character(.data[[sid]]), feature = as.character(.data[[fid]]), value = as.numeric(.data[[val]])) |>
    pivot_wider(names_from = feature, values_from = value)
}

prepare_training_features <- function(train_long, assessment_long, cfg, top_n = NULL) {
  sid <- cfg$data$sample_id; subject <- cfg$data$subject_id; fid <- cfg$data$feature_id; val <- cfg$data$abundance
  train_samples <- n_distinct(train_long[[sid]])
  stats <- train_long |>
    group_by(.data[[fid]]) |>
    summarise(
      missing_fraction = mean(is.na(.data[[val]])),
      sample_prevalence = sum(.data[[val]] > 0, na.rm = TRUE) / train_samples,
      subject_prevalence = n_distinct(.data[[subject]][.data[[val]] > 0 & !is.na(.data[[val]])]),
      mean_abundance = mean(.data[[val]], na.rm = TRUE),
      variance = var(.data[[val]], na.rm = TRUE), .groups = "drop"
    )
  eligible <- stats |>
    filter(missing_fraction <= cfg$preprocessing$max_missing_fraction,
           sample_prevalence >= cfg$preprocessing$sample_prevalence,
           subject_prevalence >= cfg$preprocessing$subject_prevalence_min,
           mean_abundance >= cfg$preprocessing$min_mean_abundance,
           is.finite(variance), variance > 0) |>
    arrange(desc(variance), .data[[fid]])
  requested_top <- top_n %||% cfg$preprocessing$top_n
  if (!is.null(requested_top)) eligible <- slice_head(eligible, n = as.integer(requested_top))
  features <- eligible[[fid]]
  if (length(features) < 2L) stop("Fewer than two training-derived eligible features.", call. = FALSE)
  make_x <- function(z, medians = NULL) {
    w <- wide_abundance(z |> filter(.data[[fid]] %in% features), cfg)
    meta_ids <- unique(as.character(z[[sid]]))
    w <- tibble(!!sid := meta_ids) |> left_join(w, by = sid)
    for (f in features) if (!f %in% names(w)) w[[f]] <- 0
    m <- as.data.frame(w[, features, drop = FALSE], check.names = FALSE)
    m[] <- lapply(m, as.numeric)
    if (is.null(medians)) medians <- vapply(m, median, numeric(1), na.rm = TRUE)
    for (f in features) m[[f]][is.na(m[[f]])] <- medians[[f]]
    if (identical(cfg$preprocessing$transform, "log10_pseudocount")) {
      pc <- cfg$preprocessing$pseudocount
      m[] <- lapply(m, function(v) log10(v + pc))
    } else if (!identical(cfg$preprocessing$transform, "none")) {
      stop("Unsupported transformation.", call. = FALSE)
    }
    list(x = m, ids = meta_ids, medians = medians)
  }
  tr <- make_x(train_long)
  te <- make_x(assessment_long, tr$medians)
  list(train_x = tr$x, assessment_x = te$x, train_ids = tr$ids, assessment_ids = te$ids,
       features = features, eligible_features = eligible[[fid]], filter_statistics = stats,
       preprocessing_fit_samples = unique(train_long[[sid]]), imputation_medians = tr$medians)
}

classification_metrics <- function(observed, predicted, levels, probabilities = NULL) {
  observed <- factor(as.character(observed), levels = levels)
  predicted <- factor(as.character(predicted), levels = levels)
  cm <- table(observed = observed, predicted = predicted)
  per <- bind_rows(lapply(levels, function(cls) {
    tp <- cm[cls, cls]; fp <- sum(cm[, cls]) - tp; fn <- sum(cm[cls, ]) - tp
    precision <- if (tp + fp == 0) 0 else tp / (tp + fp)
    recall <- if (tp + fn == 0) NA_real_ else tp / (tp + fn)
    f1 <- if (is.na(recall) || precision + recall == 0) ifelse(is.na(recall), NA_real_, 0) else 2 * precision * recall / (precision + recall)
    tibble(class = cls, precision, recall, f1, support = sum(observed == cls))
  }))
  if (anyNA(observed)) stop("Observed labels contain missing or unknown classes.", call. = FALSE)
  if (anyNA(predicted)) stop("Predicted labels contain missing or unknown classes.", call. = FALSE)
  ok <- !is.na(observed) & !is.na(predicted)
  ranks <- match(as.character(observed), levels); pred_ranks <- match(as.character(predicted), levels)
  observed_classes <- sum(per$support > 0)
  expected_weights <- seq_along(levels)
  weighted_kappa <- function(a, b) {
    n <- sum(a)
    if (n == 0) return(NA_real_)
    w <- outer(expected_weights, expected_weights, function(i, j) (i - j)^2 / (length(levels) - 1)^2)
    observed_disagreement <- sum(w * a) / n
    expected <- outer(rowSums(a), colSums(a)) / n
    expected_disagreement <- sum(w * expected) / n
    if (expected_disagreement == 0) NA_real_ else 1 - observed_disagreement / expected_disagreement
  }
  metrics <- tibble(
    accuracy = mean(observed[ok] == predicted[ok]),
    balanced_accuracy = mean(per$recall[per$support > 0], na.rm = TRUE),
    balanced_accuracy_all_classes = if (observed_classes == length(levels)) mean(per$recall) else NA_real_,
    macro_precision = mean(per$precision[per$support > 0], na.rm = TRUE),
    macro_recall = mean(per$recall[per$support > 0], na.rm = TRUE),
    macro_f1 = mean(per$f1[per$support > 0], na.rm = TRUE),
    ordinal_mae = mean(abs(ranks - pred_ranks), na.rm = TRUE),
    ordinal_rmse = sqrt(mean((ranks - pred_ranks)^2, na.rm = TRUE)),
    within_one_accuracy = mean(abs(ranks - pred_ranks) <= 1, na.rm = TRUE),
    weighted_kappa = weighted_kappa(cm),
    n_observed_classes = observed_classes
  )
  if (!is.null(probabilities)) {
    if (!all(levels %in% colnames(probabilities))) stop("Probability columns do not match configured class levels.", call. = FALSE)
    p <- as.matrix(probabilities[, levels, drop = FALSE])
    if (nrow(p) != length(observed) || any(!is.finite(p)) || any(p < 0) || any(p > 1)) stop("Probabilities must be finite values in [0, 1].", call. = FALSE)
    if (any(abs(rowSums(p) - 1) > 1e-8)) stop("Probability rows must sum to one.", call. = FALSE)
    idx <- cbind(seq_along(ranks), ranks)
    metrics$log_loss <- -mean(log(pmax(p[idx], 1e-15)), na.rm = TRUE)
  } else metrics$log_loss <- NA_real_
  list(metrics = metrics, per_class = per, confusion = as.data.frame(cm))
}

metrics_from_prediction_table <- function(predictions, levels) {
  probability_columns <- paste0("prob_", levels)
  probabilities <- NULL
  if (all(probability_columns %in% names(predictions))) {
    probabilities <- predictions[, probability_columns, drop = FALSE]
    names(probabilities) <- levels
  }
  classification_metrics(predictions$observed, predictions$predicted, levels, probabilities)
}

assert_disjoint_subjects <- function(train_subjects, assessment_subjects) {
  overlap <- intersect(as.character(train_subjects), as.character(assessment_subjects))
  if (length(overlap)) stop("Subject leakage: ", paste(overlap, collapse = ", "), call. = FALSE)
  TRUE
}
