#!/usr/bin/env Rscript
# Purpose: Fit subject-disjoint leave-one-subject-out Random Forest models with
# fold-local feature preprocessing and importance extraction.
# Inputs: Long-format model data, validation configuration, and LOSO splits.
# Outputs: LOSO predictions, metrics, confusion matrix, fold settings, feature
# eligibility, selection frequency, and importance tables.
# Dependencies: randomForest, R/common.R, and grouped-resampling helpers.
# Usage: Normally invoked by run_random_forest_subject_validation.R; set
# RF_VALIDATION_RUN_DIR for direct execution.
# Workflow: Stage 04 and shared fitter for later sensitivity stages.

source("analysis/random_forest_subject_validation/R/common.R")
source("analysis/random_forest_subject_validation/02_define_grouped_resampling.R")

derive_class_strategy <- function(train_y, levels, strategy) {
  valid <- c("unweighted", "class_weights", "balanced_sampling", "combined")
  if (!strategy %in% valid) stop("Unknown class-imbalance strategy: ", strategy, call. = FALSE)
  counts <- table(factor(as.character(train_y), levels = levels))
  if (any(counts == 0L)) stop("Every configured class must occur in an outer training fold.", call. = FALSE)
  classwt <- if (strategy %in% c("class_weights", "combined")) {
    setNames(sum(counts) / (length(levels) * as.numeric(counts)), levels)
  } else NULL
  sampsize <- if (strategy %in% c("balanced_sampling", "combined")) {
    setNames(rep(min(counts), length(levels)), levels)
  } else NULL
  list(counts = counts, classwt = classwt, sampsize = sampsize)
}

fit_loso_rf <- function(d, cfg, strategy = "unweighted", top_n = NULL, ntree = NULL) {
  if (!requireNamespace("randomForest", quietly = TRUE)) stop("randomForest is required.")
  sid <- cfg$data$sample_id; sub <- cfg$data$subject_id; outcome <- cfg$data$outcome
  levels <- unlist(cfg$data$altitude_levels); meta <- sample_metadata(d, cfg); splits <- make_loso_splits(meta, sub)
  predictions <- list(); fold_metrics <- list(); importances <- list(); selections <- list()
  eligibility <- list(); training_settings <- list()
  ntree <- ntree %||% switch(cfg$analysis$mode, smoke = cfg$model$smoke_trees,
    production = cfg$model$production_trees, cfg$model$development_trees)
  for (i in seq_along(splits)) {
    sp <- splits[[i]]; train_meta <- meta[sp$analysis, ]; test_meta <- meta[sp$assessment, ]
    assert_disjoint_subjects(train_meta[[sub]], test_meta[[sub]])
    train_long <- d[d[[sid]] %in% train_meta[[sid]], ]; test_long <- d[d[[sid]] %in% test_meta[[sid]], ]
    prep <- prepare_training_features(train_long, test_long, cfg, top_n)
    stopifnot(!any(prep$preprocessing_fit_samples %in% test_meta[[sid]]))
    train_y <- factor(as.character(train_meta[[outcome]][match(prep$train_ids, train_meta[[sid]])]), levels = levels)
    test_y <- factor(as.character(test_meta[[outcome]][match(prep$assessment_ids, test_meta[[sid]])]), levels = levels)
    class_settings <- derive_class_strategy(train_y, levels, strategy)
    counts <- class_settings$counts
    classwt <- class_settings$classwt
    sampsize <- class_settings$sampsize
    mtry <- if (identical(cfg$model$mtry, "sqrt")) max(1, floor(sqrt(ncol(prep$train_x)))) else as.integer(cfg$model$mtry)
    set.seed(cfg$analysis$seed + i + match(strategy, c("unweighted", "class_weights", "balanced_sampling", "combined"), nomatch = 0) * 1000)
    fit_args <- list(x = prep$train_x, y = train_y, ntree = ntree, mtry = mtry,
      nodesize = cfg$model$min_node_size, importance = TRUE, replace = isTRUE(cfg$model$replace))
    if (!is.null(classwt)) fit_args$classwt <- classwt
    if (!is.null(sampsize)) {
      fit_args$sampsize <- sampsize
    } else {
      fit_args$sampsize <- max(1L, floor(nrow(prep$train_x) * cfg$model$sample_fraction))
    }
    started <- proc.time()[["elapsed"]]
    model <- do.call(randomForest::randomForest, fit_args)
    elapsed <- proc.time()[["elapsed"]] - started
    pred <- predict(model, prep$assessment_x, type = "response")
    prob <- as.data.frame(predict(model, prep$assessment_x, type = "prob"))
    for (lv in setdiff(levels, names(prob))) prob[[lv]] <- 0
    prob <- prob[, levels, drop = FALSE]
    pm <- tibble(fold_id = sp$id, heldout_subject = sp$holdout_subject, strategy,
      !!sid := prep$assessment_ids, observed = as.character(test_y), predicted = as.character(pred))
    for (lv in levels) pm[[paste0("prob_", lv)]] <- prob[[lv]]
    ev <- classification_metrics(pm$observed, pm$predicted, levels, prob)
    fold_metrics[[i]] <- bind_cols(tibble(
      fold_id = sp$id, heldout_subject = sp$holdout_subject, strategy,
      n_test = nrow(pm), n_observed_classes_subject = n_distinct(pm$observed),
      n_features = length(prep$features), elapsed_seconds = elapsed
    ), ev$metrics)
    imp <- as.data.frame(randomForest::importance(model)); imp$feature <- rownames(imp)
    importances[[i]] <- as_tibble(imp) |> mutate(fold_id = sp$id, heldout_subject = sp$holdout_subject, strategy)
    selections[[i]] <- tibble(fold_id = sp$id, heldout_subject = sp$holdout_subject, strategy, feature = prep$features)
    eligibility[[i]] <- prep$filter_statistics |>
      mutate(
        fold_id = sp$id, heldout_subject = sp$holdout_subject, strategy,
        eligible = .data[[cfg$data$feature_id]] %in% prep$eligible_features,
        selected = .data[[cfg$data$feature_id]] %in% prep$features
      )
    training_settings[[i]] <- tibble(
      fold_id = sp$id, heldout_subject = sp$holdout_subject, strategy,
      train_subjects = paste(sort(unique(as.character(train_meta[[sub]]))), collapse = ";"),
      assessment_subjects = paste(sort(unique(as.character(test_meta[[sub]]))), collapse = ";"),
      n_train_samples = nrow(prep$train_x), n_assessment_samples = nrow(prep$assessment_x),
      n_features = length(prep$features), ntree = ntree, mtry = mtry,
      min_node_size = cfg$model$min_node_size, replace = isTRUE(cfg$model$replace),
      sample_fraction = cfg$model$sample_fraction,
      class_counts = paste(names(counts), as.integer(counts), sep = "=", collapse = ";"),
      class_weights = if (is.null(classwt)) NA_character_ else paste(names(classwt), signif(classwt, 8), sep = "=", collapse = ";"),
      sampling_targets = if (is.null(sampsize)) NA_character_ else paste(names(sampsize), sampsize, sep = "=", collapse = ";"),
      feature_ranking = if (is.null(top_n) && is.null(cfg$preprocessing$top_n)) "all_eligible" else "training_variance"
    )
    predictions[[i]] <- pm
  }
  list(
    predictions = bind_rows(predictions), fold_metrics = bind_rows(fold_metrics),
    importance = bind_rows(importances), selections = bind_rows(selections),
    eligibility = bind_rows(eligibility), training_settings = bind_rows(training_settings)
  )
}

if (sys.nframe() == 0L) {
  cfg <- read_validation_config(); run_dir <- Sys.getenv("RF_VALIDATION_RUN_DIR", "")
  if (!nzchar(run_dir)) run_dir <- make_run_directory(cfg)
  d <- read_long_data(cfg); fit <- fit_loso_rf(d, cfg); levels <- unlist(cfg$data$altitude_levels)
  probs <- fit$predictions |> select(all_of(paste0("prob_", levels))); names(probs) <- levels
  overall <- classification_metrics(fit$predictions$observed, fit$predictions$predicted, levels, probs)
  write.csv(fit$predictions, file.path(run_dir, "predictions", "loso_predictions.csv"), row.names = FALSE)
  write.csv(fit$fold_metrics, file.path(run_dir, "predictions", "loso_fold_metrics.csv"), row.names = FALSE)
  write.csv(overall$metrics, file.path(run_dir, "predictions", "loso_overall_metrics.csv"), row.names = FALSE)
  write.csv(overall$confusion, file.path(run_dir, "predictions", "loso_confusion_matrix.csv"), row.names = FALSE)
  write.csv(fit$importance, file.path(run_dir, "features", "feature_importance_by_fold.csv"), row.names = FALSE)
  write.csv(fit$training_settings, file.path(run_dir, "predictions", "loso_fold_training_settings.csv"), row.names = FALSE)
  write.csv(fit$eligibility, file.path(run_dir, "features", "feature_eligibility_by_fold.csv"), row.names = FALSE)
  fid <- cfg$data$feature_id
  freq <- fit$eligibility |>
    group_by(feature = .data[[fid]]) |>
    summarise(
      eligible_folds = sum(eligible), selected_folds = sum(selected),
      selection_frequency_when_eligible = ifelse(eligible_folds == 0, NA_real_, selected_folds / eligible_folds),
      .groups = "drop"
    )
  write.csv(freq, file.path(run_dir, "features", "feature_selection_frequency.csv"), row.names = FALSE)
  cat(run_dir, "\n")
}
