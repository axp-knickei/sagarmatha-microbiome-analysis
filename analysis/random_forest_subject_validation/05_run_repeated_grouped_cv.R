#!/usr/bin/env Rscript
# Purpose: Run repeated subject-grouped cross-validation as a complement to
# leave-one-subject-out validation.
# Inputs: Long-format model data, validation configuration, and grouped folds.
# Outputs: Grouped-CV predictions, fold metrics, and fold settings under
# RF_VALIDATION_RUN_DIR/predictions.
# Dependencies: randomForest, R/common.R, and grouped-resampling helpers.
# Usage: Normally invoked by run_random_forest_subject_validation.R; set
# RF_VALIDATION_RUN_DIR for direct execution.
# Workflow: Stage 05; provides repeated grouped-resampling evidence.

source("analysis/random_forest_subject_validation/R/common.R")
source("analysis/random_forest_subject_validation/02_define_grouped_resampling.R")

run_grouped_cv <- function(d, cfg, repeats = NULL) {
  sid <- cfg$data$sample_id; sub <- cfg$data$subject_id; outcome <- cfg$data$outcome; levels <- unlist(cfg$data$altitude_levels)
  meta <- sample_metadata(d, cfg); repeats <- repeats %||% ifelse(cfg$analysis$mode == "smoke", 1L, cfg$grouped_cv$repeats)
  ntree <- switch(cfg$analysis$mode, smoke = cfg$model$smoke_trees,
    production = cfg$model$production_trees, cfg$model$development_trees)
  splits <- make_grouped_folds(meta, sub, outcome, levels, cfg$grouped_cv$folds, repeats, cfg$analysis$seed)
  pred <- list(); met <- list(); settings <- list()
  for (i in seq_along(splits)) {
    sp <- splits[[i]]; train <- meta[sp$analysis, ]; test <- meta[sp$assessment, ]
    assert_disjoint_subjects(train[[sub]], test[[sub]])
    prep <- prepare_training_features(d[d[[sid]] %in% train[[sid]], ], d[d[[sid]] %in% test[[sid]], ], cfg)
    train_y <- factor(as.character(train[[outcome]][match(prep$train_ids, train[[sid]])]), levels = levels)
    test_y <- factor(as.character(test[[outcome]][match(prep$assessment_ids, test[[sid]])]), levels = levels)
    set.seed(cfg$analysis$seed + 10000 + i)
    mtry <- if (identical(cfg$model$mtry, "sqrt")) max(1, floor(sqrt(ncol(prep$train_x)))) else as.integer(cfg$model$mtry)
    fit <- randomForest::randomForest(
      prep$train_x, train_y, ntree = ntree, mtry = mtry,
      nodesize = cfg$model$min_node_size, replace = isTRUE(cfg$model$replace),
      sampsize = max(1L, floor(nrow(prep$train_x) * cfg$model$sample_fraction))
    )
    predicted <- predict(fit, prep$assessment_x)
    probability <- as.data.frame(predict(fit, prep$assessment_x, type = "prob"))
    for (lv in setdiff(levels, names(probability))) probability[[lv]] <- 0
    probability <- probability[, levels, drop = FALSE]
    pred[[i]] <- tibble(resample_id = sp$id, !!sid := prep$assessment_ids,
      heldout_subject = test[[sub]][match(prep$assessment_ids, test[[sid]])], observed = as.character(test_y), predicted = as.character(predicted))
    for (lv in levels) pred[[i]][[paste0("prob_", lv)]] <- probability[[lv]]
    met[[i]] <- bind_cols(
      tibble(resample_id = sp$id, repeat_id = sp$repeat_id, fold = sp$fold, n_test = length(test_y)),
      classification_metrics(test_y, predicted, levels, probability)$metrics
    )
    settings[[i]] <- tibble(
      resample_id = sp$id, repeat_id = sp$repeat_id, fold = sp$fold,
      train_subjects = paste(sort(unique(as.character(train[[sub]]))), collapse = ";"),
      assessment_subjects = paste(sort(unique(as.character(test[[sub]]))), collapse = ";"),
      n_train = nrow(prep$train_x), n_test = nrow(prep$assessment_x),
      n_features = length(prep$features), ntree = ntree, mtry = mtry
    )
  }
  list(predictions = bind_rows(pred), metrics = bind_rows(met), settings = bind_rows(settings))
}

if (sys.nframe() == 0L) {
  cfg <- read_validation_config(); run_dir <- Sys.getenv("RF_VALIDATION_RUN_DIR", ""); if (!nzchar(run_dir)) run_dir <- make_run_directory(cfg)
  dir.create(file.path(run_dir, "predictions"), recursive = TRUE, showWarnings = FALSE)
  z <- run_grouped_cv(read_long_data(cfg), cfg)
  write.csv(z$predictions, file.path(run_dir, "predictions", "grouped_cv_predictions.csv"), row.names = FALSE)
  write.csv(z$metrics, file.path(run_dir, "predictions", "grouped_cv_fold_metrics.csv"), row.names = FALSE)
  write.csv(z$settings, file.path(run_dir, "predictions", "grouped_cv_fold_settings.csv"), row.names = FALSE)
}
