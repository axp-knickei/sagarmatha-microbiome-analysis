repo_root <- normalizePath(file.path("..", ".."))
old_wd <- setwd(repo_root)
on.exit(setwd(old_wd), add = TRUE)
source("analysis/random_forest_subject_validation/R/common.R")
source("analysis/random_forest_subject_validation/02_define_grouped_resampling.R")
source("analysis/random_forest_subject_validation/03_run_baseline_models.R")
source("analysis/random_forest_subject_validation/04_run_loso_random_forest.R")

# Purpose: Exercise the unchanged subject-aware RF implementation using only
# fictional long-format data generated here; never load research inputs.
# Inputs: Approved config structure and this artificial fixture, not study data.
# Outputs: testthat assertions; runtime fixtures/results are temporary.
# Dependencies: testthat, randomForest and the configured manipulation packages.
# Invocation: Rscript tests/public_release/run-public-tests.R (public tree).
source("analysis/random_forest_subject_validation/05_run_repeated_grouped_cv.R")
fixture_dir <- tempfile("synthetic-rf-")
dir.create(fixture_dir)
withr::defer(unlink(fixture_dir, recursive = TRUE), envir = environment())
cfg <- read_validation_config()
# Only fixture locations, fictional labels and the test seed are overridden.
# Production configuration, preprocessing/model parameters and code are untouched.
cfg$analysis$mode <- "smoke"
cfg$analysis$seed <- 137L
cfg$data$altitude_levels <- c("11", "22", "33", "44")
cfg$data$input <- file.path(fixture_dir, "fictional_long_data.rds")
cfg$data$smoke_input_fallback <- NULL
cfg$data$existing_run <- file.path(fixture_dir, "fictional_previous_run")
cfg$results$output_root <- file.path(fixture_dir, "runs")
validate_validation_config(cfg)
# Fixture sizes below are software-test choices, unrelated to the study cohort.
fixture_subjects <- paste0("SYNTHETIC_SUBJECT_", LETTERS[1:6])
fixture_mags <- paste0("SYNTHETIC_MAG_", seq_len(8))
fixture <- expand.grid(subject_index = seq_along(fixture_subjects),
                       class_index = seq_along(cfg$data$altitude_levels),
                       visit_index = 1:2, feature_index = seq_along(fixture_mags))
artificial <- data.frame(
  sampleID = paste0("SYNTHETIC_SAMPLE_", fixture$subject_index, "_",
                    fixture$class_index, "_", fixture$visit_index),
  Subject = fixture_subjects[fixture$subject_index],
  MAG = fixture_mags[fixture$feature_index],
  Y = (1 + fixture$subject_index + fixture$visit_index +
       fixture$class_index * fixture$feature_index +
       (fixture$subject_index * fixture$feature_index) %% 7) / 1000,
  altitude = cfg$data$altitude_levels[fixture$class_index],
  Timepoint = paste0("SYNTHETIC_VISIT_", fixture$class_index, "_", fixture$visit_index),
  date_converted = as.Date("2040-01-01") + fixture$class_index * 10 + fixture$visit_index,
  days_from_baseline_subject = fixture$class_index * 10 + fixture$visit_index,
  days_from_baseline_global = fixture$class_index * 10 + fixture$visit_index)
# A deliberately constant feature tests training-derived filtering.
artificial$Y[artificial$MAG == tail(fixture_mags, 1)] <- 0.001
saveRDS(artificial, cfg$data$input)
rm(artificial, fixture)
dir.create(cfg$data$existing_run)
sentinel <- file.path(cfg$data$existing_run, "synthetic_sentinel.txt")
writeLines("Fictional pre-existing output; must not be overwritten", sentinel)
sentinel_before <- readBin(sentinel, "raw", n = file.info(sentinel)$size)
d <- read_long_data(cfg)
meta <- sample_metadata(d, cfg)
sid <- cfg$data$sample_id; sub <- cfg$data$subject_id; outcome <- cfg$data$outcome
levels <- unlist(cfg$data$altitude_levels)

testthat::test_that("LOSO holds out every subject once without overlap", {
  splits <- make_loso_splits(meta, sub)
  testthat::expect_equal(length(splits), dplyr::n_distinct(meta[[sub]]))
  testthat::expect_setequal(vapply(splits, `[[`, character(1), "holdout_subject"), unique(as.character(meta[[sub]])))
  for (sp in splits) testthat::expect_true(assert_disjoint_subjects(meta[[sub]][sp$analysis], meta[[sub]][sp$assessment]))
})

testthat::test_that("grouped resampling is deterministic and subject-disjoint", {
  a <- make_grouped_folds(meta, sub, outcome, levels, 5, 2, cfg$analysis$seed)
  b <- make_grouped_folds(meta, sub, outcome, levels, 5, 2, cfg$analysis$seed)
  testthat::expect_identical(a, b)
  for (sp in a) testthat::expect_true(assert_disjoint_subjects(meta[[sub]][sp$analysis], meta[[sub]][sp$assessment]))
})

testthat::test_that("feature fitting uses training observations only", {
  sp <- make_loso_splits(meta, sub)[[1]]
  train <- meta[sp$analysis, ]; test <- meta[sp$assessment, ]
  prep <- prepare_training_features(d[d[[sid]] %in% train[[sid]], ], d[d[[sid]] %in% test[[sid]], ], cfg, top_n = 10)
  testthat::expect_length(intersect(prep$preprocessing_fit_samples, test[[sid]]), 0)
  testthat::expect_lte(length(prep$features), 10)
  testthat::expect_false(any(c(sid, sub, outcome, cfg$data$timepoint, cfg$data$sampling_date) %in% colnames(prep$train_x)))
})

testthat::test_that("metrics retain all expected classes and totals", {
  z <- classification_metrics(c("11", "33"), c("33", "33"), levels)
  testthat::expect_equal(z$per_class$class, levels)
  testthat::expect_equal(sum(z$confusion$Freq), 2)
  testthat::expect_true(is.na(z$per_class$recall[z$per_class$class == "22"]))
  testthat::expect_true(is.na(z$metrics$balanced_accuracy_all_classes))
  testthat::expect_equal(z$metrics$n_observed_classes, 2)
})

testthat::test_that("custom metrics agree with hand calculations", {
  observed <- c("11", "11", "22", "22", "33", "33", "44", "44")
  predicted <- c("11", "22", "22", "22", "33", "44", "44", "11")
  z <- classification_metrics(observed, predicted, levels)
  testthat::expect_equal(z$metrics$accuracy, 5 / 8)
  testthat::expect_equal(z$metrics$balanced_accuracy, 0.625)
  testthat::expect_equal(z$metrics$macro_f1, mean(c(0.5, 0.8, 2 / 3, 0.5)))
  testthat::expect_equal(z$metrics$ordinal_mae, 5 / 8)
  testthat::expect_equal(z$metrics$ordinal_rmse, sqrt(11 / 8))
  testthat::expect_equal(z$metrics$within_one_accuracy, 7 / 8)
  testthat::expect_true(all(z$per_class$precision >= 0 & z$per_class$precision <= 1))
})

testthat::test_that("log loss validates probability order and normalization", {
  observed <- levels
  predicted <- levels
  probability <- matrix(0.1, 4, 4, dimnames = list(NULL, levels))
  diag(probability) <- 0.7
  z <- classification_metrics(observed, predicted, levels, probability)
  testthat::expect_equal(z$metrics$log_loss, -log(0.7))
  bad <- probability
  bad[1, ] <- 0.2
  testthat::expect_error(classification_metrics(observed, predicted, levels, bad), "sum to one")
  testthat::expect_error(classification_metrics(observed, predicted, levels, probability[, rev(levels), drop = FALSE][, -1]), "do not match")
})

testthat::test_that("baseline predictions are training-derived and reproducible", {
  a <- run_loso_baselines(meta, cfg); b <- run_loso_baselines(meta, cfg)
  testthat::expect_identical(a, b)
  repeats <- cfg$analysis$baseline_repeats[[cfg$analysis$mode]]
  testthat::expect_equal(nrow(a), (repeats + 1) * nrow(meta))
  testthat::expect_false(any(is.na(a$predicted)))
  for (fold in unique(a$fold_id)) {
    heldout <- sub("^LOSO_", "", fold)
    training <- meta[as.character(meta[[sub]]) != heldout, ]
    expected <- names(which.max(table(factor(as.character(training[[outcome]]), levels = levels))))
    actual <- unique(a$predicted[a$fold_id == fold & a$model == "majority"])
    testthat::expect_identical(actual, expected)
  }
  probability_columns <- paste0("prob_", levels)
  testthat::expect_true(all(abs(rowSums(a[, probability_columns]) - 1) < 1e-12))
})

testthat::test_that("class weights and sampling use only supplied training outcomes", {
  y <- factor(c("11", "11", "22", "33", "33", "33", "44"), levels = levels)
  weighted <- derive_class_strategy(y, levels, "class_weights")
  testthat::expect_equal(unname(weighted$classwt), 7 / (4 * c(2, 1, 3, 1)))
  balanced <- derive_class_strategy(y, levels, "balanced_sampling")
  testthat::expect_equal(unname(balanced$sampsize), rep(1, 4))
  testthat::expect_null(balanced$classwt)
})

testthat::test_that("configuration validation rejects unsupported and invalid settings", {
  invalid <- cfg
  invalid$analysis$mode <- "invalid"
  testthat::expect_error(validate_validation_config(invalid), "analysis.mode")
  invalid <- cfg
  invalid$outcomes$baseline_change <- TRUE
  testthat::expect_error(validate_validation_config(invalid), "disabled outcome")
  invalid <- cfg
  invalid$model$sample_fraction <- 0
  testthat::expect_error(validate_validation_config(invalid), "sample_fraction")
})

testthat::test_that("input validation rejects duplicate and incomplete sample-feature grids", {
  duplicate <- bind_rows(d, d[1, ])
  testthat::expect_error(validate_long_data(duplicate, cfg), "duplicate sampleID-MAG")
  incomplete <- d[-1, ]
  testthat::expect_error(validate_long_data(incomplete, cfg), "incomplete")
})

testthat::test_that("assessment values cannot change training feature fit", {
  sp <- make_loso_splits(meta, sub)[[1]]
  train <- d[d[[sid]] %in% meta[[sid]][sp$analysis], ]
  assessment <- d[d[[sid]] %in% meta[[sid]][sp$assessment], ]
  a <- prepare_training_features(train, assessment, cfg, top_n = 20)
  modified <- assessment
  modified[[cfg$data$abundance]] <- modified[[cfg$data$abundance]] * 1e6
  b <- prepare_training_features(train, modified, cfg, top_n = 20)
  testthat::expect_identical(a$features, b$features)
  testthat::expect_identical(a$imputation_medians, b$imputation_medians)
  testthat::expect_identical(a$train_x, b$train_x)
})

testthat::test_that("run directories refuse collisions", {
  temporary <- tempfile("rf-config-")
  dir.create(temporary)
  local_cfg <- cfg
  local_cfg$results$output_root <- temporary
  first <- make_run_directory(local_cfg, "collision")
  testthat::expect_true(dir.exists(first))
  testthat::expect_error(make_run_directory(local_cfg, "collision"), "Refusing to reuse")
})

testthat::test_that("synthetic previous output remains unchanged by validation", {
  testthat::expect_identical(readBin(sentinel, "raw", n = file.info(sentinel)$size),
                            sentinel_before)
})

testthat::test_that("missing files, fields and invalid values fail clearly", {
  invalid <- cfg; invalid$data$input <- file.path(fixture_dir, "absent.rds")
  testthat::expect_error(read_long_data(invalid), "unavailable")
  invalid$data$input <- file.path(fixture_dir, "missing_fields.rds")
  saveRDS(d[, setdiff(names(d), cfg$data$subject_id)], invalid$data$input)
  testthat::expect_error(read_long_data(invalid), "Input is missing")
  bad <- d; bad$Subject[1] <- ""
  testthat::expect_error(validate_long_data(bad, cfg), "blank subject")
  bad <- d; bad$Y[1] <- -1
  testthat::expect_error(validate_long_data(bad, cfg), "negative")
  bad <- d; bad$Y[1] <- NA_real_
  testthat::expect_error(validate_long_data(bad, cfg), "infinite")
  bad <- d; bad$Y[1] <- Inf
  testthat::expect_error(validate_long_data(bad, cfg), "infinite")
  bad <- d; bad$altitude[1] <- "999"
  testthat::expect_error(validate_long_data(bad, cfg), "classes do not match")
  bad <- d; bad$Timepoint[1] <- "INCONSISTENT_SYNTHETIC_VISIT"
  testthat::expect_error(validate_long_data(bad, cfg), "inconsistent")
  testthat::expect_error(assert_disjoint_subjects("SYNTHETIC_A", "SYNTHETIC_A"), "leakage")
  testthat::expect_error(make_loso_splits(meta[meta[[sub]] == fixture_subjects[1], ], sub),
                        "at least two subjects")
  testthat::expect_error(classification_metrics(levels, rep("UNKNOWN", length(levels)), levels),
                        "unknown classes")
})

testthat::test_that("RF fitting is deterministic and emits complete subject holdouts", {
  a <- fit_loso_rf(d, cfg); b <- fit_loso_rf(d, cfg)
  testthat::expect_identical(a$predictions, b$predictions)
  testthat::expect_identical(a$selections, b$selections)
  testthat::expect_identical(a$importance, b$importance)
  testthat::expect_equal(nrow(a$predictions), nrow(meta))
  testthat::expect_setequal(a$predictions[[sid]], meta[[sid]])
  testthat::expect_false(anyDuplicated(a$predictions[[sid]]) > 0)
  testthat::expect_setequal(a$predictions$heldout_subject, meta[[sub]])
  testthat::expect_true(all(a$predictions$predicted %in% levels))
  testthat::expect_equal(a$predictions$heldout_subject,
                        meta[[sub]][match(a$predictions[[sid]], meta[[sid]])])
  probabilities <- as.matrix(a$predictions[paste0("prob_", levels)])
  testthat::expect_true(all(is.finite(probabilities) & probabilities >= 0 & probabilities <= 1))
  testthat::expect_equal(rowSums(probabilities), rep(1, nrow(meta)))
  testthat::expect_false(tail(fixture_mags, 1) %in% a$selections$feature)
  for (i in seq_len(nrow(a$training_settings))) {
    z <- a$training_settings[i, ]
    testthat::expect_length(intersect(strsplit(z$train_subjects, ";")[[1]],
                                     strsplit(z$assessment_subjects, ";")[[1]]), 0)
    testthat::expect_equal(z$ntree, cfg$model$smoke_trees)
  }
  scored <- metrics_from_prediction_table(a$predictions, levels)
  testthat::expect_equal(scored$metrics$accuracy,
                        mean(a$predictions$observed == a$predictions$predicted))
  testthat::expect_equal(sum(scored$confusion$Freq), nrow(meta))
  # Exercise each existing class-balancing mode without prescribing performance.
  for (strategy in c("class_weights", "balanced_sampling", "combined")) {
    z <- fit_loso_rf(d, cfg, strategy = strategy)
    testthat::expect_equal(nrow(z$predictions), nrow(meta))
    testthat::expect_true(all(z$predictions$strategy == strategy))
    testthat::expect_true(all(z$predictions$predicted %in% levels))
  }
})

testthat::test_that("grouped RF prediction and split membership are reproducible", {
  a <- run_grouped_cv(d, cfg); b <- run_grouped_cv(d, cfg)
  testthat::expect_identical(a, b)
  testthat::expect_equal(nrow(a$predictions), nrow(meta))
  testthat::expect_setequal(a$predictions[[sid]], meta[[sid]])
  testthat::expect_equal(rowSums(a$predictions[paste0("prob_", levels)]), rep(1, nrow(meta)))
  for (i in seq_len(nrow(a$settings))) {
    z <- a$settings[i, ]
    train_subjects <- strsplit(z$train_subjects, ";")[[1]]
    assessment_subjects <- strsplit(z$assessment_subjects, ";")[[1]]
    testthat::expect_length(intersect(train_subjects, assessment_subjects), 0)
    fold_predictions <- a$predictions[a$predictions$resample_id == z$resample_id, ]
    testthat::expect_setequal(fold_predictions$heldout_subject, assessment_subjects)
  }
})
