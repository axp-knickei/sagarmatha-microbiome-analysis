#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(yaml)
  library(jsonlite)
  library(randomForest)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (is.atomic(x) && length(x) == 1 && (is.na(x) || identical(x, "")))) y else x
}

usage <- function() {
  cat(
    "Usage: Rscript scripts/run_random_forest_model.R [options]\n\n",
    "Options:\n",
    "  --config PATH                 YAML config. Default: config/random_forest_config.yml\n",
    "  --input PATH                  Processed RDS or CSV input.\n",
    "  --output-dir PATH             Base output directory. Default: results/random_forest\n",
    "  --run-name NAME               Run folder name. Default: timestamped name.\n",
    "  --response-variable COLUMN    Metadata column to predict. Default from config.\n",
    "  --sample-id-column COLUMN     Sample identifier. Default: sampleID\n",
    "  --subject-id-column COLUMN    Subject/group identifier. Default: Subject\n",
    "  --feature-id-column COLUMN    Long-format feature column. Default: MAG\n",
    "  --abundance-column COLUMN     Long-format abundance column. Default: Y\n",
    "  --task auto|classification|regression\n",
    "  --seed N\n",
    "  --n-trees N\n",
    "  --threads N                   Parallel sub-forests if doParallel/foreach are available.\n",
    "  --cv-folds N\n",
    "  --test-fraction X\n",
    "  --min-prevalence X            Keep features nonzero in at least this sample fraction.\n",
    "  --min-mean-abundance X\n",
    "  --max-missing-fraction X\n",
    "  --include-covariates          Include configured covariates as predictors.\n",
    "  --covariates a,b,c            Covariates to include if --include-covariates is set.\n",
    "  --mtry-grid a,b,c             Integers or aliases: sqrt, third, all.\n",
    "  --min-node-size-grid a,b,c     nodesize values. Default: 1,5,10\n",
    "  --no-permutation              Skip permutation importance.\n",
    "  --test-mode                   Small fast run for validation.\n",
    "  --help                        Show this help.\n",
    sep = ""
  )
}

parse_args <- function(args) {
  parsed <- list(config = "config/random_forest_config.yml", no_permutation = FALSE, test_mode = FALSE)
  key_map <- c(
    "--config" = "config",
    "--input" = "input",
    "--output-dir" = "output_dir",
    "--run-name" = "run_name",
    "--response-variable" = "response",
    "--sample-id-column" = "sample_id",
    "--subject-id-column" = "subject_id",
    "--feature-id-column" = "feature_id",
    "--abundance-column" = "abundance_value",
    "--task" = "task",
    "--seed" = "random_seed",
    "--n-trees" = "n_trees",
    "--threads" = "threads",
    "--cv-folds" = "cv_folds",
    "--test-fraction" = "test_fraction",
    "--min-prevalence" = "min_prevalence",
    "--min-mean-abundance" = "min_mean_abundance",
    "--max-missing-fraction" = "max_missing_fraction",
    "--covariates" = "covariates",
    "--mtry-grid" = "mtry_grid",
    "--min-node-size-grid" = "min_node_size_grid"
  )
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      usage()
      quit(status = 0)
    } else if (arg == "--include-covariates") {
      parsed$include_covariates <- TRUE
      i <- i + 1
    } else if (arg == "--no-permutation") {
      parsed$no_permutation <- TRUE
      i <- i + 1
    } else if (arg == "--test-mode") {
      parsed$test_mode <- TRUE
      i <- i + 1
    } else if (arg %in% names(key_map)) {
      if (i == length(args)) stop("Missing value for argument: ", arg, call. = FALSE)
      parsed[[key_map[[arg]]]] <- args[[i + 1]]
      i <- i + 2
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  parsed
}

as_bool <- function(x) isTRUE(x) || identical(tolower(as.character(x)), "true")

read_input <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(readRDS(path))
  if (ext %in% c("csv", "tsv", "txt")) {
    sep <- if (ext == "csv") "," else "\t"
    return(utils::read.table(path, sep = sep, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE))
  }
  stop("Unsupported input extension: ", ext, ". Use .rds, .csv, .tsv, or .txt.", call. = FALSE)
}

make_run_dir <- function(base_dir, run_name) {
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  label <- run_name %||% paste0("rf_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  candidate <- file.path(base_dir, label)
  if (!dir.exists(candidate)) {
    dir.create(candidate, recursive = TRUE)
    return(candidate)
  }
  for (i in seq_len(999)) {
    suffixed <- paste0(candidate, "_", sprintf("%03d", i))
    if (!dir.exists(suffixed)) {
      dir.create(suffixed, recursive = TRUE)
      return(suffixed)
    }
  }
  stop("Could not create unique output directory under: ", base_dir, call. = FALSE)
}

parse_grid <- function(value, p) {
  raw <- trimws(unlist(strsplit(as.character(value), ",", fixed = TRUE)))
  raw <- raw[nzchar(raw)]
  out <- integer(0)
  for (item in raw) {
    item_lower <- tolower(item)
    val <- switch(
      item_lower,
      sqrt = max(1L, floor(sqrt(p))),
      third = max(1L, floor(p / 3)),
      all = p,
      suppressWarnings(as.integer(item))
    )
    if (is.na(val) || val < 1) stop("Invalid grid value: ", item, call. = FALSE)
    out <- c(out, min(p, val))
  }
  sort(unique(out))
}

assign_group_folds <- function(groups, k, seed) {
  set.seed(seed)
  unique_groups <- sample(unique(as.character(groups)))
  k <- min(k, length(unique_groups))
  fold_by_group <- setNames(((seq_along(unique_groups) - 1L) %% k) + 1L, unique_groups)
  unname(fold_by_group[as.character(groups)])
}

group_train_test_split <- function(groups, y, test_fraction, seed, task) {
  set.seed(seed)
  unique_groups <- unique(as.character(groups))
  n_test <- max(1L, round(length(unique_groups) * test_fraction))
  n_test <- min(n_test, length(unique_groups) - 1L)
  best <- NULL
  for (attempt in seq_len(1000)) {
    test_groups <- sample(unique_groups, n_test)
    is_test <- as.character(groups) %in% test_groups
    if (task == "classification") {
      train_levels <- unique(as.character(y[!is_test]))
      test_levels <- unique(as.character(y[is_test]))
      if (length(setdiff(test_levels, train_levels)) == 0) return(is_test)
    } else {
      return(is_test)
    }
    best <- is_test
  }
  warning("Could not find grouped split with all test classes represented in train; using best attempted split.", call. = FALSE)
  best
}

classification_metrics <- function(obs, pred, prob = NULL) {
  obs <- factor(obs)
  pred <- factor(pred, levels = levels(obs))
  cm <- table(observed = obs, predicted = pred)
  classes <- levels(obs)
  per_class <- lapply(classes, function(cls) {
    tp <- cm[cls, cls]
    fp <- sum(cm[, cls]) - tp
    fn <- sum(cm[cls, ]) - tp
    precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    f1 <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
    data.frame(class = cls, precision = precision, recall = recall, f1 = f1, support = sum(obs == cls))
  }) |> bind_rows()
  accuracy <- sum(diag(cm)) / sum(cm)
  balanced_accuracy <- mean(per_class$recall, na.rm = TRUE)
  macro_precision <- mean(per_class$precision, na.rm = TRUE)
  macro_recall <- mean(per_class$recall, na.rm = TRUE)
  macro_f1 <- mean(per_class$f1, na.rm = TRUE)
  auc <- NA_real_
  if (!is.null(prob) && requireNamespace("pROC", quietly = TRUE)) {
    auc <- tryCatch({
      if (length(classes) == 2) {
        as.numeric(pROC::auc(pROC::roc(obs, prob[, classes[[2]]], levels = classes, quiet = TRUE)))
      } else {
        as.numeric(pROC::auc(pROC::multiclass.roc(obs, as.matrix(prob[, classes, drop = FALSE]), levels = classes, quiet = TRUE)))
      }
    }, error = function(e) NA_real_)
  }
  list(
    metrics = data.frame(
      accuracy = accuracy,
      balanced_accuracy = balanced_accuracy,
      precision_macro = macro_precision,
      recall_macro = macro_recall,
      f1_macro = macro_f1,
      roc_auc = auc
    ),
    per_class = per_class,
    confusion_matrix = as.data.frame.matrix(cm)
  )
}

regression_metrics <- function(obs, pred) {
  residual <- obs - pred
  ss_res <- sum(residual^2)
  ss_tot <- sum((obs - mean(obs))^2)
  data.frame(
    r2 = if (ss_tot == 0) NA_real_ else 1 - ss_res / ss_tot,
    rmse = sqrt(mean(residual^2)),
    mae = mean(abs(residual))
  )
}

fit_rf <- function(x, y, task, ntree, mtry, nodesize, seed) {
  set.seed(seed)
  if (task == "classification") y <- factor(y) else y <- as.numeric(y)
  randomForest(x = x, y = y, ntree = ntree, mtry = mtry, nodesize = nodesize, importance = TRUE)
}

fit_rf_parallel <- function(x, y, task, ntree, mtry, nodesize, seed, threads) {
  threads <- max(1L, min(as.integer(threads), as.integer(ntree)))
  if (threads <= 1L || !requireNamespace("doParallel", quietly = TRUE) || !requireNamespace("foreach", quietly = TRUE)) {
    return(fit_rf(x, y, task, ntree, mtry, nodesize, seed))
  }
  `%dopar%` <- foreach::`%dopar%`

  trees_per_worker <- rep(floor(ntree / threads), threads)
  trees_per_worker[seq_len(ntree %% threads)] <- trees_per_worker[seq_len(ntree %% threads)] + 1L
  trees_per_worker <- trees_per_worker[trees_per_worker > 0L]
  cluster <- tryCatch(
    parallel::makeCluster(length(trees_per_worker)),
    error = function(e) {
      warning("Parallel worker startup failed; falling back to serial Random Forest: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(cluster)) {
    return(fit_rf(x, y, task, ntree, mtry, nodesize, seed))
  }
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  doParallel::registerDoParallel(cluster)

  forests <- foreach::foreach(
    i = seq_along(trees_per_worker),
    .packages = "randomForest"
  ) %dopar% {
    set.seed(seed + i)
    if (task == "classification") {
      randomForest::randomForest(
        x = x,
        y = factor(y),
        ntree = trees_per_worker[[i]],
        mtry = mtry,
        nodesize = nodesize,
        importance = TRUE
      )
    } else {
      randomForest::randomForest(
        x = x,
        y = as.numeric(y),
        ntree = trees_per_worker[[i]],
        mtry = mtry,
        nodesize = nodesize,
        importance = TRUE
      )
    }
  }
  do.call(randomForest::combine, forests)
}

evaluate_model <- function(model, x, y, task) {
  if (task == "classification") {
    pred <- predict(model, x, type = "response")
    prob <- predict(model, x, type = "prob")
    classification_metrics(y, pred, prob)
  } else {
    pred <- predict(model, x)
    list(metrics = regression_metrics(y, pred), predictions = pred)
  }
}

permutation_importance <- function(model, x, y, task, baseline_metric, seed, max_features = 100L) {
  set.seed(seed)
  features <- colnames(x)
  if (length(features) > max_features) {
    native <- importance(model)
    score_col <- intersect(c("MeanDecreaseAccuracy", "%IncMSE", "MeanDecreaseGini", "IncNodePurity"), colnames(native))[[1]]
    ranked <- rownames(native)[order(native[, score_col], decreasing = TRUE)]
    features <- ranked[seq_len(max_features)]
  }
  bind_rows(lapply(features, function(feature) {
    x_perm <- x
    x_perm[[feature]] <- sample(x_perm[[feature]])
    if (task == "classification") {
      pred <- predict(model, x_perm, type = "response")
      metric <- classification_metrics(y, pred)$metrics$balanced_accuracy
      data.frame(feature = feature, metric = metric, importance = baseline_metric - metric)
    } else {
      pred <- predict(model, x_perm)
      metric <- regression_metrics(y, pred)$rmse
      data.frame(feature = feature, metric = metric, importance = metric - baseline_metric)
    }
  }))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
cfg <- if (file.exists(args$config)) yaml::read_yaml(args$config) else list()

input_path <- args$input %||% cfg$data$input %||% "data/processed/cleaned_data.rds"
sample_col <- args$sample_id %||% cfg$data$sample_id %||% "sampleID"
subject_col <- args$subject_id %||% cfg$data$subject_id %||% "Subject"
feature_col <- args$feature_id %||% cfg$data$feature_id %||% "MAG"
abundance_col <- args$abundance_value %||% cfg$data$abundance_value %||% "Y"
response_col <- args$response %||% cfg$data$response %||% "altitude"
task <- tolower(args$task %||% cfg$model$task %||% "auto")
seed <- as.integer(args$random_seed %||% cfg$model$random_seed %||% 20260707)
ntree <- as.integer(args$n_trees %||% cfg$model$n_trees %||% 1000)
threads <- as.integer(args$threads %||% cfg$model$threads %||% 1)
cv_folds <- as.integer(args$cv_folds %||% cfg$model$cv_folds %||% 5)
test_fraction <- as.numeric(args$test_fraction %||% cfg$model$test_fraction %||% 0.25)
min_prevalence <- as.numeric(args$min_prevalence %||% cfg$features$min_prevalence %||% 0.10)
min_mean_abundance <- as.numeric(args$min_mean_abundance %||% cfg$features$min_mean_abundance %||% 0)
max_missing_fraction <- as.numeric(args$max_missing_fraction %||% cfg$features$max_missing_fraction %||% 0.20)
include_covariates <- as_bool(args$include_covariates %||% cfg$features$include_covariates %||% FALSE)
covariates <- args$covariates %||% paste(unlist(cfg$features$covariates %||% character(0)), collapse = ",")
covariates <- trimws(unlist(strsplit(as.character(covariates), ",", fixed = TRUE)))
covariates <- covariates[nzchar(covariates)]
transform_name <- cfg$features$transform %||% "log10_pseudocount"
pseudocount <- as.numeric(cfg$features$pseudocount %||% 0)
base_output_dir <- args$output_dir %||% cfg$results$output_dir %||% "results/random_forest"
run_name <- args$run_name %||% cfg$results$run_name %||% ""

if (isTRUE(args$test_mode)) {
  ntree <- min(ntree, 100L)
  cv_folds <- min(cv_folds, 3L)
}

set.seed(seed)
run_dir <- make_run_dir(base_output_dir, run_name)
log_file <- file.path(run_dir, "run.log")
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(zz)
}, add = TRUE)

cat("=== Random Forest workflow started:", format(Sys.time()), "===\n")
cat("Input:", input_path, "\n")
cat("Output:", run_dir, "\n")
cat("Response:", response_col, "\n")
cat("Subject/group column:", subject_col, "\n")

df <- read_input(input_path)
required <- c(sample_col, subject_col, response_col)
missing_required <- setdiff(required, names(df))
if (length(missing_required) > 0) stop("Input is missing required columns: ", paste(missing_required, collapse = ", "), call. = FALSE)

if (all(c(feature_col, abundance_col) %in% names(df))) {
  meta_cols <- setdiff(names(df), c(feature_col, abundance_col))
  metadata <- df |> distinct(.data[[sample_col]], .keep_all = TRUE) |> select(all_of(meta_cols))
  abundance <- df |> select(all_of(c(sample_col, feature_col, abundance_col))) |>
    mutate(.abundance = as.numeric(.data[[abundance_col]])) |>
    select(all_of(sample_col), all_of(feature_col), .abundance)
  wide <- abundance |> tidyr::pivot_wider(names_from = all_of(feature_col), values_from = .abundance, values_fill = 0)
  model_df <- metadata |> left_join(wide, by = sample_col)
  feature_names <- setdiff(names(wide), sample_col)
} else {
  model_df <- df
  non_feature <- unique(c(sample_col, subject_col, response_col, covariates, "group", "Timepoint", "date", "date_converted"))
  feature_names <- setdiff(names(model_df), non_feature)
}

duplicate_samples <- model_df[[sample_col]][duplicated(model_df[[sample_col]])]
if (length(duplicate_samples) > 0) stop("Duplicate sample IDs after reshaping: ", paste(head(unique(duplicate_samples), 10), collapse = ", "), call. = FALSE)

feature_missing <- vapply(model_df[feature_names], function(x) mean(is.na(x)), numeric(1))
feature_mean <- vapply(model_df[feature_names], function(x) mean(as.numeric(x), na.rm = TRUE), numeric(1))
feature_prevalence <- vapply(model_df[feature_names], function(x) mean(as.numeric(x) > 0, na.rm = TRUE), numeric(1))
keep_features <- names(feature_missing)[feature_missing <= max_missing_fraction & feature_mean >= min_mean_abundance & feature_prevalence >= min_prevalence]
if (length(keep_features) < 2) stop("Fewer than two features remain after filtering; relax feature filters.", call. = FALSE)

feature_matrix <- model_df[keep_features]
feature_matrix[] <- lapply(feature_matrix, function(x) as.numeric(x))
for (name in names(feature_matrix)) {
  if (anyNA(feature_matrix[[name]])) feature_matrix[[name]][is.na(feature_matrix[[name]])] <- median(feature_matrix[[name]], na.rm = TRUE)
}
if (identical(transform_name, "log10_pseudocount")) {
  feature_matrix[] <- lapply(feature_matrix, function(x) log10(x + pseudocount))
} else if (!identical(transform_name, "none")) {
  stop("Unsupported feature transform: ", transform_name, call. = FALSE)
}

predictor_df <- feature_matrix
used_covariates <- intersect(covariates, names(model_df))
if (include_covariates && length(used_covariates) > 0) {
  cov_df <- model_df[used_covariates]
  predictor_df <- cbind(predictor_df, model.matrix(~ . - 1, data = cov_df))
}
predictor_df <- as.data.frame(predictor_df, check.names = TRUE)

model_rows <- complete.cases(model_df[, c(sample_col, subject_col, response_col), drop = FALSE])
predictor_df <- predictor_df[model_rows, , drop = FALSE]
model_df <- model_df[model_rows, , drop = FALSE]
y <- model_df[[response_col]]
if (task == "auto") {
  task <- if (is.numeric(y) && length(unique(y)) > 10) "regression" else "classification"
}
if (!task %in% c("classification", "regression")) stop("--task must be auto, classification, or regression.", call. = FALSE)
if (task == "classification") y <- factor(y)

test_idx <- group_train_test_split(model_df[[subject_col]], y, test_fraction, seed, task)
train_x <- predictor_df[!test_idx, , drop = FALSE]
test_x <- predictor_df[test_idx, , drop = FALSE]
train_y <- y[!test_idx]
test_y <- y[test_idx]
train_groups <- model_df[[subject_col]][!test_idx]

mtry_grid <- parse_grid(args$mtry_grid %||% cfg$model$mtry_grid %||% "sqrt,third,all", ncol(train_x))
nodesize_grid <- parse_grid(args$min_node_size_grid %||% cfg$model$min_node_size_grid %||% "1,5,10", ncol(train_x))
tuning_grid <- expand.grid(mtry = mtry_grid, nodesize = nodesize_grid)
fold_id <- assign_group_folds(train_groups, cv_folds, seed)

cat("Samples:", nrow(model_df), "\n")
cat("Subjects:", length(unique(model_df[[subject_col]])), "\n")
cat("Features before filtering:", length(feature_names), "\n")
cat("Features after filtering:", length(keep_features), "\n")
cat("Predictors after covariate encoding:", ncol(predictor_df), "\n")
cat("Task:", task, "\n")
cat("Train samples:", nrow(train_x), " Test samples:", nrow(test_x), "\n")
cat("Train subjects:", length(unique(model_df[[subject_col]][!test_idx])), " Test subjects:", length(unique(model_df[[subject_col]][test_idx])), "\n")

cv_results <- bind_rows(lapply(seq_len(nrow(tuning_grid)), function(i) {
  params <- tuning_grid[i, ]
  fold_metrics <- bind_rows(lapply(sort(unique(fold_id)), function(fold) {
    val <- fold_id == fold
    model <- fit_rf_parallel(
      train_x[!val, , drop = FALSE],
      train_y[!val],
      task,
      ntree,
      params$mtry,
      params$nodesize,
      seed + i + fold,
      threads
    )
    ev <- evaluate_model(model, train_x[val, , drop = FALSE], train_y[val], task)
    cbind(data.frame(fold = fold), ev$metrics)
  }))
  summary_metric <- if (task == "classification") mean(fold_metrics$balanced_accuracy, na.rm = TRUE) else mean(fold_metrics$rmse, na.rm = TRUE)
  cbind(params, data.frame(mean_metric = summary_metric), summarise(fold_metrics, across(where(is.numeric), ~ mean(.x, na.rm = TRUE))))
}))

if (task == "classification") {
  best_idx <- which.max(cv_results$mean_metric)
} else {
  best_idx <- which.min(cv_results$mean_metric)
}
best <- cv_results[best_idx, ]
cat("Best mtry:", best$mtry, " Best nodesize:", best$nodesize, "\n")

final_model <- fit_rf_parallel(train_x, train_y, task, ntree, best$mtry, best$nodesize, seed, threads)
test_eval <- evaluate_model(final_model, test_x, test_y, task)

if (task == "classification") {
  test_pred <- predict(final_model, test_x, type = "response")
  test_prob <- as.data.frame(predict(final_model, test_x, type = "prob"))
  predictions <- cbind(
    model_df[test_idx, c(sample_col, subject_col, response_col), drop = FALSE],
    predicted = as.character(test_pred),
    test_prob
  )
  utils::write.csv(test_eval$per_class, file.path(run_dir, "per_class_metrics.csv"), row.names = FALSE)
  utils::write.csv(test_eval$confusion_matrix, file.path(run_dir, "confusion_matrix.csv"), row.names = TRUE)
} else {
  test_pred <- predict(final_model, test_x)
  predictions <- cbind(
    model_df[test_idx, c(sample_col, subject_col, response_col), drop = FALSE],
    predicted = test_pred,
    residual = as.numeric(test_y) - test_pred
  )
}

importance_table <- as.data.frame(importance(final_model))
importance_table$feature <- rownames(importance_table)
importance_table <- importance_table |> relocate(feature)
utils::write.csv(importance_table, file.path(run_dir, "feature_importance.csv"), row.names = FALSE)

perm_table <- data.frame()
if (!isTRUE(args$no_permutation)) {
  baseline <- if (task == "classification") test_eval$metrics$balanced_accuracy else test_eval$metrics$rmse
  perm_table <- permutation_importance(final_model, test_x, test_y, task, baseline, seed)
  utils::write.csv(perm_table, file.path(run_dir, "permutation_importance.csv"), row.names = FALSE)
}

diagnostics <- data.frame(
  sample_id = model_df[[sample_col]],
  subject_id = model_df[[subject_col]],
  response = as.character(model_df[[response_col]]),
  split = ifelse(test_idx, "test", "train")
)
utils::write.csv(diagnostics, file.path(run_dir, "sample_split.csv"), row.names = FALSE)
utils::write.csv(cv_results, file.path(run_dir, "cv_results.csv"), row.names = FALSE)
utils::write.csv(test_eval$metrics, file.path(run_dir, "test_metrics.csv"), row.names = FALSE)
utils::write.csv(predictions, file.path(run_dir, "test_predictions.csv"), row.names = FALSE)
saveRDS(final_model, file.path(run_dir, "random_forest_model.rds"))

native_score_column <- if (task == "classification" && "MeanDecreaseAccuracy" %in% names(importance_table)) {
  "MeanDecreaseAccuracy"
} else if (task == "regression" && "%IncMSE" %in% names(importance_table)) {
  "%IncMSE"
} else if (task == "regression" && "IncNodePurity" %in% names(importance_table)) {
  "IncNodePurity"
} else {
  setdiff(names(importance_table), "feature")[[1]]
}

top_importance <- importance_table |>
  mutate(native_score = .data[[native_score_column]]) |>
  arrange(desc(native_score)) |>
  slice_head(n = 30)
p_importance <- ggplot(top_importance, aes(x = reorder(feature, native_score), y = native_score)) +
  geom_col(fill = "#2f5d50") +
  coord_flip() +
  labs(x = NULL, y = "Native RF importance", title = "Top Random Forest features") +
  theme_minimal(base_size = 11)
ggsave(file.path(run_dir, "feature_importance_top30.png"), plot = p_importance, width = 9, height = 7, dpi = 150)
ggsave(
  file.path(run_dir, "feature_importance_top30.svg"),
  plot = p_importance,
  width = 9,
  height = 7,
  device = svglite::svglite
)

if (task == "classification") {
  cm_long <- as.data.frame(as.table(as.matrix(test_eval$confusion_matrix)))
  names(cm_long) <- c("observed", "predicted", "n")
  p_confusion <- ggplot(cm_long, aes(predicted, observed, fill = n)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), size = 3) +
    scale_fill_gradient(low = "#f4efe6", high = "#7b3f2a") +
    labs(x = "Predicted", y = "Observed", title = "Test confusion matrix") +
    theme_minimal(base_size = 11)
  ggsave(file.path(run_dir, "confusion_matrix.png"), plot = p_confusion, width = 7, height = 6, dpi = 150)
  ggsave(
    file.path(run_dir, "confusion_matrix.svg"),
    plot = p_confusion,
    width = 7,
    height = 6,
    device = svglite::svglite
  )
} else {
  p_regression <- ggplot(predictions, aes(x = .data[[response_col]], y = predicted)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "Observed", y = "Predicted", title = "Test observed vs predicted") +
    theme_minimal(base_size = 11)
  ggsave(file.path(run_dir, "observed_vs_predicted.png"), plot = p_regression, width = 6, height = 5, dpi = 150)
  ggsave(
    file.path(run_dir, "observed_vs_predicted.svg"),
    plot = p_regression,
    width = 6,
    height = 5,
    device = svglite::svglite
  )
}

summary <- list(
  input = input_path,
  output_dir = run_dir,
  response = response_col,
  task = task,
  sample_id_column = sample_col,
  subject_id_column = subject_col,
  n_samples = nrow(model_df),
  n_subjects = length(unique(model_df[[subject_col]])),
  n_features_before_filtering = length(feature_names),
  n_features_after_filtering = length(keep_features),
  n_predictors = ncol(predictor_df),
  include_covariates = include_covariates,
  covariates = used_covariates,
  train_subjects = unique(as.character(model_df[[subject_col]][!test_idx])),
  test_subjects = unique(as.character(model_df[[subject_col]][test_idx])),
  best_parameters = list(mtry = best$mtry, nodesize = best$nodesize, ntree = ntree),
  test_metrics = as.list(test_eval$metrics[1, ]),
  warnings = c(
    "Subject-aware split and grouped CV were used to reduce repeated-measure leakage.",
    "Default configuration predicts altitude from MAG features only; timepoint/covariates are excluded unless requested.",
    "Small subject count means held-out performance can be unstable and should not be treated as confirmatory."
  )
)
writeLines(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), file.path(run_dir, "summary.json"))

cat("=== Random Forest workflow finished:", format(Sys.time()), "===\n")
cat("Summary:", file.path(run_dir, "summary.json"), "\n")
