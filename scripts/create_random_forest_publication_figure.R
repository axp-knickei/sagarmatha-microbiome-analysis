#!/usr/bin/env Rscript

# Build publication-formatted Random Forest figures and source data from a completed run.

required_packages <- c("dplyr", "tidyr", "ggplot2", "gridExtra", "jsonlite", "yaml", "randomForest", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && (is.na(x) || identical(x, "")))) y else x
}

usage <- function() {
  cat(
    "Usage: Rscript scripts/create_random_forest_publication_figure.R [options]\n\n",
    "Options:\n",
    "  --run-dir PATH              Completed Random Forest run. Default: newest results/random_forest/rf_*\n",
    "  --input PATH                Cleaned long-format RDS. Default: data/processed/cleaned_data.rds\n",
    "  --config PATH               Random Forest YAML. Default: config/random_forest_config.yml\n",
    "  --output-dir PATH           Default: RUN_DIR/publication\n",
    "  --bootstrap-replicates N    Subject-bootstrap replicates. Default: 1000\n",
    "  --top-features N            Importance features in panel C. Default: 12\n",
    "  --abundance-features N      Features in panel D. Default: 4\n",
    "  --seed N                    Default: 20260715\n",
    "  --force                     Replace an existing publication directory.\n",
    "  --help\n",
    sep = ""
  )
}

parse_args <- function(args) {
  out <- list(
    run_dir = NULL,
    input = "data/processed/cleaned_data.rds",
    config = "config/random_forest_config.yml",
    output_dir = NULL,
    bootstrap_replicates = 1000L,
    top_features = 12L,
    abundance_features = 4L,
    seed = 20260715L,
    force = FALSE
  )
  key_map <- c(
    "--run-dir" = "run_dir", "--input" = "input", "--config" = "config",
    "--output-dir" = "output_dir", "--bootstrap-replicates" = "bootstrap_replicates",
    "--top-features" = "top_features", "--abundance-features" = "abundance_features",
    "--seed" = "seed"
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      usage()
      quit(status = 0)
    } else if (arg == "--force") {
      out$force <- TRUE
      i <- i + 1L
    } else if (arg %in% names(key_map)) {
      if (i == length(args)) stop("Missing value for ", arg, call. = FALSE)
      out[[key_map[[arg]]]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  out$bootstrap_replicates <- as.integer(out$bootstrap_replicates)
  out$top_features <- as.integer(out$top_features)
  out$abundance_features <- as.integer(out$abundance_features)
  out$seed <- as.integer(out$seed)
  if (is.na(out$bootstrap_replicates) || out$bootstrap_replicates < 10) {
    stop("--bootstrap-replicates must be at least 10.", call. = FALSE)
  }
  if (is.na(out$top_features) || out$top_features < 1 || is.na(out$abundance_features) || out$abundance_features < 1) {
    stop("Feature counts must be positive integers.", call. = FALSE)
  }
  out
}

newest_run <- function(base = "results/random_forest") {
  candidates <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  candidates <- candidates[grepl("/rf_[0-9]{8}_[0-9]{6}($|_)", candidates)]
  candidates <- candidates[file.exists(file.path(candidates, "summary.json"))]
  if (length(candidates) == 0) stop("No completed rf_* run found under ", base, call. = FALSE)
  sort(candidates)[[length(candidates)]]
}

safe_ratio <- function(numerator, denominator) ifelse(denominator == 0, 0, numerator / denominator)

classification_by_class <- function(observed, predicted, classes) {
  observed <- factor(as.character(observed), levels = classes)
  predicted <- factor(as.character(predicted), levels = classes)
  cm <- table(observed = observed, predicted = predicted)
  bind_rows(lapply(classes, function(class_name) {
    tp <- cm[class_name, class_name]
    predicted_n <- sum(cm[, class_name])
    observed_n <- sum(cm[class_name, ])
    precision <- safe_ratio(tp, predicted_n)
    recall <- safe_ratio(tp, observed_n)
    data.frame(
      class = class_name,
      precision = precision,
      recall = recall,
      f1 = safe_ratio(2 * precision * recall, precision + recall),
      support = observed_n
    )
  }))
}

balanced_accuracy <- function(observed, predicted, classes) {
  metrics <- classification_by_class(observed, predicted, classes)
  present <- metrics$support > 0
  if (!any(present)) return(NA_real_)
  mean(metrics$recall[present])
}

interval_table <- function(values, estimate, group_columns) {
  values |>
    group_by(across(all_of(group_columns))) |>
    summarise(
      lower = stats::quantile(value, 0.025, na.rm = TRUE, names = FALSE),
      upper = stats::quantile(value, 0.975, na.rm = TRUE, names = FALSE),
      bootstrap_median = stats::median(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    left_join(estimate, by = group_columns)
}

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.title = element_text(size = base_size, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", color = "black"),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.tag = element_text(size = base_size + 3, face = "bold", color = "black"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(4, 5, 4, 5)
    )
}

save_publication_figure <- function(plot, stem, directory, width_mm = 183, height_mm = 220, dpi = 600) {
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = grDevices::cairo_pdf, bg = "white")
  ggsave(file.path(directory, paste0(stem, ".svg")), plot, width = width, height = height,
         device = svglite::svglite, bg = "white")
  ggsave(file.path(directory, paste0(stem, ".png")), plot, width = width, height = height,
         device = "png", dpi = 450, bg = "white")
  ggsave(file.path(directory, paste0(stem, ".tiff")), plot, width = width, height = height,
         device = "tiff", dpi = dpi, compression = "lzw", bg = "white")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
run_dir <- args$run_dir %||% newest_run()
output_dir <- args$output_dir %||% file.path(run_dir, "publication")
figures_dir <- file.path(output_dir, "figures")
tables_dir <- file.path(output_dir, "source_data")
text_dir <- file.path(output_dir, "text")

required_run_files <- c(
  "summary.json", "test_predictions.csv", "sample_split.csv", "feature_importance.csv",
  "random_forest_model.rds"
)
missing_run_files <- required_run_files[!file.exists(file.path(run_dir, required_run_files))]
if (length(missing_run_files) > 0) {
  stop("Run is incomplete; missing: ", paste(missing_run_files, collapse = ", "), call. = FALSE)
}
if (!file.exists(args$input)) stop("Input data not found: ", args$input, call. = FALSE)
if (!file.exists(args$config)) stop("Configuration not found: ", args$config, call. = FALSE)
if (dir.exists(output_dir) && !args$force) {
  stop("Output directory exists; use --force to replace it: ", output_dir, call. = FALSE)
}
if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)
invisible(vapply(c(figures_dir, tables_dir, text_dir), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))

summary <- jsonlite::read_json(file.path(run_dir, "summary.json"), simplifyVector = TRUE)
config <- yaml::read_yaml(args$config)
if (!identical(summary$task, "classification")) {
  stop("This publication workflow currently requires a classification Random Forest run.", call. = FALSE)
}
if (isTRUE(summary$include_covariates)) {
  stop("Publication reconstruction for models with encoded covariates is not implemented.", call. = FALSE)
}

predictions <- utils::read.csv(file.path(run_dir, "test_predictions.csv"), check.names = FALSE, stringsAsFactors = FALSE)
importance <- utils::read.csv(file.path(run_dir, "feature_importance.csv"), check.names = FALSE, stringsAsFactors = FALSE)
model <- readRDS(file.path(run_dir, "random_forest_model.rds"))
response_col <- summary$response
sample_col <- summary$sample_id_column
subject_col <- summary$subject_id_column
required_prediction_columns <- c(sample_col, subject_col, response_col, "predicted")
missing_prediction_columns <- setdiff(required_prediction_columns, names(predictions))
if (length(missing_prediction_columns) > 0) {
  stop("Prediction table is missing: ", paste(missing_prediction_columns, collapse = ", "), call. = FALSE)
}

classes <- sort(unique(as.character(c(predictions[[response_col]], predictions$predicted))))
numeric_classes <- suppressWarnings(as.numeric(classes))
if (all(!is.na(numeric_classes))) classes <- classes[order(numeric_classes)]
predictions[[response_col]] <- factor(as.character(predictions[[response_col]]), levels = classes)
predictions$predicted <- factor(as.character(predictions$predicted), levels = classes)

native_score_column <- if ("MeanDecreaseAccuracy" %in% names(importance)) {
  "MeanDecreaseAccuracy"
} else if ("MeanDecreaseGini" %in% names(importance)) {
  "MeanDecreaseGini"
} else {
  setdiff(names(importance), "feature")[[1]]
}
selected_features <- importance |>
  filter(is.finite(.data[[native_score_column]])) |>
  arrange(desc(.data[[native_score_column]]), feature) |>
  slice_head(n = args$top_features) |>
  pull(feature)
if (length(selected_features) == 0) stop("No finite feature-importance values were found.", call. = FALSE)

# Reconstruct the exact transformed predictor matrix used by the fitted model.
df <- readRDS(args$input)
feature_col <- config$data$feature_id %||% "MAG"
abundance_col <- config$data$abundance_value %||% "Y"
pseudocount <- as.numeric(config$features$pseudocount %||% 0)
transform_name <- config$features$transform %||% "log10_pseudocount"
required_input_columns <- c(sample_col, subject_col, response_col, feature_col, abundance_col)
missing_input_columns <- setdiff(required_input_columns, names(df))
if (length(missing_input_columns) > 0) {
  stop("Cleaned input is missing: ", paste(missing_input_columns, collapse = ", "), call. = FALSE)
}
metadata <- df |>
  distinct(.data[[sample_col]], .keep_all = TRUE) |>
  select(-all_of(c(feature_col, abundance_col)))
wide <- df |>
  transmute(
    !!sample_col := .data[[sample_col]],
    !!feature_col := as.character(.data[[feature_col]]),
    abundance = as.numeric(.data[[abundance_col]])
  ) |>
  pivot_wider(names_from = all_of(feature_col), values_from = abundance, values_fill = 0)
model_df <- metadata |> left_join(wide, by = sample_col)
model_features <- rownames(randomForest::importance(model))
missing_model_features <- setdiff(model_features, names(model_df))
if (length(missing_model_features) > 0) {
  stop("Could not reconstruct model features: ", paste(head(missing_model_features, 10), collapse = ", "), call. = FALSE)
}
feature_matrix <- model_df[model_features]
feature_matrix[] <- lapply(feature_matrix, as.numeric)
for (feature in names(feature_matrix)) {
  if (anyNA(feature_matrix[[feature]])) {
    feature_matrix[[feature]][is.na(feature_matrix[[feature]])] <- stats::median(feature_matrix[[feature]], na.rm = TRUE)
  }
}
if (identical(transform_name, "log10_pseudocount")) {
  feature_matrix[] <- lapply(feature_matrix, function(x) log10(x + pseudocount))
} else if (!identical(transform_name, "none")) {
  stop("Unsupported transform: ", transform_name, call. = FALSE)
}
test_match <- match(as.character(predictions[[sample_col]]), as.character(model_df[[sample_col]]))
if (anyNA(test_match)) stop("Some held-out sample IDs were not found in the cleaned input.", call. = FALSE)
test_x <- feature_matrix[test_match, , drop = FALSE]

# Panel A source data: counts and row percentages.
confusion_counts <- table(
  observed = factor(predictions[[response_col]], levels = classes),
  predicted = factor(predictions$predicted, levels = classes)
)
confusion_source <- as.data.frame(confusion_counts, stringsAsFactors = FALSE) |>
  group_by(observed) |>
  mutate(
    row_total = sum(Freq),
    row_percent = ifelse(row_total == 0, NA_real_, 100 * Freq / row_total),
    label = sprintf("%d\n%.0f%%", Freq, row_percent)
  ) |>
  ungroup()

# Subject bootstrap for precision, recall, and held-out permutation importance.
set.seed(args$seed)
subjects <- unique(as.character(predictions[[subject_col]]))
if (length(subjects) < 2) stop("At least two held-out subjects are required for bootstrap intervals.", call. = FALSE)
bootstrap_indices <- lapply(seq_len(args$bootstrap_replicates), function(i) {
  sampled <- sample(subjects, length(subjects), replace = TRUE)
  unlist(lapply(sampled, function(subject) which(as.character(predictions[[subject_col]]) == subject)), use.names = FALSE)
})

point_performance <- classification_by_class(predictions[[response_col]], predictions$predicted, classes) |>
  select(class, precision, recall) |>
  pivot_longer(c(precision, recall), names_to = "metric", values_to = "estimate")
bootstrap_performance <- bind_rows(lapply(seq_along(bootstrap_indices), function(i) {
  index <- bootstrap_indices[[i]]
  classification_by_class(predictions[[response_col]][index], predictions$predicted[index], classes) |>
    select(class, precision, recall) |>
    pivot_longer(c(precision, recall), names_to = "metric", values_to = "value") |>
    mutate(replicate = i)
}))
performance_intervals <- interval_table(bootstrap_performance, point_performance, c("class", "metric"))

point_balanced_accuracy <- balanced_accuracy(predictions[[response_col]], predictions$predicted, classes)
permutation_bootstrap <- bind_rows(lapply(seq_along(bootstrap_indices), function(i) {
  index <- bootstrap_indices[[i]]
  x_boot <- test_x[index, , drop = FALSE]
  observed <- predictions[[response_col]][index]
  baseline_prediction <- predict(model, x_boot, type = "response")
  baseline <- balanced_accuracy(observed, baseline_prediction, classes)
  bind_rows(lapply(seq_along(selected_features), function(j) {
    feature <- selected_features[[j]]
    x_permuted <- x_boot
    x_permuted[[feature]] <- sample(x_permuted[[feature]])
    permuted_prediction <- predict(model, x_permuted, type = "response")
    data.frame(
      replicate = i,
      feature = feature,
      importance = baseline - balanced_accuracy(observed, permuted_prediction, classes)
    )
  }))
}))
permutation_intervals <- permutation_bootstrap |>
  group_by(feature) |>
  summarise(
    importance = mean(importance, na.rm = TRUE),
    lower = stats::quantile(importance, 0.025, na.rm = TRUE, names = FALSE),
    upper = stats::quantile(importance, 0.975, na.rm = TRUE, names = FALSE),
    bootstrap_median = stats::median(importance, na.rm = TRUE),
    positive_fraction = mean(importance > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(importance |> select(feature, native_score = all_of(native_score_column)), by = "feature") |>
  arrange(desc(native_score), feature)

abundance_features <- head(selected_features, args$abundance_features)
abundance_source <- df |>
  filter(.data[[feature_col]] %in% abundance_features) |>
  transmute(
    sample_id = .data[[sample_col]],
    subject_id = .data[[subject_col]],
    altitude = as.character(.data[[response_col]]),
    feature = as.character(.data[[feature_col]]),
    raw_abundance = as.numeric(.data[[abundance_col]]),
    transformed_abundance = if (identical(transform_name, "log10_pseudocount")) {
      log10(raw_abundance + pseudocount)
    } else {
      raw_abundance
    },
    split = ifelse(.data[[sample_col]] %in% predictions[[sample_col]], "held-out test", "training")
  ) |>
  mutate(
    altitude = factor(altitude, levels = classes),
    feature = factor(feature, levels = abundance_features)
  )

palette_altitude <- setNames(c("#3B4CC0", "#2A9D8F", "#E9C46A", "#D1495B")[seq_along(classes)], classes)

panel_a <- ggplot(confusion_source, aes(x = predicted, y = observed, fill = row_percent)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 2.7, lineheight = 0.9) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 100), name = "Row (%)") +
  coord_equal() +
  labs(title = "Held-out classification", x = "Predicted altitude (m)", y = "Observed altitude (m)") +
  theme_nature() +
  theme(panel.grid = element_blank())

panel_b <- performance_intervals |>
  mutate(
    class = factor(class, levels = classes),
    metric = factor(metric, levels = c("recall", "precision"), labels = c("Recall", "Precision"))
  ) |>
  ggplot(aes(x = class, y = estimate, color = metric, group = metric)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, linewidth = 0.45,
                position = position_dodge(width = 0.36)) +
  geom_point(size = 2, position = position_dodge(width = 0.36)) +
  scale_color_manual(values = c(Recall = "#2166AC", Precision = "#B2182B")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.03))) +
  labs(title = "Performance by altitude", x = "Altitude (m)", y = "Metric (subject-bootstrap 95% interval)", color = NULL) +
  theme_nature() +
  theme(legend.position = c(0.78, 0.18))

feature_labels <- setNames(paste0("MAG ", seq_along(selected_features)), selected_features)
panel_c_data <- permutation_intervals |>
  mutate(feature_label = factor(feature_labels[feature], levels = rev(feature_labels[selected_features])))
panel_c <- ggplot(panel_c_data, aes(x = importance, y = feature_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "gray45") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.18,
                linewidth = 0.45, color = "#555555") +
  geom_point(aes(fill = positive_fraction), shape = 21, size = 2.2, stroke = 0.35) +
  scale_fill_gradient(
    low = "white", high = "#B2182B", limits = c(0, 1), name = "Positive fraction",
    guide = guide_colorbar(direction = "horizontal", title.position = "top",
                           barwidth = grid::unit(28, "mm"), barheight = grid::unit(2, "mm"))
  ) +
  labs(
    title = "Held-out permutation importance",
    x = "Decrease in balanced accuracy",
    y = NULL
  ) +
  theme_nature() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 6.5),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-3, 0, 0, 0)
  )

panel_d <- ggplot(abundance_source, aes(x = altitude, y = transformed_abundance, color = altitude)) +
  geom_boxplot(aes(group = altitude), width = 0.55, outlier.shape = NA, fill = "white",
               color = "black", linewidth = 0.35) +
  geom_point(aes(shape = split), position = position_jitter(width = 0.12, height = 0),
             size = 1.25, alpha = 0.75, stroke = 0.25) +
  scale_color_manual(values = palette_altitude, guide = "none") +
  scale_shape_manual(values = c("training" = 16, "held-out test" = 1), name = NULL) +
  facet_wrap(~ feature, scales = "free_y", ncol = 2, labeller = as_labeller(feature_labels)) +
  labs(
    title = "Abundance profiles of training-selected MAGs",
    x = "Altitude (m)",
    y = if (identical(transform_name, "log10_pseudocount")) "log10(relative abundance + pseudocount)" else "Abundance"
  ) +
  theme_nature(base_size = 7.5) +
  theme(legend.position = "bottom")

annotation <- paste0(
  "n = ", summary$n_samples, " samples; ", summary$n_subjects, " subjects; ",
  " held-out n = ", nrow(predictions), " samples from ", length(subjects), " subjects; ",
  " grouped CV; balanced accuracy = ", sprintf("%.2f", point_balanced_accuracy)
)
header <- gridExtra::arrangeGrob(
  grid::textGrob(
    "Subject-aware Random Forest classification of altitude from gut MAG profiles",
    x = 0, hjust = 0, gp = grid::gpar(fontfamily = "sans", fontsize = 11, fontface = "bold")
  ),
  grid::textGrob(
    annotation,
    x = 0, hjust = 0, gp = grid::gpar(fontfamily = "sans", fontsize = 8, col = "gray25")
  ),
  ncol = 1,
  heights = c(1.2, 1)
)
combined <- gridExtra::arrangeGrob(
  header,
  panel_a + labs(tag = "A"),
  panel_b + labs(tag = "B"),
  panel_c + labs(tag = "C"),
  panel_d + labs(tag = "D"),
  layout_matrix = rbind(c(1, 1), c(2, 3), c(4, 5)),
  heights = c(0.12, 0.78, 1.1)
)

figure_stem <- "Figure_RF1_subject_aware_altitude_classification"
save_publication_figure(combined, figure_stem, figures_dir)

utils::write.csv(confusion_source |> select(-label), file.path(tables_dir, "Figure_RF1a_confusion_matrix.csv"), row.names = FALSE)
utils::write.csv(performance_intervals, file.path(tables_dir, "Figure_RF1b_class_performance.csv"), row.names = FALSE)
utils::write.csv(permutation_intervals, file.path(tables_dir, "Figure_RF1c_permutation_importance.csv"), row.names = FALSE)
utils::write.csv(abundance_source, file.path(tables_dir, "Figure_RF1d_MAG_abundance.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(
    samples = summary$n_samples,
    subjects = summary$n_subjects,
    held_out_samples = nrow(predictions),
    held_out_subjects = length(subjects),
    accuracy = mean(predictions[[response_col]] == predictions$predicted),
    balanced_accuracy = point_balanced_accuracy,
    majority_class_accuracy = max(prop.table(table(predictions[[response_col]]))),
    trees = summary$best_parameters$ntree,
    mtry = summary$best_parameters$mtry,
    minimum_node_size = summary$best_parameters$nodesize
  ),
  file.path(tables_dir, "Figure_RF1_run_summary.csv"), row.names = FALSE
)
utils::write.csv(
  data.frame(figure_label = unname(feature_labels), feature = names(feature_labels)),
  file.path(tables_dir, "Figure_RF1_MAG_key.csv"), row.names = FALSE
)

caption <- paste0(
  "# Figure RF1 | Subject-aware Random Forest classification of altitude from gut MAG profiles\n\n",
  "**a,** Row-normalized confusion matrix for the held-out test set; cells show sample count and row percentage. ",
  "**b,** Per-class recall and precision with percentile 95% intervals from ", args$bootstrap_replicates,
  " subject-level bootstrap replicates of the held-out subjects. ",
  "**c,** Change in held-out balanced accuracy after permutation of features preselected by training-set mean decrease in accuracy. ",
  "Points show mean importance, horizontal bars show subject-bootstrap percentile 95% intervals, and fill shows the fraction of bootstrap replicates with positive importance. ",
  "**d,** Log10-pseudocount-transformed abundance distributions for the four highest-ranked training-selected MAGs; open points denote held-out samples. ",
  "The model used ", summary$n_samples, " samples from ", summary$n_subjects, " subjects and ", summary$n_features_after_filtering,
  " MAG predictors. Training/test separation and cross-validation were grouped by subject; the held-out set contained ", nrow(predictions),
  " samples from ", length(subjects), " subjects. Feature importance describes prediction, not direction, causality, or statistical association. ",
  "Intervals condition on the fitted forest and reflect resampling of only the held-out subjects; with ", length(subjects),
  " held-out subjects, uncertainty is coarse and the analysis remains exploratory.\n"
)
writeLines(caption, file.path(text_dir, "Figure_RF1_caption.md"))

methods <- paste0(
  "# Random Forest figure methods\n\n",
  "MAG relative abundances were transformed as log10(x + ", format(pseudocount, scientific = TRUE), "). ",
  "A Random Forest classifier containing ", summary$best_parameters$ntree, " trees was tuned by ",
  config$model$cv_folds, "-fold subject-grouped cross-validation. The selected mtry was ",
  summary$best_parameters$mtry, " and the selected minimum terminal node size was ", summary$best_parameters$nodesize,
  ". Subjects were separated between training and testing. Features displayed in the importance panel were selected using training-derived mean decrease in accuracy. ",
  "Held-out permutation importance was defined as the decrease in balanced accuracy after shuffling one feature. ",
  "Intervals were estimated with ", args$bootstrap_replicates, " subject-level bootstrap replicates of the held-out set.\n"
)
writeLines(methods, file.path(text_dir, "Figure_RF1_methods.md"))

manifest <- list(
  generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  source_run = normalizePath(run_dir),
  input = normalizePath(args$input),
  bootstrap_replicates = args$bootstrap_replicates,
  seed = args$seed,
  held_out_subjects = subjects,
  caution = paste0(
    "Publication-formatted output does not make the analysis confirmatory. Intervals condition on the fitted model and only ",
    length(subjects), " held-out subjects."
  ),
  figures = paste0("figures/", figure_stem, c(".pdf", ".svg", ".png", ".tiff")),
  source_data = basename(list.files(tables_dir, full.names = TRUE))
)
yaml::write_yaml(manifest, file.path(output_dir, "publication_output_manifest.yml"))

cat("Random Forest publication outputs completed: ", normalizePath(output_dir), "\n", sep = "")
cat("Main figure: ", file.path(figures_dir, paste0(figure_stem, ".pdf")), "\n", sep = "")
