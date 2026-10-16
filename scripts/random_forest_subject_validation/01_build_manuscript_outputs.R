#!/usr/bin/env Rscript

# Purpose: Build reader-first manuscript figures and tables from a completed
# subject-aware Random Forest validation run.
# Inputs: Only CSV/YAML outputs produced by the canonical subject-validation
# workflow under results/random_forest_subject_validation/<run-id>/.
# Outputs: One three-panel main figure, two main tables, exact panel data,
# validation records, and rendering/session metadata under <run-id>/manuscript/.
# Dependencies: R packages declared in environment.yml. No model is refitted,
# no raw or processed data are read, and canonical result files are unchanged.
# Usage: Rscript scripts/random_forest_subject_validation/01_build_manuscript_outputs.R \
#          --run-dir results/random_forest_subject_validation/<run-id>
# Relationship: Publication post-processing for the validation workflow only;
# the completed run remains the scientific source of truth.

required_packages <- c(
  "ggplot2", "gridExtra", "grid", "scales", "svglite", "ragg", "yaml"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    ". Update the glmm-hpc environment from environment.yml.",
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage: Rscript scripts/random_forest_subject_validation/01_build_manuscript_outputs.R [options]\n\n",
    "Options:\n",
    "  --run-dir PATH     Completed subject-validation run. Default: newest complete non-smoke run.\n",
    "  --output-dir PATH  Output directory. Default: RUN_DIR/manuscript.\n",
    "  --allow-smoke      Permit a smoke-mode run when --run-dir is explicit.\n",
    "  --force            Replace only the selected manuscript output directory.\n",
    "  --help             Show this message.\n",
    sep = ""
  )
}

arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop(flag, " requires a value.", call. = FALSE)
  args[[idx + 1L]]
}

has_flag <- function(flag) flag %in% args

known_flags <- c("--run-dir", "--output-dir", "--allow-smoke", "--force", "--help")
unknown_flags <- args[grepl("^--", args) & !(args %in% known_flags)]
if (length(unknown_flags) > 0L) {
  stop("Unknown option(s): ", paste(unique(unknown_flags), collapse = ", "), call. = FALSE)
}
if (has_flag("--help")) {
  usage()
  quit(status = 0L)
}

repo_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(repo_root, "analysis", "random_forest_subject_validation", "run_random_forest_subject_validation.R"))) {
  stop("Run this script from the repository root.", call. = FALSE)
}

validation_root <- file.path(repo_root, "results", "random_forest_subject_validation")

read_csv <- function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

write_csv <- function(data, path) {
  utils::write.csv(data, path, row.names = FALSE, na = "")
}

write_tsv <- function(data, path) {
  utils::write.table(
    data,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

require_columns <- function(data, columns, source_name) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stop(
      source_name, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

select_default_run <- function(root) {
  candidates <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  if (length(candidates) == 0L) {
    stop("No subject-validation run directories were found under ", root, call. = FALSE)
  }

  records <- lapply(candidates, function(candidate) {
    manifest_path <- file.path(candidate, "reports", "run_manifest.yml")
    if (!file.exists(manifest_path)) return(NULL)
    manifest <- tryCatch(yaml::read_yaml(manifest_path), error = function(e) NULL)
    if (is.null(manifest) || !identical(as.character(manifest$status), "complete")) return(NULL)
    mode <- as.character(manifest$mode)
    if (!(mode %in% c("production", "development"))) return(NULL)
    data.frame(
      path = candidate,
      mode = mode,
      priority = ifelse(mode == "production", 2L, 1L),
      completed = as.character(manifest$completed),
      stringsAsFactors = FALSE
    )
  })
  records <- records[!vapply(records, is.null, logical(1))]
  if (length(records) == 0L) {
    stop(
      "No complete production/development run was found. Supply --run-dir explicitly; smoke runs also require --allow-smoke.",
      call. = FALSE
    )
  }
  index <- do.call(rbind, records)
  index <- index[order(index$priority, index$completed, index$path), , drop = FALSE]
  index$path[[nrow(index)]]
}

run_dir_input <- arg_value("--run-dir")
run_dir <- if (is.null(run_dir_input)) {
  select_default_run(validation_root)
} else {
  normalizePath(run_dir_input, mustWork = TRUE)
}

normalized_validation_root <- normalizePath(validation_root, mustWork = TRUE)
normalized_run_dir <- normalizePath(run_dir, mustWork = TRUE)
if (!startsWith(normalized_run_dir, paste0(normalized_validation_root, .Platform$file.sep))) {
  stop(
    "The selected run must be below results/random_forest_subject_validation/.",
    call. = FALSE
  )
}

manifest_path <- file.path(run_dir, "reports", "run_manifest.yml")
if (!file.exists(manifest_path)) {
  stop("Run manifest not found: ", manifest_path, call. = FALSE)
}
run_manifest <- yaml::read_yaml(manifest_path)
if (!identical(as.character(run_manifest$status), "complete")) {
  stop("The selected run is not complete according to its manifest.", call. = FALSE)
}
run_mode <- as.character(run_manifest$mode)
if (identical(run_mode, "smoke") && !has_flag("--allow-smoke")) {
  stop("Smoke-mode results are not suitable for manuscript output. Use --allow-smoke only for layout testing.", call. = FALSE)
}
if (!(run_mode %in% c("production", "development", "smoke"))) {
  stop("Unsupported run mode: ", run_mode, call. = FALSE)
}

output_dir_input <- arg_value("--output-dir")
output_dir <- if (is.null(output_dir_input)) {
  file.path(run_dir, "manuscript")
} else {
  normalizePath(output_dir_input, mustWork = FALSE)
}

normalized_output_parent <- normalizePath(dirname(output_dir), mustWork = TRUE)
if (!identical(normalized_output_parent, normalized_run_dir) || basename(output_dir) != "manuscript") {
  stop(
    "For result safety, --output-dir must be RUN_DIR/manuscript.",
    call. = FALSE
  )
}
if (dir.exists(output_dir) && !has_flag("--force")) {
  stop(
    "Manuscript outputs already exist. Review them first, then use --force to replace only: ",
    output_dir,
    call. = FALSE
  )
}
if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)

paths <- list(
  dataset_summary = file.path(run_dir, "audit", "dataset_summary.csv"),
  confusion = file.path(run_dir, "predictions", "loso_confusion_matrix.csv"),
  overall = file.path(run_dir, "predictions", "loso_overall_metrics.csv"),
  per_class = file.path(run_dir, "predictions", "loso_per_class_metrics.csv"),
  subject_macro = file.path(run_dir, "predictions", "loso_subject_macro_metrics.csv"),
  baseline = file.path(run_dir, "comparisons", "baseline_model_comparison.csv"),
  class_balance = file.path(run_dir, "comparisons", "class_balance_strategy_comparison.csv"),
  tree_stability = file.path(run_dir, "comparisons", "tree_count_stability.csv"),
  output_validation = file.path(run_dir, "reports", "output_validation.csv"),
  run_manifest = manifest_path,
  resolved_config = file.path(run_dir, "reports", "resolved_config.yml")
)
missing_inputs <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_inputs) > 0L) {
  stop(
    "Selected run is missing required result files: ",
    paste(paste0(missing_inputs, "=", unlist(paths[missing_inputs])), collapse = "; "),
    call. = FALSE
  )
}

source_validation <- read_csv(paths$output_validation)
require_columns(source_validation, c("check", "passed", "detail"), "Source output validation")
passed_values <- tolower(as.character(source_validation$passed))
if (any(!(passed_values %in% c("true", "t", "1")))) {
  failed_checks <- source_validation$check[!(passed_values %in% c("true", "t", "1"))]
  stop("Source validation contains failed checks: ", paste(failed_checks, collapse = ", "), call. = FALSE)
}

dataset_summary <- read_csv(paths$dataset_summary)
confusion <- read_csv(paths$confusion)
overall <- read_csv(paths$overall)
per_class <- read_csv(paths$per_class)
subject_macro <- read_csv(paths$subject_macro)
baseline <- read_csv(paths$baseline)
class_balance <- read_csv(paths$class_balance)
tree_stability <- read_csv(paths$tree_stability)

require_columns(dataset_summary, c("metric", "value"), "Dataset summary")
require_columns(confusion, c("observed", "predicted", "Freq"), "LOSO confusion matrix")
require_columns(
  overall,
  c(
    "accuracy", "balanced_accuracy", "macro_precision", "macro_recall",
    "macro_f1", "ordinal_mae", "ordinal_rmse", "within_one_accuracy",
    "weighted_kappa", "log_loss"
  ),
  "LOSO overall metrics"
)
require_columns(per_class, c("class", "precision", "recall", "f1", "support"), "LOSO per-class metrics")
require_columns(subject_macro, c("metric", "estimate", "rule"), "LOSO subject-macro metrics")
require_columns(
  baseline,
  c(
    "model", "repeats", "accuracy_mean", "accuracy_sd",
    "balanced_accuracy_mean", "balanced_accuracy_sd", "macro_f1_mean",
    "macro_f1_sd", "ordinal_mae_mean", "ordinal_mae_sd",
    "within_one_accuracy_mean", "within_one_accuracy_sd",
    "weighted_kappa_mean", "weighted_kappa_sd", "log_loss_mean", "log_loss_sd"
  ),
  "Baseline comparison"
)
require_columns(class_balance, c("strategy", "balanced_accuracy", "macro_f1"), "Class-balance comparison")
require_columns(tree_stability, c("trees", "balanced_accuracy", "macro_f1"), "Tree-count stability")

if (nrow(overall) != 1L) stop("LOSO overall metrics must contain exactly one row.", call. = FALSE)
if (any(!is.finite(confusion$Freq)) || any(confusion$Freq < 0) || any(confusion$Freq != round(confusion$Freq))) {
  stop("Confusion-matrix counts must be finite non-negative integers.", call. = FALSE)
}

class_levels <- c("0", "3650", "5200", "6500")
grid <- expand.grid(
  observed = class_levels,
  predicted = class_levels,
  stringsAsFactors = FALSE
)
confusion$observed <- as.character(confusion$observed)
confusion$predicted <- as.character(confusion$predicted)
confusion <- merge(grid, confusion, by = c("observed", "predicted"), all.x = TRUE, sort = FALSE)
confusion$Freq[is.na(confusion$Freq)] <- 0
if (nrow(confusion) != 16L || any(!(confusion$observed %in% class_levels)) || any(!(confusion$predicted %in% class_levels))) {
  stop("The confusion matrix must represent the four configured outcome classes.", call. = FALSE)
}

safe_ratio <- function(numerator, denominator) {
  ifelse(denominator == 0, 0, numerator / denominator)
}

recomputed_class <- do.call(rbind, lapply(class_levels, function(class_name) {
  observed_rows <- confusion$observed == class_name
  predicted_rows <- confusion$predicted == class_name
  support <- sum(confusion$Freq[observed_rows])
  predicted_count <- sum(confusion$Freq[predicted_rows])
  correct_count <- sum(confusion$Freq[observed_rows & predicted_rows])
  precision <- safe_ratio(correct_count, predicted_count)
  recall <- safe_ratio(correct_count, support)
  f1 <- safe_ratio(2 * precision * recall, precision + recall)
  data.frame(
    class = class_name,
    support = support,
    predicted_count = predicted_count,
    correct_count = correct_count,
    precision = precision,
    recall = recall,
    f1 = f1,
    stringsAsFactors = FALSE
  )
}))

n_samples <- sum(confusion$Freq)
recomputed_overall <- list(
  accuracy = sum(recomputed_class$correct_count) / n_samples,
  balanced_accuracy = mean(recomputed_class$recall),
  macro_precision = mean(recomputed_class$precision),
  macro_recall = mean(recomputed_class$recall),
  macro_f1 = mean(recomputed_class$f1)
)

assert_close <- function(actual, expected, label, tolerance = 1e-10) {
  if (length(actual) != length(expected) || any(!is.finite(actual)) || any(!is.finite(expected)) || any(abs(actual - expected) > tolerance)) {
    stop(
      label, " failed independent recomputation. Actual=", paste(signif(actual, 12), collapse = ","),
      "; expected=", paste(signif(expected, 12), collapse = ","),
      call. = FALSE
    )
  }
}

per_class <- per_class[match(class_levels, as.character(per_class$class)), , drop = FALSE]
if (any(is.na(per_class$class))) stop("Per-class metrics do not contain all four configured classes.", call. = FALSE)
assert_close(recomputed_class$support, per_class$support, "Per-class support")
assert_close(recomputed_class$precision, per_class$precision, "Per-class precision")
assert_close(recomputed_class$recall, per_class$recall, "Per-class recall")
assert_close(recomputed_class$f1, per_class$f1, "Per-class F1")
for (metric_name in names(recomputed_overall)) {
  assert_close(recomputed_overall[[metric_name]], overall[[metric_name]], paste("Overall", metric_name))
}

dataset_lookup <- stats::setNames(as.character(dataset_summary$value), dataset_summary$metric)
expected_samples <- as.integer(dataset_lookup[["samples"]])
expected_subjects <- as.integer(dataset_lookup[["subjects"]])
expected_features <- as.integer(dataset_lookup[["MAG_predictors"]])
if (!isTRUE(n_samples == expected_samples)) {
  stop("Confusion-matrix total does not match dataset_summary samples.", call. = FALSE)
}

majority <- baseline[baseline$model == "majority", , drop = FALSE]
stratified <- baseline[baseline$model == "stratified_random", , drop = FALSE]
if (nrow(majority) != 1L || nrow(stratified) != 1L) {
  stop("Baseline comparison must contain one majority and one stratified_random row.", call. = FALSE)
}

subject_macro_accuracy <- subject_macro$estimate[subject_macro$metric == "subject_macro_accuracy"]
if (length(subject_macro_accuracy) != 1L) stop("Subject-macro accuracy is missing.", call. = FALSE)

figures_dir <- file.path(output_dir, "figures")
tables_dir <- file.path(output_dir, "tables")
panel_data_dir <- file.path(output_dir, "panel_data")
metadata_dir <- file.path(output_dir, "metadata")
invisible(lapply(
  c(output_dir, figures_dir, tables_dir, panel_data_dir, metadata_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

# Exact palette values from Ordnance Survey GeoDataViz Toolkit colour palettes
# v0.9.1 and its accompanying GDV-colour-palettes.json file. Labels, ordering,
# shapes, and printed values ensure that interpretation does not depend on hue.
gdv <- list(
  ink = "#252525",
  text = "#4D4D4D",
  neutral = "#767676",
  light_grey = "#D9D9D9",
  faint = "#F4F4F4",
  white = "#FFFFFF",
  purple_dark = "#443F90",
  blue = "#009ADE",
  purple = "#AF58BA",
  pink = "#FF1F5B",
  gold = "#FFC61E",
  orange = "#FD8D3C",
  sequential_blue = c(
    "#E4F1F7", "#C5E1EF", "#9EC9E2", "#6CB0D6",
    "#3C93C2", "#226E9C", "#0D4A70"
  )
)

class_labels <- c(
  "0" = "0 m",
  "3650" = "3,650 m",
  "5200" = "5,200 m",
  "6500" = "6,500 m"
)
class_palette <- c(
  "0" = gdv$pink,
  "3650" = gdv$blue,
  "5200" = gdv$purple,
  "6500" = gdv$gold
)
model_labels <- c(
  "Random Forest" = "Random Forest",
  "Majority" = "Majority baseline",
  "Stratified random" = "Stratified-random baseline"
)
model_colors <- c(
  "Random Forest" = gdv$purple_dark,
  "Majority" = gdv$neutral,
  "Stratified random" = gdv$orange
)
model_shapes <- c("Random Forest" = 21, "Majority" = 22, "Stratified random" = 24)

nature_theme <- function(base_size = 7.4) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = base_size + 1.0, color = gdv$ink,
        margin = ggplot2::margin(b = 2.5, l = 18)
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 0.25, color = gdv$text, lineheight = 1.05,
        margin = ggplot2::margin(b = 4.5, l = 18)
      ),
      plot.tag = ggplot2::element_text(
        face = "bold", size = 11, color = gdv$ink,
        margin = ggplot2::margin(r = 4)
      ),
      plot.tag.position = c(0.005, 0.995),
      plot.title.position = "plot",
      axis.title = ggplot2::element_text(size = base_size, color = gdv$ink),
      axis.text = ggplot2::element_text(size = base_size - 0.35, color = gdv$text),
      axis.line = ggplot2::element_line(linewidth = 0.35, color = gdv$ink),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, color = gdv$ink),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      legend.key.height = grid::unit(0.12, "in"),
      legend.spacing.x = grid::unit(0.08, "in"),
      plot.margin = ggplot2::margin(8, 8, 6, 8)
    )
}

# Panel A: the complete confusion matrix, normalized within observed class.
confusion$support <- ave(confusion$Freq, confusion$observed, FUN = sum)
confusion$row_fraction <- safe_ratio(confusion$Freq, confusion$support)
confusion$observed_label <- factor(
  class_labels[confusion$observed],
  levels = unname(class_labels)
)
confusion$predicted_label <- factor(
  class_labels[confusion$predicted],
  levels = unname(class_labels)
)
confusion$tile_label <- sprintf("%d\n(%.1f%%)", confusion$Freq, 100 * confusion$row_fraction)
confusion$label_color <- ifelse(confusion$row_fraction >= 0.48, gdv$white, gdv$ink)

panel_a <- ggplot2::ggplot(
  confusion,
  ggplot2::aes(x = predicted_label, y = observed_label, fill = row_fraction)
) +
  ggplot2::geom_tile(color = gdv$white, linewidth = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = tile_label, color = label_color),
    size = 2.65,
    lineheight = 0.92,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_fill_gradientn(
    colours = gdv$sequential_blue,
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_percent(accuracy = 1),
    name = "Within-row\nfraction"
  ) +
  ggplot2::scale_x_discrete(position = "top", drop = FALSE) +
  ggplot2::scale_y_discrete(limits = rev(unname(class_labels)), drop = FALSE) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Out-of-subject classification",
    subtitle = paste0(
      "Observed categories are rows; cells show sample count and row percentage (n = ",
      n_samples, ")"
    ),
    x = "Predicted expedition-linked category",
    y = "Observed expedition-linked category",
    tag = "A"
  ) +
  nature_theme() +
  ggplot2::theme(
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
    legend.position = "right",
    legend.key.height = grid::unit(0.28, "in"),
    legend.key.width = grid::unit(0.10, "in"),
    plot.margin = ggplot2::margin(8, 12, 7, 10)
  )

# Panel B: per-category sensitivity (recall) with observed support.
class_panel <- recomputed_class
class_panel$altitude_label <- class_labels[class_panel$class]
class_panel$axis_label <- sprintf("%s  (n = %d)", class_panel$altitude_label, class_panel$support)
class_panel$axis_label <- factor(
  class_panel$axis_label,
  levels = rev(sprintf("%s  (n = %d)", unname(class_labels), class_panel$support))
)
class_panel$label_x <- ifelse(class_panel$recall > 0.84, class_panel$recall - 0.035, class_panel$recall + 0.028)
class_panel$label_hjust <- ifelse(class_panel$recall > 0.84, 1, 0)
class_panel$label_color <- ifelse(class_panel$recall > 0.84, gdv$white, gdv$ink)

panel_b <- ggplot2::ggplot(
  class_panel,
  ggplot2::aes(x = recall, y = axis_label, fill = class)
) +
  ggplot2::geom_vline(
    xintercept = 0.25,
    linewidth = 0.45,
    linetype = "22",
    color = gdv$neutral
  ) +
  ggplot2::geom_col(width = 0.63, color = gdv$ink, linewidth = 0.35) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = label_x,
      label = scales::percent(recall, accuracy = 0.1),
      hjust = label_hjust,
      color = label_color
    ),
    size = 2.55,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::annotate(
    "text",
    x = 0.25,
    y = Inf,
    label = "25% reference",
    vjust = 1.25,
    hjust = -0.08,
    size = 2.1,
    color = gdv$neutral
  ) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_fill_manual(values = class_palette, guide = "none") +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "Recall by category",
    subtitle = "Percentage of observed samples correctly identified\nwithin each expedition-linked category",
    x = "Recall",
    y = NULL,
    tag = "B"
  ) +
  nature_theme() +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(linewidth = 0.28, color = gdv$light_grey),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(9, 8, 8, 10)
  )

# Panel C: three complementary performance measures against simple baselines.
benchmark_metrics <- data.frame(
  metric = rep(c("accuracy", "balanced_accuracy", "macro_f1"), each = 3L),
  model = rep(c("Random Forest", "Majority", "Stratified random"), times = 3L),
  estimate = c(
    overall$accuracy, majority$accuracy_mean, stratified$accuracy_mean,
    overall$balanced_accuracy, majority$balanced_accuracy_mean, stratified$balanced_accuracy_mean,
    overall$macro_f1, majority$macro_f1_mean, stratified$macro_f1_mean
  ),
  sd = c(
    NA, NA, stratified$accuracy_sd,
    NA, NA, stratified$balanced_accuracy_sd,
    NA, NA, stratified$macro_f1_sd
  ),
  stringsAsFactors = FALSE
)
metric_labels <- c(
  "accuracy" = "Overall accuracy",
  "balanced_accuracy" = "Balanced accuracy",
  "macro_f1" = "Macro F1"
)
benchmark_metrics$metric_label <- factor(
  metric_labels[benchmark_metrics$metric],
  levels = rev(unname(metric_labels))
)
benchmark_metrics$model <- factor(
  benchmark_metrics$model,
  levels = c("Random Forest", "Majority", "Stratified random")
)
benchmark_metrics$lower <- pmax(0, benchmark_metrics$estimate - benchmark_metrics$sd)
benchmark_metrics$upper <- pmin(1, benchmark_metrics$estimate + benchmark_metrics$sd)

dodge <- ggplot2::position_dodge(width = 0.58)
panel_c <- ggplot2::ggplot(
  benchmark_metrics,
  ggplot2::aes(x = metric_label, y = estimate, color = model, fill = model, shape = model)
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lower, ymax = upper),
    width = 0.16,
    linewidth = 0.5,
    position = dodge,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    size = 2.7,
    stroke = 0.65,
    position = dodge
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_color_manual(values = model_colors, labels = model_labels, drop = FALSE) +
  ggplot2::scale_fill_manual(values = model_colors, labels = model_labels, drop = FALSE) +
  ggplot2::scale_shape_manual(values = model_shapes, labels = model_labels, drop = FALSE) +
  ggplot2::scale_y_continuous(
    limits = c(0, 0.70),
    breaks = seq(0, 0.7, 0.1),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "Performance against baselines",
    subtitle = paste0(
      "Pooled leave-one-subject-out estimates;\nrandom baseline: mean +/-1 SD (",
      stratified$repeats, " repeats)"
    ),
    x = NULL,
    y = "Metric value",
    tag = "C"
  ) +
  nature_theme() +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(
        size = 2.4,
        shape = unname(model_shapes),
        fill = unname(model_colors),
        color = unname(model_colors),
        linetype = 0
      )
    ),
    fill = "none",
    shape = "none"
  ) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(linewidth = 0.28, color = gdv$light_grey),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "left",
    plot.margin = ggplot2::margin(9, 8, 8, 10)
  )

bottom_row <- gridExtra::arrangeGrob(panel_b, panel_c, ncol = 2, widths = c(0.95, 1.05))
main_figure <- gridExtra::arrangeGrob(
  panel_a,
  bottom_row,
  ncol = 1,
  heights = c(1.08, 1.00)
)

save_grob_bundle <- function(grob, stem, width, height) {
  ggplot2::ggsave(
    paste0(stem, ".pdf"), grob,
    width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".svg"), grob,
    width = width, height = height, units = "in",
    device = svglite::svglite, bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".png"), grob,
    width = width, height = height, units = "in", dpi = 600,
    device = ragg::agg_png, bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".tiff"), grob,
    width = width, height = height, units = "in", dpi = 600,
    device = ragg::agg_tiff, compression = "lzw", bg = "white"
  )
}

figure_stem <- file.path(figures_dir, "main_figure_random_forest_subject_validation")
save_grob_bundle(main_figure, figure_stem, width = 7.2, height = 7.1)

panel_a_data <- confusion[
  , c("observed", "predicted", "Freq", "support", "row_fraction")
]
names(panel_a_data) <- c("observed_class", "predicted_class", "sample_count", "observed_class_support", "within_observed_class_fraction")
panel_b_data <- class_panel[
  , c("class", "altitude_label", "support", "predicted_count", "correct_count", "precision", "recall", "f1")
]
panel_c_data <- benchmark_metrics[
  , c("metric", "model", "estimate", "sd", "lower", "upper")
]
write_csv(panel_a_data, file.path(panel_data_dir, "figure_panel_a_confusion_matrix.csv"))
write_csv(panel_b_data, file.path(panel_data_dir, "figure_panel_b_class_recall.csv"))
write_csv(panel_c_data, file.path(panel_data_dir, "figure_panel_c_baseline_comparison.csv"))

metric_table <- data.frame(
  metric = c(
    "accuracy", "balanced_accuracy", "macro_f1", "ordinal_mae",
    "within_one_accuracy", "weighted_kappa", "log_loss"
  ),
  metric_label = c(
    "Overall accuracy", "Balanced accuracy", "Macro F1", "Ordinal MAE",
    "Within-one-category accuracy", "Quadratic weighted kappa", "Multiclass log loss"
  ),
  random_forest = c(
    overall$accuracy, overall$balanced_accuracy, overall$macro_f1,
    overall$ordinal_mae, overall$within_one_accuracy,
    overall$weighted_kappa, overall$log_loss
  ),
  majority_baseline = c(
    majority$accuracy_mean, majority$balanced_accuracy_mean, majority$macro_f1_mean,
    majority$ordinal_mae_mean, majority$within_one_accuracy_mean,
    majority$weighted_kappa_mean, majority$log_loss_mean
  ),
  stratified_random_mean = c(
    stratified$accuracy_mean, stratified$balanced_accuracy_mean, stratified$macro_f1_mean,
    stratified$ordinal_mae_mean, stratified$within_one_accuracy_mean,
    stratified$weighted_kappa_mean, stratified$log_loss_mean
  ),
  stratified_random_sd = c(
    stratified$accuracy_sd, stratified$balanced_accuracy_sd, stratified$macro_f1_sd,
    stratified$ordinal_mae_sd, stratified$within_one_accuracy_sd,
    stratified$weighted_kappa_sd, stratified$log_loss_sd
  ),
  reader_interpretation = c(
    "Correct samples; sensitive to the dominant 5,200 m category.",
    "Mean recall across all four categories; 0.25 is the four-class equal-recall reference.",
    "Class-balanced summary of precision and recall; low values indicate weak discrimination.",
    "Mean absolute error in ordered category steps; lower values are better.",
    "Samples predicted no more than one ordered category from the observed category.",
    "Agreement beyond chance with larger penalties for more distant categories; 0 indicates chance-level agreement.",
    "Quality of predicted probabilities; lower values are better."
  ),
  stringsAsFactors = FALSE
)

format_metric_value <- function(metric, value) {
  if (metric %in% c("accuracy", "balanced_accuracy", "macro_f1", "within_one_accuracy")) {
    return(scales::percent(value, accuracy = 0.1))
  }
  sprintf("%.3f", value)
}

format_random <- function(metric, mean_value, sd_value) {
  if (metric %in% c("accuracy", "balanced_accuracy", "macro_f1", "within_one_accuracy")) {
    return(sprintf(
      "%s ± %s",
      scales::percent(mean_value, accuracy = 0.1),
      scales::percent(sd_value, accuracy = 0.1)
    ))
  }
  sprintf("%.3f ± %.3f", mean_value, sd_value)
}

metric_display <- data.frame(
  Metric = vapply(
    metric_table$metric_label,
    function(x) paste(strwrap(x, width = 22), collapse = "\n"),
    character(1)
  ),
  RF = mapply(format_metric_value, metric_table$metric, metric_table$random_forest),
  Majority = mapply(format_metric_value, metric_table$metric, metric_table$majority_baseline),
  `Stratified random\nmean ± s.d.` = mapply(
    format_random,
    metric_table$metric,
    metric_table$stratified_random_mean,
    metric_table$stratified_random_sd
  ),
  Interpretation = vapply(
    metric_table$reader_interpretation,
    function(x) paste(strwrap(x, width = 42), collapse = "\n"),
    character(1)
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

class_table <- recomputed_class
class_table$altitude_m <- as.integer(class_table$class)
class_table <- class_table[
  , c("class", "altitude_m", "support", "predicted_count", "correct_count", "precision", "recall", "f1")
]
class_display <- data.frame(
  `Category` = unname(class_labels[class_table$class]),
  `Observed, n` = class_table$support,
  `Predicted, n` = class_table$predicted_count,
  `Correct, n` = class_table$correct_count,
  Precision = scales::percent(class_table$precision, accuracy = 0.1),
  Recall = scales::percent(class_table$recall, accuracy = 0.1),
  `F1 score` = sprintf("%.3f", class_table$f1),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

markdown_escape <- function(x) {
  x <- gsub("\n", "<br>", as.character(x), fixed = TRUE)
  gsub("|", "\\|", x, fixed = TRUE)
}

write_markdown_table <- function(display_table, path, title, subtitle, footnote) {
  escaped <- as.data.frame(lapply(display_table, markdown_escape), check.names = FALSE)
  header <- paste0("| ", paste(names(escaped), collapse = " | "), " |")
  divider <- paste0("| ", paste(rep("---", ncol(escaped)), collapse = " | "), " |")
  rows <- apply(escaped, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  writeLines(
    c(
      paste0("# ", title), "", subtitle, "", header, divider, rows, "",
      paste0("Note: ", footnote)
    ),
    path
  )
}

wrap_text <- function(text, width) paste(strwrap(text, width = width), collapse = "\n")

build_table_grob <- function(display_table, title, subtitle, footnote, column_widths, base_size) {
  table_theme <- gridExtra::ttheme_minimal(
    base_size = base_size,
    base_family = "sans",
    core = list(
      fg_params = list(col = gdv$ink, hjust = 0, x = 0.03, lineheight = 0.95),
      bg_params = list(
        fill = rep(c(gdv$white, gdv$faint), length.out = nrow(display_table)),
        col = NA
      )
    ),
    colhead = list(
      fg_params = list(col = gdv$white, fontface = "bold", hjust = 0, x = 0.03),
      bg_params = list(fill = gdv$purple_dark, col = gdv$purple_dark)
    )
  )
  table_grob <- gridExtra::tableGrob(display_table, rows = NULL, theme = table_theme)
  table_grob$widths <- grid::unit(column_widths / sum(column_widths), "npc")

  title_grob <- grid::textGrob(
    wrap_text(title, 82), x = 0, hjust = 0,
    gp = grid::gpar(fontfamily = "sans", fontsize = 12, fontface = "bold", col = gdv$ink)
  )
  subtitle_grob <- grid::textGrob(
    wrap_text(subtitle, 108), x = 0, hjust = 0,
    gp = grid::gpar(fontfamily = "sans", fontsize = 8, lineheight = 1.05, col = gdv$text)
  )
  note_grob <- grid::textGrob(
    wrap_text(footnote, 126), x = 0, hjust = 0,
    gp = grid::gpar(fontfamily = "sans", fontsize = 6.5, lineheight = 1.07, col = gdv$text)
  )

  gridExtra::arrangeGrob(
    title_grob,
    subtitle_grob,
    table_grob,
    note_grob,
    ncol = 1,
    heights = grid::unit.c(
      grid::unit(0.38, "in"),
      grid::unit(0.46, "in"),
      grid::unit(1, "null"),
      grid::unit(0.72, "in")
    ),
    padding = grid::unit(0.05, "in")
  )
}

write_table_bundle <- function(
    raw_table,
    display_table,
    stem,
    title,
    subtitle,
    footnote,
    width,
    height,
    column_widths,
    base_size = 7.1) {
  write_csv(raw_table, paste0(stem, ".csv"))
  write_tsv(raw_table, paste0(stem, ".tsv"))
  write_markdown_table(
    display_table,
    paste0(stem, ".md"),
    title,
    subtitle,
    footnote
  )
  table_grob <- build_table_grob(
    display_table,
    title,
    subtitle,
    footnote,
    column_widths,
    base_size
  )
  save_grob_bundle(table_grob, stem, width = width, height = height)
}

table_1_title <- "Main Table 1 | Out-of-subject Random Forest performance"
table_1_subtitle <- paste0(
  "Development-mode leave-one-subject-out validation across ", n_samples,
  " samples from ", expected_subjects, " participants; stratified-random values are mean ± s.d. across ",
  stratified$repeats, " deterministic repeats."
)
table_1_note <- paste0(
  "The majority baseline always predicts the most common training category. Balanced accuracy and macro F1 weight categories equally. ",
  "Altitude was perfectly determined by sampling timepoint/date, so these are expedition-linked prediction categories rather than isolated altitude effects. ",
  "Current source run: ", basename(run_dir), " (", run_mode, "; manifest status complete)."
)
write_table_bundle(
  metric_table,
  metric_display,
  file.path(tables_dir, "main_table_1_model_performance"),
  table_1_title,
  table_1_subtitle,
  table_1_note,
  width = 7.2,
  height = 5.4,
  column_widths = c(0.19, 0.09, 0.10, 0.20, 0.42),
  base_size = 6.7
)

table_2_title <- "Main Table 2 | Performance within each expedition-linked category"
table_2_subtitle <- paste0(
  "Exact class-level results from pooled leave-one-subject-out predictions (n = ",
  n_samples, ")."
)
table_2_note <- paste0(
  "Precision is the fraction of predictions for a category that were correct; recall is the fraction of observed samples in that category that were recovered. ",
  "F1 is their harmonic mean. The categories are labelled by altitude for continuity, but altitude, timepoint and sampling date cannot be separated in this dataset."
)
write_table_bundle(
  class_table,
  class_display,
  file.path(tables_dir, "main_table_2_class_performance"),
  table_2_title,
  table_2_subtitle,
  table_2_note,
  width = 7.2,
  height = 3.6,
  column_widths = c(0.16, 0.14, 0.14, 0.13, 0.14, 0.14, 0.15),
  base_size = 7.1
)

class_balance_range <- range(class_balance$balanced_accuracy, na.rm = TRUE)
tree_balance_range <- range(tree_stability$balanced_accuracy, na.rm = TRUE)
run_summary <- data.frame(
  field = c(
    "run_id", "mode", "manifest_status", "samples", "subjects", "MAG_predictors",
    "loso_accuracy", "loso_balanced_accuracy", "loso_macro_f1",
    "subject_macro_accuracy", "majority_accuracy", "majority_balanced_accuracy",
    "stratified_random_balanced_accuracy_mean", "class_balance_balanced_accuracy_min",
    "class_balance_balanced_accuracy_max", "tree_stability_balanced_accuracy_min",
    "tree_stability_balanced_accuracy_max"
  ),
  value = c(
    basename(run_dir), run_mode, as.character(run_manifest$status), n_samples,
    expected_subjects, expected_features, overall$accuracy,
    overall$balanced_accuracy, overall$macro_f1, subject_macro_accuracy,
    majority$accuracy_mean, majority$balanced_accuracy_mean,
    stratified$balanced_accuracy_mean, class_balance_range[[1]],
    class_balance_range[[2]], tree_balance_range[[1]], tree_balance_range[[2]]
  ),
  stringsAsFactors = FALSE
)
write_csv(run_summary, file.path(metadata_dir, "manuscript_run_summary.csv"))

source_md5 <- unname(tools::md5sum(unlist(paths)))
names(source_md5) <- names(paths)
render_manifest <- list(
  status = "complete",
  generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  source_run = basename(run_dir),
  source_mode = run_mode,
  source_manifest_status = as.character(run_manifest$status),
  source_git_commit = as.character(run_manifest$git_commit),
  source_git_dirty = as.character(run_manifest$git_dirty),
  rendering_git_commit = tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[[1]],
    error = function(e) "unknown"
  ),
  palette = list(
    name = "Ordnance Survey GeoDataViz Toolkit colour palettes",
    version = "0.9.1",
    repository = "https://github.com/OrdnanceSurvey/GeoDataViz-Toolkit",
    definitions = "Colours/GDV-colour-palettes.json"
  ),
  source_files_md5 = as.list(source_md5),
  figure = "figures/main_figure_random_forest_subject_validation",
  tables = c(
    "tables/main_table_1_model_performance",
    "tables/main_table_2_class_performance"
  ),
  interpretation_boundary = paste(
    "Exploratory prediction among altitude/time-confounded expedition categories;",
    "not evidence of isolated altitude biology, association, mechanism, or causality."
  )
)
yaml::write_yaml(render_manifest, file.path(metadata_dir, "manuscript_render_manifest.yml"))
writeLines(capture.output(sessionInfo()), file.path(metadata_dir, "session_info.txt"))

expected_outputs <- c(
  paste0(figure_stem, c(".pdf", ".svg", ".png", ".tiff")),
  unlist(lapply(
    c(
      file.path(tables_dir, "main_table_1_model_performance"),
      file.path(tables_dir, "main_table_2_class_performance")
    ),
    function(stem) paste0(stem, c(".csv", ".tsv", ".md", ".pdf", ".svg", ".png", ".tiff"))
  )),
  file.path(panel_data_dir, c(
    "figure_panel_a_confusion_matrix.csv",
    "figure_panel_b_class_recall.csv",
    "figure_panel_c_baseline_comparison.csv"
  )),
  file.path(metadata_dir, c(
    "manuscript_run_summary.csv",
    "manuscript_render_manifest.yml",
    "session_info.txt"
  ))
)
output_checks <- data.frame(
  check = c(
    "source_output_validation_passed",
    "confusion_total_matches_dataset",
    "overall_metrics_recomputed",
    "per_class_metrics_recomputed",
    "all_expected_manuscript_outputs_exist",
    "all_expected_manuscript_outputs_nonempty"
  ),
  passed = c(
    TRUE,
    isTRUE(n_samples == expected_samples),
    TRUE,
    TRUE,
    all(file.exists(expected_outputs)),
    all(file.exists(expected_outputs)) && all(file.info(expected_outputs)$size > 0, na.rm = TRUE)
  ),
  detail = c(
    paste(nrow(source_validation), "canonical source checks passed"),
    paste(n_samples, "confusion observations for", expected_samples, "dataset samples"),
    "accuracy, balanced accuracy, macro precision/recall/F1 independently reproduced",
    "support, precision, recall and F1 independently reproduced for four classes",
    paste(length(expected_outputs), "expected files"),
    "all expected files have positive byte size"
  ),
  stringsAsFactors = FALSE
)
write_csv(output_checks, file.path(metadata_dir, "manuscript_output_validation.csv"))
if (!all(output_checks$passed)) {
  stop(
    "Manuscript output validation failed: ",
    paste(output_checks$check[!output_checks$passed], collapse = ", "),
    call. = FALSE
  )
}

cat("Random Forest manuscript post-processing completed.\n")
cat("Source run: ", run_dir, "\n", sep = "")
cat("Mode: ", run_mode, "\n", sep = "")
cat("Output: ", output_dir, "\n", sep = "")
cat("Main figure: ", paste0(figure_stem, ".pdf"), "\n", sep = "")
cat("Validation: all checks passed.\n")
