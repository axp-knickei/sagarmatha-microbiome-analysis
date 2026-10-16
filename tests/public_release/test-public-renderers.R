#!/usr/bin/env Rscript
# Synthetic regression checks; no study inputs or output paths are used.
# Run from repository root: Rscript tests/public_release/test-public-renderers.R
repo <- normalizePath(getwd())
root <- tempfile("public-renderers-")
dir.create(root)
set.seed(47)
ids <- paste0("synthetic_MAG_", seq_len(12))
d <- expand.grid(MAG = ids, Subject = paste0("synthetic_subject_", seq_len(9)), Timepoint = paste0("T", seq_len(10)), stringsAsFactors = FALSE)
d$sampleID <- paste(d$Subject, d$Timepoint, sep = "_")
d$days_from_baseline_subject <- as.integer(sub("T", "", d$Timepoint)) - 1L
d$altitude <- 100 + 20 * d$days_from_baseline_subject
d$age <- 20 + match(d$Subject, unique(d$Subject)); d$gender <- "synthetic"
effects <- rep(c(-0.025, 0.025), 6)
subject_shift <- rnorm(9, sd = 0.3)
d$Y <- 10^(-2 + effects[match(d$MAG, ids)] * d$days_from_baseline_subject + subject_shift[match(d$Subject, unique(d$Subject))] + rnorm(nrow(d), sd = 0.08))
input <- file.path(root, "cleaned.rds"); saveRDS(d, input)
tax <- data.frame(Genome = ids, phylum = rep(paste0("p__Synthetic", 1:3), 4), label = "synthetic taxonomy", domain = "d__Synthetic", class = "c__Synthetic", ordeer = "o__Synthetic", family = rep(paste0("f__Synthetic", 1:4), 3), genus = paste0("g__Synthetic", 1:12), species = "s__Synthetic")
tax_path <- file.path(root, "taxonomy.csv"); write.csv(tax, tax_path, row.names = FALSE)
r <- data.frame(MAG = ids, status = "lmm_ok", primary_analysis = TRUE, fit_succeeded = TRUE, primary_significant = FALSE, time_estimate = effects, time_std_error = 0.005, time_p_value = seq(0.001, 0.024, length.out = 12), time_p_adj_bh_primary = 0.2, zero_fraction = 0, n_subjects = 9, n_timepoints = 10, singular = FALSE, response_variance = 0.1, n_rows = 90, nonzero_sample_count = 90, subjects_with_nonzero = 9, warnings = "", errors = "")
scripts <- c("create_table_1_hpc_per_mag_lmertest_results", "create_table_2_hpc_per_mag_lmertest_results", "create_table_2_top_longitudinal_mags", paste0("plot_figure_", c(4:8, 10), "_hpc_per_mag_lmertest_results"))
for (discoveries in c(FALSE, TRUE)) {
  r$primary_significant <- discoveries
  r$time_p_adj_bh_primary <- if (discoveries) 0.02 else 0.2
  results_path <- file.path(root, paste0("results-", discoveries, ".csv")); write.csv(r, results_path, row.names = FALSE)
  for (stem in scripts) {
    out <- file.path(root, paste0(stem, "-", discoveries))
    Sys.setenv(PUBLIC_CLEANED_DATA = input, PUBLIC_LMM_RESULTS = results_path, PUBLIC_TAXONOMY = tax_path, PUBLIC_RENDER_OUTPUT_DIR = out)
    e <- new.env(parent = globalenv())
    # Reduce raster resolution only to keep smoke checks small; rendering is real.
    e$png <- function(...) { a <- list(...); a$res <- 60; do.call(grDevices::png, a) }
    sys.source(file.path(repo, "scripts", paste0(stem, "_public.R")), envir = e)
    artifacts <- list.files(out, full.names = TRUE)
    stopifnot(length(artifacts) >= 3L, all(file.info(artifacts)$size > 0))
    if (grepl("table_1", stem)) {
      stopifnot(e$n_subjects == 9, e$fdr_mags == if (discoveries) 12 else 0)
      stopifnot(!any(grepl("because no MAG", e$table_1$notes)))
    }
    if (grepl("table_2", stem)) {
      ranked <- if (exists("ranked", e, inherits = FALSE)) e$ranked else e$ranked_table
      stopifnot(nrow(ranked) == 12, setequal(ranked$MAG, ids))
    }
    if (grepl("figure_4_", stem)) stopifnot(e$fdr_mags == if (discoveries) 12 else 0, !grepl("no MAG", e$panel_c$labels$subtitle))
    if (grepl("figure_5_", stem)) stopifnot(nrow(e$selected) == 8, !grepl("no MAG", e$selection_text))
    if (grepl("figure_6_", stem)) stopifnot(e$fdr_mags == if (discoveries) 12 else 0, !grepl("no BH", e$panel_a$labels$subtitle))
    if (grepl("figure_7_", stem)) stopifnot(grepl("Subject counts range 9-9", e$summary_text, fixed = TRUE))
    if (grepl("figure_8_", stem)) stopifnot(nrow(e$selected) == 4, !grepl("no MAG", e$summary_text))
    if (grepl("figure_10_", stem)) stopifnot(nrow(e$selected) == 12, if (discoveries) e$selection_label == "BH-FDR-significant MAGs" else grepl("Exploratory", e$selection_label))
    cat("PASS", stem, "discoveries =", discoveries, "\n")
  }
}
Sys.unsetenv(c("PUBLIC_CLEANED_DATA", "PUBLIC_LMM_RESULTS", "PUBLIC_TAXONOMY", "PUBLIC_RENDER_OUTPUT_DIR"))
cat("Synthetic renderer artifacts:", root, "\n")

# Composite diversity renderer: all values below are artificial and independent
# of study inputs. Its canonical path is verified before use in assembled trees.
diversity_relative <- "scripts/manuscript_diversity_per_mag/01_build_main_figures_tables.R"
diversity_source <- sub("[.]R$", "_public.R", file.path(repo, diversity_relative))
if (!file.exists(diversity_source)) {
  diversity_source <- file.path(repo, diversity_relative)
  replacement_map <- read.csv(file.path(repo, "docs/releases/v1.0.0-public-replacements.csv"))
  mapping <- replacement_map[replacement_map$public_destination_path == diversity_relative, ]
  stopifnot(nrow(mapping) == 1L, identical(digest::digest(file = diversity_source, algo = "sha256"), mapping$replacement_sha256))
}
fixture <- file.path(root, "diversity-inputs")
write_fixture <- function(table, relative, tsv = FALSE) {
  path <- file.path(fixture, relative)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (tsv) write.table(table, path, sep = "\t", row.names = FALSE, quote = FALSE) else write.csv(table, path, row.names = FALSE)
}
alpha_fixture <- unique(d[c("sampleID", "Subject", "Timepoint", "days_from_baseline_subject", "altitude")])
names(alpha_fixture)[names(alpha_fixture) == "Subject"] <- "subjectID"
names(alpha_fixture)[names(alpha_fixture) == "days_from_baseline_subject"] <- "days_from_baseline"
alpha_fixture$Shannon <- 1 + alpha_fixture$days_from_baseline * 0.01 + seq_len(nrow(alpha_fixture)) * 0.0001
alpha_fixture$Simpson <- 0.4 + alpha_fixture$days_from_baseline * 0.001 + seq_len(nrow(alpha_fixture)) * 0.00001
alpha_fixture$Observed_richness <- 4 + alpha_fixture$days_from_baseline * 0.1 + seq_len(nrow(alpha_fixture)) * 0.01
alpha_dir <- "results/alpha_diversity_lmm/tables"
write_fixture(alpha_fixture, file.path(alpha_dir, "alpha_lmm_input.csv"))
alpha_effect_fixture <- data.frame(metric = c("Shannon", "Simpson", "Observed_richness"), metric_label = c("Shannon", "Simpson", "Observed richness"), model_type = "main", term = "days_from_baseline", estimate = c(0.01, 0.001, 0.1), standard_error = c(0.002, 0.0002, 0.02), p_value = c(0.01, 0.02, 0.03), confidence_interval_low = c(0.006, 0.0006, 0.06), confidence_interval_high = c(0.014, 0.0014, 0.14), n_samples = nrow(alpha_fixture), n_subjects = length(unique(alpha_fixture$subjectID)))
write_fixture(alpha_effect_fixture, file.path(alpha_dir, "alpha_lmm_fixed_effects.csv"))
alpha_prediction_fixture <- expand.grid(metric = alpha_effect_fixture$metric, time_value = 0:9)
alpha_prediction_fixture$metric_label <- alpha_prediction_fixture$metric
alpha_prediction_fixture$model_type <- "main"
alpha_prediction_fixture$estimate <- 1 + alpha_prediction_fixture$time_value * 0.01
alpha_prediction_fixture$confidence_interval_low <- alpha_prediction_fixture$estimate - 0.1
alpha_prediction_fixture$confidence_interval_high <- alpha_prediction_fixture$estimate + 0.1
write_fixture(alpha_prediction_fixture, file.path(alpha_dir, "alpha_lmm_predictions.csv"))
write_fixture(data.frame(metric = alpha_effect_fixture$metric, diagnostic_status = "synthetic"), file.path(alpha_dir, "alpha_lmm_diagnostics_summary.csv"))
metrics <- c("bray_curtis", "jaccard_binary", "aitchison")
beta_fixture <- data.frame(metric = metrics, hypothesis = "timepoint", r_squared = c(0.1, 0.2, 0.3), statistic = c(2, 3, 4), p_value = c(0.1, 0.2, 0.3), p_adjusted_bh = 0.3, sample_count = nrow(alpha_fixture), permutations = 99)
beta_dir <- "results/beta_diversity"
write_fixture(beta_fixture, file.path(beta_dir, "statistics/permanova_results.tsv"), TRUE)
dispersion_fixture <- beta_fixture; dispersion_fixture$hypothesis <- "Groups"
write_fixture(dispersion_fixture, file.path(beta_dir, "statistics/dispersion_results.tsv"), TRUE)
write_fixture(data.frame(metric = metrics, estimate = 0.02, std_error = 0.01, conf_low = 0, conf_high = 0.04, sample_count = nrow(alpha_fixture), subject_count = 9, singular = FALSE), file.path(beta_dir, "statistics/longitudinal_lmm_results.tsv"), TRUE)
baseline_fixture <- do.call(rbind, lapply(metrics, function(metric) data.frame(metric = metric, sample_id = alpha_fixture$sampleID, subject_id = alpha_fixture$subjectID, timepoint = alpha_fixture$Timepoint, days_from_subject_baseline = alpha_fixture$days_from_baseline, distance_to_baseline = alpha_fixture$days_from_baseline / 12)))
write_fixture(baseline_fixture, file.path(beta_dir, "statistics/distance_to_baseline.tsv"), TRUE)
coordinates_fixture <- data.frame(sample_id = alpha_fixture$sampleID, subject_id = alpha_fixture$subjectID, timepoint = alpha_fixture$Timepoint, sample_date = as.Date("2020-01-01") + alpha_fixture$days_from_baseline, PCoA1 = sin(seq_len(nrow(alpha_fixture))), PCoA2 = cos(seq_len(nrow(alpha_fixture))))
for (metric in c("bray_curtis", "aitchison")) {
  write_fixture(coordinates_fixture, file.path(beta_dir, "ordination", paste0("pcoa_", metric, "_coordinates.tsv")), TRUE)
  write_fixture(data.frame(axis = 1:3, variance_percent = c(40, 30, 30)), file.path(beta_dir, "ordination", paste0("pcoa_", metric, "_eigenvalues.tsv")), TRUE)
}
write_fixture(tax, "data/raw/GDTB.csv")
dir.create(file.path(fixture, "data/processed"), recursive = TRUE)
saveRDS(d, file.path(fixture, "data/processed/cleaned_data.rds"))
dir.create(file.path(fixture, "config"))
yaml::write_yaml(list(model_decision_rules = list(pseudocount = list(value = 0.01))), file.path(fixture, "config/model_decision_rules.yml"))
mag_relative <- "results/hpc_per_mag_lmertest/tables/per_mag_lmertest_hpc_results.csv"
# Avoid inheriting an external research-data override in a synthetic test.
vars <- c("PUBLIC_DIVERSITY_INPUT_ROOT", "PUBLIC_DIVERSITY_OUTPUT_DIR", "PUBLIC_DIVERSITY_RASTER_DPI", "BLACKICE_RAW_DATA_DIR", "BLACKICE_PROCESSED_DATA_DIR")
old_values <- Sys.getenv(vars, unset = NA_character_)
Sys.unsetenv(vars)
run_composite <- function(scenario, required_error = NULL) {
  out <- file.path(root, paste0("composite-", scenario))
  Sys.setenv(PUBLIC_DIVERSITY_INPUT_ROOT = fixture, PUBLIC_DIVERSITY_OUTPUT_DIR = out, PUBLIC_DIVERSITY_RASTER_DPI = "60")
  env <- new.env(parent = globalenv())
  error <- tryCatch({sys.source(diversity_source, envir = env); NULL}, error = conditionMessage)
  if (!is.null(required_error)) {
    stopifnot(!is.null(error), grepl(required_error, error, fixed = TRUE), !dir.exists(out))
    return(invisible(NULL))
  }
  if (!is.null(error)) stop(scenario, ": ", error)
  for (branch in c("alpha_diversity", "beta_diversity", "per_mag_lmm")) {
    for (extension in c("pdf", "svg", "png", "tiff")) stopifnot(file.info(file.path(out, branch, paste0("main_figure_", branch, ".", extension)))$size > 0)
    for (extension in c("csv", "tsv", "md", "pdf", "svg", "png", "tiff")) stopifnot(file.info(file.path(out, branch, paste0("main_table_", branch, ".", extension)))$size > 0)
  }
  stopifnot(env$alpha_sample_count == nrow(alpha_fixture), env$alpha_subject_count == 9, env$primary_count == sum(composite_results$primary_analysis))
  stopifnot(env$fdr_count == sum(composite_results$primary_significant %in% TRUE))
  notes <- paste(readLines(file.path(out, "per_mag_lmm/main_table_per_mag_lmm.md")), collapse = " ")
  stopifnot(grepl(paste0(env$fdr_count, " supplied primary BH-FDR discovery flags"), notes, fixed = TRUE))
  stopifnot(grepl(paste0(nrow(env$top_mags), " smallest finite raw P values among ", env$primary_count, " primary valid models"), notes, fixed = TRUE))
  stopifnot(grepl(env$minimum_q_text, env$mag_panel_forest$labels$subtitle, fixed = TRUE))
  stopifnot(all(env$beta_table$permutations == 99), nrow(env$mag_table) == min(10L, sum(composite_results$primary_analysis)))
  status <- read.csv(file.path(out, "supplementary/supplementary_inventory_status.csv"))
  expected_status <- if (dir.exists(file.path(fixture, "results/alpha_diversity"))) "empty_optional_directory" else "missing_optional_directory"
  stopifnot(status$status[status$source_dir == "results/alpha_diversity"] == expected_status)
  cat("PASS composite", scenario, "\n")
  env
}
composite_results <- r
composite_results$nominal_p_lt_0_05 <- composite_results$time_p_value < 0.05
for (discoveries in c(TRUE, FALSE)) {
  composite_results$primary_significant <- discoveries
  composite_results$time_p_adj_bh_primary <- if (discoveries) 0.02 else 0.2
  write_fixture(composite_results, mag_relative)
  run_composite(paste0("discoveries-", discoveries))
}
# Fewer candidates, no nominal association and missing q-values must not fabricate evidence.
composite_results <- composite_results[1:4, ]
composite_results$time_p_value <- 0.5
composite_results$nominal_p_lt_0_05 <- FALSE
composite_results$primary_significant <- NA
composite_results$time_p_adj_bh_primary <- NA_real_
write_fixture(composite_results, mag_relative)
dir.create(file.path(fixture, "results/alpha_diversity"))
empty_case <- run_composite("no-nominal-missing-q")
stopifnot(nrow(empty_case$trajectory_selection) == 0, nrow(empty_case$trajectory_predictions) == 0, empty_case$minimum_q_text == "unavailable", grepl("unavailable", empty_case$screening_text))
composite_results$primary_analysis <- FALSE
write_fixture(composite_results, mag_relative)
zero_case <- run_composite("zero-primary")
stopifnot(nrow(zero_case$mag_table) == 0, nrow(zero_case$trajectory_predictions) == 0)
required_path <- file.path(fixture, alpha_dir, "alpha_lmm_fixed_effects.csv")
file.rename(required_path, paste0(required_path, ".held"))
run_composite("missing-required-file", "Missing required manuscript input(s)")
file.rename(paste0(required_path, ".held"), required_path)
beta_bad <- beta_fixture; beta_bad$r_squared <- NULL
write_fixture(beta_bad, file.path(beta_dir, "statistics/permanova_results.tsv"), TRUE)
run_composite("missing-required-column", "missing required columns: r_squared")
write_fixture(beta_fixture, file.path(beta_dir, "statistics/permanova_results.tsv"), TRUE)
Sys.unsetenv(vars)
for (i in seq_along(vars)) if (!is.na(old_values[i])) do.call(Sys.setenv, setNames(list(old_values[i]), vars[i]))

# RF feature stability: artificial values satisfy the renderer's existing
# 441-feature/13-fold input contract; no cohort values or identities are copied.
rf_relative <- "scripts/plot_random_forest_feature_stability.R"
rf_source <- sub("[.]R$", "_public.R", file.path(repo, rf_relative))
if (!file.exists(rf_source)) {
  rf_source <- file.path(repo, rf_relative)
  map <- read.csv(file.path(repo, "docs/releases/v1.0.0-public-replacements.csv"))
  row <- map[map$public_destination_path == rf_relative, ]
  stopifnot(nrow(row) == 1L, identical(digest::digest(file = rf_source, algo = "sha256"), row$replacement_sha256))
}
# Reading a numeric summary must preserve its full precision, including NA.
probability_env <- new.env(parent = globalenv())
for (expression in parse(rf_source)) {
  if (is.call(expression) && identical(expression[[1]], as.name("<-")) && identical(expression[[2]], as.name("validate_probability"))) eval(expression, probability_env)
}
precision_fixture <- c(0.12345678901234567, 0.9876543210987654, NA_real_)
stopifnot(identical(precision_fixture, probability_env$validate_probability(precision_fixture, "artificial probability")))
rf_root <- file.path(root, "artificial-rf-stability")
rf_features <- file.path(rf_root, "features")
dir.create(rf_features, recursive = TRUE)
dir.create(file.path(rf_root, "reports"))
dir.create(file.path(rf_root, "predictions"))
feature_ids_test <- sprintf("ARTIFICIAL_MAG_%04d", seq_len(441))
subject_ids_test <- sprintf("ARTIFICIAL_SUBJECT_%02d", seq_len(13))
egrid <- expand.grid(MAG = feature_ids_test, heldout_subject = subject_ids_test, stringsAsFactors = FALSE)
fi <- match(egrid$MAG, feature_ids_test); si <- match(egrid$heldout_subject, subject_ids_test)
egrid$fold_id <- paste0("ARTIFICIAL_FOLD_", si)
egrid$eligible <- !(fi == 1 & si %in% c(1, 2))
egrid$selected <- egrid$eligible
egrid$missing_fraction <- 0
egrid$sample_prevalence <- ifelse(egrid$eligible, 0.8, 0.01)
egrid$subject_prevalence <- ifelse(egrid$eligible, 9, 1)
egrid$mean_abundance <- 0.1
egrid$variance <- 0.02
importance_test <- data.frame(feature = egrid$MAG[egrid$eligible], fold_id = egrid$fold_id[egrid$eligible], heldout_subject = egrid$heldout_subject[egrid$eligible], MeanDecreaseAccuracy = (sin(fi * 0.17) + cos(si * 0.41 + fi * 0.03) * (1 + fi %% 7))[egrid$eligible], MeanDecreaseGini = (1 + (fi %% 11) / 10)[egrid$eligible])
importance_test$MeanDecreaseAccuracy[importance_test$feature == feature_ids_test[1]] <- 50
# High performance does not require a particular sign of feature importance;
# both negative and entirely positive supplied values must be preserved.
write_rf_inputs <- function(positive_only = FALSE) {
  imp <- importance_test
  if (positive_only) imp$MeanDecreaseAccuracy <- imp$MeanDecreaseAccuracy + 10
  imp$rank <- ave(imp$MeanDecreaseAccuracy, imp$fold_id, FUN = function(x) rank(-x, ties.method = "average"))
  stab <- do.call(rbind, lapply(split(imp, imp$feature), function(z) data.frame(feature = z$feature[1], importance_folds = length(unique(z$fold_id)), median_native_importance = median(z$MeanDecreaseAccuracy), native_importance_iqr = IQR(z$MeanDecreaseAccuracy, type = 7), median_importance_rank = median(z$rank), importance_rank_iqr = IQR(z$rank, type = 7))))
  rownames(stab) <- NULL
  freq <- aggregate(cbind(eligible_folds = as.integer(egrid$eligible), selected_folds = as.integer(egrid$selected)), list(feature = egrid$MAG), sum)
  freq$selection_frequency_when_eligible <- freq$selected_folds / freq$eligible_folds
  candidates <- merge(stab, freq, by = "feature", sort = FALSE)
  candidates <- candidates[order(candidates$median_importance_rank, -candidates$median_native_importance, candidates$feature), ][1:20, ]
  imp$rank <- NULL
  for (name in c("candidate_MAG_summary", "feature_eligibility_by_fold", "feature_importance_by_fold", "feature_importance_stability", "feature_selection_frequency")) {
    obj <- switch(name, candidate_MAG_summary = candidates, feature_eligibility_by_fold = egrid, feature_importance_by_fold = imp, feature_importance_stability = stab, feature_selection_frequency = freq)
    write.csv(obj, file.path(rf_features, paste0(name, ".csv")), row.names = FALSE)
  }
  invisible(imp)
}
tax_rf <- data.frame(Genome = feature_ids_test, label = sprintf("ARTIFICIAL_LABEL_%04d", seq_along(feature_ids_test)), domain = "d__Artificial", phylum = "p__Artificial", class = "c__Artificial", ordeer = "o__Artificial", family = "f__Artificial", genus = "g__Artificial", species = "s__Artificial")
rf_taxonomy <- file.path(rf_root, "artificial-taxonomy.csv"); write.csv(tax_rf, rf_taxonomy, row.names = FALSE)
yaml::write_yaml(list(preprocessing = list(top_n = NULL, max_missing_fraction = 0.2, sample_prevalence = 0.1, subject_prevalence_min = 2, min_mean_abundance = 0)), file.path(rf_root, "reports/resolved_config.yml"))
rf_overall <- file.path(rf_root, "predictions/loso_overall_metrics.csv")
rf_class <- file.path(rf_root, "predictions/loso_per_class_metrics.csv")
class_names <- paste0("ARTIFICIAL_CLASS_", LETTERS[1:4])
run_rf_renderer <- function(name, expected_error = NULL) {
  out <- file.path(root, paste0("rf-render-", name)); log <- file.path(root, paste0("rf-render-", name, ".log"))
  old_dpi <- Sys.getenv("PUBLIC_RF_RASTER_DPI", unset = NA_character_)
  Sys.setenv(PUBLIC_RF_RASTER_DPI = "60")
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(rf_source, "--feature-results-dir", rf_features, "--taxonomy-file", rf_taxonomy, "--output-dir", out)), stdout = log, stderr = log)
  if (is.na(old_dpi)) Sys.unsetenv("PUBLIC_RF_RASTER_DPI") else Sys.setenv(PUBLIC_RF_RASTER_DPI = old_dpi)
  if (!is.null(expected_error)) {
    stopifnot(status != 0, any(grepl(expected_error, readLines(log), fixed = TRUE)))
    cat("PASS RF rejected", name, "\n"); return(invisible(NULL))
  }
  if (status != 0) stop(paste(readLines(log), collapse = "\n"))
  text <- paste(readLines(file.path(out, "figure_RF_MAG_feature_stability_caption.md")), readLines(file.path(out, "figure_RF_MAG_feature_stability_methods.md")), collapse = " ")
  stopifnot(!grepl("\\b(weak|strong|poor|excellent)\\b", text, ignore.case = TRUE))
  validation <- read.csv(file.path(out, "figure_generation_validation.csv"));stopifnot(all(validation$passed))
  # Exact input identities and importance values must survive presentation joins.
  plotted <- read.csv(file.path(out, "feature_importance_by_fold_with_taxonomy.csv"))
  input_imp <- read.csv(file.path(rf_features, "feature_importance_by_fold.csv"))
  stopifnot(identical(plotted$feature, input_imp$feature), identical(plotted$MeanDecreaseAccuracy, input_imp$MeanDecreaseAccuracy))
  candidates <- read.csv(file.path(rf_features, "candidate_MAG_summary.csv"))
  presented <- read.csv(file.path(out, "candidate_MAG_summary_with_taxonomy.csv"))
  stopifnot(identical(candidates$feature, presented$feature))
  cat("PASS RF renderer", name, "\n")
  text
}
write_rf_inputs()
write.csv(data.frame(balanced_accuracy = 0.2, macro_f1 = 0.18), rf_overall, row.names = FALSE)
write.csv(data.frame(class = class_names, recall = c(0, 0.2, 0.2, 0.4)), rf_class, row.names = FALSE)
low_text <- run_rf_renderer("low");stopifnot(grepl("2 displayed candidate-fold cells were ineligible", low_text, fixed = TRUE), grepl("0.200", low_text), grepl("0.180", low_text), grepl("1/4 classes with available recall", low_text, fixed = TRUE))
write_rf_inputs(TRUE)
write.csv(data.frame(balanced_accuracy = 0.93, macro_f1 = 0.92), rf_overall, row.names = FALSE)
write.csv(data.frame(class = class_names, recall = c(0.9, 0.92, 0.94, 0.96)), rf_class, row.names = FALSE)
high_text <- run_rf_renderer("high-positive-importance");stopifnot(grepl("0.930", high_text), grepl("0.920", high_text), grepl("0/4 classes with available recall", high_text, fixed = TRUE))
file.rename(rf_overall, paste0(rf_overall, ".held"));file.rename(rf_class, paste0(rf_class, ".held"))
missing_text <- run_rf_renderer("missing");stopifnot(grepl("balanced accuracy = unavailable", missing_text, fixed = TRUE), grepl("Class recall unavailable", missing_text, fixed = TRUE))
file.rename(paste0(rf_class, ".held"), rf_class)
write.csv(data.frame(balanced_accuracy = 0.93), rf_overall, row.names = FALSE)
write.csv(data.frame(class = class_names, recall = c(0, NA, 0.5, NA)), rf_class, row.names = FALSE)
partial_text <- run_rf_renderer("partial");stopifnot(grepl("macro F1 = unavailable", partial_text, fixed = TRUE), grepl("1/2 classes with available recall", partial_text, fixed = TRUE), grepl("unavailable for 2 additional classes", partial_text, fixed = TRUE))
write.csv(data.frame(balanced_accuracy = 1.5, macro_f1 = 0.5), rf_overall, row.names = FALSE)
run_rf_renderer("invalid-probability", "must contain probabilities in [0, 1]")
write.csv(data.frame(balanced_accuracy = 0.93, macro_f1 = 0.92), rf_overall, row.names = FALSE)
required_rf <- file.path(rf_features, "feature_importance_stability.csv")
file.rename(required_rf, paste0(required_rf, ".held"));run_rf_renderer("missing-required", "Missing feature input(s)");file.rename(paste0(required_rf, ".held"), required_rf)
bad_eligibility <- egrid;bad_eligibility$eligible <- NULL
write.csv(bad_eligibility, file.path(rf_features, "feature_eligibility_by_fold.csv"), row.names = FALSE)
run_rf_renderer("missing-required-column", "missing required column(s): eligible")
write.csv(egrid[-1, ], file.path(rf_features, "feature_eligibility_by_fold.csv"), row.names = FALSE)
run_rf_renderer("incomplete-grid", "Eligibility table is not the complete")
