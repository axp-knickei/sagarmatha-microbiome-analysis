#!/usr/bin/env Rscript

# Purpose: Construct sample-level predicted functional potential from the
# unnormalized gut MAG table, fit participant-aware longitudinal models, and
# test prespecified equal-recorded-altitude paired contrasts.
# Inputs: Stage-02 MAG-function tables, gut.csv, and sample_info.csv.
# Outputs: Functional-potential matrices, longitudinal model tables, timepoint
# summaries, and hysteresis-like paired-contrast tables.
# Dependencies: R, data.table, Matrix, lme4, yaml, and common.R.
# Usage: Rscript analysis/mag_functional_annotation/03_longitudinal_analyses.R --config config/mag_functional_annotation.yml
# Workflow: Implements priority questions 4-5; gut.csv rows are never renormalized.

source("analysis/mag_functional_annotation/R/common.R")
args <- parse_stage_args()
cfg <- read_magfunc_config(args$config)
assert_packages("lme4")
inputs <- resolve_input_paths(cfg)
paths <- output_paths(cfg)
if (!file.exists(file.path(paths$qc, "02_catalogue_analyses.complete"))) stop("Stage 02 completion marker is missing.", call. = FALSE)
guard_stage(file.path(paths$qc, "03_longitudinal_analyses.complete"), args, cfg)
set.seed(as.integer(cfg$analysis$random_seed))

gut <- data.table::fread(inputs$abundance_file)
sample_col <- cfg$columns$sample_id
gut_mags <- setdiff(names(gut), sample_col)
metadata_all <- data.table::fread(inputs$metadata_file)
missing_metadata <- setdiff(as.character(gut[[sample_col]]), as.character(metadata_all[[sample_col]]))
if (length(missing_metadata) > 0) stop(length(missing_metadata), " gut sample IDs have no metadata row.", call. = FALSE)
metadata <- prepare_metadata(inputs$metadata_file, cfg)
metadata_ids <- as.character(metadata[[sample_col]])

eligible_ids <- intersect(as.character(gut[[sample_col]]), metadata_ids)
gut <- gut[get(sample_col) %in% eligible_ids]
data.table::setorderv(gut, sample_col)
metadata <- metadata[get(sample_col) %in% eligible_ids]
data.table::setorderv(metadata, sample_col)
if (!identical(as.character(gut[[sample_col]]), as.character(metadata[[sample_col]]))) stop("Gut and metadata sample order could not be aligned.", call. = FALSE)

abundance <- as.matrix(gut[, ..gut_mags])
storage.mode(abundance) <- "double"
if (!identical(cfg$analysis$abundance_row_normalization, "none")) stop("This workflow is approved only for abundance_row_normalization: none.", call. = FALSE)
row_sums <- rowSums(abundance)
sample_ids <- as.character(gut[[sample_col]])
longitudinal_ontologies <- unlist(cfg$analysis$longitudinal_ontologies)
threads <- max(1L, as.integer(Sys.getenv("MAGFUNC_THREADS", unset = cfg$analysis$threads)))

fit_one_function <- function(target_term, dt, carrier_count) {
  d <- dt[term == target_term]
  n_samples <- nrow(d)
  n_subjects <- data.table::uniqueN(d$Subject)
  n_nonzero <- sum(d$potential > 0)
  response_variance <- stats::var(d$potential)
  if (n_nonzero < as.integer(cfg$analysis$minimum_nonzero_samples) || n_subjects < as.integer(cfg$analysis$minimum_subjects) || !is.finite(response_variance) || response_variance == 0) {
    return(data.table(
      term = target_term, carrier_mags = carrier_count, n_samples = n_samples, n_subjects = n_subjects,
      n_nonzero_samples = n_nonzero, pseudocount = NA_real_, omnibus_chisq = NA_real_, omnibus_df = NA_real_,
      omnibus_p_value = NA_real_, day_estimate = NA_real_, day_std_error = NA_real_, day_p_value = NA_real_,
      singular = NA, status = "insufficient_variation_or_support"
    ))
  }
  positive <- d$potential[d$potential > 0]
  pseudocount <- min(positive) / 2
  d[, response := log10(potential + pseudocount)]
  d[, Timepoint_unordered := factor(as.character(Timepoint), levels = unlist(cfg$analysis$ordered_timepoints))]
  result <- tryCatch({
    control <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
    null_model <- lme4::lmer(response ~ 1 + (1 | Subject), data = d, REML = FALSE, control = control)
    full_model <- lme4::lmer(response ~ Timepoint_unordered + (1 | Subject), data = d, REML = FALSE, control = control)
    comparison <- stats::anova(null_model, full_model)
    day_null_model <- lme4::lmer(response ~ 1 + (1 | Subject), data = d, REML = FALSE, control = control)
    day_model <- lme4::lmer(response ~ days_from_baseline_subject + (1 | Subject), data = d, REML = FALSE, control = control)
    day_comparison <- stats::anova(day_null_model, day_model)
    day_coef <- summary(day_model)$coefficients["days_from_baseline_subject", ]
    data.table(
      term = target_term,
      carrier_mags = carrier_count,
      n_samples = n_samples,
      n_subjects = n_subjects,
      n_nonzero_samples = n_nonzero,
      pseudocount = pseudocount,
      omnibus_chisq = comparison$Chisq[[2]],
      omnibus_df = as.numeric(comparison$Df[[2]]),
      omnibus_p_value = comparison$`Pr(>Chisq)`[[2]],
      day_estimate = unname(day_coef[["Estimate"]]),
      day_std_error = unname(day_coef[["Std. Error"]]),
      day_p_value = day_comparison$`Pr(>Chisq)`[[2]],
      singular = lme4::isSingular(full_model, tol = 1e-4) || lme4::isSingular(day_model, tol = 1e-4),
      status = "fitted"
    )
  }, error = function(e) {
    data.table(
      term = target_term, carrier_mags = carrier_count, n_samples = n_samples, n_subjects = n_subjects,
      n_nonzero_samples = n_nonzero, pseudocount = pseudocount, omnibus_chisq = NA_real_, omnibus_df = NA_real_,
      omnibus_p_value = NA_real_, day_estimate = NA_real_, day_std_error = NA_real_, day_p_value = NA_real_,
      singular = NA, status = paste0("model_error: ", conditionMessage(e))
    )
  })
  result
}

run_hysteresis <- function(potential_long, model_results, ontology) {
  contrast_results <- lapply(cfg$hysteresis_contrasts, function(contrast) {
    before <- as.character(contrast$before)
    after <- as.character(contrast$after)
    subset <- potential_long[as.character(Timepoint) %in% c(before, after)]
    wide <- data.table::dcast(subset, Subject + term ~ Timepoint, value.var = "potential")
    if (!all(c(before, after) %in% names(wide))) return(data.table())
    pseudocounts <- potential_long[, {
      positive <- potential[potential > 0]
      list(pseudocount = if (length(positive) > 0) min(positive) / 2 else NA_real_)
    }, by = term]
    wide <- merge(wide, pseudocounts, by = "term", all.x = TRUE)
    wide <- wide[is.finite(get(before)) & is.finite(get(after)) & is.finite(pseudocount)]
    if (nrow(wide) == 0) return(data.table())
    wide[, delta_log10 := log10(get(after) + pseudocount) - log10(get(before) + pseudocount)]
    results <- wide[, {
      n_pairs <- .N
      if (n_pairs < as.integer(cfg$analysis$minimum_paired_subjects) || stats::sd(delta_log10) == 0) {
        list(n_pairs = n_pairs, mean_log10_difference = mean(delta_log10), standard_error = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_, paired_ratio = 10^mean(delta_log10), p_value = NA_real_, wilcoxon_p_value = NA_real_, status = "insufficient_pairs_or_variation")
      } else {
        test <- stats::t.test(delta_log10, mu = 0)
        wilcox <- tryCatch(stats::wilcox.test(delta_log10, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
        list(n_pairs = n_pairs, mean_log10_difference = mean(delta_log10), standard_error = stats::sd(delta_log10) / sqrt(n_pairs), lower_95 = unname(test$conf.int[[1]]), upper_95 = unname(test$conf.int[[2]]), paired_ratio = 10^mean(delta_log10), p_value = test$p.value, wilcoxon_p_value = wilcox, status = "fitted")
      }
    }, by = term]
    results[, `:=`(
      ontology = ontology,
      contrast_id = as.character(contrast$id),
      before_timepoint = before,
      after_timepoint = after,
      recorded_altitude_m = as.numeric(contrast$recorded_altitude_m),
      primary_contrast = isTRUE(contrast$primary),
      evidence = as.character(contrast$evidence)
    )]
    results[, p_adj_bh := safe_bh(p_value)]
    results[, wilcoxon_p_adj_bh := safe_bh(wilcoxon_p_value)]
    results
  })
  data.table::rbindlist(contrast_results, fill = TRUE)
}

all_models <- list()
all_timepoint_summaries <- list()
all_hysteresis <- list()
potential_dir <- file.path(paths$intermediate, "functional_potential")
dir.create(potential_dir, recursive = TRUE, showWarnings = FALSE)

for (ontology in longitudinal_ontologies) {
  message("Constructing unnormalized predicted functional potential: ", ontology)
  mag_function_path <- file.path(paths$intermediate, "mag_function_tables", paste0(ontology, ".rds"))
  mag_function <- read_object(mag_function_path)
  carriers <- mag_function[, .(carrier_mags = data.table::uniqueN(MAG)), by = term]
  eligible_terms <- carriers[carrier_mags >= as.integer(cfg$analysis$minimum_mag_carriers), term]
  mag_function <- unique(mag_function[term %in% eligible_terms, .(MAG, term)])
  terms <- sort(unique(mag_function$term))
  if (length(terms) == 0) {
    warning("No eligible terms for ontology ", ontology, call. = FALSE)
    next
  }
  sparse_presence <- Matrix::sparseMatrix(
    i = match(mag_function$MAG, gut_mags),
    j = match(mag_function$term, terms),
    x = 1,
    dims = c(length(gut_mags), length(terms)),
    dimnames = list(gut_mags, terms)
  )
  potential_matrix <- as.matrix(abundance %*% sparse_presence)
  potential_long <- data.table(
    sampleID = rep(sample_ids, times = length(terms)),
    term = rep(terms, each = length(sample_ids)),
    potential = as.vector(potential_matrix)
  )
  potential_long[, ontology := ontology]
  metadata_for_join <- metadata[, .(
    sampleID = get(sample_col), Subject, Timepoint, analysis_date,
    altitude = altitude_numeric, days_from_baseline_subject
  )]
  potential_long <- merge(potential_long, metadata_for_join, by = "sampleID", all.x = TRUE)
  if (anyNA(potential_long$Subject)) stop("Functional-potential rows failed the sample metadata join.", call. = FALSE)
  save_object(potential_long, file.path(potential_dir, paste0(ontology, ".rds")))
  write_table(potential_long, file.path(paths$longitudinal, paste0("predicted_functional_potential_", ontology, ".csv.gz")))

  timepoint_summary <- potential_long[, .(
    n_samples = .N,
    n_subjects = data.table::uniqueN(Subject),
    mean_potential = mean(potential),
    standard_error = stats::sd(potential) / sqrt(.N),
    median_potential = stats::median(potential),
    nonzero_fraction = mean(potential > 0),
    abundance_sum_min = min(row_sums[match(sampleID, sample_ids)]),
    abundance_sum_max = max(row_sums[match(sampleID, sample_ids)])
  ), by = .(ontology, term, Timepoint)]
  all_timepoint_summaries[[ontology]] <- timepoint_summary

  carrier_lookup <- carriers[term %in% terms]
  fit_terms <- terms
  fitted <- parallel::mclapply(fit_terms, function(target_term) {
    fit_one_function(target_term, potential_long, carrier_lookup[term == target_term, carrier_mags][[1]])
  }, mc.cores = threads, mc.preschedule = TRUE)
  model_results <- data.table::rbindlist(fitted, fill = TRUE)
  model_results[, ontology := ontology]
  model_results[, omnibus_p_adj_bh := safe_bh(omnibus_p_value)]
  model_results[, day_p_adj_bh := safe_bh(day_p_value)]
  model_results[, significant_omnibus := !is.na(omnibus_p_adj_bh) & omnibus_p_adj_bh <= as.numeric(cfg$analysis$fdr_threshold)]
  model_results[, significant_day := !is.na(day_p_adj_bh) & day_p_adj_bh <= as.numeric(cfg$analysis$fdr_threshold)]
  all_models[[ontology]] <- model_results
  all_hysteresis[[ontology]] <- run_hysteresis(potential_long, model_results, ontology)
  rm(potential_matrix, potential_long, sparse_presence, mag_function); invisible(gc())
}

timepoint_summaries <- data.table::rbindlist(all_timepoint_summaries, fill = TRUE)
model_results <- data.table::rbindlist(all_models, fill = TRUE)
hysteresis_results <- data.table::rbindlist(all_hysteresis, fill = TRUE)
data.table::setorder(model_results, omnibus_p_adj_bh, omnibus_p_value, ontology, term, na.last = TRUE)
if (nrow(hysteresis_results) > 0) data.table::setorder(hysteresis_results, -primary_contrast, p_adj_bh, p_value, na.last = TRUE)

longitudinal_path <- file.path(paths$longitudinal, "Table_06_longitudinal_functional_models.csv.gz")
hysteresis_path <- file.path(paths$hysteresis, "Table_07_equal_altitude_paired_contrasts.csv.gz")
write_table(timepoint_summaries, file.path(paths$longitudinal, "functional_potential_timepoint_summaries.csv.gz"))
write_table(model_results, longitudinal_path)
write_table(hysteresis_results, hysteresis_path)
write_table(data.table(sampleID = sample_ids, abundance_sum_as_provided = row_sums), file.path(paths$longitudinal, "gut_abundance_sums_used_without_normalization.csv"))
writeLines(c(
  "Stage 03 complete",
  "Functional potential is the sum of unnormalized gut MAG values among MAGs carrying a function.",
  "It is predicted genomic potential, not direct gene abundance, expression, activity, or metabolite production.",
  "Equal-altitude paired contrasts are hysteresis-like tests and do not establish a continuous hysteresis loop."
), file.path(paths$qc, "03_longitudinal_analyses.complete"))

register_results(cfg, data.table::data.table(
  result_number = c("Table 06", "Table 07"),
  category = c("longitudinal_functional_potential", "ascent_descent_hysteresis"),
  question = c(
    "Does predicted functional potential vary across T1-T10 after accounting for repeated participants?",
    "Do equal-recorded-altitude functional profiles differ across prespecified expedition phases?"
  ),
  artifact_type = "table",
  path = vapply(c(longitudinal_path, hysteresis_path), relative_to_root, character(1)),
  description = c(
    "Omnibus timepoint and continuous-day mixed-model results using gut.csv values as provided; Ta/Tb are excluded.",
    "Within-participant log10 paired differences for prespecified 5200 m contrasts, with t-test and Wilcoxon sensitivity p-values."
  ),
  stage = "03_longitudinal_analyses"
))
write_session_info(cfg, "03_longitudinal_analyses")
message("Stage 03 completed.")
