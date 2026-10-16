#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(yaml) })
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 2L && args[[1]] == "--output-dir") args[[2]] else tempfile("blackice_funcdyn_smoke_")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_dir <- file.path(output_dir, "raw")
magfunc_dir <- file.path(output_dir, "mag_functional_annotation", "04_longitudinal_functional_potential")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(magfunc_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260812)
subjects <- paste0("S", 1:6)
timepoints <- paste0("T", 1:10)
meta <- rbindlist(lapply(subjects, function(subject) {
  available <- if (subject %in% c("S4", "S5", "S6")) timepoints[-1L] else timepoints
  index <- match(available, timepoints)
  data.table(
    sampleID = paste(subject, available, sep = "_"), group = subject, Timepoint = available,
    date = format(as.IDate("2026-01-01") + index * 7L + match(subject, subjects), "%Y%m%d")
  )
}))
fwrite(meta, file.path(raw_dir, "sample_info.csv"))

time_index <- match(meta$Timepoint, timepoints)
subject_index <- match(meta$group, subjects)
mag_names <- paste0("MAG_", sprintf("%02d", 1:12))
tax <- sapply(seq_along(mag_names), function(j) {
  phase <- sin((time_index - 1) / 9 * 2 * pi + j / 2)
  shift <- ifelse(j <= 4, (time_index - 1)^1.4 / 20, ifelse(j <= 8, (10 - time_index)^1.3 / 18, 0.25))
  pmax(0.0001, exp(-2 + 1.4 * phase + shift + subject_index * 0.025 + rnorm(nrow(meta), 0, 0.08)))
})
tax <- tax / rowSums(tax)
abundance <- cbind(meta[, .(sampleID)], as.data.table(tax))
setnames(abundance, c("sampleID", mag_names))
fwrite(abundance, file.path(raw_dir, "gut.csv"))

make_functional <- function(layer, terms, amplitude) {
  rows <- rbindlist(lapply(seq_along(terms), function(j) {
    broad <- 1 + amplitude * sin((time_index - 2) / 9 * 2 * pi + j / 3)
    recovery <- 1 + amplitude * 0.7 * exp(-((time_index - 6) / 2)^2) * ifelse(j %% 2L, 1, -0.4)
    subject <- 1 + subject_index * 0.015
    value <- pmax(1e-5, (0.03 + j / 300) * broad * recovery * subject + rnorm(nrow(meta), 0, 0.0006))
    data.table(sampleID = meta$sampleID, term = terms[[j]], potential = value, ontology = layer,
               Subject = meta$group, Timepoint = meta$Timepoint)
  }))
  fwrite(rows, file.path(magfunc_dir, paste0("predicted_functional_potential_", layer, ".csv")))
}
make_functional("KO", paste0("K", sprintf("%05d", 1:10)), 0.55)
make_functional("KEGG_module", paste0("M", sprintf("%05d", 1:7)), 0.42)
make_functional("EC", paste0("1.2.3.", 1:8), 0.36)
make_functional("COG", LETTERS[1:8], 0.12)

config <- list(
  project = list(name = "Functional community dynamics smoke test", decision_date = "2026-08-12"),
  inputs = list(
    abundance_file = file.path(raw_dir, "gut.csv"), metadata_file = file.path(raw_dir, "sample_info.csv"),
    mag_functional_results_dir = file.path(output_dir, "mag_functional_annotation")
  ),
  input_environment_overrides = list(abundance_file = "", metadata_file = "", mag_functional_results_dir = ""),
  columns = list(sample_id = "sampleID", subject_id = "group", timepoint = "Timepoint", date = "date"),
  analysis = list(
    excluded_timepoints = c("Ta", "Tb"), ordered_timepoints = timepoints,
    expected_sample_count = nrow(meta), enforce_expected_sample_count = TRUE,
    minimum_feature_prevalence = 0.05, minimum_features_per_layer = 2L,
    permutation_count = 99L, smoke_permutation_count = 99L, random_seed = 20260812L,
    fdr_method = "BH", fdr_threshold = 0.05, primary_layer = "KO", primary_metric = "hellinger",
    functional_layers = list(
      list(id = "KO", role = "primary", file = "04_longitudinal_functional_potential/predicted_functional_potential_KO.csv"),
      list(id = "KEGG_module", role = "sensitivity", file = "04_longitudinal_functional_potential/predicted_functional_potential_KEGG_module.csv"),
      list(id = "EC", role = "sensitivity", file = "04_longitudinal_functional_potential/predicted_functional_potential_EC.csv"),
      list(id = "COG", role = "descriptive", file = "04_longitudinal_functional_potential/predicted_functional_potential_COG.csv")
    ),
    metrics = list(list(id = "hellinger", role = "primary"), list(id = "bray_curtis", role = "sensitivity"), list(id = "aitchison", role = "sensitivity"))
  ),
  outputs = list(output_dir = file.path(output_dir, "results"), overwrite = FALSE, png_dpi = 100L, figure_width_in = 7.4, figure_height_in = 5.4)
)
config_path <- file.path(output_dir, "functional_community_dynamics_smoke.yml")
write_yaml(config, config_path)
cat(normalizePath(config_path), "\n", sep = "")
