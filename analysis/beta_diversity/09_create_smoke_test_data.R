#!/usr/bin/env Rscript
# Purpose: Generate deterministic, wholly synthetic beta-diversity test inputs.
# Inputs: Only the seed and column names in the supplied YAML; no study data.
# Outputs: Three CSV fixtures and README under data/test_beta_diversity/.
# Dependencies: R and the shared YAML/argument helpers (environment.yml).
# Usage: Rscript analysis/beta_diversity/09_create_smoke_test_data.R --config config/beta_diversity_test.yml
# Workflow: Test-data generation for make beta-test; never a scientific cohort.

source("analysis/beta_diversity/utils/io_helpers.R")
args <- parse_beta_args()
config <- load_beta_config(args$config)
set.seed(config$analysis$random_seed)

# Fictional IDs and dates prevent test fixtures from disclosing study records.
subjects <- sprintf("SYNTHETIC_SUBJECT_%02d", 1:4)
timepoints <- c("T2", "T3", "T7", "T10")
features <- sprintf("SYNTHETIC_MAG_%03d", 1:60)
metadata <- data.frame(
  sampleID = paste(rep(subjects, each = 4), rep(timepoints, 4), sep = "_"),
  group = rep(subjects, each = 4),
  Timepoint = rep(timepoints, 4),
  date = rep(c("20000101", "20000111", "20000121", "20000131"), 4),
  gender = rep(c("female", "male", "female", "male"), each = 4),
  age = rep(c(25, 35, 45, 55), each = 4),
  altitude = rep(c(100, 200, 300, 100), 4),
  stringsAsFactors = FALSE
)
names(metadata)[1:2] <- c(config$columns$sample_id, config$columns$subject_id %||% "group")

# Positive variation and structural zeros exercise abundance transformations;
# no rankings, values, dates or taxonomy are sampled from production inputs.
values <- matrix(stats::rgamma(16 * 60, shape = 2), nrow = 16, ncol = 60)
values[matrix(stats::runif(16 * 60) < 0.20, nrow = 16)] <- 0
values[, 1:2] <- values[, 1:2] + 1
values <- values / rowSums(values)
colnames(values) <- features
abundance <- data.frame(metadata[1], values, check.names = FALSE)
taxonomy <- data.frame(
  Genome = features, label = features, domain = "d__Synthetic",
  phylum = "p__Synthetic", class = "c__Synthetic", ordeer = "o__Synthetic",
  family = paste0("f__Synthetic_", rep(1:10, each = 6)),
  genus = paste0("g__Synthetic_", 1:60), species = paste0("s__Synthetic_", 1:60),
  species_note = "Fictional test taxonomy", stringsAsFactors = FALSE
)
names(taxonomy)[1] <- config$columns$taxonomy_feature_id
stopifnot(nrow(abundance) == 16, ncol(abundance) == 61,
          all(colSums(values > 0) >= 4), all(apply(values, 2, var) > 0))

out <- "data/test_beta_diversity"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(abundance, file.path(out, "gut_test.csv"), row.names = FALSE)
utils::write.csv(metadata, file.path(out, "sample_info_test.csv"), row.names = FALSE)
utils::write.csv(taxonomy, file.path(out, "GDTB_test.csv"), row.names = FALSE)
writeLines(c(
  "# Synthetic beta-diversity smoke-test data", "",
  sprintf("Generated with seed %s using only R random draws and fictional metadata.", config$analysis$random_seed),
  "Contains 16 samples, four fictional subjects, four timepoints and 60 fictional MAGs.",
  "No study input, participant record, real MAG identifier or taxonomy is copied.",
  "The synthetic dates, ages, sex labels and altitudes have no scientific meaning.",
  "No tree is supplied, allowing explicit UniFrac-skip behavior to be tested.", "",
  "Regenerate from the repository root:", "", "```bash",
  "Rscript analysis/beta_diversity/09_create_smoke_test_data.R --config config/beta_diversity_test.yml",
  "```", "",
  "Earlier Git history contains real-data-derived fixtures. Replacing current files",
  "does not remove those historical disclosures; review history before publication."
), file.path(out, "README.md"))
cat("COMPLETED_SYNTHETIC_FIXTURES: 16 samples x 60 MAGs\n")
