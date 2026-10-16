setwd(Sys.getenv("BETA_PROJECT_ROOT"))
source("analysis/beta_diversity/utils/io_helpers.R")

test_that("ordination, statistics, figures, and skips are complete", {
  config <- load_beta_config("config/beta_diversity_test.yml")
  coords <- list.files(beta_output_dir(config, "ordination"), "pcoa_.*_coordinates.tsv$", full.names = TRUE)
  expect_gte(length(coords), 5)
  for (file in coords) {
    table <- read.delim(file)
    expect_gte(ncol(table), 3)
    expect_equal(nrow(table), 16)
  }
  permanova <- read.delim(beta_output_dir(config, "statistics", "permanova_results.tsv"))
  expect_true(all(c("metric", "model", "hypothesis", "permutation_scheme", "statistic", "r_squared", "p_value") %in% names(permanova)))
  expect_true(all(grepl("within-subject restricted permutations", permanova$permutation_scheme)))
  figures <- list.files(beta_output_dir(config, "figures"), "\\.(pdf|png)$", full.names = TRUE)
  expect_gte(length(figures), 10)
  expect_true(all(file.info(figures)$size > 0))
  manifest <- yaml::read_yaml(beta_output_dir(config, "distances", "distance_unweighted_unifrac_manifest.yml"))
  expect_identical(manifest$status, "SKIPPED_UNIFRAC_NO_TREE")
})
