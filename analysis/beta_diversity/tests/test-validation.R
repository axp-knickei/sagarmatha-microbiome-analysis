setwd(Sys.getenv("BETA_PROJECT_ROOT"))
source("analysis/beta_diversity/utils/io_helpers.R")
source("analysis/beta_diversity/utils/validation_helpers.R")
source("analysis/beta_diversity/utils/distance_helpers.R")

test_that("community and metadata invariants hold", {
  config <- load_beta_config("config/beta_diversity_test.yml")
  data <- readRDS(beta_output_dir(config, "prepared_data", "prepared_data.rds"))
  expect_false(anyDuplicated(rownames(data$community_matrix)) > 0)
  expect_false(anyDuplicated(colnames(data$community_matrix)) > 0)
  expect_true(is.numeric(data$community_matrix))
  expect_true(all(data$community_matrix >= 0))
  expect_true(all(rowSums(data$community_matrix) > 0))
  expect_identical(rownames(data$community_matrix), data$metadata$sample_id)
})

test_that("distance and transformation invariants hold", {
  config <- load_beta_config("config/beta_diversity_test.yml")
  files <- list.files(beta_output_dir(config, "distances"), "^distance_.*\\.rds$", full.names = TRUE)
  files <- files[!grepl("_transformed", files)]
  expect_gte(length(files), 5)
  for (file in files) {
    d <- readRDS(file); m <- as.matrix(d)
    expect_equal(m, t(m), tolerance = 1e-12)
    expect_equal(unname(diag(m)), rep(0, nrow(m)), tolerance = 1e-12)
    expect_true(all(is.finite(m)))
    expect_identical(rownames(m), colnames(m))
  }
  binary <- readRDS(beta_output_dir(config, "distances", "distance_jaccard_binary_transformed.rds"))
  expect_true(all(binary %in% c(0, 1)))
  clr <- readRDS(beta_output_dir(config, "distances", "distance_aitchison_transformed.rds"))
  expect_equal(unname(rowMeans(clr)), rep(0, nrow(clr)), tolerance = 1e-10)
  bounded <- c("bray_curtis", "jaccard_binary", "quantitative_jaccard_ruzička")
  for (metric in bounded) expect_true(all(as.numeric(readRDS(beta_output_dir(config, "distances", paste0("distance_", metric, ".rds")))) >= 0 & as.numeric(readRDS(beta_output_dir(config, "distances", paste0("distance_", metric, ".rds")))) <= 1))
})

test_that("toy distances match hand calculations", {
  toy <- rbind(a = c(1, 1, 0), b = c(1, 0, 1))
  expect_equal(as.numeric(vegan::vegdist(toy, "bray")), 0.5)
  expect_equal(as.numeric(vegan::vegdist(toy > 0, "jaccard", binary = TRUE)), 2 / 3)
  expect_equal(as.numeric(quantitative_jaccard_ruzička(toy)), 2 / 3)
})

test_that("fixed seed and permutation blocks are reproducible and aligned", {
  config <- load_beta_config("config/beta_diversity_test.yml")
  data <- readRDS(beta_output_dir(config, "prepared_data", "prepared_data.rds"))
  set.seed(config$analysis$random_seed); first <- vegan::vegdist(sqrt(total_sum_scale(data$community_matrix)), "bray")
  set.seed(config$analysis$random_seed); second <- vegan::vegdist(sqrt(total_sum_scale(data$community_matrix)), "bray")
  expect_identical(first, second)
  expect_identical(attr(first, "Labels"), data$metadata$sample_id)
  expect_length(data$metadata$subject_id, attr(first, "Size"))
})

test_that("model summaries are rejected as community matrices", {
  config <- load_beta_config("config/beta_diversity_test.yml")
  fake <- data.frame(sampleID = "MAG1", estimate = 1, p_value = 0.1)
  metadata <- data.frame(sampleID = "MAG1")
  taxonomy <- data.frame(Genome = "estimate")
  expect_error(validate_community(fake, metadata, taxonomy, config), "model-summary")
})
