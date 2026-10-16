# Purpose: Validate resampling and configuration with hand-constructed data only.
# Inputs: Public YAML configuration and synthetic matrices; no study files.
# Dependencies: testthat, yaml, and the altitude-stability helper dependencies.
# Run from repo root: NETMOSS_PROJECT_ROOT="$PWD" Rscript scripts/netmoss/tests/testthat-altitude-stability.R
# Outputs: Test assertions to stdout; real-cohort regression checks remain private.
project_root <- Sys.getenv("NETMOSS_PROJECT_ROOT")
source(file.path(project_root, "scripts/netmoss/altitude_stability_helpers.R"))

test_that("paired bootstrap and leave-one-out specifications are deterministic", {
  taxa <- paste0("x", 1:3)
  subjects <- paste0("S", 1:4)
  group_a <- matrix(seq_len(12), nrow = 3, dimnames = list(taxa, subjects))
  group_b <- group_a + 1
  first <- make_resample_specifications(group_a, group_b, 5L, 100L)
  second <- make_resample_specifications(group_a, group_b, 5L, 100L)
  expect_identical(first, second)
  expect_length(first, 9L)
  expect_equal(first[[1]]$subjects_a, first[[1]]$subjects_b)
  expect_equal(vapply(first[6:9], `[[`, character(1), "omitted"), subjects)
})

test_that("stability rules apply all predeclared gates", {
  config <- yaml::read_yaml(file.path(project_root, "config/netmoss_altitude_stability.yml"))
  expect_true(operational_network_pass(0.95, 0.75, 0.72, config))
  expect_false(operational_network_pass(0.89, 0.90, 0.90, config))
  expect_equal(
    operational_selection(c(TRUE, TRUE), c(0.9, 0.7), c(9, 11), 11, TRUE, config),
    c(TRUE, FALSE)
  )
})

test_that("configuration rejects overlapping biological groups", {
  config <- yaml::read_yaml(file.path(project_root, "config/netmoss_altitude_stability.yml"))
  expect_silent(validate_stability_config(config))
  config$design$group_b_timepoints <- c(unlist(config$design$group_b_timepoints), "T1")
  expect_error(validate_stability_config(config), "groups overlap")
})

test_that("module and failed-resample diagnostics are explicit", {
  expect_true(module_is_degenerate(c(1, 1, 1, 1), 0.8))
  expect_false(module_is_degenerate(c(1, 1, 2, 2), 0.8))
  expect_equal(adjusted_rand_index(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1)
  fake_reference <- list(
    network_a = diag(2), network_b = diag(2), scores = c(a = 1, b = 0),
    module_a = c(a = 1, b = 2), module_b = c(a = 1, b = 2)
  )
  config <- yaml::read_yaml(file.path(project_root, "config/netmoss_altitude_stability.yml"))
  invalid <- list(type = "bootstrap", id = "bad", omitted = NA_character_, subjects_a = "missing", subjects_b = "missing", seed = 1L)
  matrix <- matrix(1:4, nrow = 2, dimnames = list(c("a", "b"), c("S1", "S2")))
  result <- analyze_resample(invalid, matrix, matrix, fake_reference, config, "spearman")
  expect_false(result$success)
  expect_match(result$error, "subscript out of bounds")
})
