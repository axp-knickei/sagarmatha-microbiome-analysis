library(testthat)
Sys.setenv(BETA_PROJECT_ROOT = normalizePath("."))
test_dir("analysis/beta_diversity/tests", filter = "validation|outputs", reporter = "summary")
