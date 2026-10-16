library(testthat)
Sys.setenv(NETMOSS_PROJECT_ROOT = normalizePath("."))
test_file("scripts/netmoss/tests/test-matched-timepoint.R", reporter = "summary")
