library(testthat)
Sys.setenv(NETMOSS_PROJECT_ROOT = normalizePath("."))
test_dir("scripts/netmoss/tests", filter = "hpc|publication", reporter = "summary")
