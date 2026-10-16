library(testthat)
library(yaml)
Sys.setenv(NETMOSS_PROJECT_ROOT = normalizePath("."))
test_file("scripts/netmoss/tests/test-altitude-stability-public.R", reporter = "summary")
