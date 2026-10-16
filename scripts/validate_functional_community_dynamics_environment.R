#!/usr/bin/env Rscript
source("analysis/functional_community_dynamics/R/common.R")
packages <- assert_dynamics_packages()
cat("Functional community dynamics environment validated\n")
cat("R:", R.version.string, "\n")
for (package in packages) cat(package, as.character(utils::packageVersion(package)), "\n")
