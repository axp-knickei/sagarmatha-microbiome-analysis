#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(testthat))
testthat::test_dir("analysis/mag_functional_annotation/tests", reporter = "summary")
