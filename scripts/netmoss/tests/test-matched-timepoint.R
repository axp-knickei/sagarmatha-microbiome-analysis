project_root <- Sys.getenv("NETMOSS_PROJECT_ROOT")
source(file.path(project_root, "scripts/netmoss/matched_timepoint_helpers.R"))

test_that("ten approved timepoints enumerate all 45 unique pairs", {
  timepoints <- paste0("T", 1:10)
  pairs <- make_timepoint_pairs(timepoints)
  expect_equal(nrow(pairs), 45L)
  expect_equal(length(unique(pairs$pair_id)), 45L)
  expect_equal(pairs$pair_id[[1]], "T1_vs_T2")
  expect_equal(pairs$pair_id[[45]], "T9_vs_T10")
  expect_false(any(grepl("Ta|Tb", pairs$pair_id)))
})

test_that("exact sign-flip p-values include both extreme directions", {
  changes <- rbind(all_positive = c(1, 1, 1), mixed = c(1, -1, 1))
  p_values <- exact_sign_flip_pvalues(changes)
  expect_equal(unname(p_values[[1]]), 0.25)
  expect_true(p_values[[2]] >= p_values[[1]])
})

test_that("leave-one-out stability detects direction robustness", {
  changes <- rbind(stable = c(2, 3, 4, 5), fragile = c(10, -1, -1, -1))
  stability <- leave_one_out_sign_stability(changes)
  expect_equal(unname(stability[[1]]), 1)
  expect_true(stability[[2]] < 1)
})
