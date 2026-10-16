project_root <- Sys.getenv("NETMOSS_PROJECT_ROOT")

test_that("publication smoke mode produces source data and three figure formats", {
  output <- tempfile("netmoss_publication_test_")
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      file.path(project_root, "scripts/netmoss/create_publication_outputs.R"),
      "--config", file.path(project_root, "config/netmoss_analysis.yml"),
      "--output-dir", output,
      "--test-mode", "--force"
    ),
    stdout = FALSE,
    stderr = FALSE
  )
  expect_equal(status, 0L)
  expected <- c(
    "figures/Figure_1_NetMoss_network_rewiring.pdf",
    "figures/Figure_1_NetMoss_network_rewiring.svg",
    "figures/Figure_1_NetMoss_network_rewiring.tiff",
    "tables/Supplementary_Data_1_NetMoss_taxon_results.csv",
    "tables/Supplementary_Data_2_NetMoss_network_edges.csv",
    "tables/Source_Data_Figure_1.csv",
    "text/Figure_1_legend.txt",
    "text/NetMoss_methods_draft.txt",
    "publication_output_manifest.yml",
    "sessionInfo.txt"
  )
  paths <- file.path(output, expected)
  expect_true(all(file.exists(paths)))
  expect_true(all(file.info(paths)$size > 0))

  table <- read.csv(file.path(output, "tables/Supplementary_Data_1_NetMoss_taxon_results.csv"))
  expect_true(all(c("taxon_names", "netmoss_score", "q_value", "marker_status") %in% names(table)))
  expect_equal(nrow(table), 24L)
})

test_that("publication workflow refuses to relabel historical analysis outputs", {
  script <- readLines(
    file.path(project_root, "scripts/netmoss/create_publication_outputs.R"),
    warn = FALSE
  )
  expect_true(any(grepl("netmoss_manifest.yml", script, fixed = TRUE)))
  expect_true(any(grepl("validate_analysis_manifest", script, fixed = TRUE)))
  expect_true(any(grepl("Historical results will not be relabeled", script, fixed = TRUE)))
})
