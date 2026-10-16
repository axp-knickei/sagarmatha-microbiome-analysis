source(testthat::test_path("..", "R", "common.R"))

testthat::test_that("MAG identifiers parse only through approved terminal suffixes", {
  pattern <- "^(.+\\.(?:strict|permissive|orig))_.*$"
  ids <- c(
    "QML1234-000101-91B_bin.9.permissive_NODE_177_length_2881_cov_2.186837_2",
    "QML1234-000102-92C_bin.39.orig_FP400003218_L01_141_1_ctg184469_3",
    "malformed_protein"
  )
  expected <- c("QML1234-000101-91B_bin.9.permissive", "QML1234-000102-92C_bin.39.orig", NA_character_)
  testthat::expect_equal(parse_mag_id(ids, pattern), expected)
})

testthat::test_that("annotation term splitting preserves unique query-term relations", {
  annotation <- data.table::data.table(query = c("q1", "q2"), KEGG_ko = c("ko:K00001,ko:K00002", "ko:K00001"))
  spec <- list(id = "KO", column = "KEGG_ko", split_mode = "comma", strip_prefix = "ko:")
  terms <- split_annotation_terms(annotation, "query", spec)
  testthat::expect_equal(nrow(terms), 3L)
  testthat::expect_setequal(terms$term, c("K00001", "K00002"))
})

testthat::test_that("COG combined categories split into individual letters", {
  annotation <- data.table::data.table(query = "q1", COG_category = "GM")
  spec <- list(id = "COG", column = "COG_category", split_mode = "characters", strip_prefix = "")
  terms <- split_annotation_terms(annotation, "query", spec)
  testthat::expect_setequal(terms$term, c("G", "M"))
})

testthat::test_that("production decisions exclude Ta/Tb and forbid row normalization", {
  cfg <- read_magfunc_config(testthat::test_path("..", "..", "..", "config", "mag_functional_annotation.yml"))
  testthat::expect_setequal(unlist(cfg$analysis$excluded_timepoints), c("Ta", "Tb"))
  testthat::expect_identical(cfg$analysis$abundance_row_normalization, "none")
})

testthat::test_that("eggNOG annotation reader preserves tab fields with long descriptions", {
  cfg <- read_magfunc_config(testthat::test_path("..", "..", "..", "config", "mag_functional_annotation.yml"))
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)

  columns <- c(
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
    "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs",
    "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction",
    "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"
  )
  long_description <- paste(rep("ATP synthase description with many space-separated words", 12L), collapse = " ")
  row_one <- c(
    "QML1234-000105-95C_bin.20.orig_NODE_1_1", "1234.FICTIONAL_1", "1e-20",
    "340.0", "COG0712@1|root", "976|Bacteroidetes", "C", long_description,
    "atpH", "-", "-", "ko:K02113", "ko00190,map00190", "M00157", "-", "-",
    "ko00000", "3.A.2.1", "-", "-", "OSCP"
  )
  row_two <- replace(row_one, c(1L, 8L), c("QML1234-000105-95C_bin.20.orig_NODE_2_1", "short description"))
  writeLines(c(
    paste(columns, collapse = "\t"),
    paste(row_one, collapse = "\t"),
    paste(row_two, collapse = "\t")
  ), path)

  selected <- c("query", "Description", "KEGG_ko", "PFAMs")
  testthat::expect_warning(
    parsed <- read_annotation(path, cfg, select_columns = selected, nrows = Inf),
    NA
  )
  testthat::expect_identical(names(parsed), selected)
  testthat::expect_equal(nrow(parsed), 2L)
  testthat::expect_identical(parsed$Description[[1]], long_description)
  testthat::expect_identical(parsed$KEGG_ko[[1]], "ko:K02113")
})

testthat::test_that("grouped annotation-field medians use a consistent numeric type", {
  annotation_status <- data.table::data.table(
    size_group = c("singleton", "larger", "larger"),
    annotation_field_count = c(4L, 4L, 5L)
  )
  summary <- annotation_status[, .(
    median_annotation_fields = as.numeric(stats::median(annotation_field_count, na.rm = TRUE))
  ), by = size_group]

  testthat::expect_type(summary$median_annotation_fields, "double")
  testthat::expect_equal(summary[size_group == "singleton", median_annotation_fields], 4)
  testthat::expect_equal(summary[size_group == "larger", median_annotation_fields], 4.5)
})

testthat::test_that("compressed result tables read without the optional R.utils package", {
  testthat::skip_if(Sys.which("gzip") == "", "gzip command is unavailable")
  path <- tempfile(fileext = ".csv.gz")
  on.exit(unlink(path), add = TRUE)
  expected <- data.table::data.table(term = c("K00001", "K00002"), score = c(1.5, 2.5))
  data.table::fwrite(expected, path)

  observed <- read_table(path)
  testthat::expect_equal(observed, expected)
})

testthat::test_that("forced stages remove stale completion markers before work begins", {
  marker <- tempfile(fileext = ".complete")
  writeLines("stale completion", marker)
  cfg <- list(outputs = list(overwrite = FALSE))

  testthat::expect_invisible(guard_stage(marker, list(force = TRUE), cfg))
  testthat::expect_false(file.exists(marker))
})

testthat::test_that("follow-up candidates are explicit and unique", {
  source(testthat::test_path("..", "R", "followup.R"), local = TRUE)
  cfg <- read_magfunc_config(testthat::test_path("..", "..", "..", "config", "mag_functional_annotation.yml"))
  candidates <- configured_followup_candidates(cfg)

  testthat::expect_setequal(candidates$candidate_id, c("KO:K04028", "EC:1.2.1.87"))
  testthat::expect_false(anyDuplicated(candidates[, .(ontology, term)]) > 0L)
})

testthat::test_that("candidate curation is complete and bounded", {
  source(testthat::test_path("..", "R", "followup.R"), local = TRUE)
  cfg <- read_magfunc_config(testthat::test_path("..", "..", "..", "config", "mag_functional_annotation.yml"))
  curation <- configured_candidate_curation(cfg)

  testthat::expect_setequal(curation$candidate_id, c("KO:K04028", "EC:1.2.1.87"))
  testthat::expect_true(all(grepl("^https://", curation$authoritative_source_urls)))
  testthat::expect_true(all(nzchar(curation$candidate_claim_boundary)))
  testthat::expect_true(all(grepl("does not establish", curation$relationship_claim_boundary, fixed = TRUE)))
})

testthat::test_that("MAG construction provenance is parsed without becoming sample metadata", {
  source(testthat::test_path("..", "R", "followup.R"), local = TRUE)
  parsed <- parse_mag_construction_provenance(c(
    "QML1234-000101-91B_bin.9.permissive",
    "QML1234-000103-93C_bin.5.strict",
    "QML1234-000104-94D_bin.11.orig"
  ))

  testthat::expect_identical(parsed$construction_participant, c("91", "93", "94"))
  testthat::expect_identical(parsed$construction_specimen, c("fecal", "saliva", "skin"))
  testthat::expect_identical(parsed$construction_tier, c("permissive", "strict", "orig"))
})

testthat::test_that("within-subject correlation removes subject-level offsets", {
  source(testthat::test_path("..", "R", "followup.R"), local = TRUE)
  subject <- rep(c("S01", "S02"), each = 3L)
  x <- c(1, 2, 3, 101, 102, 103)
  y <- c(3, 2, 1, 103, 102, 101)

  testthat::expect_equal(within_subject_correlation(x, y, subject), -1)
})

testthat::test_that("generic timepoint response models retain finite omnibus results", {
  testthat::skip_if_not_installed("lme4")
  source(testthat::test_path("..", "R", "followup.R"), local = TRUE)
  cfg <- read_magfunc_config(testthat::test_path("..", "..", "..", "config", "mag_functional_annotation.yml"))
  d <- data.table::CJ(
    Subject = sprintf("S%02d", 1:6),
    Timepoint = unlist(cfg$analysis$ordered_timepoints)
  )
  set.seed(20260810)
  d[, response_value := as.numeric(factor(Timepoint)) +
    as.numeric(factor(Subject)) / 10 + stats::rnorm(.N, sd = 0.05)]

  fit <- fit_timepoint_response_model(d, cfg, "response_value")

  testthat::expect_identical(fit$status, "fitted")
  testthat::expect_true(is.finite(fit$p_value))
  testthat::expect_equal(fit$df, 9)
})

testthat::test_that("paired driver changes add exactly to the total paired change", {
  source(testthat::test_path("..", "R", "followup.R"), local = TRUE)
  d <- data.table::data.table(
    Subject = rep(c("S01", "S02"), each = 4L),
    Timepoint = rep(rep(c("T1", "T2"), each = 2L), 2L),
    feature = rep(c("A", "B"), 4L),
    value = c(1, 2, 3, 1, 2, 1, 4, 3)
  )

  drivers <- paired_change_drivers(
    d,
    before_timepoint = "T1",
    after_timepoint = "T2",
    id_columns = "feature",
    value_column = "value"
  )

  testthat::expect_equal(sum(drivers$mean_delta), 2.5)
  testthat::expect_equal(unique(drivers$net_change), 2.5)
  testthat::expect_equal(sum(drivers$absolute_change_share), 1)
  testthat::expect_true(all(drivers$n_pairs == 2L))
})
