#!/usr/bin/env Rscript
# Validate production or smoke outputs from the manuscript synthesis workflow.

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/manuscript_functional_synthesis.yml"
if (length(args) == 2L && args[[1]] == "--config") config_path <- args[[2]]
if (!file.exists(config_path)) stop("Configuration not found: ", config_path, call. = FALSE)
cfg <- yaml::read_yaml(config_path)
root <- Sys.getenv("MANSYN_OUTPUT_DIR", unset = cfg$outputs$output_dir)
main <- file.path(root, "main")

required <- c(
  file.path(root, "00_quality_control", "00_build_main_figure_and_tables.complete"),
  file.path(root, "RESULTS_INDEX.md"),
  file.path(main, "Main_Figure_01_functional_catalogue_and_dynamics.pdf"),
  file.path(main, "Main_Figure_01_functional_catalogue_and_dynamics.png"),
  file.path(main, "Main_Figure_01_functional_catalogue_and_dynamics.svg"),
  file.path(main, "Main_Table_01_study_and_catalogue_summary.csv"),
  file.path(main, "Main_Table_02_key_statistical_results.csv")
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing synthesis outputs: ", paste(missing, collapse = ", "), call. = FALSE)
if (any(file.info(required)$size <= 0)) stop("One or more synthesis outputs are empty.", call. = FALSE)

table1 <- fread(file.path(main, "Main_Table_01_study_and_catalogue_summary.csv"))
table2 <- fread(file.path(main, "Main_Table_02_key_statistical_results.csv"))
if (nrow(table1) != 10L) stop("Main Table 1 must contain 10 rows.", call. = FALSE)
if (!all(c("K04028 timepoint screen", "EC 1.2.1.87 timepoint screen", "Taxonomic-KO transition coupling") %in% table2$Analysis)) {
  stop("Main Table 2 is missing required manuscript results.", call. = FALSE)
}
source("analysis/manuscript_functional_synthesis/R/common.R")
inputs <- resolve_synthesis_inputs(cfg)
mag_root <- inputs$mag_functional_results_dir
dyn_root <- inputs$functional_dynamics_results_dir
panel_c <- fread(file.path(root, "reproducibility", "Main_Figure_01_panel_c.csv"))
panel_e <- fread(file.path(root, "reproducibility", "Main_Figure_01_panel_e.csv"))
overlap <- read_result_table(file.path(mag_root, "reproducibility", "Figure_11_plot_data.csv"))
concordance <- read_result_table(file.path(mag_root, "10_candidate_annotation_audit", "Table_24_candidate_carrier_MAG_concordance.csv"))
classes <- c("KO:K04028 only", "Shared by both candidates", "EC:1.2.1.87 only")
expected <- c(sum(concordance$candidate_1_carried & !concordance$candidate_2_carried), sum(concordance$candidate_1_carried & concordance$candidate_2_carried), sum(concordance$candidate_2_carried & !concordance$candidate_1_carried))
if (anyNA(expected) || anyDuplicated(concordance$MAG) || sum(expected) != nrow(concordance) || !sum(expected) ||
    nrow(panel_c) != 3L || anyDuplicated(panel_c$display_class) || !setequal(panel_c$display_class, classes) ||
    nrow(overlap) != 3L || anyDuplicated(overlap$display_class) || !setequal(overlap$display_class, classes)) stop("Invalid carrier-set structure.", call. = FALSE)
setorder(panel_c, display_class)
expected <- expected[match(panel_c$display_class, classes)]
if (!all(panel_c$carrier_mags == expected) || !all(overlap$carrier_mags[match(panel_c$display_class, overlap$display_class)] == expected)) stop("Carrier counts disagree with source membership.", call. = FALSE)
union_n <- sum(expected)
shared_n <- panel_c[display_class == "Shared by both candidates", carrier_mags]
labels <- as.character(panel_c$carrier_mags)
labels[panel_c$display_class == "Shared by both candidates"] <- sprintf("%d MAGs\n%.1f%% of union", shared_n, 100 * shared_n / union_n)
if (!identical(panel_c$segment_label, labels) || any(abs(panel_c$union_fraction - expected / union_n) > 1e-12) ||
    any(panel_c$union_total != union_n) || any(abs(panel_c$shared_union_fraction - shared_n / union_n) > 1e-12) ||
    !identical(unique(panel_c$diagram_geometry), "schematic_equal_circles") || !identical(unique(panel_c$areas_proportional), FALSE)) stop("Carrier diagram labels, fractions or geometry are inconsistent.", call. = FALSE)
source_e <- read_result_table(file.path(dyn_root, "reproducibility", "Figure_18_plot_data.csv"))[layer == cfg$analysis$primary_functional_layer]
if (!isTRUE(all.equal(as.data.frame(panel_e), as.data.frame(source_e), check.attributes = FALSE))) stop("Transition panel differs from supplied primary-layer transitions.", call. = FALSE)
trajectories <- read_result_table(file.path(mag_root, "06_targeted_candidate_validation", "candidate_participant_trajectories.csv.gz"))
if (table1[Measure == "Included samples", Value] != as.character(uniqueN(trajectories$sampleID)) ||
    table1[Measure == "Participants", Value] != as.character(uniqueN(trajectories$Subject))) stop("Cohort counts disagree with supplied trajectories.", call. = FALSE)
singleton <- read_result_table(file.path(mag_root, "02_unknown_by_cluster_size", "singleton_unknown_enrichment_model.csv"))
expected_evidence <- sprintf("Descriptive catalogue model; p = %s; status = %s", format(singleton$p_value, scientific = TRUE, digits = 4), singleton$status)
if (!identical(table2[Analysis == "Singleton unknown enrichment", `Statistical evidence`], expected_evidence)) stop("Singleton evidence differs from supplied model output.", call. = FALSE)
svg_text <- toupper(paste(readLines(file.path(main, "Main_Figure_01_functional_catalogue_and_dynamics.svg"), warn = FALSE), collapse = "\n"))
for (colour in c("#487DB0", "#9D96BC", "#EBADC8")) {
  if (!grepl(toupper(colour), svg_text, fixed = TRUE)) {
    stop("Main Figure SVG is missing required GDV colour: ", colour, call. = FALSE)
  }
}
cat("Manuscript functional synthesis outputs validated.\n")
