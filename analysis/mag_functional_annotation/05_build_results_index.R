#!/usr/bin/env Rscript

# Purpose: Build a human- and AI-readable index of every non-intermediate
# functional-analysis result produced on Nebula.
# Inputs: Results registry and generated files from stages 00-04.
# Outputs: RESULTS_INDEX.md and RESULTS_INDEX.csv in the workflow output root.
# Dependencies: R, data.table, yaml, and common.R.
# Usage: Rscript analysis/mag_functional_annotation/05_build_results_index.R --config config/mag_functional_annotation.yml
# Workflow: Final documentation stage; does not alter scientific results.

source("analysis/mag_functional_annotation/R/common.R")
args <- parse_stage_args()
cfg <- read_magfunc_config(args$config)
assert_packages()
paths <- output_paths(cfg)
if (!file.exists(file.path(paths$qc, "04_make_publication_outputs.complete"))) stop("Stage 04 completion marker is missing.", call. = FALSE)
guard_stage(file.path(paths$root, "RESULTS_INDEX.md"), args, cfg)

session_info_path <- write_session_info(cfg, "05_build_results_index")
source("scripts/utils/write_hpc_execution_record.R")
inputs <- resolve_input_paths(cfg)
model_table <- file.path(paths$longitudinal, "Table_06_longitudinal_functional_models.csv.gz")
number_of_tests <- if (file.exists(model_table)) nrow(read_table(model_table, select = "term")) else NA_integer_
write_hpc_execution_record(
  analysis_name = "MAG functional annotation and longitudinal predicted functional potential",
  script_path = "scripts/run_mag_functional_annotation.sh",
  result_dir = paths$root,
  input_data_file = unlist(inputs, use.names = FALSE),
  config_files = cfg$.config_path,
  number_of_tests = number_of_tests,
  fdr_family = "BH within each ontology and model family or paired contrast",
  main_results_file = c(model_table, file.path(paths$hysteresis, "Table_07_equal_altitude_paired_contrasts.csv.gz")),
  session_info_file = session_info_path,
  notes = paste(
    "Ta/Tb excluded; gut.csv values not row-renormalized;",
    "protein member counts are not abundance;",
    "equal-altitude contrasts are hysteresis-like and do not establish adaptation."
  )
)

registry_path <- file.path(paths$intermediate, "results_registry.csv")
registry <- data.table::fread(registry_path)
all_files <- list.files(paths$root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
all_files <- all_files[file.info(all_files)$isdir %in% FALSE]
all_files <- all_files[!grepl("(^|/)intermediate/|RESULTS_INDEX\\.(md|csv)$", all_files)]
relative_files <- vapply(all_files, relative_to_root, character(1))
inventory <- data.table(path = relative_files)
inventory[, filename := basename(path)]
inventory[, extension := tolower(tools::file_ext(filename))]
inventory[, size_bytes := file.info(all_files)$size]
inventory[, artifact_type := data.table::fcase(
  extension %in% c("pdf", "png", "svg", "tiff"), "figure",
  extension %in% c("csv", "tsv", "gz"), "table_or_data",
  extension %in% c("md", "txt"), "documentation",
  extension == "rds", "model_ready_intermediate",
  default = "other"
)]

index <- merge(inventory, registry, by = "path", all.x = TRUE, suffixes = c("_inventory", "_registry"))
index[is.na(result_number), result_number := "Supporting"]
output_root_relative <- relative_to_root(paths$root)
supporting_relative <- sub(paste0("^", gsub("([][{}()+*^$|\\?.])", "\\\\\\1", output_root_relative), "/?"), "", index$path)
missing_category <- which(is.na(index$category))
index[missing_category, category := sub("/.*$", "", supporting_relative[missing_category])]
index[is.na(question), question := "Supporting output for the registered analysis stage"]
index[is.na(description), description := "Generated supporting, diagnostic, or reproducibility artifact."]
index[is.na(stage), stage := "supporting"]
index[is.na(artifact_type_registry), artifact_type_registry := artifact_type_inventory]
index[, artifact_type := artifact_type_registry]
index[, c("artifact_type_inventory", "artifact_type_registry") := NULL]
index[, supporting_sort := result_number == "Supporting"]
data.table::setorder(index, supporting_sort, result_number, category, path)
index[, supporting_sort := NULL]

csv_path <- file.path(paths$root, "RESULTS_INDEX.csv")
md_path <- file.path(paths$root, "RESULTS_INDEX.md")
write_table(index, csv_path)

escape_md <- function(x) gsub("\\|", "\\\\|", as.character(x))
lines <- c(
  "# MAG Functional Annotation Results Index",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Analysis contract",
  "",
  "- Complete raw files and generated outputs remain on Nebula and are not Git evidence.",
  "- `Ta` and `Tb` are excluded; longitudinal analyses use `T1` through `T10`.",
  "- `gut.csv` values are used exactly as provided and are not row-renormalized.",
  "- Protein-cluster member counts are catalogue membership, not abundance.",
  "- Sample-level values are predicted functional potential, not direct gene abundance, expression, activity, or metabolite production.",
  "- Equal-altitude paired contrasts are hysteresis-like tests, not proof of a continuous hysteresis loop or adaptation.",
  "",
  "## Indexed results",
  "",
  "Number | Category | Type | File | Path | Purpose",
  "--- | --- | --- | --- | --- | ---"
)
for (i in seq_len(nrow(index))) {
  lines <- c(lines, paste(
    escape_md(index$result_number[[i]]),
    escape_md(index$category[[i]]),
    escape_md(index$artifact_type[[i]]),
    paste0("`", escape_md(index$filename[[i]]), "`"),
    paste0("`", escape_md(index$path[[i]]), "`"),
    escape_md(index$description[[i]]),
    sep = " | "
  ))
}
lines <- c(
  lines,
  "",
  "## Interpretation order",
  "",
  "1. Review Tables 01-02 and all files under `00_quality_control/` before biological interpretation.",
  "2. Review catalogue questions in Tables 03-05 and Figures 01-04.",
  "3. Review longitudinal model support in Table 06 before using Figure 05.",
  "4. Review paired-subject counts and evidence labels in Table 07 before using Figure 06.",
  "5. Consult `reproducibility/` for session information and the chart contract map."
)
writeLines(lines, md_path)
message("Stage 05 completed: ", md_path)
