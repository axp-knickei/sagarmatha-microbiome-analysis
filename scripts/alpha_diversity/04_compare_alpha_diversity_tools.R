# Compare vegan and scikit-bio alpha-diversity outputs and attach metadata.

required_packages <- c("dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

library(dplyr)

output_dir <- Sys.getenv("BLACKICE_ALPHA_OUTPUT_DIR", unset = "results/alpha_diversity")
tables_dir <- file.path(output_dir, "tables")
logs_dir <- file.path(output_dir, "logs")
vegan_path <- file.path(tables_dir, "alpha_diversity_vegan.csv")
skbio_path <- file.path(tables_dir, "alpha_diversity_scikitbio.csv")
metadata_path <- file.path(tables_dir, "alpha_diversity_metadata_prepared.csv")

if (!file.exists(vegan_path)) stop("Missing vegan output: ", vegan_path, call. = FALSE)
if (!file.exists(skbio_path)) stop("Missing scikit-bio output: ", skbio_path, call. = FALSE)
if (!file.exists(metadata_path)) stop("Missing prepared metadata: ", metadata_path, call. = FALSE)

vegan <- utils::read.csv(vegan_path, check.names = FALSE, stringsAsFactors = FALSE)
skbio <- utils::read.csv(skbio_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- utils::read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)

metrics <- c("Shannon", "Simpson", "Observed_richness", "Chao1", "ACE")
comparison <- vegan %>%
  select(sampleID, all_of(metrics)) %>%
  rename_with(~ paste0(.x, "_vegan"), all_of(metrics)) %>%
  inner_join(
    skbio %>%
      select(sampleID, all_of(metrics)) %>%
      rename_with(~ paste0(.x, "_scikitbio"), all_of(metrics)),
    by = "sampleID"
  )

for (metric in metrics) {
  lhs <- comparison[[paste0(metric, "_vegan")]]
  rhs <- comparison[[paste0(metric, "_scikitbio")]]
  comparison[[paste0(metric, "_absolute_difference")]] <- abs(lhs - rhs)
  comparison[[paste0(metric, "_consistent")]] <- ifelse(
    is.na(lhs) & is.na(rhs),
    TRUE,
    abs(lhs - rhs) <= 1e-8
  )
}

utils::write.csv(comparison, file.path(tables_dir, "alpha_diversity_tool_comparison.csv"), row.names = FALSE)

metadata_cols <- intersect(
  c("sampleID", "subjectID", "Timepoint", "date", "date_converted", "days_from_baseline", "altitude", "age", "gender", "group"),
  names(metadata)
)
with_metadata <- vegan %>%
  select(sampleID, all_of(metrics), tool, abundance_scale_decision, chao1_ace_decision) %>%
  left_join(metadata[, metadata_cols, drop = FALSE], by = "sampleID") %>%
  relocate(any_of(c("sampleID", "subjectID", "Timepoint", "date", "date_converted", "days_from_baseline", "altitude", "age", "gender")))

utils::write.csv(with_metadata, file.path(tables_dir, "alpha_diversity_with_metadata.csv"), row.names = FALSE)

summary_lines <- c(
  "# Alpha-Diversity Tool Comparison",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("- Samples in vegan output: ", nrow(vegan)),
  paste0("- Samples in scikit-bio output: ", nrow(skbio)),
  paste0("- Samples compared: ", nrow(comparison)),
  paste0("- Metrics compared: ", paste(metrics, collapse = ", ")),
  paste0("- Maximum Shannon absolute difference: ", max(comparison$Shannon_absolute_difference, na.rm = TRUE)),
  paste0("- Maximum Simpson absolute difference: ", max(comparison$Simpson_absolute_difference, na.rm = TRUE)),
  paste0("- Maximum Observed richness absolute difference: ", max(comparison$Observed_richness_absolute_difference, na.rm = TRUE))
)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(summary_lines, file.path(logs_dir, "alpha_diversity_tool_comparison_summary.md"))
cat(paste(summary_lines, collapse = "\n"), "\n")
