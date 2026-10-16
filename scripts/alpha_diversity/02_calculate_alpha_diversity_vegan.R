# Calculate alpha-diversity metrics with vegan.

required_packages <- c("vegan")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

output_dir <- Sys.getenv("BLACKICE_ALPHA_OUTPUT_DIR", unset = "results/alpha_diversity")
tables_dir <- file.path(output_dir, "tables")
logs_dir <- file.path(output_dir, "logs")
prepared_path <- file.path(tables_dir, "alpha_diversity_prepared_inputs.rds")
if (!file.exists(prepared_path)) stop("Prepared inputs not found. Run 01_prepare_alpha_diversity_inputs.R first.", call. = FALSE)

prepared <- readRDS(prepared_path)
x <- prepared$abundance_matrix
scale_decision <- prepared$diagnostics$abundance_scale_decision
count_like <- identical(scale_decision, "count_like")

shannon <- vegan::diversity(x, index = "shannon")
simpson <- vegan::diversity(x, index = "simpson")
observed <- vegan::specnumber(x)

if (count_like) {
  estimates <- vegan::estimateR(t(round(x)))
  chao1 <- as.numeric(estimates["S.chao1", ])
  ace <- as.numeric(estimates["S.ACE", ])
} else {
  chao1 <- rep(NA_real_, nrow(x))
  ace <- rep(NA_real_, nrow(x))
  warning("Chao1 and ACE were not calculated because the abundance table is not count-like.", call. = FALSE)
}

result <- data.frame(
  sampleID = rownames(x),
  Shannon = as.numeric(shannon),
  Simpson = as.numeric(simpson),
  Observed_richness = as.numeric(observed),
  Chao1 = chao1,
  ACE = ace,
  tool = "vegan",
  abundance_scale_decision = scale_decision,
  chao1_ace_decision = if (count_like) "calculated" else "not_calculated_relative_or_non_count_input",
  stringsAsFactors = FALSE
)

utils::write.csv(result, file.path(tables_dir, "alpha_diversity_vegan.csv"), row.names = FALSE)

session_file <- file.path(logs_dir, "alpha_diversity_vegan_session_info.txt")
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
capture.output(utils::sessionInfo(), file = session_file)
cat("Wrote vegan alpha-diversity table: ", file.path(tables_dir, "alpha_diversity_vegan.csv"), "\n", sep = "")
