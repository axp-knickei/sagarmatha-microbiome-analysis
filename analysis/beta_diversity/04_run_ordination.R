#!/usr/bin/env Rscript
# Purpose: Run PCoA and configured NMDS ordinations for available
# beta-diversity distance matrices.
# Inputs: Prepared metadata and distance RDS files from stages 02 and 03.
# Outputs: Ordination coordinates, eigenvalues, diagnostics, RDS objects, and
# Shepard diagnostics under the configured results directory.
# Dependencies: Shared beta-diversity helpers and ordination packages declared
# in environment.yml.
# Usage: Rscript analysis/beta_diversity/04_run_ordination.R --config config/beta_diversity.yml
# Workflow: Stage 04; supplies coordinates to 07_make_figures.R.

source("analysis/beta_diversity/utils/io_helpers.R")
source("analysis/beta_diversity/utils/provenance_helpers.R")
args <- parse_beta_args(); config <- load_beta_config(args$config); ensure_beta_dirs(config)
require_beta_packages(c("vegan", "ggplot2")); source("analysis/beta_diversity/utils/plotting_theme.R")

with_stage_log(config, "04_run_ordination", {
  data <- readRDS(beta_output_dir(config, "prepared_data", "prepared_data.rds")); metadata <- data$metadata
  files <- list.files(beta_output_dir(config, "distances"), pattern = "^distance_.*\\.rds$", full.names = TRUE)
  files <- files[!grepl("_transformed\\.rds$", files)]
  for (file in files) {
    metric <- sub("\\.rds$", "", sub("^distance_", "", basename(file))); distance <- readRDS(file)
    n <- attr(distance, "Size"); k <- min(config$analysis$ordination_dimensions, n - 1)
    raw <- vegan::wcmdscale(distance, k = k, eig = TRUE, add = FALSE)
    negative <- raw$eig[raw$eig < 0]; positive <- raw$eig[raw$eig > 0]
    negative_ratio <- if (length(negative)) sum(abs(negative)) / sum(abs(raw$eig)) else 0
    correction <- if (negative_ratio > 0.01) "lingoes" else "none"
    fit <- if (correction == "lingoes") vegan::wcmdscale(distance, k = k, eig = TRUE, add = "lingoes") else raw
    coords <- as.data.frame(fit$points); names(coords) <- paste0("PCoA", seq_len(ncol(coords)))
    coords$sample_id <- rownames(coords); coords <- merge(metadata, coords, by = "sample_id", sort = FALSE)
    coords <- coords[match(attr(distance, "Labels"), coords$sample_id), ]
    eig <- data.frame(axis = seq_along(fit$eig), eigenvalue = fit$eig,
                      variance_percent = 100 * pmax(fit$eig, 0) / sum(pmax(fit$eig, 0)))
    write_tsv(coords, beta_output_dir(config, "ordination", paste0("pcoa_", metric, "_coordinates.tsv")))
    write_tsv(eig, beta_output_dir(config, "ordination", paste0("pcoa_", metric, "_eigenvalues.tsv")))
    write_tsv(data.frame(metric = metric, negative_eigenvalues = length(negative),
                         negative_eigenvalue_sum = sum(negative), negative_absolute_ratio = negative_ratio,
                         correction = correction), beta_output_dir(config, "ordination", paste0("pcoa_", metric, "_diagnostics.tsv")))
    saveRDS(fit, beta_output_dir(config, "ordination", paste0("pcoa_", metric, ".rds")))
    if (correction != "none") saveRDS(raw, beta_output_dir(config, "ordination", paste0("pcoa_", metric, "_uncorrected.rds")))
    if (metric %in% c("bray_curtis", "jaccard_binary") && n >= 4) {
      set.seed(config$analysis$random_seed)
      nmds <- vegan::metaMDS(distance, k = 2, trymax = config$analysis$nmds_trymax, autotransform = FALSE, trace = FALSE)
      scores <- data.frame(sample_id = rownames(nmds$points), NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2])
      scores <- merge(metadata, scores, by = "sample_id", sort = FALSE)
      write_tsv(scores, beta_output_dir(config, "ordination", paste0("nmds_", metric, "_coordinates.tsv")))
      write_tsv(data.frame(metric = metric, stress = nmds$stress, converged = nmds$converged,
                           tries = nmds$tries), beta_output_dir(config, "ordination", paste0("nmds_", metric, "_diagnostics.tsv")))
      saveRDS(nmds, beta_output_dir(config, "ordination", paste0("nmds_", metric, ".rds")))
      png <- beta_output_dir(config, "figures", paste0("nmds_", metric, "_shepard.png"))
      grDevices::png(png, width = 1800, height = 1600, res = 250); vegan::stressplot(nmds); grDevices::dev.off()
    }
    beta_log(config, "04_run_ordination", paste0("COMPLETED_PCOA_", toupper(metric)), paste("correction=", correction))
  }
})
