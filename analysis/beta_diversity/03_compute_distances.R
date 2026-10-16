#!/usr/bin/env Rscript
# Purpose: Compute the configured ecological distance matrices from the
# filtered sample-by-MAG community matrix.
# Inputs: prepared_data/prepared_data.rds, optional phylogenetic tree, and
# distance settings from the selected beta-diversity YAML file.
# Outputs: Distance RDS/TSV files, summaries, transformation objects, manifests,
# and related QC files under the configured results directory.
# Dependencies: Shared beta-diversity helpers plus vegan, ape, GUniFrac, and
# other distance packages declared in environment.yml.
# Usage: Rscript analysis/beta_diversity/03_compute_distances.R --config config/beta_diversity.yml
# Workflow: Stage 03; consumed by ordination, inference, figures, and reporting.

source("analysis/beta_diversity/utils/io_helpers.R")
source("analysis/beta_diversity/utils/provenance_helpers.R")
source("analysis/beta_diversity/utils/validation_helpers.R")
source("analysis/beta_diversity/utils/distance_helpers.R")
args <- parse_beta_args(); config <- load_beta_config(args$config); ensure_beta_dirs(config)
require_beta_packages(c("vegan", "yaml"))

compute_unifrac <- function(metric, x) {
  tree_path <- config$inputs$tree_file
  if (is.null(tree_path) || is.na(tree_path) || !nzchar(tree_path)) return(NULL)
  require_beta_packages(c("ape", "phyloseq"))
  tree <- ape::read.tree(tree_path)
  if (is.null(tree$edge.length) || any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0)) {
    stop("FAILED_INVALID_TREE_BRANCH_LENGTHS", call. = FALSE)
  }
  if (!ape::is.rooted(tree)) stop("FAILED_UNROOTED_PHYLOGENETIC_TREE", call. = FALSE)
  shared <- intersect(colnames(x), tree$tip.label)
  mismatch <- rbind(data.frame(identifier = setdiff(colnames(x), tree$tip.label), source = "community_only"),
                    data.frame(identifier = setdiff(tree$tip.label, colnames(x)), source = "tree_only"))
  write_tsv(mismatch, beta_output_dir(config, "qc", "tree_feature_matching_report.tsv"))
  if (length(shared) < 2) stop("FAILED_INSUFFICIENT_TREE_MAG_MATCH", call. = FALSE)
  x <- x[, shared, drop = FALSE]; x <- x[rowSums(x) > 0, , drop = FALSE]
  tree <- ape::keep.tip(tree, shared)
  ps <- phyloseq::phyloseq(phyloseq::otu_table(x, taxa_are_rows = FALSE), phyloseq::phy_tree(tree))
  if (metric == "unweighted_unifrac") return(phyloseq::UniFrac(ps, weighted = FALSE, normalized = TRUE, parallel = FALSE, fast = TRUE))
  if (metric == "weighted_unifrac") return(phyloseq::UniFrac(ps, weighted = TRUE, normalized = TRUE, parallel = FALSE, fast = TRUE))
  require_beta_packages("GUniFrac")
  alpha <- as.numeric(sub("generalized_unifrac_alpha_", "", metric))
  result <- GUniFrac::GUniFrac(x, tree, alpha = alpha)$unifracs
  as.dist(result[, , paste0("d_", alpha)])
}

save_distance <- function(metric, distance, preprocessing, transformed = NULL) {
  validate_distance(distance, attr(distance, "Labels"))
  stem <- file.path(beta_output_dir(config, "distances"), paste0("distance_", metric))
  assert_writable_output(paste0(stem, ".rds"), args$force)
  saveRDS(distance, paste0(stem, ".rds"))
  write_tsv(distance_long(distance), paste0(stem, "_long.tsv"))
  if (length(distance) <= 2000000) write_matrix_tsv(as.matrix(distance), paste0(stem, "_matrix.tsv"))
  write_tsv(distance_summary(distance), paste0(stem, "_summary.tsv"))
  write_yaml_manifest(config, paste0(stem, "_manifest.yml"), list(
    metric = metric, status = "COMPLETED", preprocessing = preprocessing,
    sample_order = attr(distance, "Labels"), checksum = sha256_file(paste0(stem, ".rds"))))
  if (!is.null(transformed)) saveRDS(transformed, paste0(stem, "_transformed.rds"))
  beta_log(config, "03_compute_distances", paste0("COMPLETED_", toupper(metric)), preprocessing)
}

with_stage_log(config, "03_compute_distances", {
  data <- readRDS(beta_output_dir(config, "prepared_data", "prepared_data.rds")); x <- data$community_matrix
  metrics <- unlist(config$analysis$metrics); if (!is.null(args$metric)) metrics <- args$metric
  relative <- total_sum_scale(x)
  for (metric in metrics) {
    if (metric == "bray_curtis") {
      transformed <- if (isTRUE(config$data$bray_square_root)) sqrt(relative) else relative
      save_distance(metric, vegan::vegdist(transformed, method = "bray"),
                    if (isTRUE(config$data$bray_square_root)) "total-sum scaling; square root" else "total-sum scaling", transformed)
      if (isTRUE(config$data$bray_square_root)) {
        save_distance("bray_curtis_untransformed_sensitivity", vegan::vegdist(relative, method = "bray"),
                      "total-sum scaling; no square-root transformation", relative)
      }
    } else if (metric == "jaccard_binary") {
      binary <- (x > config$data$detection_threshold) * 1
      if (any(!binary %in% c(0, 1))) stop("Jaccard branch is not binary.")
      save_distance(metric, vegan::vegdist(binary, method = "jaccard", binary = TRUE),
                    paste("presence threshold >", config$data$detection_threshold), binary)
      thresholds <- setdiff(unlist(config$data$detection_threshold_sensitivity), config$data$detection_threshold)
      for (threshold in thresholds) {
        sensitivity_binary <- (x > threshold) * 1
        sensitivity_metric <- paste0("jaccard_binary_detection_", format(threshold, scientific = TRUE))
        save_distance(sensitivity_metric, vegan::vegdist(sensitivity_binary, method = "jaccard", binary = TRUE),
                      paste("presence threshold >", threshold), sensitivity_binary)
      }
    } else if (metric == "quantitative_jaccard_ruzička") {
      save_distance(metric, quantitative_jaccard_ruzička(relative), "total-sum scaling; quantitative Jaccard/Ružička", relative)
    } else if (metric == "euclidean_hellinger") {
      transformed <- hellinger_transform(x)
      save_distance(metric, stats::dist(transformed, method = "euclidean"), "total-sum scaling; Hellinger transform", transformed)
    } else if (metric == "aitchison") {
      replaced <- multiplicative_replace(x, config$data$zero_replacement_fraction)
      clr <- clr_transform(replaced)
      if (max(abs(rowMeans(clr))) > 1e-10) stop("CLR row means are not approximately zero.")
      save_distance(metric, stats::dist(clr, method = "euclidean"),
                    paste("multiplicative zero replacement; fraction", config$data$zero_replacement_fraction, "; CLR"), clr)
      pseudocount <- config$data$pseudocount %||% (min(x[x > 0]) / 2)
      pseudo_clr <- clr_transform(total_sum_scale(x + pseudocount))
      save_distance("aitchison_pseudocount_sensitivity", stats::dist(pseudo_clr),
                    paste("pseudocount", format(pseudocount, scientific = TRUE), "; total-sum scaling; CLR"), pseudo_clr)
      write_tsv(data.frame(zero_count = sum(x == 0), zero_fraction = mean(x == 0),
                           replacement_method = "multiplicative_replacement", replacement_fraction = config$data$zero_replacement_fraction),
                beta_output_dir(config, "qc", "aitchison_zero_replacement.tsv"))
    } else if (grepl("unifrac", metric)) {
      tree <- config$inputs$tree_file
      if (is.null(tree) || is.na(tree) || !nzchar(tree)) {
        status <- "SKIPPED_UNIFRAC_NO_TREE"
        write_yaml_manifest(config, file.path(beta_output_dir(config, "distances"), paste0("distance_", metric, "_manifest.yml")),
                            list(metric = metric, status = status, reason = "A rooted MAG tree with non-negative branch lengths is required."))
        beta_log(config, "03_compute_distances", status, metric)
      } else {
        distance <- compute_unifrac(metric, x)
        save_distance(metric, distance, "matched rooted MAG phylogeny; unmatched tips/features audited and pruned")
      }
    } else stop("Unknown metric: ", metric, call. = FALSE)
  }
})
