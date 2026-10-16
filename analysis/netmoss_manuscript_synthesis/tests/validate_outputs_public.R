#!/usr/bin/env Rscript

# Validate exact reconciliation and required deliverables for the manuscript
# synthesis. This test reads generated tables; it does not rerun NetMoss.

suppressPackageStartupMessages(library(yaml))

args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/netmoss_manuscript_synthesis.yml"
if (length(args) == 2L && args[[1]] == "--config") config_path <- args[[2]]
if (!file.exists(config_path)) stop("Configuration not found: ", config_path, call. = FALSE)
cfg <- yaml::read_yaml(config_path)
root <- Sys.getenv("NETMOSS_MANUSCRIPT_OUTPUT_DIR", unset = cfg$outputs$output_dir)
main <- file.path(root, "main")
repro <- file.path(root, "reproducibility")

required <- c(
  file.path(root, "00_quality_control", "00_build_main_figure_and_tables.complete"),
  file.path(root, "RESULTS_INDEX.md"),
  file.path(main, "Main_Figure_01_exploratory_NetMoss_landscape.pdf"),
  file.path(main, "Main_Figure_01_exploratory_NetMoss_landscape.png"),
  file.path(main, "Main_Figure_01_exploratory_NetMoss_landscape.svg"),
  file.path(main, "Main_Figure_01_caption.txt"),
  file.path(main, "Main_Table_01_analysis_scope_and_evidence.csv"),
  file.path(main, "Main_Table_02_top_ranked_genera.csv"),
  file.path(repro, "Main_Figure_01_panels_b_c.csv"),
  file.path(repro, "Main_Figure_01_panel_d_edges.csv"),
  file.path(repro, "manuscript_synthesis_manifest.yml")
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing synthesis outputs: ", paste(missing, collapse = ", "), call. = FALSE)
if (any(file.info(required)$size <= 0)) stop("One or more required outputs are empty.", call. = FALSE)

taxa <- utils::read.csv(file.path(repro, "Main_Figure_01_panels_b_c.csv"), check.names = FALSE)
edges <- utils::read.csv(file.path(repro, "Main_Figure_01_panel_d_edges.csv"), check.names = FALSE)
table1 <- utils::read.csv(file.path(main, "Main_Table_01_analysis_scope_and_evidence.csv"), check.names = FALSE)
table2 <- utils::read.csv(file.path(main, "Main_Table_02_top_ranked_genera.csv"), check.names = FALSE)
manifest <- yaml::read_yaml(file.path(repro, "manuscript_synthesis_manifest.yml"))

input_dir <- cfg$inputs$netmoss_results_dir
key <- cfg$input_environment_overrides$netmoss_results_dir
if (!is.null(key) && nzchar(Sys.getenv(key))) input_dir <- Sys.getenv(key)
scores <- utils::read.csv(file.path(input_dir, "tables", "netmoss_scores.csv"), check.names = FALSE)
source_manifest <- yaml::read_yaml(file.path(input_dir, "reproducibility", "netmoss_manifest.yml"))
if (anyDuplicated(scores$taxon_names) || anyDuplicated(taxa$taxon_names) || !setequal(taxa$taxon_names, scores$taxon_names)) stop("Scored taxa differ from supplied scores.", call. = FALSE)
source_rows <- scores[match(taxa$taxon_names, scores$taxon_names), ]
for (pair in list(c("netmoss_score", "NetMoss_Score"), c("p_value", "p.val"), c("q_value", "p.adj"))) {
  if (!isTRUE(all.equal(taxa[[pair[1]]], source_rows[[pair[2]]], check.attributes = FALSE))) stop("Score values differ from supplied scores.", call. = FALSE)
}
if (any(!is.na(taxa$q_value) & (taxa$q_value < 0 | taxa$q_value > 1)) ||
    !identical(taxa$fdr_significant, is.finite(taxa$q_value) & taxa$q_value < cfg$figure$fdr_threshold) ||
    !identical(taxa$passes_score_threshold, taxa$netmoss_score >= cfg$figure$score_threshold)) stop("Invalid probability or threshold classification.", call. = FALSE)
if (nrow(table1) != 9L || nrow(table2) != min(nrow(taxa), as.integer(cfg$figure$top_taxa_table))) stop("Main table row counts do not match the configured contract.", call. = FALSE)
ranked <- taxa[order(-taxa$netmoss_score, taxa$q_value), ]
if (!identical(table2$Genus, sub("^g__", "", head(ranked$taxon_names, nrow(table2))))) stop("Top table ranking disagrees with supplied scores.", call. = FALSE)
expected_edges <- lapply(c("a", "b"), function(group) {
  network <- utils::read.csv(file.path(input_dir, "networks", paste0("group_", group, "_network.csv")), check.names = FALSE)
  ids <- as.character(network[[1]])
  mat <- as.matrix(network[-1]); storage.mode(mat) <- "double"; rownames(mat) <- ids
  if (!identical(rownames(mat), colnames(mat)) || anyNA(mat) || any(abs(mat) > 1) || max(abs(mat - t(mat))) > 1e-8) stop("Invalid source network.", call. = FALSE)
  mat <- mat[scores$taxon_names, scores$taxon_names, drop = FALSE]
  index <- which(upper.tri(mat), arr.ind = TRUE)
  label <- if (group == "a") paste0("Early expedition\nT1/T2/Ta; n = ", source_manifest$group_a$independent_profiles) else paste0("Late/recovery\nT8/T9/T10; n = ", source_manifest$group_b$independent_profiles)
  out <- data.frame(taxon_1 = rownames(mat)[index[, 1]], taxon_2 = colnames(mat)[index[, 2]], correlation = mat[index], condition = label)
  out[abs(out$correlation) >= cfg$figure$context_edge_threshold, ]
})
expected_edges <- do.call(rbind, expected_edges)
edge_key <- function(x) paste(x$condition, x$taxon_1, x$taxon_2, sep = "|")
if (anyDuplicated(edge_key(edges)) || !setequal(edge_key(edges), edge_key(expected_edges))) stop("Displayed edge identities disagree with supplied networks.", call. = FALSE)
expected_edges <- expected_edges[match(edge_key(edges), edge_key(expected_edges)), ]
if (!isTRUE(all.equal(as.numeric(edges$correlation), as.numeric(expected_edges$correlation), check.attributes = FALSE)) ||
    !identical(as.logical(edges$retained), abs(as.numeric(edges$correlation)) >= cfg$figure$retained_edge_threshold)) stop("Edge values or retention flags disagree with supplied networks.", call. = FALSE)
for (pair in list(c("scored_taxa", nrow(taxa)), c("score_threshold_taxa", sum(taxa$passes_score_threshold)), c("fdr_taxa", sum(taxa$fdr_significant)), c("both_threshold_taxa", sum(taxa$passes_score_threshold & taxa$fdr_significant)))) {
 if (as.numeric(manifest[[pair[1]]]) != as.numeric(pair[2])) stop("Synthesis manifest counts disagree with supplied scores.", call. = FALSE)
}

svg <- toupper(paste(readLines(file.path(main, "Main_Figure_01_exploratory_NetMoss_landscape.svg"), warn = FALSE), collapse = "\n"))
used_colors <- "#009ADE" # participant support is always blue
if (any(taxa$fdr_significant & taxa$passes_score_threshold)) used_colors <- c(used_colors, "#FF1F5B")
if (any(taxa$fdr_significant & !taxa$passes_score_threshold)) used_colors <- c(used_colors, "#F28522")
if (length(unique(taxa$phylum)) >= 2L) used_colors <- c(used_colors, "#FFC61E")
if (length(unique(taxa$phylum)) >= 3L) used_colors <- c(used_colors, "#F28522")
if (length(unique(taxa$phylum)) >= 5L) used_colors <- c(used_colors, "#AF58BA")
for (color in unique(used_colors)) {
  if (!grepl(color, svg, fixed = TRUE)) stop("Main figure SVG lacks required GeoDataViz color: ", color, call. = FALSE)
}

if (!identical(manifest$scientific_status, "historical_exploratory_supporting_only")) {
  stop("Scientific status boundary is missing or changed.", call. = FALSE)
}
caption <- paste(readLines(file.path(main, "Main_Figure_01_caption.txt"), warn = FALSE), collapse = " ")
for (phrase in c("faint dashed edges", "not retained associations", "are not direct ecological interactions", "exploratory")) {
  if (!grepl(phrase, caption, fixed = TRUE)) stop("Caption lacks safeguard: ", phrase, call. = FALSE)
}

cat("NetMoss manuscript synthesis outputs validated.\n")
