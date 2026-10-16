#!/usr/bin/env Rscript

# Create journal-ready tables, vector figures, high-resolution TIFF, captions,
# and reporting text from completed NetMoss2 outputs. This script visualizes and
# documents results; it does not upgrade exploratory evidence into validation.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(gridExtra)
  library(yaml)
})

`%||%` <- function(x, y) if (is.null(x) || !length(x) || all(is.na(x)) || identical(x, "")) y else x

parse_args <- function(args) {
  out <- list(config = "config/netmoss_analysis.yml", input_dir = NULL, output_dir = NULL,
              test_mode = FALSE, force = FALSE)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--test-mode", "--force")) {
      out[[if (arg == "--test-mode") "test_mode" else "force"]] <- TRUE
      i <- i + 1L
    } else if (arg %in% c("--config", "--input-dir", "--output-dir")) {
      if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
      key <- switch(arg, `--config` = "config", `--input-dir` = "input_dir", `--output-dir` = "output_dir")
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    } else if (arg == "--help") {
      cat("Usage: Rscript scripts/netmoss/create_publication_outputs.R [--input-dir PATH] [--output-dir PATH] [--test-mode] [--force]\n")
      quit(status = 0)
    } else stop("Unknown argument: ", arg, call. = FALSE)
  }
  out
}

read_matrix_csv <- function(path) {
  if (!file.exists(path)) stop("Matrix file not found: ", path, call. = FALSE)
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 2) stop("Matrix file has fewer than two columns: ", path, call. = FALSE)
  ids <- as.character(x[[1]])
  matrix <- as.matrix(x[-1])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- ids
  matrix
}

normalize_score_table <- function(scores) {
  clean <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  names_clean <- clean(names(scores))
  pick <- function(aliases, required = TRUE) {
    hit <- match(clean(aliases), names_clean, nomatch = 0L)
    hit <- hit[hit > 0]
    if (length(hit)) names(scores)[hit[[1]]] else if (required) stop("Missing score column: ", paste(aliases, collapse = "/"), call. = FALSE) else NA_character_
  }
  taxon_col <- pick(c("taxon_names", "taxon", "name"), required = FALSE)
  if (is.na(taxon_col)) scores$taxon_names <- rownames(scores) else scores$taxon_names <- as.character(scores[[taxon_col]])
  score_col <- pick(c("NetMoss_score", "NetMoss_Score", "score"))
  p_col <- pick(c("p.val", "p_value", "pvalue"), required = FALSE)
  q_col <- pick(c("p.adj", "q_value", "padj", "fdr"), required = FALSE)
  data.frame(
    taxon_names = scores$taxon_names,
    netmoss_score = suppressWarnings(as.numeric(scores[[score_col]])),
    p_value = if (is.na(p_col)) NA_real_ else suppressWarnings(as.numeric(scores[[p_col]])),
    q_value = if (is.na(q_col)) NA_real_ else suppressWarnings(as.numeric(scores[[q_col]])),
    stringsAsFactors = FALSE
  )
}

validate_inputs <- function(scores, abundance_a, abundance_b, network_a, network_b) {
  taxa <- scores$taxon_names
  if (anyDuplicated(taxa) || anyNA(taxa) || any(!nzchar(taxa))) stop("Score table has invalid taxon names.", call. = FALSE)
  if (any(!is.finite(scores$netmoss_score))) stop("NetMoss scores must be finite.", call. = FALSE)
  matrices <- list(abundance_a, abundance_b, network_a, network_b)
  labels <- c("group A abundance", "group B abundance", "group A network", "group B network")
  for (i in seq_along(matrices)) {
    if (length(setdiff(taxa, rownames(matrices[[i]])))) stop(labels[[i]], " lacks scored taxa.", call. = FALSE)
  }
  if (nrow(network_a) != ncol(network_a) || nrow(network_b) != ncol(network_b)) stop("Networks must be square.", call. = FALSE)
  if (max(abs(network_a - t(network_a)), na.rm = TRUE) > 1e-8 || max(abs(network_b - t(network_b)), na.rm = TRUE) > 1e-8) stop("Networks must be symmetric.", call. = FALSE)
}

validate_analysis_manifest <- function(manifest, config, manifest_path) {
  expected <- list(
    scenario = as.character(config$scenario$name),
    group_a = as.character(unlist(config$scenario$group_a_values)),
    group_b = as.character(unlist(config$scenario$group_b_values)),
    excluded = as.character(unlist(config$scenario$excluded_values %||% character(0)))
  )
  observed <- list(
    scenario = as.character(manifest$scenario %||% character(0)),
    group_a = as.character(unlist(manifest$group_a$timepoints %||% character(0))),
    group_b = as.character(unlist(manifest$group_b$timepoints %||% character(0))),
    excluded = as.character(unlist(manifest$excluded_timepoints %||% character(0)))
  )
  mismatched <- names(expected)[!vapply(
    names(expected),
    function(field) identical(observed[[field]], expected[[field]]),
    logical(1)
  )]
  if (length(mismatched)) {
    stop(
      "Existing NetMoss outputs do not match the active configuration (mismatch: ",
      paste(mismatched, collapse = ", "), "). Rerun the NetMoss analysis before ",
      "creating publication outputs. Historical results will not be relabeled. Manifest: ",
      manifest_path,
      call. = FALSE
    )
  }
}

synthetic_inputs <- function() {
  set.seed(20260714)
  taxa <- paste0("Genus_", sprintf("%02d", 1:24))
  abundance_a <- matrix(rexp(24 * 12, 30), nrow = 24, dimnames = list(taxa, paste0("A", 1:12)))
  abundance_b <- matrix(rexp(24 * 12, 28), nrow = 24, dimnames = list(taxa, paste0("B", 1:12)))
  abundance_a <- sweep(abundance_a, 2, colSums(abundance_a), "/")
  abundance_b <- sweep(abundance_b, 2, colSums(abundance_b), "/")
  network_a <- stats::cor(t(abundance_a), method = "spearman")
  network_b <- stats::cor(t(abundance_b), method = "spearman")
  diag(network_a) <- diag(network_b) <- 1
  scores <- data.frame(taxon_names = taxa, netmoss_score = sort(runif(24), decreasing = TRUE),
                       p_value = sort(runif(24, 0.0001, 0.2)), stringsAsFactors = FALSE)
  scores$q_value <- p.adjust(scores$p_value, method = "BH")
  list(scores = scores, abundance_a = abundance_a, abundance_b = abundance_b,
       network_a = network_a, network_b = network_b)
}

network_edges <- function(network_a, network_b, threshold) {
  taxa <- rownames(network_a)
  index <- which(upper.tri(network_a), arr.ind = TRUE)
  data.frame(
    taxon_1 = taxa[index[, 1]], taxon_2 = taxa[index[, 2]],
    group_a_correlation = network_a[index], group_b_correlation = network_b[index],
    delta_correlation = network_b[index] - network_a[index],
    retained_at_threshold = abs(network_a[index]) >= threshold | abs(network_b[index]) >= threshold,
    stringsAsFactors = FALSE
  )
}

theme_nature <- function() {
  theme_classic(base_size = 8, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0),
      axis.title = element_text(size = 8), axis.text = element_text(size = 7, color = "black"),
      legend.title = element_text(size = 7), legend.text = element_text(size = 7),
      strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8),
      plot.margin = margin(4, 4, 4, 4)
    )
}

save_publication_figure <- function(plot, stem, directory, width_mm, height_mm, dpi) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot, width = width_in, height = height_in, device = grDevices::cairo_pdf)
  ggsave(file.path(directory, paste0(stem, ".svg")), plot, width = width_in, height = height_in, device = svglite::svglite)
  ggsave(file.path(directory, paste0(stem, ".tiff")), plot, width = width_in, height = height_in,
         device = "tiff", dpi = dpi, compression = "lzw", bg = "white")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!file.exists(args$config)) stop("Configuration not found: ", args$config, call. = FALSE)
config <- yaml::read_yaml(args$config)
input_dir <- args$input_dir %||% config$outputs$output_dir
output_dir <- args$output_dir %||% file.path(input_dir, "publication")
figures_dir <- file.path(output_dir, "figures")
tables_dir <- file.path(output_dir, "tables")
text_dir <- file.path(output_dir, "text")
invisible(vapply(c(figures_dir, tables_dir, text_dir), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))

sentinel <- file.path(output_dir, "publication_output_manifest.yml")
if (file.exists(sentinel) && !args$force) stop("Refusing to overwrite publication outputs. Use --force.", call. = FALSE)

if (args$test_mode) {
  inputs <- synthetic_inputs()
} else {
  analysis_manifest_path <- file.path(input_dir, "reproducibility", "netmoss_manifest.yml")
  if (!file.exists(analysis_manifest_path)) {
    stop("NetMoss analysis manifest not found: ", analysis_manifest_path, call. = FALSE)
  }
  analysis_manifest <- yaml::read_yaml(analysis_manifest_path)
  validate_analysis_manifest(analysis_manifest, config, analysis_manifest_path)
  score_path <- file.path(input_dir, "tables", "netmoss_scores.csv")
  if (!file.exists(score_path)) stop("NetMoss score table not found: ", score_path, call. = FALSE)
  raw_scores <- utils::read.csv(score_path, check.names = FALSE, stringsAsFactors = FALSE)
  inputs <- list(
    scores = normalize_score_table(raw_scores),
    abundance_a = read_matrix_csv(file.path(input_dir, "prepared", "group_a_abundance.csv")),
    abundance_b = read_matrix_csv(file.path(input_dir, "prepared", "group_b_abundance.csv")),
    network_a = read_matrix_csv(file.path(input_dir, "networks", "group_a_network.csv")),
    network_b = read_matrix_csv(file.path(input_dir, "networks", "group_b_network.csv"))
  )
}

scores <- normalize_score_table(inputs$scores)
validate_inputs(scores, inputs$abundance_a, inputs$abundance_b, inputs$network_a, inputs$network_b)
taxa <- scores$taxon_names
abundance_a <- inputs$abundance_a[taxa, , drop = FALSE]
abundance_b <- inputs$abundance_b[taxa, , drop = FALSE]
network_a <- inputs$network_a[taxa, taxa, drop = FALSE]
network_b <- inputs$network_b[taxa, taxa, drop = FALSE]

fdr_threshold <- as.numeric(config$publication$fdr_threshold)
score_threshold <- as.numeric(config$publication$score_threshold)
edge_threshold <- as.numeric(config$publication$edge_threshold)
pseudocount <- min(c(abundance_a[abundance_a > 0], abundance_b[abundance_b > 0]), na.rm = TRUE) / 2

degree_a <- rowSums(abs(network_a) >= edge_threshold) - 1
degree_b <- rowSums(abs(network_b) >= edge_threshold) - 1
publication_table <- scores %>%
  mutate(
    mean_abundance_group_a = rowMeans(abundance_a[taxon_names, , drop = FALSE]),
    mean_abundance_group_b = rowMeans(abundance_b[taxon_names, , drop = FALSE]),
    prevalence_group_a = rowMeans(abundance_a[taxon_names, , drop = FALSE] > 0),
    prevalence_group_b = rowMeans(abundance_b[taxon_names, , drop = FALSE] > 0),
    log2_fold_change_b_vs_a = log2((mean_abundance_group_b + pseudocount) / (mean_abundance_group_a + pseudocount)),
    degree_group_a = degree_a[taxon_names], degree_group_b = degree_b[taxon_names],
    degree_change_b_minus_a = degree_group_b - degree_group_a,
    plot_mean_group_a = mean_abundance_group_a + pseudocount,
    plot_mean_group_b = mean_abundance_group_b + pseudocount,
    fdr_significant = is.finite(q_value) & q_value < fdr_threshold,
    passes_score_threshold = netmoss_score >= score_threshold,
    marker_status = case_when(
      fdr_significant & passes_score_threshold ~ "FDR and score threshold",
      fdr_significant ~ "FDR only", passes_score_threshold ~ "Score threshold only",
      TRUE ~ "Neither threshold"
    )
  ) %>% arrange(desc(netmoss_score), q_value)

edges <- network_edges(network_a, network_b, edge_threshold)
utils::write.csv(publication_table, file.path(tables_dir, "Supplementary_Data_1_NetMoss_taxon_results.csv"), row.names = FALSE)
utils::write.csv(edges, file.path(tables_dir, "Supplementary_Data_2_NetMoss_network_edges.csv"), row.names = FALSE)
utils::write.csv(publication_table %>% slice_head(n = min(as.integer(config$publication$top_taxa), nrow(publication_table))),
                 file.path(tables_dir, "Source_Data_Figure_1.csv"), row.names = FALSE)

top_n <- min(as.integer(config$publication$top_taxa), nrow(publication_table))
top <- publication_table %>% slice_head(n = top_n) %>% mutate(taxon_names = factor(taxon_names, levels = rev(taxon_names)))
status_colors <- c("FDR and score threshold" = "#B2182B", "FDR only" = "#EF8A62",
                   "Score threshold only" = "#2166AC", "Neither threshold" = "#BDBDBD")

panel_a <- ggplot(top, aes(netmoss_score, taxon_names, fill = marker_status)) +
  geom_col(width = 0.72) + geom_vline(xintercept = score_threshold, linetype = 2, linewidth = 0.35) +
  scale_fill_manual(values = status_colors, drop = FALSE) +
  labs(title = "a  Ranked NetMoss scores", x = "NetMoss score", y = NULL, fill = NULL) + theme_nature()

panel_b <- ggplot(publication_table, aes(plot_mean_group_a, plot_mean_group_b,
                                         color = netmoss_score, shape = fdr_significant)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55", linewidth = 0.35) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_x_log10(labels = label_scientific()) + scale_y_log10(labels = label_scientific()) +
  scale_color_viridis_c(option = "C") + scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16)) +
  labs(title = "b  Mean relative abundance", x = config$scenario$group_a_label,
       y = config$scenario$group_b_label, color = "NetMoss", shape = "FDR < 0.05") + theme_nature()

panel_c <- ggplot(publication_table, aes(netmoss_score, degree_change_b_minus_a,
                                         color = marker_status)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey55", linewidth = 0.35) +
  geom_point(size = 1.8, alpha = 0.85) + scale_color_manual(values = status_colors, drop = FALSE) +
  labs(title = "c  Network degree change", x = "NetMoss score", y = "Degree change (group B − group A)", color = NULL) + theme_nature()

highlight <- publication_table %>% slice_head(n = min(as.integer(config$publication$highlight_taxa), nrow(publication_table)))
node_positions <- data.frame(
  taxon_names = highlight$taxon_names,
  angle = seq(0, 2 * pi, length.out = nrow(highlight) + 1)[-(nrow(highlight) + 1)]
) %>% mutate(x = cos(angle), y = sin(angle))
network_panel_data <- bind_rows(
  edges %>% transmute(taxon_1, taxon_2, correlation = group_a_correlation, condition = config$scenario$group_a_label),
  edges %>% transmute(taxon_1, taxon_2, correlation = group_b_correlation, condition = config$scenario$group_b_label)
) %>%
  filter(taxon_1 %in% node_positions$taxon_names, taxon_2 %in% node_positions$taxon_names,
         abs(correlation) >= edge_threshold) %>%
  left_join(node_positions %>% select(taxon_1 = taxon_names, x1 = x, y1 = y), by = "taxon_1") %>%
  left_join(node_positions %>% select(taxon_2 = taxon_names, x2 = x, y2 = y), by = "taxon_2")
node_panel_data <- tidyr::crossing(node_positions, condition = c(config$scenario$group_a_label, config$scenario$group_b_label)) %>%
  left_join(publication_table %>% select(taxon_names, netmoss_score), by = "taxon_names")

panel_d <- ggplot() +
  geom_segment(data = network_panel_data, aes(x = x1, y = y1, xend = x2, yend = y2,
                                               color = correlation, linewidth = abs(correlation)), alpha = 0.65) +
  geom_point(data = node_panel_data, aes(x, y, size = netmoss_score), shape = 21, fill = "white", color = "black") +
  geom_text(data = node_panel_data, aes(x * 1.18, y * 1.18, label = taxon_names), size = 2, check_overlap = TRUE) +
  facet_wrap(~condition) + coord_equal(clip = "off") +
  scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0, limits = c(-1, 1)) +
  scale_linewidth(range = c(0.2, 1), guide = "none") + scale_size(range = c(2, 5)) +
  labs(title = "d  Rewiring among top-ranked taxa", color = "Correlation", size = "NetMoss") +
  theme_void(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), strip.text = element_text(face = "bold"), legend.position = "bottom")

figure_top <- gridExtra::arrangeGrob(
  panel_a,
  gridExtra::arrangeGrob(panel_b, panel_c, ncol = 1),
  ncol = 2,
  widths = c(1.15, 1)
)
figure <- gridExtra::arrangeGrob(figure_top, panel_d, ncol = 1, heights = c(1.25, 1))
save_publication_figure(
  figure, "Figure_1_NetMoss_network_rewiring", figures_dir,
  as.numeric(config$publication$figure_width_mm), as.numeric(config$publication$figure_height_mm),
  as.integer(config$publication$raster_dpi)
)

caption <- paste0(
  "Figure 1 | Exploratory NetMoss analysis of longitudinal gut microbiome profiles. ",
  "a, Ranked NetMoss scores for the top ", top_n, " taxa; the dashed line denotes the prespecified score threshold of ", score_threshold, ". ",
  "Colours distinguish taxa meeting the NetMoss score and/or Benjamini-Hochberg FDR criteria. ",
  "b, Mean relative abundance in ", config$scenario$group_a_label, " and ", config$scenario$group_b_label, ". ",
  "c, Association-network degree change at |r| ≥ ", edge_threshold, ". ",
  "d, Paired association networks among the highest-ranked taxa. Edges denote associations rather than direct ecological interactions. ",
  "Profiles were aggregated within participant before network construction. Owing to the small number of independent participants, all results are exploratory."
)
writeLines(caption, file.path(text_dir, "Figure_1_legend.txt"))

methods <- c(
  "NetMoss exploratory analysis and visualization",
  "===============================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
  paste0("Microbial features were aggregated to the ", config$columns$taxonomy_rank, " level and filtered at a prevalence threshold of ", config$preprocessing$prevalence_threshold, "."),
  paste0("Repeated samples within each participant and analysis phase were summarized using the ", config$preprocessing$subject_aggregation, "."),
  paste0("Separate ", config$network$method, " association networks were constructed for ", config$scenario$group_a_label, " and ", config$scenario$group_b_label, "."),
  paste0("NetMoss2 was pinned to Git commit ", config$netmoss$package_commit, "; scores were calculated with deepSplit=", config$netmoss$deep_split,
         ", minimum module size=", config$netmoss$minimum_module_size, ", and soft power=", config$netmoss$soft_power, "."),
  paste0("P values supplied by NetMoss2 were interpreted using Benjamini-Hochberg adjusted q < ", fdr_threshold,
         " together with a prespecified score threshold of ", score_threshold, "."),
  "All input matrices, network edges, source data, configuration, package versions, and figure-generation code were retained.",
  "Because the study has few independent participants, analyses were designated exploratory/supporting and were not interpreted as validated biomarkers."
)
writeLines(methods, file.path(text_dir, "NetMoss_methods_draft.txt"))

checklist <- data.frame(
  reporting_item = c("Analysis status", "Biological contrast", "Independent unit", "Taxonomic rank", "Prevalence threshold",
                     "Network estimator", "Network stability", "Multiple testing", "Software pinning", "Source data", "Vector figure", "High-resolution raster"),
  value = c("Exploratory/supporting only", config$scenario$name, "Participant after within-phase aggregation", config$columns$taxonomy_rank,
            config$preprocessing$prevalence_threshold, config$network$method, "Required before biological claims", paste0("BH FDR < ", fdr_threshold),
            config$netmoss$package_commit, "Source_Data_Figure_1.csv", "PDF and SVG", paste0(config$publication$raster_dpi, " dpi TIFF")),
  stringsAsFactors = FALSE
)
utils::write.csv(checklist, file.path(tables_dir, "NetMoss_reporting_checklist.csv"), row.names = FALSE)

manifest <- list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  test_mode = args$test_mode,
  input_dir = input_dir,
  taxa = nrow(publication_table),
  fdr_significant_taxa = sum(publication_table$fdr_significant),
  score_threshold_taxa = sum(publication_table$passes_score_threshold),
  figure_formats = c("PDF", "SVG", "TIFF"),
  raster_dpi = config$publication$raster_dpi,
  caution = "Publication-ready formatting does not change the exploratory evidential status."
)
yaml::write_yaml(manifest, sentinel)
old_tz <- Sys.getenv("TZ", unset = "")
Sys.setenv(TZ = "UTC")
suppressWarnings(utils::capture.output(utils::sessionInfo(), file = file.path(output_dir, "sessionInfo.txt")))
if (nzchar(old_tz)) Sys.setenv(TZ = old_tz) else Sys.unsetenv("TZ")
cat("Publication outputs completed: ", normalizePath(output_dir, mustWork = TRUE), "\n", sep = "")
