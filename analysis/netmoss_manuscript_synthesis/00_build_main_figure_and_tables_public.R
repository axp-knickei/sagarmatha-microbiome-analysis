#!/usr/bin/env Rscript

# Build a four-panel manuscript figure and two main tables from a completed,
# preserved NetMoss analysis. This script is post-processing only: it never
# recalculates condition networks, NetMoss scores, or inferential statistics.
#
# Usage:
#   Rscript analysis/netmoss_manuscript_synthesis/00_build_main_figure_and_tables_public.R \
#     --config config/netmoss_manuscript_synthesis.yml [--force]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(igraph)
  library(scales)
  library(gridExtra)
  library(grid)
  library(svglite)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x)) || identical(x, "")) y else x
}

parse_args <- function(args) {
  out <- list(config = "config/netmoss_manuscript_synthesis.yml", force = FALSE)
  i <- 1L
  while (i <= length(args)) {
    if (args[[i]] == "--config" && i < length(args)) {
      out$config <- args[[i + 1L]]
      i <- i + 2L
    } else if (args[[i]] == "--force") {
      out$force <- TRUE
      i <- i + 1L
    } else {
      stop("Unknown argument: ", args[[i]], call. = FALSE)
    }
  }
  out
}

resolve_inputs <- function(cfg) {
  paths <- cfg$inputs
  overrides <- cfg$input_environment_overrides %||% list()
  for (key in names(paths)) {
    env_name <- overrides[[key]] %||% ""
    env_value <- if (nzchar(env_name)) Sys.getenv(env_name, unset = "") else ""
    if (nzchar(env_value)) paths[[key]] <- env_value
  }
  paths
}

read_matrix_csv <- function(path) {
  if (!file.exists(path)) stop("Required matrix not found: ", path, call. = FALSE)
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- as.character(x[[1]])
  matrix <- as.matrix(x[-1])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- ids
  if (nrow(matrix) != ncol(matrix) || !identical(rownames(matrix), colnames(matrix))) {
    stop("Network row and column identifiers do not form a matching square matrix: ", path, call. = FALSE)
  }
  if (max(abs(matrix - t(matrix)), na.rm = TRUE) > 1e-8) {
    stop("Network is not symmetric: ", path, call. = FALSE)
  }
  matrix
}

read_feature_matrix_csv <- function(path) {
  if (!file.exists(path)) stop("Required feature matrix not found: ", path, call. = FALSE)
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- as.character(x[[1]])
  matrix <- as.matrix(x[-1])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- ids
  matrix
}

normalize_scores <- function(path) {
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("taxon_names", "NetMoss_Score", "p.val", "p.adj")
  absent <- setdiff(required, names(x))
  if (length(absent)) stop("Score table lacks: ", paste(absent, collapse = ", "), call. = FALSE)
  x %>%
    transmute(
      taxon_names = as.character(taxon_names),
      netmoss_score = as.numeric(NetMoss_Score),
      p_value = as.numeric(p.val),
      q_value = as.numeric(p.adj)
    )
}

build_taxonomy <- function(abundance_path, taxonomy_path, taxa) {
  abundance <- utils::read.csv(abundance_path, check.names = FALSE, stringsAsFactors = FALSE)
  taxonomy <- utils::read.csv(taxonomy_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"sampleID" %in% names(abundance)) stop("gut.csv lacks sampleID.", call. = FALSE)
  required <- c("Genome", "genus", "phylum")
  absent <- setdiff(required, names(taxonomy))
  if (length(absent)) stop("Taxonomy lacks: ", paste(absent, collapse = ", "), call. = FALSE)
  mag_ids <- setdiff(names(abundance), "sampleID")
  selected <- data.frame(Genome = mag_ids, stringsAsFactors = FALSE) %>%
    left_join(taxonomy %>% select(Genome, genus, phylum), by = "Genome")
  if (anyNA(selected$genus) || anyNA(selected$phylum)) {
    stop("Every gut MAG must join exactly to taxonomy through Genome.", call. = FALSE)
  }
  genus_phylum <- selected %>%
    count(genus, phylum, name = "mag_count") %>%
    arrange(genus, desc(mag_count), phylum) %>%
    group_by(genus) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(
      taxon_names = genus,
      phylum = if_else(phylum == "p__", "Unclassified", sub("^p__", "", phylum))
    )
  data.frame(taxon_names = taxa, stringsAsFactors = FALSE) %>%
    left_join(genus_phylum, by = "taxon_names") %>%
    mutate(phylum = replace_na(phylum, "Unclassified"))
}

network_edges <- function(network, condition, context_threshold, retained_threshold) {
  index <- which(upper.tri(network), arr.ind = TRUE)
  data.frame(
    taxon_1 = rownames(network)[index[, 1]],
    taxon_2 = colnames(network)[index[, 2]],
    correlation = network[index],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      condition = condition,
      abs_correlation = abs(correlation),
      retained = abs_correlation >= retained_threshold,
      evidence_class = if_else(retained, "Retained association", "Below-threshold context"),
      correlation_sign = if_else(correlation >= 0, "Positive", "Negative")
    ) %>%
    filter(abs_correlation >= context_threshold)
}

curve_edges <- function(edges, positions, points = 24L) {
  if (!nrow(edges)) return(data.frame(edge_id = character(), condition = character(), correlation = numeric(), abs_correlation = numeric(), retained = logical(), evidence_class = character(), correlation_sign = character(), x = numeric(), y = numeric()))
  joined <- edges %>%
    mutate(edge_id = paste(condition, taxon_1, taxon_2, sep = "|")) %>%
    left_join(positions %>% select(taxon_1 = taxon_names, x1 = x, y1 = y), by = "taxon_1") %>%
    left_join(positions %>% select(taxon_2 = taxon_names, x2 = x, y2 = y), by = "taxon_2")
  bind_rows(lapply(seq_len(nrow(joined)), function(i) {
    edge <- joined[i, ]
    dx <- edge$x2 - edge$x1
    dy <- edge$y2 - edge$y1
    distance <- sqrt(dx^2 + dy^2)
    direction <- if ((match(edge$taxon_1, positions$taxon_names) +
                       match(edge$taxon_2, positions$taxon_names)) %% 2L) 1 else -1
    curvature <- direction * 0.13 * distance
    control_x <- (edge$x1 + edge$x2) / 2 - curvature * dy / max(distance, .Machine$double.eps)
    control_y <- (edge$y1 + edge$y2) / 2 + curvature * dx / max(distance, .Machine$double.eps)
    t <- seq(0, 1, length.out = points)
    data.frame(
      edge_id = edge$edge_id,
      condition = edge$condition,
      correlation = edge$correlation,
      abs_correlation = edge$abs_correlation,
      retained = edge$retained,
      evidence_class = edge$evidence_class,
      correlation_sign = edge$correlation_sign,
      x = (1 - t)^2 * edge$x1 + 2 * (1 - t) * t * control_x + t^2 * edge$x2,
      y = (1 - t)^2 * edge$y1 + 2 * (1 - t) * t * control_y + t^2 * edge$y2,
      stringsAsFactors = FALSE
    )
  }))
}

theme_nature <- function(base_size = 7.3) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(color = "#202124"),
      axis.text = element_text(color = "#35373A"),
      axis.title = element_text(color = "#202124"),
      axis.line = element_line(linewidth = 0.32, color = "#202124"),
      axis.ticks = element_line(linewidth = 0.28, color = "#202124"),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(2.8, "mm"),
      legend.key.width = unit(4.8, "mm"),
      plot.title = element_text(face = "bold", size = base_size + 0.7, margin = margin(b = 2, l = 9)),
      plot.subtitle = element_text(size = base_size - 0.3, color = "#55585C", margin = margin(b = 3, l = 9)),
      plot.tag = element_text(face = "bold", size = base_size + 3),
      plot.tag.position = c(0, 1),
      plot.margin = margin(5, 5, 5, 5)
    )
}

wrap_table_text <- function(x, width) {
  vapply(as.character(x), function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

save_table_preview <- function(x, stem, output_dir, widths, height, dpi) {
  preview <- data.frame(lapply(seq_along(x), function(i) wrap_table_text(x[[i]], widths[[i]])), check.names = FALSE)
  names(preview) <- names(x)
  table_theme <- ttheme_minimal(
    base_size = 7.2,
    base_family = "sans",
    core = list(
      fg_params = list(hjust = 0, x = 0.025, col = "#202124"),
      bg_params = list(fill = rep(c("#FFFFFF", "#F4F6F8"), length.out = nrow(preview)), col = "#D9DDE2", lwd = 0.35)
    ),
    colhead = list(
      fg_params = list(col = "white", fontface = "bold", hjust = 0, x = 0.025),
      bg_params = list(fill = "#443F90", col = "#443F90")
    )
  )
  grob <- tableGrob(preview, rows = NULL, theme = table_theme)
  grob$widths <- unit(widths / sum(widths), "npc")
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), grob, width = 183 / 25.4, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(output_dir, paste0(stem, ".png")), grob, width = 183 / 25.4, height = height, dpi = dpi, bg = "white")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!file.exists(args$config)) stop("Configuration not found: ", args$config, call. = FALSE)
cfg <- yaml::read_yaml(args$config)
inputs <- resolve_inputs(cfg)
input_dir <- inputs$netmoss_results_dir
output_dir <- Sys.getenv("NETMOSS_MANUSCRIPT_OUTPUT_DIR", unset = cfg$outputs$output_dir)
paths <- list(
  root = output_dir,
  main = file.path(output_dir, "main"),
  supplementary = file.path(output_dir, "supplementary"),
  reproducibility = file.path(output_dir, "reproducibility"),
  qc = file.path(output_dir, "00_quality_control")
)
invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
marker <- file.path(paths$qc, "00_build_main_figure_and_tables.complete")
overwrite <- isTRUE(args$force) || isTRUE(cfg$outputs$overwrite)
if (file.exists(marker) && !overwrite) stop("Completion marker exists; use --force intentionally.", call. = FALSE)
if (file.exists(marker) && overwrite) unlink(marker)

required <- c(
  file.path(input_dir, "reproducibility", "netmoss_manifest.yml"),
  file.path(input_dir, "tables", "netmoss_scores.csv"),
  file.path(input_dir, "prepared", "group_a_abundance.csv"),
  file.path(input_dir, "prepared", "group_b_abundance.csv"),
  file.path(input_dir, "networks", "group_a_network.csv"),
  file.path(input_dir, "networks", "group_b_network.csv"),
  inputs$gut_abundance_file,
  inputs$taxonomy_file
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Required inputs missing:\n", paste("-", missing, collapse = "\n"), call. = FALSE)

manifest <- yaml::read_yaml(required[[1]])
if (!identical(as.character(manifest$status), "COMPLETED_EXPLORATORY")) {
  stop("The source manifest must record COMPLETED_EXPLORATORY.", call. = FALSE)
}
historical_a <- c("T1", "T2", "Ta")
historical_b <- c("T8", "T9", "T10")
if (!identical(as.character(unlist(manifest$group_a$timepoints)), historical_a) ||
    !identical(as.character(unlist(manifest$group_b$timepoints)), historical_b)) {
  stop("This manuscript draft is explicitly bound to the preserved historical T1/T2/Ta versus T8/T9/T10 result.", call. = FALSE)
}

scores <- normalize_scores(file.path(input_dir, "tables", "netmoss_scores.csv"))
abundance_a <- read_feature_matrix_csv(file.path(input_dir, "prepared", "group_a_abundance.csv"))
abundance_b <- read_feature_matrix_csv(file.path(input_dir, "prepared", "group_b_abundance.csv"))
network_a <- read_matrix_csv(file.path(input_dir, "networks", "group_a_network.csv"))
network_b <- read_matrix_csv(file.path(input_dir, "networks", "group_b_network.csv"))
taxa <- scores$taxon_names
for (x in list(abundance_a, abundance_b, network_a, network_b)) {
  if (length(setdiff(taxa, rownames(x)))) stop("A source matrix lacks scored taxa.", call. = FALSE)
}
abundance_a <- abundance_a[taxa, , drop = FALSE]
abundance_b <- abundance_b[taxa, , drop = FALSE]
network_a <- network_a[taxa, taxa, drop = FALSE]
network_b <- network_b[taxa, taxa, drop = FALSE]

score_threshold <- as.numeric(cfg$figure$score_threshold)
fdr_threshold <- as.numeric(cfg$figure$fdr_threshold)
retained_threshold <- as.numeric(cfg$figure$retained_edge_threshold)
context_threshold <- as.numeric(cfg$figure$context_edge_threshold)
pseudocount <- min(c(abundance_a[abundance_a > 0], abundance_b[abundance_b > 0]), na.rm = TRUE) / 2
degree_a <- rowSums(abs(network_a) >= retained_threshold) - 1
degree_b <- rowSums(abs(network_b) >= retained_threshold) - 1

results <- scores %>%
  mutate(
    mean_abundance_early = rowMeans(abundance_a[taxon_names, , drop = FALSE]),
    mean_abundance_late = rowMeans(abundance_b[taxon_names, , drop = FALSE]),
    log2_fold_change_late_vs_early = log2((mean_abundance_late + pseudocount) / (mean_abundance_early + pseudocount)),
    degree_early = unname(degree_a[taxon_names]),
    degree_late = unname(degree_b[taxon_names]),
    fdr_significant = is.finite(q_value) & q_value < fdr_threshold,
    passes_score_threshold = netmoss_score >= score_threshold,
    evidence_group = case_when(
      fdr_significant & passes_score_threshold ~ "Score and BH q thresholds",
      fdr_significant ~ "BH q threshold only",
      passes_score_threshold ~ "Score threshold only",
      TRUE ~ "Neither threshold"
    )
  ) %>%
  arrange(desc(netmoss_score), q_value) %>%
  mutate(rank = row_number(), genus = sub("^g__", "", taxon_names))

taxonomy <- build_taxonomy(inputs$gut_abundance_file, inputs$taxonomy_file, taxa)
results <- results %>% left_join(taxonomy, by = "taxon_names")

condition_labels <- c(
  early = paste0("Early expedition\nT1/T2/Ta; n = ", manifest$group_a$independent_profiles),
  late = paste0("Late/recovery\nT8/T9/T10; n = ", manifest$group_b$independent_profiles)
)
edges <- bind_rows(
  network_edges(network_a, condition_labels[["early"]], context_threshold, retained_threshold),
  network_edges(network_b, condition_labels[["late"]], context_threshold, retained_threshold)
) %>%
  mutate(condition = factor(condition, levels = unname(condition_labels)))

# Exact GeoDataViz v0.9.1 colors.
gdv_blue <- "#009ADE"
gdv_pink <- "#FF1F5B"
gdv_purple <- "#AF58BA"
gdv_gold <- "#FFC61E"
gdv_orange <- "#F28522"
gdv_green <- "#00CD6C"
gdv_ink <- "#202124"
gdv_grey <- "#A0B1BA"
evidence_palette <- c(
  "Score and BH q thresholds" = gdv_pink,
  "BH q threshold only" = gdv_orange,
  "Score threshold only" = gdv_blue,
  "Neither threshold" = "#A0B1BA"
)

# Panel a: support before and after within-participant aggregation.
support <- data.frame(
  condition = factor(c("Early\nexpedition", "Early\nexpedition", "Late /\nrecovery", "Late /\nrecovery"),
                     levels = c("Early\nexpedition", "Late /\nrecovery")),
  measure = factor(c("Raw samples", "Independent profiles", "Raw samples", "Independent profiles"),
                   levels = c("Raw samples", "Independent profiles")),
  count = c(manifest$group_a$raw_samples, manifest$group_a$independent_profiles,
            manifest$group_b$raw_samples, manifest$group_b$independent_profiles)
)
panel_a <- ggplot(support, aes(condition, count, fill = measure)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, color = "white", linewidth = 0.3) +
  geom_text(aes(label = count), position = position_dodge(width = 0.72), vjust = -0.35, size = 2.35, color = gdv_ink) +
  annotate("label", x = 1.5, y = max(support$count) * 1.15, label = paste0(manifest$retained_taxa, " retained genera\n", nrow(results), " received NetMoss scores"),
           size = 2.25, linewidth = 0, fill = "#F4F6F8", color = gdv_ink) +
  scale_fill_manual(values = c("Raw samples" = "#99D7F0", "Independent profiles" = gdv_blue),
                    labels = c("Raw samples" = "Samples", "Independent profiles" = "Profiles")) +
  scale_y_continuous(limits = c(0, max(support$count) * 1.3), expand = expansion(mult = c(0, 0.02))) +
  labs(tag = "a", title = "Participant support", subtitle = "Median aggregation within participant", x = NULL, y = "Number entering each step", fill = NULL) +
  theme_nature(base_size = 6.6) + theme(legend.position = "top", axis.text.x = element_text(size = 6.2))

# Panel b: ranked score landscape.
top_rank <- results %>% slice_head(n = as.integer(cfg$figure$top_taxa_ranked_panel))
panel_b <- ggplot(results, aes(rank, netmoss_score, color = evidence_group)) +
  geom_hline(yintercept = score_threshold, linetype = 2, linewidth = 0.4, color = "#59636E") +
  geom_segment(aes(xend = rank, y = 0, yend = netmoss_score), color = "#D9DDE2", linewidth = 0.25) +
  geom_point(size = 1.65, alpha = 0.9) +
  geom_text_repel(
    data = top_rank %>% slice_head(n = 4), aes(label = genus),
    size = 1.9, color = gdv_ink, box.padding = 0.18, point.padding = 0.12,
    min.segment.length = 0, segment.size = 0.18, seed = as.integer(cfg$figure$random_seed), max.overlaps = Inf
  ) +
  scale_color_manual(values = evidence_palette, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1.08), breaks = c(0, 0.3, 0.6, 0.9), expand = expansion(mult = c(0, 0.01))) +
  scale_x_continuous(breaks = unique(c(1, pretty(c(1, nrow(results))), nrow(results)))) +
  labs(tag = "b", title = "Ranked NetMoss scores", subtitle = paste0("Dashed line: score threshold = ", score_threshold), x = "Genus rank", y = "NetMoss score", color = NULL) +
  theme_nature(base_size = 6.6) +
  theme(legend.position = "none")

# Panel c: marginal abundance change is distinct from a network score.
rho <- suppressWarnings(cor(results$netmoss_score, abs(results$log2_fold_change_late_vs_early), method = "spearman"))
label_c <- bind_rows(
  results %>% slice_head(n = 5),
  results %>% slice_max(abs(log2_fold_change_late_vs_early), n = 2, with_ties = FALSE)
) %>% distinct(taxon_names, .keep_all = TRUE)
panel_c <- ggplot(results, aes(netmoss_score, log2_fold_change_late_vs_early)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "#59636E") +
  geom_vline(xintercept = score_threshold, linetype = 2, linewidth = 0.35, color = "#59636E") +
  geom_point(aes(fill = log2_fold_change_late_vs_early >= 0, shape = fdr_significant),
             size = 1.75, alpha = 0.82, color = "white", stroke = 0.3) +
  geom_text_repel(
    data = label_c, aes(label = genus), size = 1.85, color = gdv_ink,
    box.padding = 0.16, point.padding = 0.1, min.segment.length = 0,
    segment.size = 0.18, seed = as.integer(cfg$figure$random_seed), max.overlaps = Inf
  ) +
  annotate("label", x = 0.03, y = Inf, hjust = 0, vjust = 1.12,
           label = sprintf("Spearman rho = %.2f", rho), size = 2.1,
           linewidth = 0, fill = "#F4F6F8", color = gdv_ink) +
  scale_fill_manual(values = c(`TRUE` = gdv_orange, `FALSE` = gdv_blue),
                    labels = c(`TRUE` = "Higher late", `FALSE` = "Higher early")) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 21), guide = "none") +
  scale_x_continuous(limits = c(0, 1.03), breaks = c(0, 0.3, 0.6, 0.9)) +
  labs(tag = "c", title = "Score versus abundance shift", subtitle = "Mean relative abundance: late versus early", x = "NetMoss score", y = expression(log[2]~"fold change"), fill = NULL) +
  theme_nature(base_size = 6.6) + theme(legend.position = "top")

# Panel d: shared-layout enhanced networks. Color denotes taxonomy, never an
# inferred NetMoss module. Only retained associations receive solid emphasis.
layout_edges <- edges %>%
  transmute(taxon_1, taxon_2, weight = abs_correlation) %>%
  group_by(taxon_1, taxon_2) %>%
  summarise(weight = max(weight), .groups = "drop")
graph <- graph_from_data_frame(layout_edges, directed = FALSE, vertices = data.frame(name = taxa))
set.seed(as.integer(cfg$figure$random_seed))
layout <- layout_with_fr(graph, weights = E(graph)$weight, niter = 3000, grid = "nogrid")
positions <- data.frame(taxon_names = V(graph)$name, x = rescale(layout[, 1], c(-1, 1)), y = rescale(layout[, 2], c(-1, 1)))
top_phyla <- results %>% count(phylum, sort = TRUE) %>% slice_head(n = 5) %>% pull(phylum)
nodes <- positions %>%
  left_join(results, by = "taxon_names") %>%
  mutate(phylum_group = if_else(phylum %in% top_phyla, phylum, "Other / unclassified"))
node_panel <- tidyr::crossing(nodes, condition = factor(unname(condition_labels), levels = unname(condition_labels)))
retained_taxa <- unique(c(edges$taxon_1[edges$retained], edges$taxon_2[edges$retained]))
label_taxa <- unique(c(results$taxon_names[seq_len(min(as.integer(cfg$figure$top_taxa_network_labels), nrow(results)))], retained_taxa))
label_panel <- node_panel %>% filter(taxon_names %in% label_taxa)
curve_data <- curve_edges(edges, positions)
curve_data$condition <- factor(curve_data$condition, levels = unname(condition_labels))
edge_counts <- edges %>%
  group_by(condition) %>%
  summarise(retained_edges = sum(retained), context_edges = sum(!retained), .groups = "drop") %>%
  tidyr::complete(condition = factor(unname(condition_labels), levels = unname(condition_labels)), fill = list(retained_edges = 0L, context_edges = 0L)) %>%
  mutate(x = -Inf, y = Inf, note = paste0("Retained: ", retained_edges, " | context: ", context_edges))
phylum_palette <- c(c(gdv_blue, gdv_gold, gdv_orange, "#6FA253", gdv_purple)[seq_along(top_phyla)], gdv_grey)
names(phylum_palette) <- c(top_phyla, "Other / unclassified")
panel_d <- ggplot() +
  geom_path(
    data = curve_data,
    aes(x, y, group = edge_id, color = correlation_sign, linewidth = abs_correlation,
        alpha = evidence_class, linetype = evidence_class), lineend = "round"
  ) +
  geom_point(data = node_panel, aes(x, y, size = netmoss_score, fill = phylum_group),
             shape = 21, color = "white", stroke = 0.45) +
  geom_point(data = filter(node_panel, fdr_significant), aes(x, y, size = netmoss_score),
             shape = 21, fill = NA, color = gdv_ink, stroke = 0.55, show.legend = FALSE) +
  geom_text_repel(
    data = label_panel, aes(x, y, label = genus), size = 1.65, color = "#303740",
    box.padding = 0.12, point.padding = 0.08, min.segment.length = 0,
    segment.color = alpha("#68717B", 0.35), segment.size = 0.15,
    seed = as.integer(cfg$figure$random_seed), max.overlaps = Inf
  ) +
  geom_label(data = edge_counts, aes(x, y, label = note), hjust = -0.02, vjust = 1.05,
             size = 1.9, linewidth = 0, fill = alpha("white", 0.85), color = "#4E5965") +
  facet_wrap(~condition, nrow = 1) +
  coord_equal(clip = "off") +
  scale_fill_manual(values = phylum_palette, drop = FALSE) +
  scale_color_manual(values = c("Positive" = gdv_orange, "Negative" = gdv_blue)) +
  scale_size_continuous(range = c(1.2, 4.5), limits = c(0, 1), breaks = c(0.3, 0.6, 0.9)) +
  scale_linewidth_continuous(range = c(0.18, 0.95), limits = c(context_threshold, 1), guide = "none") +
  scale_alpha_manual(values = c("Retained association" = 0.88, "Below-threshold context" = 0.14)) +
  scale_linetype_manual(values = c("Retained association" = "solid", "Below-threshold context" = "22")) +
  guides(
    fill = guide_legend(title = "Phylum", override.aes = list(size = 3), order = 4, nrow = 1),
    color = guide_legend(title = "Correlation", override.aes = list(linewidth = 0.9, alpha = 0.9), order = 3),
    size = guide_legend(title = "NetMoss score", override.aes = list(fill = "#D7E6F4", color = "#5B6570"), order = 1),
    alpha = guide_legend(title = "Edge status", override.aes = list(linewidth = 0.9), order = 2),
    linetype = "none"
  ) +
  labs(
    tag = "d",
    title = "Condition-specific genus association networks",
    subtitle = paste0("Solid edges meet |r| >= ", sprintf("%.2f", retained_threshold),
                      "; faint dashed edges (|r| >= ", sprintf("%.2f", context_threshold), ") are context only")
  ) +
  theme_void(base_size = 6.8, base_family = "sans") +
  theme(
    text = element_text(color = gdv_ink),
    plot.title = element_text(face = "bold", size = 8, margin = margin(b = 2, l = 9)),
    plot.subtitle = element_text(size = 6.8, color = "#55585C", margin = margin(b = 3, l = 9)),
    plot.tag = element_text(face = "bold", size = 10.3),
    plot.tag.position = c(0, 1),
    strip.text = element_text(face = "bold", size = 7.5, lineheight = 1.0),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 6.3),
    legend.text = element_text(size = 6.0),
    legend.key.width = unit(4, "mm"),
    legend.key.height = unit(2.6, "mm"),
    panel.spacing = unit(5, "mm"),
    plot.margin = margin(5, 7, 4, 7)
  )

figure_top <- arrangeGrob(panel_a, panel_b, panel_c, ncol = 3, widths = c(0.9, 1.05, 1.35))
main_figure <- arrangeGrob(figure_top, panel_d, ncol = 1, heights = c(0.92, 1.55))
figure_stem <- file.path(paths$main, "Main_Figure_01_exploratory_NetMoss_landscape")
width_in <- as.numeric(cfg$figure$width_mm) / 25.4
height_in <- as.numeric(cfg$figure$height_mm) / 25.4
ggsave(paste0(figure_stem, ".pdf"), main_figure, width = width_in, height = height_in, device = cairo_pdf, bg = "white")
ggsave(paste0(figure_stem, ".svg"), main_figure, width = width_in, height = height_in, device = svglite, bg = "white")
ggsave(paste0(figure_stem, ".png"), main_figure, width = width_in, height = height_in, dpi = as.integer(cfg$figure$dpi), bg = "white")

# Main Table 1: analytical scope and evidential boundary.
table1 <- data.frame(
  Component = c(
    "Biological contrast", "Raw abundance samples", "Independent network profiles",
    "Taxonomic resolution", "Feature screen", "Network inference",
    "NetMoss output", "Retained network structure", "Scientific status"
  ),
  Result = c(
    "T1/T2/Ta versus T8/T9/T10",
    paste0(manifest$group_a$raw_samples, " early; ", manifest$group_b$raw_samples, " late"),
    paste0(manifest$group_a$independent_profiles, " early; ", manifest$group_b$independent_profiles, " late"),
    "Genus",
    paste0(manifest$retained_taxa, " retained genera at 50% prevalence; ", nrow(results), " scored"),
    "Separate SparCC matrices",
    paste0(sum(results$passes_score_threshold), " score >= ", score_threshold, "; ", sum(results$fdr_significant), " BH q < ", fdr_threshold, "; ",
           sum(results$passes_score_threshold & results$fdr_significant), " met both"),
    paste0(sum(edges$retained & edges$condition == condition_labels[["early"]]), " early edges; ",
           sum(edges$retained & edges$condition == condition_labels[["late"]]), " late edges at |r| >= ", retained_threshold),
    "Historical exploratory/supporting result; not validated biomarker evidence"
  ),
  Reader_interpretation = c(
    "Early expedition was compared with a heterogeneous late/post-expedition period.",
    "Repeated longitudinal observations were available before aggregation.",
    "Participant aggregation reduced pseudoreplication but left very small networks.",
    "MAG abundances were summed to genera before network inference.",
    paste0("Retained taxa = ", manifest$retained_taxa, "; profiles = ", manifest$group_a$independent_profiles, " early and ", manifest$group_b$independent_profiles, " late."),
    "Edges are compositional associations, not direct ecological interactions.",
    "Nominal score and q-value patterns require stability assessment before biological claims.",
    paste0("Retained edges = ", sum(edges$retained), " across both conditions; inspect condition-specific counts."),
    "No stable network rewiring or biomarker conclusion is supported."
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
names(table1)[[3]] <- "Reader interpretation"

# Main Table 2: exact top-ranked genera with abundance context.
table2 <- results %>%
  slice_head(n = as.integer(cfg$figure$top_taxa_table)) %>%
  transmute(
    Genus = genus,
    `NetMoss score` = sprintf("%.3f", netmoss_score),
    `BH q` = format(q_value, digits = 3, scientific = TRUE),
    `Early mean (%)` = sprintf("%.3f", 100 * mean_abundance_early),
    `Late mean (%)` = sprintf("%.3f", 100 * mean_abundance_late),
    `log2FC late/early` = sprintf("%+.2f", log2_fold_change_late_vs_early),
    `Degree early/late` = paste0(degree_early, "/", degree_late)
  )

for (item in list(
  list(x = table1, stem = "Main_Table_01_analysis_scope_and_evidence", widths = c(21, 27, 52), height = 4.8),
  list(x = table2, stem = "Main_Table_02_top_ranked_genera", widths = c(22, 13, 12, 14, 14, 16, 16), height = 3.8)
)) {
  utils::write.csv(item$x, file.path(paths$main, paste0(item$stem, ".csv")), row.names = FALSE)
  utils::write.table(item$x, file.path(paths$main, paste0(item$stem, ".tsv")), row.names = FALSE, sep = "\t", quote = FALSE)
  save_table_preview(item$x, item$stem, paths$main, item$widths, item$height, min(400L, as.integer(cfg$figure$dpi)))
}

# Exact panel and supplementary source data.
utils::write.csv(support, file.path(paths$reproducibility, "Main_Figure_01_panel_a.csv"), row.names = FALSE)
utils::write.csv(results, file.path(paths$reproducibility, "Main_Figure_01_panels_b_c.csv"), row.names = FALSE)
utils::write.csv(nodes, file.path(paths$reproducibility, "Main_Figure_01_panel_d_nodes.csv"), row.names = FALSE)
utils::write.csv(edges %>% mutate(condition = as.character(condition)),
                 file.path(paths$reproducibility, "Main_Figure_01_panel_d_edges.csv"), row.names = FALSE)
utils::write.csv(results, file.path(paths$supplementary, "Supplementary_Table_01_all_taxon_results.csv"), row.names = FALSE)
utils::write.csv(edges %>% mutate(condition = as.character(condition)),
                 file.path(paths$supplementary, "Supplementary_Table_02_displayed_network_edges.csv"), row.names = FALSE)

caption <- paste0(
  "Main Figure 1 | Exploratory NetMoss analysis of participant-aggregated longitudinal gut microbiome profiles. ",
  "a, Raw samples and independent participant profiles contributing to the historical early-expedition (T1/T2/Ta) and late/recovery (T8/T9/T10) periods. Repeated samples within each participant and period were summarized by the median before network construction. ",
  paste0("b, Ranked NetMoss scores for ", nrow(results), " scored genera. The dashed horizontal line denotes the prespecified score threshold of ", score_threshold, "; colors distinguish the score and Benjamini-Hochberg q-value criteria. "),
  "c, NetMoss score plotted against the log2 ratio of mean late to early relative abundance. Positive values indicate higher mean abundance late and negative values indicate higher mean abundance early. The displayed Spearman coefficient relates NetMoss score to the absolute abundance shift and is descriptive. ",
  paste0("d, Shared-layout SparCC association networks for the two periods. Node color denotes taxonomic phylum, node size denotes NetMoss score, and a dark outline denotes BH q < ", fdr_threshold, ". Orange and blue edges denote positive and negative correlations, respectively. Solid edges meet the prespecified |r| >= ", retained_threshold, " threshold; faint dashed edges between |r| >= ", context_threshold, " and < ", retained_threshold, " are shown only as visual context and are not retained associations. Curvature and node placement are visual encodings. "),
  paste0("Associations are not direct ecological interactions. The networks contain ", manifest$group_a$independent_profiles, " and ", manifest$group_b$independent_profiles, " independent profiles for ", manifest$retained_taxa, " retained genera. All findings remain exploratory and do not constitute validated biomarkers or stable network-rewiring evidence.")
)
writeLines(caption, file.path(paths$main, "Main_Figure_01_caption.txt"))

chart_map <- c(
  "# NetMoss manuscript synthesis chart map", "",
  "| Panel | Reader question | Visual form | Fields | Supported statement |",
  "|---|---|---|---|---|",
  "| a | How much independent support entered each network? | Grouped bar | raw samples, participant profiles | Profile counts are calculated from the supplied source manifest. |",
  "| b | How were NetMoss scores distributed? | Ranked stem-and-dot | rank, score, score/FDR status | Score and q classifications use configured thresholds. |",
  "| c | Are network scores equivalent to marginal abundance shifts? | Scatter with zero and score references | score, log2 abundance ratio | The displayed Spearman coefficient describes the supplied score and abundance-shift relationship. |",
  "| d | What association structure remained at the declared threshold? | Faceted force-directed network | SparCC r, score, q, phylum | Retained associations are calculated from the supplied networks at the configured threshold. |",
  "", "Exact GeoDataViz v0.9.1 palette values are declared in the generating script."
)
writeLines(chart_map, file.path(paths$reproducibility, "chart_map.md"))

index <- c(
  "# NetMoss Manuscript Synthesis Results Index", "",
  paste0("Source analysis: `", manifest$scenario, "`, generated ", manifest$generated_at, "."), "",
  "- Main Figure 1: `main/Main_Figure_01_exploratory_NetMoss_landscape.{pdf,png,svg}`",
  "- Main Figure 1 caption: `main/Main_Figure_01_caption.txt`",
  "- Main Table 1: `main/Main_Table_01_analysis_scope_and_evidence.{csv,tsv,pdf,png}`",
  "- Main Table 2: `main/Main_Table_02_top_ranked_genera.{csv,tsv,pdf,png}`",
  "- Exact panel data: `reproducibility/Main_Figure_01_panel_*.csv`",
  "- Supplementary all-taxon table: `supplementary/Supplementary_Table_01_all_taxon_results.csv`",
  "- Supplementary displayed-edge table: `supplementary/Supplementary_Table_02_displayed_network_edges.csv`", "",
  "The preserved result is historical, exploratory, and Ta-inclusive. It must not be relabeled as the current Ta/Tb-excluded configuration."
)
writeLines(index, file.path(paths$root, "RESULTS_INDEX.md"))
capture.output(sessionInfo(), file = file.path(paths$reproducibility, "sessionInfo.txt"))
yaml::write_yaml(list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  source_manifest = normalizePath(required[[1]], mustWork = TRUE),
  source_git_commit = manifest$git_commit,
  source_analysis_generated_at = manifest$generated_at,
  source_scenario = manifest$scenario,
  scored_taxa = nrow(results),
  score_threshold_taxa = sum(results$passes_score_threshold),
  fdr_taxa = sum(results$fdr_significant),
  both_threshold_taxa = sum(results$passes_score_threshold & results$fdr_significant),
  retained_early_edges = sum(edges$retained & edges$condition == condition_labels[["early"]]),
  retained_late_edges = sum(edges$retained & edges$condition == condition_labels[["late"]]),
  scientific_status = "historical_exploratory_supporting_only"
), file.path(paths$reproducibility, "manuscript_synthesis_manifest.yml"))
writeLines(c(
  "NetMoss manuscript synthesis complete.",
  paste0("Source Git commit: ", manifest$git_commit),
  "Source results were read without modification.",
  "No SparCC correlation, NetMoss score, p-value, q-value, or threshold was recalculated.",
  "Scientific status: historical exploratory/supporting only."
), marker)

cat("NetMoss manuscript synthesis completed: ", normalizePath(paths$root, mustWork = TRUE), "\n", sep = "")
