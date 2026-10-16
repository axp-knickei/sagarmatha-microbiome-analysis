#!/usr/bin/env Rscript

# Purpose: Create a visually enhanced, scientifically conservative rendering of
# the preserved exploratory NetMoss condition networks. This is post-processing:
# it never recalculates SparCC correlations, NetMoss scores, or significance.
#
# Inputs:
#   results/netmoss/networks/group_a_network.csv
#   results/netmoss/networks/group_b_network.csv
#   results/netmoss/tables/netmoss_scores.csv
#   results/netmoss/reproducibility/netmoss_manifest.yml
#   data/raw/gut.csv
#   data/raw/GDTB.csv
#
# Outputs:
#   results/netmoss/enhanced_network/figures/NetMoss_enhanced_network.{pdf,svg,png}
#   results/netmoss/enhanced_network/source_data/{nodes,edges}.csv
#   results/netmoss/enhanced_network/figure_notes.txt
#
# Dependencies: dplyr, tidyr, ggplot2, ggrepel, igraph, scales, svglite, yaml.
# Run from the repository root:
#   Rscript scripts/netmoss/create_enhanced_network_figure.R --force
#
# Relationship to adjacent stages: this consumes the historical NetMoss output
# and writes a separate visualization. It does not replace Figure 1 or change the
# exploratory scientific status of the analysis.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(igraph)
  library(scales)
  library(svglite)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x)) || identical(x, "")) y else x
}

parse_args <- function(args) {
  out <- list(
    input_dir = "results/netmoss",
    output_dir = "results/netmoss/enhanced_network",
    abundance = "data/raw/gut.csv",
    taxonomy = "data/raw/GDTB.csv",
    edge_threshold = 0.40,
    context_threshold = 0.25,
    force = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--force") {
      out$force <- TRUE
      i <- i + 1L
    } else if (arg %in% c("--input-dir", "--output-dir", "--abundance", "--taxonomy",
                          "--edge-threshold", "--context-threshold")) {
      if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
      key <- switch(
        arg,
        `--input-dir` = "input_dir",
        `--output-dir` = "output_dir",
        `--abundance` = "abundance",
        `--taxonomy` = "taxonomy",
        `--edge-threshold` = "edge_threshold",
        `--context-threshold` = "context_threshold"
      )
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    } else if (arg == "--help") {
      cat(paste0(
        "Usage: Rscript scripts/netmoss/create_enhanced_network_figure.R [options]\n",
        "  --input-dir PATH          Default: results/netmoss\n",
        "  --output-dir PATH         Default: results/netmoss/enhanced_network\n",
        "  --abundance PATH          Default: data/raw/gut.csv\n",
        "  --taxonomy PATH           Default: data/raw/GDTB.csv\n",
        "  --edge-threshold NUMBER   Preserved analysis threshold; default: 0.40\n",
        "  --context-threshold NUM   Faint descriptive-edge floor; default: 0.25\n",
        "  --force                   Replace this script's enhanced outputs only\n"
      ))
      quit(status = 0)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  out$edge_threshold <- suppressWarnings(as.numeric(out$edge_threshold))
  out$context_threshold <- suppressWarnings(as.numeric(out$context_threshold))
  if (!is.finite(out$edge_threshold) || out$edge_threshold <= 0 || out$edge_threshold > 1) {
    stop("--edge-threshold must be in (0, 1].", call. = FALSE)
  }
  if (!is.finite(out$context_threshold) || out$context_threshold < 0 ||
      out$context_threshold >= out$edge_threshold) {
    stop("--context-threshold must be >= 0 and below --edge-threshold.", call. = FALSE)
  }
  out
}

read_matrix_csv <- function(path) {
  if (!file.exists(path)) stop("Required matrix not found: ", path, call. = FALSE)
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- as.character(x[[1]])
  matrix <- as.matrix(x[-1])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- ids
  if (nrow(matrix) != ncol(matrix)) stop("Network is not square: ", path, call. = FALSE)
  if (!identical(rownames(matrix), colnames(matrix))) {
    stop("Network row and column identifiers differ: ", path, call. = FALSE)
  }
  if (max(abs(matrix - t(matrix)), na.rm = TRUE) > 1e-8) {
    stop("Network is not symmetric: ", path, call. = FALSE)
  }
  matrix
}

normalize_scores <- function(path) {
  if (!file.exists(path)) stop("NetMoss score table not found: ", path, call. = FALSE)
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("taxon_names", "NetMoss_Score", "p.adj")
  if (length(setdiff(required, names(x)))) {
    stop("Score table lacks required columns: ", paste(setdiff(required, names(x)), collapse = ", "), call. = FALSE)
  }
  x %>%
    transmute(
      taxon_names = as.character(taxon_names),
      netmoss_score = as.numeric(NetMoss_Score),
      q_value = as.numeric(p.adj),
      fdr_significant = is.finite(q_value) & q_value < 0.05
    )
}

build_taxonomy <- function(abundance_path, taxonomy_path, taxa) {
  if (!file.exists(abundance_path)) stop("Abundance file not found: ", abundance_path, call. = FALSE)
  if (!file.exists(taxonomy_path)) stop("Taxonomy file not found: ", taxonomy_path, call. = FALSE)
  abundance <- utils::read.csv(abundance_path, check.names = FALSE, stringsAsFactors = FALSE)
  taxonomy <- utils::read.csv(taxonomy_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("Genome", "genus", "phylum")
  if (length(setdiff(required, names(taxonomy)))) {
    stop("Taxonomy file lacks Genome/genus/phylum columns.", call. = FALSE)
  }
  if (!"sampleID" %in% names(abundance)) stop("Abundance file lacks sampleID.", call. = FALSE)
  mag_ids <- setdiff(names(abundance), "sampleID")
  selected_taxonomy <- data.frame(Genome = mag_ids, stringsAsFactors = FALSE) %>%
    left_join(taxonomy %>% select(Genome, genus, phylum), by = "Genome")
  if (anyNA(selected_taxonomy$genus) || anyNA(selected_taxonomy$phylum)) {
    stop("One or more gut MAG columns did not join exactly to GDTB.csv Genome.", call. = FALSE)
  }
  genus_phylum <- selected_taxonomy %>%
    filter(!is.na(genus), nzchar(genus), genus != "g__") %>%
    count(genus, phylum, name = "mag_count") %>%
    arrange(genus, desc(mag_count), phylum) %>%
    group_by(genus) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(taxon_names = genus, phylum = if_else(is.na(phylum) | phylum == "p__", "Unclassified", sub("^p__", "", phylum)))

  data.frame(taxon_names = taxa, stringsAsFactors = FALSE) %>%
    left_join(genus_phylum, by = "taxon_names") %>%
    mutate(phylum = replace_na(phylum, "Unclassified"))
}

matrix_edges <- function(network, condition, context_threshold, edge_threshold) {
  index <- which(upper.tri(network), arr.ind = TRUE)
  data.frame(
    taxon_1 = rownames(network)[index[, 1]],
    taxon_2 = colnames(network)[index[, 2]],
    correlation = network[index],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      condition = condition,
      abs_correlation = abs(correlation)
    ) %>%
    filter(abs_correlation >= context_threshold) %>%
    mutate(
      evidence_class = if_else(
        abs_correlation >= edge_threshold,
        paste0("Retained: |r| ≥ ", sprintf("%.2f", edge_threshold)),
        paste0("Context only: ", sprintf("%.2f", context_threshold), " ≤ |r| < ", sprintf("%.2f", edge_threshold))
      ),
      correlation_sign = if_else(correlation >= 0, "Positive correlation", "Negative correlation")
    )
}

curve_edges <- function(edges, positions, points = 30L) {
  if (!nrow(edges)) return(data.frame())
  edges <- edges %>%
    mutate(edge_id = paste(condition, taxon_1, taxon_2, sep = "|")) %>%
    left_join(positions %>% select(taxon_1 = taxon_names, x1 = x, y1 = y), by = "taxon_1") %>%
    left_join(positions %>% select(taxon_2 = taxon_names, x2 = x, y2 = y), by = "taxon_2")
  pieces <- lapply(seq_len(nrow(edges)), function(i) {
    edge <- edges[i, ]
    dx <- edge$x2 - edge$x1
    dy <- edge$y2 - edge$y1
    distance <- sqrt(dx^2 + dy^2)
    direction <- if ((match(edge$taxon_1, positions$taxon_names) +
                       match(edge$taxon_2, positions$taxon_names)) %% 2L) 1 else -1
    curvature <- direction * 0.14 * distance
    midpoint_x <- (edge$x1 + edge$x2) / 2
    midpoint_y <- (edge$y1 + edge$y2) / 2
    control_x <- midpoint_x - curvature * dy / max(distance, .Machine$double.eps)
    control_y <- midpoint_y + curvature * dx / max(distance, .Machine$double.eps)
    t <- seq(0, 1, length.out = points)
    data.frame(
      edge_id = edge$edge_id,
      condition = edge$condition,
      correlation = edge$correlation,
      abs_correlation = edge$abs_correlation,
      evidence_class = edge$evidence_class,
      correlation_sign = edge$correlation_sign,
      x = (1 - t)^2 * edge$x1 + 2 * (1 - t) * t * control_x + t^2 * edge$x2,
      y = (1 - t)^2 * edge$y1 + 2 * (1 - t) * t * control_y + t^2 * edge$y2,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(pieces)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
manifest_path <- file.path(args$input_dir, "reproducibility", "netmoss_manifest.yml")
if (!file.exists(manifest_path)) stop("NetMoss manifest not found: ", manifest_path, call. = FALSE)
manifest <- yaml::read_yaml(manifest_path)
if (!identical(as.character(manifest$status), "COMPLETED_EXPLORATORY")) {
  stop("Expected a COMPLETED_EXPLORATORY NetMoss manifest.", call. = FALSE)
}

network_a <- read_matrix_csv(file.path(args$input_dir, "networks", "group_a_network.csv"))
network_b <- read_matrix_csv(file.path(args$input_dir, "networks", "group_b_network.csv"))
scores <- normalize_scores(file.path(args$input_dir, "tables", "netmoss_scores.csv"))
taxa <- scores$taxon_names
if (length(setdiff(taxa, rownames(network_a))) || length(setdiff(taxa, rownames(network_b)))) {
  stop("One or both network matrices lack scored taxa.", call. = FALSE)
}
network_a <- network_a[taxa, taxa, drop = FALSE]
network_b <- network_b[taxa, taxa, drop = FALSE]

group_a_label <- as.character(manifest$group_a$label %||% "Group A")
group_b_label <- as.character(manifest$group_b$label %||% "Group B")
display_condition <- c(
  setNames(paste0("Early expedition\n", paste(unlist(manifest$group_a$timepoints), collapse = "/"),
                  "; n = ", manifest$group_a$independent_profiles), group_a_label),
  setNames(paste0("Late/recovery\n", paste(unlist(manifest$group_b$timepoints), collapse = "/"),
                  "; n = ", manifest$group_b$independent_profiles), group_b_label)
)

edges <- bind_rows(
  matrix_edges(network_a, group_a_label, args$context_threshold, args$edge_threshold),
  matrix_edges(network_b, group_b_label, args$context_threshold, args$edge_threshold)
) %>%
  mutate(condition = factor(condition, levels = c(group_a_label, group_b_label), labels = unname(display_condition)))

# A shared force-directed layout supports direct comparison between conditions.
# All context-level edges affect placement, but placement itself is descriptive
# and is never interpreted as an additional inferred community or test result.
layout_edges <- edges %>%
  transmute(taxon_1, taxon_2, weight = abs_correlation) %>%
  group_by(taxon_1, taxon_2) %>%
  summarise(weight = max(weight), .groups = "drop")
graph <- graph_from_data_frame(
  layout_edges,
  directed = FALSE,
  vertices = data.frame(name = taxa, stringsAsFactors = FALSE)
)
set.seed(20260714)
layout <- layout_with_fr(graph, weights = E(graph)$weight, niter = 3000, grid = "nogrid")
positions <- data.frame(taxon_names = V(graph)$name, x = layout[, 1], y = layout[, 2], stringsAsFactors = FALSE)
positions$x <- rescale(positions$x, to = c(-1, 1))
positions$y <- rescale(positions$y, to = c(-1, 1))

taxonomy <- build_taxonomy(args$abundance, args$taxonomy, taxa)
top_phyla <- taxonomy %>% count(phylum, sort = TRUE) %>% slice_head(n = 5) %>% pull(phylum)
nodes <- positions %>%
  left_join(scores, by = "taxon_names") %>%
  left_join(taxonomy, by = "taxon_names") %>%
  mutate(
    phylum_group = if_else(phylum %in% top_phyla, phylum, "Other / unclassified"),
    label = sub("^g__", "", taxon_names)
  )
node_panel <- tidyr::crossing(
  nodes,
  condition = factor(unname(display_condition), levels = unname(display_condition))
)

curve_data <- curve_edges(edges, positions)
if (nrow(curve_data)) {
  curve_data$condition <- factor(curve_data$condition, levels = unname(display_condition))
}

counts <- edges %>%
  count(condition, evidence_class, name = "edges") %>%
  tidyr::complete(condition, evidence_class, fill = list(edges = 0L))
strong_label <- paste0("Retained: |r| ≥ ", sprintf("%.2f", args$edge_threshold))
context_label <- paste0("Context only: ", sprintf("%.2f", args$context_threshold),
                        " ≤ |r| < ", sprintf("%.2f", args$edge_threshold))
panel_notes <- counts %>%
  pivot_wider(names_from = evidence_class, values_from = edges) %>%
  mutate(
    retained = .data[[strong_label]] %||% 0L,
    context = .data[[context_label]] %||% 0L,
    note = paste0("Retained edges: ", retained, "  |  below-threshold context: ", context),
    x = -Inf,
    y = Inf
  )

phylum_palette <- c(
  "#2F6B9A", "#E0A11A", "#D96C3F", "#718C42", "#C05A7A",
  "#9AA1A8"
)
names(phylum_palette) <- c(top_phyla, "Other / unclassified")
edge_palette <- c("Positive correlation" = "#D96C3F", "Negative correlation" = "#2F6B9A")
edge_alpha <- setNames(c(0.82, 0.16), c(strong_label, context_label))
edge_linetype <- setNames(c("solid", "22"), c(strong_label, context_label))

plot <- ggplot() +
  geom_path(
    data = curve_data,
    aes(x, y, group = edge_id, color = correlation_sign,
        linewidth = abs_correlation, alpha = evidence_class, linetype = evidence_class),
    lineend = "round"
  ) +
  geom_point(
    data = node_panel,
    aes(x, y, size = netmoss_score, fill = phylum_group),
    shape = 21, color = "white", stroke = 0.7
  ) +
  geom_point(
    data = filter(node_panel, fdr_significant),
    aes(x, y, size = netmoss_score),
    shape = 21, fill = NA, color = "#20262E", stroke = 0.9, show.legend = FALSE
  ) +
  geom_text_repel(
    data = node_panel,
    aes(x, y, label = label),
    size = 2.15,
    color = "#303740",
    family = "sans",
    box.padding = 0.18,
    point.padding = 0.10,
    min.segment.length = 0,
    segment.color = alpha("#68717B", 0.32),
    segment.size = 0.18,
    max.overlaps = Inf,
    seed = 20260714,
    force = 1.1
  ) +
  geom_label(
    data = panel_notes,
    aes(x, y, label = note),
    hjust = -0.02,
    vjust = 1.05,
    size = 2.5,
    linewidth = 0,
    fill = alpha("white", 0.86),
    color = "#4E5965",
    inherit.aes = FALSE
  ) +
  facet_wrap(~condition, nrow = 1) +
  coord_equal(clip = "off") +
  scale_fill_manual(values = phylum_palette, drop = FALSE) +
  scale_color_manual(values = edge_palette, drop = FALSE) +
  scale_size_continuous(range = c(2.4, 8.2), limits = c(0, 1), breaks = c(0.25, 0.50, 0.75, 1.00)) +
  scale_linewidth_continuous(range = c(0.25, 1.35), limits = c(args$context_threshold, 1)) +
  scale_alpha_manual(values = edge_alpha, drop = FALSE) +
  scale_linetype_manual(values = edge_linetype, drop = FALSE) +
  guides(
    fill = guide_legend(title = "Taxonomic phylum", override.aes = list(size = 4)),
    color = guide_legend(title = "Association sign", override.aes = list(linewidth = 1.1, alpha = 0.9)),
    size = guide_legend(
      title = "NetMoss score",
      override.aes = list(fill = "#D7E6F4", color = "#5B6570", alpha = 1)
    ),
    linewidth = "none",
    alpha = guide_legend(title = "Edge status", override.aes = list(linewidth = 1.1)),
    linetype = "none"
  ) +
  labs(
    title = "Exploratory NetMoss association networks",
    subtitle = paste0(
      "Historical early-versus-late run; shared layout across panels. Node color = phylum, size = NetMoss score, dark ring = BH q < 0.05.\n",
      "Solid edges meet the prespecified |r| ≥ ", sprintf("%.2f", args$edge_threshold),
      " threshold; faint dashed edges are descriptive context only (|r| ≥ ",
      sprintf("%.2f", args$context_threshold), ")."
    ),
    caption = paste0(
      "SparCC associations are not direct ecological interactions. The layout and curved paths are visual encodings only. ",
      "This 8-versus-12-participant analysis remains exploratory and is not validated biomarker evidence."
    )
  ) +
  theme_void(base_size = 10, base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "#FBFCFD", color = NA),
    panel.background = element_rect(fill = "#FBFCFD", color = NA),
    plot.title = element_text(size = 17, face = "bold", color = "#20262E", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 9.5, color = "#4E5965", lineheight = 1.15, margin = margin(b = 10)),
    plot.caption = element_text(size = 8.5, color = "#59636E", hjust = 0, lineheight = 1.1, margin = margin(t = 9)),
    strip.text = element_text(size = 11.5, face = "bold", color = "#20262E", lineheight = 1.05, margin = margin(b = 4)),
    panel.spacing = unit(1.2, "lines"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 3),
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(14, 20, 12, 20)
  )

figures_dir <- file.path(args$output_dir, "figures")
source_dir <- file.path(args$output_dir, "source_data")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
sentinel <- file.path(args$output_dir, "figure_notes.txt")
if (file.exists(sentinel) && !args$force) {
  stop("Refusing to replace enhanced outputs. Use --force.", call. = FALSE)
}

stem <- file.path(figures_dir, "NetMoss_enhanced_network")
ggsave(paste0(stem, ".pdf"), plot, width = 18, height = 10.5, device = grDevices::cairo_pdf)
ggsave(paste0(stem, ".svg"), plot, width = 18, height = 10.5, device = svglite::svglite)
ggsave(paste0(stem, ".png"), plot, width = 18, height = 10.5, dpi = 300, bg = "#FBFCFD")

utils::write.csv(
  nodes %>% arrange(desc(netmoss_score), taxon_names),
  file.path(source_dir, "NetMoss_enhanced_network_nodes.csv"),
  row.names = FALSE
)
utils::write.csv(
  edges %>% mutate(condition = as.character(condition)) %>% arrange(condition, desc(abs_correlation), taxon_1, taxon_2),
  file.path(source_dir, "NetMoss_enhanced_network_edges.csv"),
  row.names = FALSE
)

strong_counts <- edges %>%
  filter(evidence_class == strong_label) %>%
  count(condition, name = "retained_edges") %>%
  tidyr::complete(condition = factor(unname(display_condition), levels = unname(display_condition)), fill = list(retained_edges = 0L))
notes <- c(
  "Enhanced exploratory NetMoss network visualization",
  "=================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Source analysis generated: ", manifest$generated_at),
  paste0("Source analysis Git commit: ", manifest$git_commit),
  paste0("Source scenario: ", manifest$scenario),
  paste0("Retained-edge threshold: |r| >= ", sprintf("%.2f", args$edge_threshold)),
  paste0("Descriptive context floor: |r| >= ", sprintf("%.2f", args$context_threshold)),
  paste0(as.character(strong_counts$condition), ": ", strong_counts$retained_edges, " retained edges"),
  "Node colors encode taxonomic phylum; they do not represent inferred NetMoss modules.",
  "The shared force-directed layout uses context-level edge magnitudes only for placement.",
  "Faint dashed edges are below the analysis threshold and must not be reported as retained associations.",
  "No network values, NetMoss scores, p-values, q-values, or scientific thresholds were recalculated.",
  "The source analysis remains exploratory/supporting only because it contains 8 versus 12 independent profiles."
)
writeLines(notes, sentinel)

cat("Enhanced NetMoss network outputs completed: ", normalizePath(args$output_dir, mustWork = TRUE), "\n", sep = "")
