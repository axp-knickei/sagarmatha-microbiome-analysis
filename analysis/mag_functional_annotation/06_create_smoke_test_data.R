#!/usr/bin/env Rscript

# Purpose: Create deterministic synthetic MMseqs, eggNOG, abundance, metadata,
# and taxonomy fixtures for end-to-end workflow validation.
# Inputs: Production configuration for schema defaults only; no scientific data.
# Outputs: Synthetic files and a smoke-test YAML under a caller-selected directory.
# Dependencies: R, data.table, yaml.
# Usage: Rscript analysis/mag_functional_annotation/06_create_smoke_test_data.R --output-dir /tmp/magfunc_smoke
# Workflow: Test-only fixture generator; outputs must never be interpreted biologically.

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
output_dir <- "/tmp/blackice_magfunc_smoke"
config_path <- "config/mag_functional_annotation.yml"
i <- 1L
while (i <= length(args)) {
  if (args[[i]] == "--output-dir" && i < length(args)) {
    output_dir <- args[[i + 1L]]; i <- i + 2L
  } else if (args[[i]] == "--config" && i < length(args)) {
    config_path <- args[[i + 1L]]; i <- i + 2L
  } else {
    stop("Unknown argument: ", args[[i]], call. = FALSE)
  }
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data_dir <- file.path(output_dir, "data")
result_dir <- file.path(output_dir, "results")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

mags <- c(
  "QML1234-000101-91B_bin.1.strict",
  "QML1234-000102-92B_bin.2.permissive",
  "QML1234-000103-93C_bin.3.orig",
  "QML1234-000104-94D_bin.4.strict",
  "QML1234-000105-95B_bin.5.permissive"
)
protein <- function(mag, node, gene) paste0(mag, "_NODE_", node, "_length_1000_cov_2.0_", gene)
representatives <- vapply(seq_len(8L), function(i) protein(mags[((i - 1L) %% 5L) + 1L], 100L + i, 1L), character(1))
cluster <- rbindlist(lapply(seq_along(representatives), function(i) {
  if (i %in% c(5L, 7L, 8L)) {
    return(data.table(representative = representatives[[i]], member = representatives[[i]]))
  }
  member_mags <- if (i %in% c(1L, 3L, 6L)) c(mags[((i - 1L) %% 5L) + 1L], mags[(i %% 5L) + 1L]) else mags[((i - 1L) %% 5L) + 1L]
  members <- c(representatives[[i]], vapply(seq_along(member_mags), function(j) protein(member_mags[[j]], 500L + i * 10L + j, j + 1L), character(1)))
  data.table(representative = representatives[[i]], member = unique(members))
}))
fwrite(cluster, file.path(data_dir, "cluster.tsv"), sep = "\t", col.names = FALSE)

annotation <- data.table(
  query = representatives[1:7],
  seed_ortholog = paste0("1234.TEST_", 1:7),
  evalue = c(0, 1e-100, 1e-80, 1e-40, 1e-20, 1e-10, 1e-5),
  score = c(900, 700, 600, 400, 250, 150, 80),
  eggNOG_OGs = paste0("COG000", 1:7, "@2|Bacteria"),
  max_annot_lvl = "2|Bacteria",
  COG_category = c("G", "C", "GM", "P", "S", "E", "O"),
  Description = c("carbohydrate enzyme", "energy enzyme", "membrane protein", "transport protein", "uncharacterized protein", "amino acid enzyme", "oxidative stress protein"),
  Preferred_name = c("geneA", "geneB", "geneC", "geneD", "-", "geneF", "geneG"),
  GOs = "-",
  EC = c("1.1.1.1", "2.2.2.2", "3.3.3.3", "-", "1.1.1.1", "4.4.4.4", "1.11.1.6"),
  KEGG_ko = c("ko:K00001", "ko:K00002", "ko:K00001", "ko:K00003", "-", "ko:K00004", "ko:K00005"),
  KEGG_Pathway = c("ko00010,map00010", "ko00020", "ko00010", "ko02010", "-", "ko00250", "ko00480"),
  KEGG_Module = c("M00001", "M00002", "M00001", "M00240", "-", "M00003", "M00118"),
  KEGG_Reaction = "-",
  KEGG_rclass = "-",
  BRITE = "ko00000",
  KEGG_TC = "-",
  CAZy = c("GH13", "-", "GH13,CBM48", "-", "-", "-", "-"),
  BiGG_Reaction = "-",
  PFAMs = c("PFAM_A", "PFAM_B", "PFAM_A,PFAM_C", "PFAM_D", "DUF999", "PFAM_E", "PFAM_F")
)
fwrite(annotation, file.path(data_dir, "annotations.tsv"), sep = "\t", quote = FALSE)

taxonomy <- data.table(
  Genome = mags,
  label = paste0("MAG", seq_along(mags)),
  domain = "d__Bacteria",
  phylum = c("p__Bacillota", "p__Bacillota", "p__Bacteroidota", "p__Actinomycetota", "p__Bacillota"),
  class = c("c__Clostridia", "c__Clostridia", "c__Bacteroidia", "c__Actinomycetia", "c__Clostridia"),
  ordeer = paste0("o__Order", 1:5),
  family = paste0("f__Family", 1:5),
  genus = paste0("g__Genus", 1:5),
  species = paste0("s__Species", 1:5),
  species_note = paste0("s__Species", 1:5)
)
fwrite(taxonomy, file.path(data_dir, "GDTB.csv"))

timepoints <- c("T1", "T2", "Ta", "T3", "Tb", "T4", "T5", "T6", "T7", "T8", "T9", "T10")
dates <- as.Date("2000-01-01") + seq(0, by = 7, length.out = length(timepoints))
# Artificial dates: synthetic fixture metadata are not copied from study records.
altitudes <- c(0, 3650, 3800, 5200, 5800, 5200, 6500, 5200, 5200, 5200, 0, 0)
metadata <- CJ(group = sprintf("S%02d", 1:6), timepoint_index = seq_along(timepoints))
metadata[, `:=`(
  Timepoint = timepoints[timepoint_index],
  date = format(dates[timepoint_index], "%Y%m%d"),
  gender = ifelse(as.integer(sub("S", "", group)) %% 2L == 0L, "female", "male"),
  age = 30L + as.integer(sub("S", "", group)),
  altitude = altitudes[timepoint_index],
  sampleID = paste(group, timepoints[timepoint_index], sep = "_")
)]
setcolorder(metadata, c("sampleID", "group", "Timepoint", "date", "gender", "age", "altitude", "timepoint_index"))
metadata[, timepoint_index := NULL]
fwrite(metadata, file.path(data_dir, "sample_info.csv"))

gut <- metadata[, .(sampleID, group, Timepoint)]
for (j in seq_len(5L)) {
  subject_number <- as.integer(sub("S", "", gut$group))
  time_number <- match(gut$Timepoint, timepoints)
  value <- 0.08 + 0.01 * j + 0.004 * subject_number + 0.015 * sin((time_number + j) / 2)
  if (j == 1L) value <- value + ifelse(gut$Timepoint %in% c("T5", "T6"), 0.08, 0)
  if (j == 2L) value <- value + ifelse(gut$Timepoint %in% c("T7", "T8"), 0.06, 0)
  gut[, (mags[[j]]) := pmax(0, value)]
}
gut[, c("group", "Timepoint") := NULL]
fwrite(gut, file.path(data_dir, "gut.csv"))

cfg <- yaml::read_yaml(config_path)
cfg$inputs$cluster_file <- file.path(data_dir, "cluster.tsv")
cfg$inputs$annotation_file <- file.path(data_dir, "annotations.tsv")
cfg$inputs$abundance_file <- file.path(data_dir, "gut.csv")
cfg$inputs$metadata_file <- file.path(data_dir, "sample_info.csv")
cfg$inputs$taxonomy_file <- file.path(data_dir, "GDTB.csv")
cfg$quality_control$reported_cluster_lines <- nrow(cluster)
cfg$quality_control$reported_annotation_lines <- nrow(annotation) + 1L
cfg$quality_control$enforce_reported_line_counts <- TRUE
cfg$analysis$minimum_mag_carriers <- 1L
cfg$analysis$minimum_nonzero_samples <- 5L
cfg$analysis$minimum_subjects <- 3L
cfg$analysis$minimum_paired_subjects <- 3L
cfg$analysis$threads <- 1L
cfg$followup$candidates[[1]]$term <- "K00001"
cfg$followup$candidates[[1]]$selection_basis <- "synthetic smoke-test candidate"
cfg$followup$candidates[[1]]$curation <- list(
  accepted_name = "Synthetic microcompartment shell candidate",
  gene_symbol = "synA",
  entity_type = "synthetic protein",
  reaction = "",
  authoritative_sources = list(list(
    name = "Synthetic fixture source A",
    url = "https://example.invalid/synthetic-a"
  )),
  claim_boundary = "Synthetic fixture identity only; no biological claim."
)
cfg$followup$candidates[[2]]$term <- "1.1.1.1"
cfg$followup$candidates[[2]]$selection_basis <- "synthetic smoke-test candidate"
cfg$followup$candidates[[2]]$curation <- list(
  accepted_name = "Synthetic enzyme candidate",
  gene_symbol = "synB",
  entity_type = "synthetic enzyme activity",
  reaction = "synthetic substrate = synthetic product",
  authoritative_sources = list(list(
    name = "Synthetic fixture source B",
    url = "https://example.invalid/synthetic-b"
  )),
  claim_boundary = "Synthetic fixture identity only; no biological claim."
)
cfg$followup$candidate_curation$review_date <- "2026-08-12"
cfg$followup$candidate_curation$relationship_statement <-
  "Synthetic co-carriage used only to test workflow contracts"
cfg$followup$candidate_curation$relationship_claim_boundary <-
  "Synthetic fixture overlap does not support biological interpretation."
cfg$followup$minimum_co_carriage_mags <- 1L
cfg$followup$minimum_unknown_group_clusters <- 1L
cfg$outputs$output_dir <- result_dir
cfg$outputs$overwrite <- FALSE
smoke_config <- file.path(output_dir, "mag_functional_annotation_smoke.yml")
yaml::write_yaml(cfg, smoke_config)
cat(smoke_config, "\n", sep = "")
