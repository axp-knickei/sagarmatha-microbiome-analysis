# Public renderer: values come from external inputs; see docs/public-rendering.md.
# Create manuscript Table 2: top longitudinal MAGs from per-MAG LMM results.
#
# Ranking rule:
#   1. primary_significant first
#   2. lower BH-FDR q-value
#   3. larger absolute longitudinal slope
#
# The script does not refit models. It formats the existing per-MAG lmerTest
# result table and joins GTDB taxonomy for biological interpretation.

required_packages <- c("dplyr", "grid", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

library(dplyr)

repo_root <- getwd()
results_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "tables",
  "per_mag_lmertest_hpc_results.csv"
)
taxonomy_path <- file.path(repo_root, "data", "raw", "GDTB.csv")
output_dir <- file.path(repo_root, "results", "manuscript_tables")
output_stem <- "table_2_top_longitudinal_mags"
display_n <- as.integer(Sys.getenv("TABLE2_TOP_N", unset = "25"))

# External inputs and isolated output destination (run from repository root).
results_path <- Sys.getenv("PUBLIC_LMM_RESULTS", unset = results_path)
taxonomy_path <- Sys.getenv("PUBLIC_TAXONOMY", unset = taxonomy_path)
output_dir <- Sys.getenv("PUBLIC_RENDER_OUTPUT_DIR", unset = output_dir)

required_files <- c(results_path, taxonomy_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
}
if (is.na(display_n) || display_n < 1) {
  stop("TABLE2_TOP_N must be a positive integer.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
taxonomy <- utils::read.csv(taxonomy_path, check.names = FALSE, stringsAsFactors = FALSE)

required_result_columns <- c(
  "MAG", "time_estimate", "time_std_error", "time_p_value",
  "time_p_adj_bh_primary", "primary_significant", "zero_fraction",
  "n_subjects", "n_timepoints", "status", "primary_analysis"
)
required_taxonomy_columns <- c("Genome", "phylum", "class", "ordeer", "family", "genus", "species")
missing_result_columns <- setdiff(required_result_columns, names(results))
missing_taxonomy_columns <- setdiff(required_taxonomy_columns, names(taxonomy))
if (length(missing_result_columns) > 0) {
  stop("Missing required result columns: ", paste(missing_result_columns, collapse = ", "), call. = FALSE)
}
if (length(missing_taxonomy_columns) > 0) {
  stop("Missing required taxonomy columns: ", paste(missing_taxonomy_columns, collapse = ", "), call. = FALSE)
}

clean_taxon <- function(x, prefix) {
  x <- as.character(x)
  empty_token <- paste0(prefix, "__")
  x[is.na(x) | x == "" | x == empty_token] <- "Unclassified"
  x <- sub(paste0("^", prefix, "__"), "", x)
  x <- gsub("_", " ", x)
  x
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits, trim = TRUE))
}

fmt_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 2), sprintf("%.3f", x))
  )
}

fmt_pct <- function(x, digits = 0) {
  ifelse(is.na(x), "NA", paste0(format(round(100 * x, digits), nsmall = digits, trim = TRUE), "%"))
}

wrap_text <- function(x, width) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

shorten <- function(x, n = 42) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

taxonomy_clean <- taxonomy %>%
  transmute(
    MAG = Genome,
    phylum = clean_taxon(phylum, "p"),
    class = clean_taxon(class, "c"),
    order = clean_taxon(ordeer, "o"),
    family = clean_taxon(family, "f"),
    genus = clean_taxon(genus, "g"),
    species = clean_taxon(species, "s")
  ) %>%
  mutate(
    GTDB_taxonomy = paste0(
      "p__", phylum,
      "; c__", class,
      "; o__", order,
      "; f__", family,
      "; g__", genus,
      "; s__", species
    ),
    visual_taxonomy = case_when(
      genus != "Unclassified" & species != "Unclassified" ~ paste0(genus, "; ", species),
      genus != "Unclassified" & family != "Unclassified" ~ paste0(genus, " (", family, ")"),
      family != "Unclassified" ~ family,
      order != "Unclassified" ~ order,
      class != "Unclassified" ~ class,
      phylum != "Unclassified" ~ phylum,
      TRUE ~ "Unclassified"
    )
  )

ranked <- results %>%
  filter(primary_analysis %in% TRUE) %>%
  left_join(taxonomy_clean, by = "MAG") %>%
  mutate(
    primary_significant = primary_significant %in% TRUE,
    abs_time_estimate = abs(time_estimate),
    GTDB_taxonomy = ifelse(is.na(GTDB_taxonomy), "Unclassified", GTDB_taxonomy),
    visual_taxonomy = ifelse(is.na(visual_taxonomy), "Unclassified", visual_taxonomy)
  ) %>%
  arrange(
    desc(primary_significant),
    time_p_adj_bh_primary,
    desc(abs_time_estimate)
  ) %>%
  mutate(rank = row_number())

if (nrow(ranked) == 0) {
  stop("No primary-analysis MAG rows were available in the per-MAG LMM result table.", call. = FALSE)
}

table_2 <- ranked %>%
  slice_head(n = display_n) %>%
  transmute(
    rank,
    MAG,
    GTDB_taxonomy,
    time_estimate,
    time_std_error,
    time_p_value,
    time_p_adj_bh_primary,
    primary_significant,
    zero_fraction,
    n_subjects,
    n_timepoints,
    status
  )

csv_path <- file.path(output_dir, paste0(output_stem, ".csv"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))
png_path <- file.path(output_dir, paste0(output_stem, ".png"))
utils::write.csv(table_2, csv_path, row.names = FALSE)

fdr_count <- sum(ranked$primary_significant, na.rm = TRUE)
nominal_count <- sum(ranked$time_p_value < 0.05, na.rm = TRUE)
best_q <- min(ranked$time_p_adj_bh_primary, na.rm = TRUE)

display_table <- ranked %>%
  slice_head(n = display_n) %>%
  transmute(
    rank = as.character(rank),
    MAG = wrap_text(shorten(MAG, 34), 24),
    taxonomy = wrap_text(shorten(visual_taxonomy, 38), 28),
    slope = fmt_num(time_estimate, 3),
    SE = fmt_num(time_std_error, 3),
    raw_p = fmt_p(time_p_value),
    q_value = fmt_p(time_p_adj_bh_primary),
    FDR = ifelse(primary_significant, "Yes", "No"),
    zero_fraction = fmt_pct(zero_fraction, 0),
    support = paste0(n_subjects, " subj; ", n_timepoints, " TP"),
    status = status
  )

col_labels <- c("Rank", "MAG", "GTDB taxonomy", "Slope", "SE", "raw p", "q", "FDR", "Zeros", "Support", "Status")
col_widths <- c(0.045, 0.195, 0.190, 0.070, 0.060, 0.070, 0.060, 0.045, 0.055, 0.090, 0.120)
col_x <- cumsum(c(0, col_widths))

pal <- list(
  ink = "#111827",
  text = "#1F2937",
  muted = "#4B5563",
  hair = "#D1D5DB",
  header = "#E8F1F8",
  stripe = "#F8FAFC",
  blue = "#0072B2",
  orange = "#D55E00"
)

draw_table <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  grid::grid.text(
    "Table 2 | Top longitudinal MAGs from per-MAG LMM",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.965, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 17, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    paste0(
      "Top ", min(display_n, nrow(ranked)),
      " primary-analysis MAGs ranked by FDR significance, q-value, and absolute slope. ",
      "FDR-significant MAGs: ", fdr_count,
      "; nominal raw p < 0.05: ", nominal_count,
      "; minimum q = ", fmt_p(best_q), "."
    ),
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.925, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9.2, col = pal$muted)
  )

  table_top <- 0.865
  table_bottom <- 0.080
  header_h <- 0.052
  row_h <- (table_top - table_bottom - header_h) / nrow(display_table)

  grid::grid.rect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(table_top - header_h / 2, "npc"),
    width = grid::unit(0.95, "npc"),
    height = grid::unit(header_h, "npc"),
    gp = grid::gpar(fill = pal$header, col = pal$hair, lwd = 0.8)
  )

  for (j in seq_along(col_labels)) {
    x_left <- 0.025 + 0.95 * col_x[j]
    grid::grid.text(
      col_labels[j],
      x = grid::unit(x_left + 0.0045, "npc"),
      y = grid::unit(table_top - header_h / 2, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 7.6, fontface = "bold", col = pal$ink)
    )
  }

  for (i in seq_len(nrow(display_table))) {
    y_top <- table_top - header_h - (i - 1) * row_h
    y_mid <- y_top - row_h / 2
    fill <- if (i %% 2 == 0) "#FFFFFF" else pal$stripe
    is_fdr <- display_table$FDR[[i]] == "Yes"

    grid::grid.rect(
      x = grid::unit(0.5, "npc"),
      y = grid::unit(y_mid, "npc"),
      width = grid::unit(0.95, "npc"),
      height = grid::unit(row_h, "npc"),
      gp = grid::gpar(fill = fill, col = "#E5E7EB", lwd = 0.42)
    )

    row_values <- as.character(display_table[i, ])
    for (j in seq_along(row_values)) {
      x_left <- 0.025 + 0.95 * col_x[j]
      color <- if (is_fdr && j %in% c(1, 7, 8)) pal$orange else pal$text
      fontface <- if (is_fdr && j %in% c(1, 7, 8)) "bold" else "plain"
      grid::grid.text(
        row_values[j],
        x = grid::unit(x_left + 0.0045, "npc"),
        y = grid::unit(y_mid, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(fontsize = 6.25, col = color, fontface = fontface, lineheight = 0.86)
      )
    }
  }

  for (j in seq_along(col_x)) {
    x_pos <- 0.025 + 0.95 * col_x[j]
    grid::grid.lines(
      x = grid::unit(c(x_pos, x_pos), "npc"),
      y = grid::unit(c(table_bottom, table_top), "npc"),
      gp = grid::gpar(col = "#E5E7EB", lwd = 0.4)
    )
  }

  grid::grid.text(
    paste0(
      "Notes: Slope is the fixed effect for days from subject baseline from the per-MAG LMM. ",
      "q is Benjamini-Hochberg FDR within the primary MAG family. ",
      "Zeros is the sample-level zero fraction. No rows are highlighted when no MAG passes q < 0.05."
    ),
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.045, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.1, col = pal$muted)
  )
  grid::grid.text(
    "Abbreviations: FDR, false discovery rate; GTDB, Genome Taxonomy Database; LMM, linear mixed model; MAG, metagenome-assembled genome; SE, standard error; TP, timepoints.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.025, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.1, col = pal$muted)
  )
}

grDevices::png(filename = png_path, width = 13.8, height = 9.0, units = "in", res = 450, bg = "white", type = "cairo")
draw_table()
grDevices::dev.off()

svglite::svglite(file = svg_path, width = 13.8, height = 9.0, bg = "white")
draw_table()
grDevices::dev.off()

message("Table 2 saved:")
message("- ", csv_path)
message("- ", svg_path)
message("- ", png_path)
