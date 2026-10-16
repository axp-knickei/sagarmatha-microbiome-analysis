# Public renderer: values come from external inputs; see docs/public-rendering.md.
# ================================================
# create_table_2_hpc_per_mag_lmertest_results.R
# ================================================
# Table 2: top longitudinal MAGs from the full T1-T10 per-MAG lmerTest LMM.
#
# This table ranks MAGs using the existing primary result table. It does not
# refit models or change the original FDR family. The rendered manuscript table
# shows the top 25 ranked MAGs; the CSV/Markdown outputs contain the full ranked
# primary-analysis table.

required_packages <- c("dplyr", "grid")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
results_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "tables",
  "per_mag_lmertest_hpc_results.csv"
)
taxonomy_path <- file.path(repo_root, "data", "GDTB.csv")
output_dir <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "tables")
output_stem <- "table_2_hpc_per_mag_lmertest_results"

# External inputs and isolated output destination (run from repository root).
results_path <- Sys.getenv("PUBLIC_LMM_RESULTS", unset = results_path)
taxonomy_path <- Sys.getenv("PUBLIC_TAXONOMY", unset = taxonomy_path)
output_dir <- Sys.getenv("PUBLIC_RENDER_OUTPUT_DIR", unset = output_dir)

if (!file.exists(results_path)) {
  stop("Missing per-MAG LMM result table: ", results_path, call. = FALSE)
}
if (!file.exists(taxonomy_path)) {
  stop("Missing GTDB-Tk taxonomy table: ", taxonomy_path, call. = FALSE)
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
  x[is.na(x) | x == "" | x == paste0(prefix, "__")] <- "Unclassified"
  x <- sub(paste0("^", prefix, "__"), "", x)
  x <- gsub("_", " ", x)
  x
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, formatC(x, format = "e", digits = 2), sprintf("%.3f", x)))
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits, trim = TRUE))
}

wrap_text <- function(x, width) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

shorten <- function(x, n = 36) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

taxonomy_clean <- taxonomy |>
  dplyr::transmute(
    MAG = Genome,
    phylum = clean_taxon(phylum, "p"),
    class = clean_taxon(class, "c"),
    order = clean_taxon(ordeer, "o"),
    family = clean_taxon(family, "f"),
    genus = clean_taxon(genus, "g"),
    species = clean_taxon(species, "s")
  ) |>
  dplyr::mutate(
    gtdb_taxonomy = paste0(
      "p__", phylum,
      "; c__", class,
      "; o__", order,
      "; f__", family,
      "; g__", genus,
      "; s__", species
    ),
    display_taxonomy = dplyr::case_when(
      genus != "Unclassified" & family != "Unclassified" ~ paste0(genus, " (", family, ")"),
      family != "Unclassified" ~ family,
      order != "Unclassified" ~ order,
      class != "Unclassified" ~ class,
      phylum != "Unclassified" ~ phylum,
      TRUE ~ "Unclassified"
    )
  )

ranked_table <- results |>
  dplyr::filter(primary_analysis %in% TRUE) |>
  dplyr::left_join(taxonomy_clean, by = "MAG") |>
  dplyr::mutate(
    primary_significant = primary_significant %in% TRUE,
    significance = ifelse(primary_significant, "BH-FDR q < 0.05", "Not FDR-significant"),
    direction = dplyr::case_when(
      time_estimate > 0 ~ "Increasing",
      time_estimate < 0 ~ "Decreasing",
      TRUE ~ "No directional estimate"
    ),
    abs_time_estimate = abs(time_estimate),
    display_taxonomy = ifelse(is.na(display_taxonomy), "Unclassified", display_taxonomy),
    gtdb_taxonomy = ifelse(is.na(gtdb_taxonomy), "Unclassified", gtdb_taxonomy)
  ) |>
  dplyr::arrange(
    dplyr::desc(primary_significant),
    time_p_adj_bh_primary,
    dplyr::desc(abs_time_estimate)
  ) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::select(
    rank,
    MAG,
    gtdb_taxonomy,
    display_taxonomy,
    direction,
    time_estimate,
    time_std_error,
    time_p_value,
    time_p_adj_bh_primary,
    primary_significant,
    significance,
    zero_fraction,
    n_subjects,
    n_timepoints,
    status
  )

if (nrow(ranked_table) == 0) {
  stop("No primary-analysis MAG rows were available for Table 2.", call. = FALSE)
}

csv_path <- file.path(output_dir, paste0(output_stem, ".csv"))
md_path <- file.path(output_dir, paste0(output_stem, ".md"))
png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

utils::write.csv(ranked_table, csv_path, row.names = FALSE)

markdown_lines <- c(
  "# Table 2. Top Longitudinal MAGs from Per-MAG LMM",
  "",
  "Rows are sorted by primary significance first, then BH-FDR q-value, then absolute longitudinal slope. This full Markdown table contains all primary-analysis MAGs; the rendered manuscript table shows the top 25 rows.",
  "",
  paste(names(ranked_table), collapse = " | "),
  paste(rep("---", ncol(ranked_table)), collapse = " | ")
)
for (i in seq_len(nrow(ranked_table))) {
  row_values <- vapply(ranked_table[i, ], as.character, character(1))
  row_values <- gsub("\\|", "/", row_values)
  row_values[is.na(row_values)] <- ""
  markdown_lines <- c(markdown_lines, paste(row_values, collapse = " | "))
}
writeLines(markdown_lines, md_path)

display_n <- min(25L, nrow(ranked_table))
display_table <- ranked_table |>
  dplyr::slice_head(n = display_n) |>
  dplyr::transmute(
    rank = as.character(rank),
    MAG = wrap_text(shorten(MAG, 35), 25),
    taxonomy = wrap_text(shorten(display_taxonomy, 31), 24),
    slope = fmt_num(time_estimate, 3),
    SE = fmt_num(time_std_error, 3),
    raw_p = fmt_p(time_p_value),
    q_value = fmt_p(time_p_adj_bh_primary),
    FDR = ifelse(primary_significant, "Yes", "No"),
    zeros = paste0(round(100 * zero_fraction), "%"),
    support = paste0(n_subjects, " subj; ", n_timepoints, " TP"),
    status = status
  )

total_primary <- nrow(ranked_table)
fdr_mags <- sum(ranked_table$primary_significant, na.rm = TRUE)
nominal_mags <- sum(ranked_table$time_p_value < 0.05, na.rm = TRUE)
best_q <- min(ranked_table$time_p_adj_bh_primary, na.rm = TRUE)

col_labels <- c("Rank", "MAG", "GTDB taxonomy", "Slope", "SE", "raw p", "q", "FDR", "Zeros", "Support", "Status")
col_widths <- c(0.045, 0.205, 0.18, 0.07, 0.06, 0.075, 0.065, 0.045, 0.055, 0.09, 0.11)
col_x <- cumsum(c(0, col_widths))

pal <- list(
  ink = "#111827",
  text = "#1F2937",
  muted = "#4B5563",
  hair = "#D1D5DB",
  header = "#E8F1F8",
  stripe = "#F8FAFC",
  blue = "#0072B2",
  orange = "#D55E00",
  red = "#CC3311"
)

draw_table <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  grid::grid.text(
    "Table 2 | Top longitudinal MAGs from per-MAG LMM",
    x = grid::unit(0.018, "npc"),
    y = grid::unit(0.975, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 17, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    paste0(
      "Ranked by primary significance, BH-FDR q-value, and absolute slope. ",
      "Primary MAGs: ", total_primary, "; nominal raw p < 0.05: ", nominal_mags,
      "; BH-FDR q < 0.05: ", fdr_mags, "; best q = ", sprintf("%.3f", best_q), "."
    ),
    x = grid::unit(0.018, "npc"),
    y = grid::unit(0.943, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9.2, col = pal$muted)
  )

  table_top <- 0.895
  table_bottom <- 0.055
  header_h <- 0.044
  row_h <- (table_top - table_bottom - header_h) / nrow(display_table)

  grid::grid.rect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(table_top - header_h / 2, "npc"),
    width = grid::unit(0.964, "npc"),
    height = grid::unit(header_h, "npc"),
    gp = grid::gpar(fill = pal$header, col = pal$hair, lwd = 0.75)
  )

  for (j in seq_along(col_labels)) {
    x_left <- 0.018 + 0.964 * col_x[j]
    grid::grid.text(
      col_labels[j],
      x = grid::unit(x_left + 0.0035, "npc"),
      y = grid::unit(table_top - header_h / 2, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 7.65, fontface = "bold", col = pal$ink)
    )
  }

  for (i in seq_len(nrow(display_table))) {
    y_top <- table_top - header_h - (i - 1) * row_h
    y_mid <- y_top - row_h / 2
    fill <- if (i %% 2 == 0) "#FFFFFF" else pal$stripe
    grid::grid.rect(
      x = grid::unit(0.5, "npc"),
      y = grid::unit(y_mid, "npc"),
      width = grid::unit(0.964, "npc"),
      height = grid::unit(row_h, "npc"),
      gp = grid::gpar(fill = fill, col = "#E5E7EB", lwd = 0.35)
    )

    row_values <- as.character(display_table[i, ])
    for (j in seq_along(row_values)) {
      x_left <- 0.018 + 0.964 * col_x[j]
      color <- pal$text
      if (col_labels[j] == "FDR" && row_values[j] == "Yes") color <- pal$blue
      if (col_labels[j] == "FDR" && row_values[j] == "No") color <- pal$muted
      if (col_labels[j] == "Slope" && suppressWarnings(as.numeric(row_values[j])) < 0) color <- pal$red
      if (col_labels[j] == "Slope" && suppressWarnings(as.numeric(row_values[j])) > 0) color <- pal$blue
      grid::grid.text(
        row_values[j],
        x = grid::unit(x_left + 0.0035, "npc"),
        y = grid::unit(y_mid, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(fontsize = 5.95, col = color, lineheight = 0.82)
      )
    }
  }

  for (j in seq_along(col_x)) {
    x_pos <- 0.018 + 0.964 * col_x[j]
    grid::grid.lines(
      x = grid::unit(c(x_pos, x_pos), "npc"),
      y = grid::unit(c(table_bottom, table_top), "npc"),
      gp = grid::gpar(col = "#E5E7EB", lwd = 0.35)
    )
  }

  grid::grid.text(
    paste0("Slope is the fixed effect for days_from_baseline_subject on log10(Y + pseudocount). Nominal rows are exploratory; BH-FDR discoveries = ", fdr_mags, "."),
    x = grid::unit(0.018, "npc"),
    y = grid::unit(0.022, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.3, col = pal$muted)
  )
}

grDevices::png(filename = png_path, width = 14.6, height = 9.2, units = "in", res = 450, bg = "white", type = "cairo")
draw_table()
grDevices::dev.off()

grDevices::pdf(file = pdf_path, width = 14.6, height = 9.2, bg = "white", useDingbats = FALSE)
draw_table()
grDevices::dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = svg_path, width = 14.6, height = 9.2, bg = "white")
  draw_table()
  grDevices::dev.off()
} else {
  warning("Package 'svglite' is not installed; SVG output was skipped.", call. = FALSE)
}

message("Table 2 saved:")
message("- ", csv_path)
message("- ", md_path)
message("- ", png_path)
message("- ", pdf_path)
if (file.exists(svg_path)) {
  message("- ", svg_path)
}
