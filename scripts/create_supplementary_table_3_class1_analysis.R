# ============================================
# create_supplementary_table_3_class1_analysis.R
# ============================================
# Supplementary Table 3: Top nominal MAG-window candidates, not FDR-significant.
# Ranks nominal candidates by global Class 1 q-value, then raw p-value, then
# non-zero sample support to avoid a raw-p-value-only presentation.

required_packages <- c("readr", "dplyr", "grid")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
results_path <- file.path(
  repo_root,
  "results",
  "class1_window_lmertest",
  "tables",
  "class1_window_nominal_candidates_raw_p_lt_0_05.csv"
)
output_dir <- file.path(repo_root, "results", "class1_window_lmertest", "tables")
figure_dir <- file.path(repo_root, "results", "class1_window_lmertest", "figures")

if (!file.exists(results_path)) {
  stop("Missing nominal candidate table: ", results_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

class1_nominal <- readr::read_csv(results_path, show_col_types = FALSE)

required_columns <- c(
  "MAG",
  "window_label",
  "forward_timepoint_path",
  "estimate_window_order",
  "direction",
  "raw_p_value",
  "q_global_class1_bh",
  "zero_fraction_window",
  "nonzero_samples_window",
  "subjects_with_nonzero_window",
  "status",
  "singular_fit"
)
missing_columns <- setdiff(required_columns, names(class1_nominal))
if (length(missing_columns) > 0) {
  stop("Nominal candidate table is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

# Taxonomy is not yet present in the current Class 1 outputs, so retain a
# transparent placeholder field rather than inventing an annotation.
taxonomy_columns <- grep(
  "taxonomy|taxon|phylum|family|genus|species|^class$|^order$",
  names(class1_nominal),
  ignore.case = TRUE,
  value = TRUE
)
has_taxonomy <- length(taxonomy_columns) > 0

pick_taxonomy <- function(df) {
  if (!has_taxonomy) {
    return(rep("not_available_in_current_output", nrow(df)))
  }
  taxonomy_values <- as.character(df[[taxonomy_columns[[1]]]])
  taxonomy_values[is.na(taxonomy_values) | !nzchar(taxonomy_values)] <- "not_available_in_current_output"
  taxonomy_values
}

top_n <- 30L

supplementary_table <- class1_nominal |>
  dplyr::mutate(taxonomy = pick_taxonomy(class1_nominal)) |>
  dplyr::arrange(
    q_global_class1_bh,
    raw_p_value,
    dplyr::desc(nonzero_samples_window),
    dplyr::desc(subjects_with_nonzero_window)
  ) |>
  dplyr::transmute(
    MAG = MAG,
    taxonomy = taxonomy,
    window_label = window_label,
    timepoint_path = forward_timepoint_path,
    estimate_window_order = estimate_window_order,
    direction = direction,
    raw_p_value = raw_p_value,
    q_global_class1_bh = q_global_class1_bh,
    zero_fraction_window = zero_fraction_window,
    nonzero_samples_window = nonzero_samples_window,
    subjects_with_nonzero_window = subjects_with_nonzero_window,
    status = status,
    singular_fit = singular_fit
  ) |>
  dplyr::slice_head(n = top_n)

fmt_num <- function(x, digits = 4) ifelse(is.na(x), NA_character_, sprintf(paste0("%.", digits, "f"), x))
fmt_p <- function(x) ifelse(is.na(x), NA_character_, formatC(x, format = "e", digits = 3))
fmt_bool <- function(x) ifelse(is.na(x), "", ifelse(x, "TRUE", "FALSE"))

display_table <- supplementary_table |>
  dplyr::mutate(
    estimate_window_order = fmt_num(estimate_window_order, 4),
    raw_p_value = fmt_p(raw_p_value),
    q_global_class1_bh = fmt_num(q_global_class1_bh, 4),
    zero_fraction_window = fmt_num(zero_fraction_window, 3),
    singular_fit = fmt_bool(singular_fit)
  )

csv_path <- file.path(output_dir, "supplementary_table_3_top_nominal_mag_window_candidates.csv")
md_path <- file.path(output_dir, "supplementary_table_3_top_nominal_mag_window_candidates.md")
png_path <- file.path(figure_dir, "supplementary_table_3_top_nominal_mag_window_candidates.png")
pdf_path <- file.path(figure_dir, "supplementary_table_3_top_nominal_mag_window_candidates.pdf")
svg_path <- file.path(figure_dir, "supplementary_table_3_top_nominal_mag_window_candidates.svg")

readr::write_csv(supplementary_table, csv_path)

markdown_lines <- c(
  "# Supplementary Table 3. Top nominal MAG-window candidates, not FDR-significant",
  "",
  "Nominal candidates are ranked by q_global_class1_bh, then raw_p_value, then nonzero_samples_window. No candidate in this table passed global Class 1 BH-FDR.",
  "",
  paste(names(display_table), collapse = " | "),
  paste(rep("---", ncol(display_table)), collapse = " | ")
)
for (i in seq_len(nrow(display_table))) {
  row_values <- vapply(display_table[i, ], as.character, character(1))
  row_values <- gsub("\\|", "/", row_values)
  row_values[is.na(row_values)] <- ""
  markdown_lines <- c(markdown_lines, paste(row_values, collapse = " | "))
}
writeLines(markdown_lines, md_path)

wrap_text <- function(x, width) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

pretty_table <- display_table |>
  dplyr::mutate(
    MAG = wrap_text(MAG, 23),
    taxonomy = wrap_text(taxonomy, 24),
    window_label = wrap_text(window_label, 19),
    timepoint_path = wrap_text(timepoint_path, 12),
    direction = wrap_text(direction, 11),
    status = wrap_text(status, 11)
  )

col_labels <- c(
  "MAG",
  "taxonomy",
  "window_label",
  "timepoint\npath",
  "estimate\nwindow_order",
  "direction",
  "raw\np-value",
  "global\nq",
  "zero\nfraction",
  "nonzero\nsamples",
  "subjects with\nnonzero",
  "status",
  "singular\nfit"
)
col_widths <- c(0.16, 0.12, 0.11, 0.07, 0.08, 0.06, 0.07, 0.06, 0.06, 0.06, 0.07, 0.05, 0.03)
col_x <- cumsum(c(0, col_widths / sum(col_widths)))

draw_table <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text(
    "Supplementary Table 3 | Top nominal MAG-window candidates, not FDR-significant",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.975, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 16, fontface = "bold", col = "#111827")
  )
  grid::grid.text(
    paste0(
      "Rows are ranked by global q-value, then raw p-value, then non-zero support. ",
      "Taxonomy is shown only if present in current outputs; otherwise a transparent placeholder is retained."
    ),
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.94, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9.1, col = "#4b5563")
  )

  table_top <- 0.90
  table_bottom <- 0.055
  header_h <- 0.06
  row_h <- (table_top - table_bottom - header_h) / nrow(pretty_table)

  grid::grid.rect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(table_top - header_h / 2, "npc"),
    width = grid::unit(0.95, "npc"),
    height = grid::unit(header_h, "npc"),
    gp = grid::gpar(fill = "#f3f4f6", col = "#d1d5db", lwd = 0.75)
  )

  for (j in seq_along(col_labels)) {
    x_left <- 0.025 + 0.95 * col_x[j]
    x_right <- 0.025 + 0.95 * col_x[j + 1]
    grid::grid.text(
      col_labels[j],
      x = grid::unit((x_left + x_right) / 2, "npc"),
      y = grid::unit(table_top - header_h / 2, "npc"),
      just = "center",
      gp = grid::gpar(fontsize = 7.1, fontface = "bold", col = "#111827", lineheight = 0.84)
    )
  }

  for (i in seq_len(nrow(pretty_table))) {
    y_top <- table_top - header_h - (i - 1) * row_h
    y_mid <- y_top - row_h / 2
    fill <- if (i %% 2 == 0) "#ffffff" else "#f8fafc"
    grid::grid.rect(
      x = grid::unit(0.5, "npc"),
      y = grid::unit(y_mid, "npc"),
      width = grid::unit(0.95, "npc"),
      height = grid::unit(row_h, "npc"),
      gp = grid::gpar(fill = fill, col = "#e5e7eb", lwd = 0.4)
    )
    for (j in seq_along(pretty_table)) {
      x_left <- 0.025 + 0.95 * col_x[j]
      x_right <- 0.025 + 0.95 * col_x[j + 1]
      value <- as.character(pretty_table[i, j, drop = TRUE])
      left_justified <- j %in% c(1, 2, 3, 4, 6, 12)
      just <- if (left_justified) c("left", "center") else "center"
      x_pos <- if (left_justified) x_left + 0.003 else (x_left + x_right) / 2
      grid::grid.text(
        value,
        x = grid::unit(x_pos, "npc"),
        y = grid::unit(y_mid, "npc"),
        just = just,
        gp = grid::gpar(fontsize = 5.5, col = "#1f2937", lineheight = 0.84)
      )
    }
  }

  for (j in seq_along(col_x)) {
    x_pos <- 0.025 + 0.95 * col_x[j]
    grid::grid.lines(
      x = grid::unit(c(x_pos, x_pos), "npc"),
      y = grid::unit(c(table_bottom, table_top), "npc"),
      gp = grid::gpar(col = "#e5e7eb", lwd = 0.4)
    )
  }

  grid::grid.text(
    "Interpretation: these are nominal stage-local candidate signals only; none passed global Class 1 BH-FDR.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.02, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.8, col = "#4b5563")
  )
}

png(filename = png_path, width = 16, height = 12, units = "in", res = 360, bg = "white")
draw_table()
dev.off()

pdf(file = pdf_path, width = 16, height = 12, bg = "white")
draw_table()
dev.off()

svglite::svglite(file = svg_path, width = 16, height = 12, bg = "white")
draw_table()
dev.off()

message("Supplementary Table 3 nominal candidate outputs created:")
message(" - ", csv_path)
message(" - ", md_path)
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Input:")
message(" - ", results_path)
message("Rows: ", nrow(supplementary_table))
message("Taxonomy available: ", has_taxonomy)
message("Best q_global_class1_bh in table: ", min(supplementary_table$q_global_class1_bh, na.rm = TRUE))
