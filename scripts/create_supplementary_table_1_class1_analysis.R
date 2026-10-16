# ============================================
# create_supplementary_table_1_class1_analysis.R
# ============================================
# Supplementary Table 1: Class 1 window definitions and evidence type.
# This table documents the biological justification for predefined windows.

required_packages <- c("readr", "dplyr", "grid")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
windows_path <- file.path(
  repo_root,
  "results",
  "metadata_expedition_windows",
  "proposed_biological_windows.csv"
)
output_dir <- file.path(repo_root, "results", "class1_window_lmertest", "tables")
figure_dir <- file.path(repo_root, "results", "class1_window_lmertest", "figures")

if (!file.exists(windows_path)) {
  stop("Missing window definition table: ", windows_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

windows <- readr::read_csv(windows_path, show_col_types = FALSE)

required_columns <- c(
  "window_label",
  "forward_timepoint_path",
  "biological_interpretation",
  "sampling_altitude_span",
  "figure_informed_altitude_exposure_interpretation",
  "evidence_type",
  "n_complete_subjects_all_three_timepoints",
  "caveats"
)
missing_columns <- setdiff(required_columns, names(windows))
if (length(missing_columns) > 0) {
  stop("Window table is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

main_caveat <- function(caveats) {
  caveats <- as.character(caveats)
  caveats[is.na(caveats) | !nzchar(caveats)] <- "standard longitudinal caveats"
  vapply(strsplit(caveats, ";", fixed = TRUE), function(parts) trimws(parts[[1]]), character(1))
}

clean_stage <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  x
}

supplementary_table <- windows |>
  dplyr::transmute(
    window_label = window_label,
    timepoint_path = forward_timepoint_path,
    expedition_stage = biological_interpretation,
    sampling_altitude_span = sampling_altitude_span,
    figure_informed_exposure_context = figure_informed_altitude_exposure_interpretation,
    evidence_type = evidence_type,
    complete_subjects = n_complete_subjects_all_three_timepoints,
    main_caveat = main_caveat(caveats)
  )

csv_path <- file.path(output_dir, "supplementary_table_1_class1_window_definitions.csv")
md_path <- file.path(output_dir, "supplementary_table_1_class1_window_definitions.md")
png_path <- file.path(figure_dir, "supplementary_table_1_class1_window_definitions.png")
pdf_path <- file.path(figure_dir, "supplementary_table_1_class1_window_definitions.pdf")
svg_path <- file.path(figure_dir, "supplementary_table_1_class1_window_definitions.svg")

readr::write_csv(supplementary_table, csv_path)

markdown_lines <- c(
  "# Supplementary Table 1. Class 1 Window Definitions and Evidence Type",
  "",
  "This table documents the predefined biological windows used for the Class 1 secondary LMM analysis. Windows were defined before MAG p-value interpretation using sampling metadata and figure-informed expedition exposure context.",
  "",
  paste(names(supplementary_table), collapse = " | "),
  paste(rep("---", ncol(supplementary_table)), collapse = " | ")
)
for (i in seq_len(nrow(supplementary_table))) {
  row_values <- vapply(supplementary_table[i, ], as.character, character(1))
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

display_table <- supplementary_table |>
  dplyr::mutate(
    window_label = wrap_text(window_label, 18),
    timepoint_path = wrap_text(timepoint_path, 14),
    expedition_stage = wrap_text(expedition_stage, 27),
    sampling_altitude_span = wrap_text(sampling_altitude_span, 12),
    figure_informed_exposure_context = wrap_text(figure_informed_exposure_context, 32),
    evidence_type = wrap_text(evidence_type, 22),
    complete_subjects = as.character(complete_subjects),
    main_caveat = wrap_text(main_caveat, 25)
  )

col_labels <- c(
  "Window",
  "Timepoint\npath",
  "Expedition\nstage",
  "Sampling\naltitude span",
  "Figure-informed exposure context",
  "Evidence type",
  "Complete\nsubjects",
  "Main caveat"
)
col_widths <- c(0.13, 0.09, 0.17, 0.09, 0.23, 0.13, 0.06, 0.10)
col_x <- cumsum(c(0, col_widths))

draw_table <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text(
    "Supplementary Table 1 | Class 1 window definitions and evidence type",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.965, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 16, fontface = "bold", col = "#111827")
  )
  grid::grid.text(
    "Predefined windows were justified using expedition chronology, sampling altitude, and figure-informed exposure trajectory before MAG p-value interpretation.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.925, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9.4, col = "#4b5563")
  )

  table_top <- 0.865
  table_bottom <- 0.045
  header_h <- 0.075
  row_h <- (table_top - table_bottom - header_h) / nrow(display_table)

  grid::grid.rect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(table_top - header_h / 2, "npc"),
    width = grid::unit(0.95, "npc"),
    height = grid::unit(header_h, "npc"),
    gp = grid::gpar(fill = "#f3f4f6", col = "#d1d5db", lwd = 0.8)
  )

  for (j in seq_along(col_labels)) {
    x_left <- 0.025 + 0.95 * col_x[j]
    x_right <- 0.025 + 0.95 * col_x[j + 1]
    grid::grid.text(
      col_labels[j],
      x = grid::unit((x_left + x_right) / 2, "npc"),
      y = grid::unit(table_top - header_h / 2, "npc"),
      just = "center",
      gp = grid::gpar(fontsize = 7.9, fontface = "bold", col = "#111827", lineheight = 0.88)
    )
  }

  for (i in seq_len(nrow(display_table))) {
    y_top <- table_top - header_h - (i - 1) * row_h
    y_mid <- y_top - row_h / 2
    fill <- if (i %% 2 == 0) "#ffffff" else "#f8fafc"
    grid::grid.rect(
      x = grid::unit(0.5, "npc"),
      y = grid::unit(y_mid, "npc"),
      width = grid::unit(0.95, "npc"),
      height = grid::unit(row_h, "npc"),
      gp = grid::gpar(fill = fill, col = "#e5e7eb", lwd = 0.45)
    )
    for (j in seq_along(display_table)) {
      x_left <- 0.025 + 0.95 * col_x[j]
      x_right <- 0.025 + 0.95 * col_x[j + 1]
      value <- as.character(display_table[i, j, drop = TRUE])
      grid::grid.text(
        value,
        x = grid::unit(x_left + 0.004, "npc"),
        y = grid::unit(y_mid, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(fontsize = 6.45, col = "#1f2937", lineheight = 0.84)
      )
    }
  }

  for (j in seq_along(col_x)) {
    x_pos <- 0.025 + 0.95 * col_x[j]
    grid::grid.lines(
      x = grid::unit(c(x_pos, x_pos), "npc"),
      y = grid::unit(c(table_bottom, table_top), "npc"),
      gp = grid::gpar(col = "#e5e7eb", lwd = 0.45)
    )
  }
  grid::grid.text(
    "Critical evidence labels: T6-T7-T8 = figure_informed_critical_interval; T7-T8-T9 = figure_informed_plus_sampling_descent.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.018, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.6, col = "#4b5563")
  )
}

png(filename = png_path, width = 13.6, height = 8.0, units = "in", res = 360, bg = "white")
draw_table()
dev.off()

pdf(file = pdf_path, width = 13.6, height = 8.0, bg = "white")
draw_table()
dev.off()

svglite::svglite(file = svg_path, width = 13.6, height = 8.0, bg = "white")
draw_table()
dev.off()

message("Supplementary Table 1 Class 1 window definition outputs created:")
message(" - ", csv_path)
message(" - ", md_path)
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Input:")
message(" - ", windows_path)
message("Rows: ", nrow(supplementary_table))
message("Critical evidence types:")
message(" - T6 -> T7 -> T8: ", supplementary_table$evidence_type[supplementary_table$timepoint_path == "T6 -> T7 -> T8"])
message(" - T7 -> T8 -> T9: ", supplementary_table$evidence_type[supplementary_table$timepoint_path == "T7 -> T8 -> T9"])
