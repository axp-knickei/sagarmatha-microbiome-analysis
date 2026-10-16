# ============================================
# create_supplementary_table_2_class1_analysis.R
# ============================================
# Supplementary Table 2: per-window Class 1 model-status summary.
# This table supports the conclusion that the workflow was usable, while no
# MAG-window test passed global Class 1 BH-FDR correction.

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
  "class1_window_lmertest_results.csv"
)
windows_path <- file.path(
  repo_root,
  "results",
  "metadata_expedition_windows",
  "proposed_biological_windows.csv"
)
output_dir <- file.path(repo_root, "results", "class1_window_lmertest", "tables")
figure_dir <- file.path(repo_root, "results", "class1_window_lmertest", "figures")

if (!file.exists(results_path)) {
  stop("Missing Class 1 results table: ", results_path, call. = FALSE)
}
if (!file.exists(windows_path)) {
  stop("Missing window definition table: ", windows_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

class1_results <- readr::read_csv(results_path, show_col_types = FALSE)
windows <- readr::read_csv(windows_path, show_col_types = FALSE)

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
fmt_q <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.4f", x))

window_order <- windows$window_label
status_levels <- c("lmm_ok", "lmm_singular", "too_sparse", "low_variance", "model_failed")

status_counts <- class1_results |>
  dplyr::mutate(
    window_label = factor(window_label, levels = window_order),
    status = factor(status, levels = status_levels)
  ) |>
  dplyr::count(window_label, status, name = "n")

all_statuses <- expand.grid(
  window_label = factor(window_order, levels = window_order),
  status = factor(status_levels, levels = status_levels),
  KEEP.OUT.ATTRS = FALSE
)

status_counts <- all_statuses |>
  dplyr::left_join(status_counts, by = c("window_label", "status")) |>
  dplyr::mutate(n = dplyr::coalesce(n, 0L))

status_wide <- status_counts |>
  dplyr::group_by(window_label) |>
  dplyr::summarise(
    lmm_ok = sum(n[status == "lmm_ok"]),
    lmm_singular = sum(n[status == "lmm_singular"]),
    too_sparse = sum(n[status == "too_sparse"]),
    low_variance = sum(n[status == "low_variance"]),
    model_failed = sum(n[status == "model_failed"]),
    .groups = "drop"
  )

significance_summary <- class1_results |>
  dplyr::mutate(window_label = factor(window_label, levels = window_order)) |>
  dplyr::group_by(window_label) |>
  dplyr::summarise(
    total_MAGs = dplyr::n_distinct(MAG),
    nominal_raw_p_lt_0_05 = sum(primary_analysis %in% TRUE & raw_p_value < 0.05, na.rm = TRUE),
    q_global_class1_bh_lt_0_05 = sum(primary_analysis %in% TRUE & q_global_class1_bh < 0.05, na.rm = TRUE),
    best_q_global_class1_bh = suppressWarnings(min(q_global_class1_bh[primary_analysis %in% TRUE], na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    best_q_global_class1_bh = ifelse(is.infinite(best_q_global_class1_bh), NA_real_, best_q_global_class1_bh)
  )

supplementary_table <- windows |>
  dplyr::select(window_label, forward_timepoint_path) |>
  dplyr::mutate(window_label = factor(window_label, levels = window_order)) |>
  dplyr::left_join(significance_summary, by = "window_label") |>
  dplyr::left_join(status_wide, by = "window_label") |>
  dplyr::mutate(
    percent_lmm_ok = lmm_ok / total_MAGs,
    window_label = as.character(window_label)
  ) |>
  dplyr::select(
    window_label,
    total_MAGs,
    lmm_ok,
    lmm_singular,
    too_sparse,
    low_variance,
    model_failed,
    percent_lmm_ok,
    nominal_raw_p_lt_0_05,
    q_global_class1_bh_lt_0_05,
    best_q_global_class1_bh
  )

csv_path <- file.path(output_dir, "supplementary_table_2_class1_per_window_model_status.csv")
md_path <- file.path(output_dir, "supplementary_table_2_class1_per_window_model_status.md")
png_path <- file.path(figure_dir, "supplementary_table_2_class1_per_window_model_status.png")
pdf_path <- file.path(figure_dir, "supplementary_table_2_class1_per_window_model_status.pdf")
svg_path <- file.path(figure_dir, "supplementary_table_2_class1_per_window_model_status.svg")

readr::write_csv(supplementary_table, csv_path)

display_table <- supplementary_table |>
  dplyr::mutate(
    percent_lmm_ok = fmt_pct(percent_lmm_ok),
    best_q_global_class1_bh = fmt_q(best_q_global_class1_bh)
  )

markdown_lines <- c(
  "# Supplementary Table 2. Per-Window Class 1 Model Status Summary",
  "",
  "This table summarizes model usability and global Class 1 BH-FDR results for each predefined biological window.",
  "",
  paste(names(display_table), collapse = " | "),
  paste(rep("---", ncol(display_table)), collapse = " | ")
)
for (i in seq_len(nrow(display_table))) {
  row_values <- vapply(display_table[i, ], as.character, character(1))
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
  dplyr::mutate(window_label = wrap_text(window_label, 20))

col_labels <- c(
  "window_label",
  "total\nMAGs",
  "lmm_ok",
  "lmm\nsingular",
  "too\nsparse",
  "low\nvariance",
  "model\nfailed",
  "%\nlmm_ok",
  "nominal\nraw p<0.05",
  "global\nq<0.05",
  "best\nglobal q"
)
col_widths <- c(0.22, 0.075, 0.075, 0.08, 0.08, 0.08, 0.08, 0.085, 0.095, 0.075, 0.075)
col_x <- cumsum(c(0, col_widths / sum(col_widths)))

draw_table <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text(
    "Supplementary Table 2 | Per-window Class 1 model status summary",
    x = grid::unit(0.035, "npc"),
    y = grid::unit(0.94, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 16, fontface = "bold", col = "#111827")
  )
  grid::grid.text(
    "The workflow was statistically usable across predefined windows, but no MAG-window test passed global Class 1 BH-FDR.",
    x = grid::unit(0.035, "npc"),
    y = grid::unit(0.895, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9.4, col = "#4b5563")
  )

  table_top <- 0.815
  table_bottom <- 0.14
  header_h <- 0.11
  row_h <- (table_top - table_bottom - header_h) / nrow(pretty_table)

  grid::grid.rect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(table_top - header_h / 2, "npc"),
    width = grid::unit(0.93, "npc"),
    height = grid::unit(header_h, "npc"),
    gp = grid::gpar(fill = "#f3f4f6", col = "#d1d5db", lwd = 0.8)
  )

  for (j in seq_along(col_labels)) {
    x_left <- 0.035 + 0.93 * col_x[j]
    x_right <- 0.035 + 0.93 * col_x[j + 1]
    grid::grid.text(
      col_labels[j],
      x = grid::unit((x_left + x_right) / 2, "npc"),
      y = grid::unit(table_top - header_h / 2, "npc"),
      just = "center",
      gp = grid::gpar(fontsize = 7.6, fontface = "bold", col = "#111827", lineheight = 0.86)
    )
  }

  for (i in seq_len(nrow(pretty_table))) {
    y_top <- table_top - header_h - (i - 1) * row_h
    y_mid <- y_top - row_h / 2
    fill <- if (i %% 2 == 0) "#ffffff" else "#f8fafc"
    grid::grid.rect(
      x = grid::unit(0.5, "npc"),
      y = grid::unit(y_mid, "npc"),
      width = grid::unit(0.93, "npc"),
      height = grid::unit(row_h, "npc"),
      gp = grid::gpar(fill = fill, col = "#e5e7eb", lwd = 0.45)
    )
    for (j in seq_along(pretty_table)) {
      x_left <- 0.035 + 0.93 * col_x[j]
      x_right <- 0.035 + 0.93 * col_x[j + 1]
      value <- as.character(pretty_table[i, j, drop = TRUE])
      just <- if (j == 1) c("left", "center") else "center"
      x_pos <- if (j == 1) x_left + 0.004 else (x_left + x_right) / 2
      grid::grid.text(
        value,
        x = grid::unit(x_pos, "npc"),
        y = grid::unit(y_mid, "npc"),
        just = just,
        gp = grid::gpar(fontsize = 7.3, col = "#1f2937", lineheight = 0.86)
      )
    }
  }

  for (j in seq_along(col_x)) {
    x_pos <- 0.035 + 0.93 * col_x[j]
    grid::grid.lines(
      x = grid::unit(c(x_pos, x_pos), "npc"),
      y = grid::unit(c(table_bottom, table_top), "npc"),
      gp = grid::gpar(col = "#e5e7eb", lwd = 0.45)
    )
  }

  grid::grid.text(
    "Interpretation: lmm_ok rows are eligible for primary Class 1 claims; singular/sparse/failed statuses are retained for auditability but excluded from primary significance claims.",
    x = grid::unit(0.035, "npc"),
    y = grid::unit(0.075, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 8.2, col = "#4b5563")
  )
}

png(filename = png_path, width = 12.6, height = 6.6, units = "in", res = 360, bg = "white")
draw_table()
dev.off()

pdf(file = pdf_path, width = 12.6, height = 6.6, bg = "white")
draw_table()
dev.off()

svglite::svglite(file = svg_path, width = 12.6, height = 6.6, bg = "white")
draw_table()
dev.off()

message("Supplementary Table 2 Class 1 per-window status outputs created:")
message(" - ", csv_path)
message(" - ", md_path)
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Input:")
message(" - ", results_path)
message("Rows: ", nrow(supplementary_table))
message("Total lmm_ok: ", sum(supplementary_table$lmm_ok))
message("Total q_global_class1_bh_lt_0_05: ", sum(supplementary_table$q_global_class1_bh_lt_0_05))
