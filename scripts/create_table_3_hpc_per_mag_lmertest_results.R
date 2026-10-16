# ================================================
# create_table_3_hpc_per_mag_lmertest_results.R
# ================================================
# Table 3: model status summary for the full T1-T10 per-MAG lmerTest LMM.
#
# This table summarizes model quality-control status from the existing per-MAG
# result table. It does not refit models or alter statistical decisions.

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
output_dir <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "tables")
output_stem <- "table_3_hpc_per_mag_lmertest_results_20260708_024546"

if (!file.exists(results_path)) {
  stop("Missing per-MAG LMM result table: ", results_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)

required_result_columns <- c(
  "MAG", "status", "fit_succeeded", "singular", "too_sparse",
  "low_variance", "primary_analysis", "zero_fraction", "response_variance",
  "n_subjects", "n_timepoints"
)
missing_result_columns <- setdiff(required_result_columns, names(results))
if (length(missing_result_columns) > 0) {
  stop("Missing required result columns: ", paste(missing_result_columns, collapse = ", "), call. = FALSE)
}

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_pct <- function(x) paste0(format(round(100 * x, 1), nsmall = 1, trim = TRUE), "%")

wrap_text <- function(x, width) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

expected_status <- data.frame(
  status = c("lmm_ok", "lmm_singular", "too_sparse", "low_variance", "model_failed"),
  interpretation = c(
    "Usable model included in the primary BH-FDR testing family.",
    "Model fit succeeded, but random-effect variance was effectively zero or singular.",
    "MAG had insufficient non-zero abundance support for reliable LMM fitting.",
    "Transformed response was too constant for reliable longitudinal modeling.",
    "Model fitting or required coefficient extraction failed."
  ),
  stringsAsFactors = FALSE
)

observed_counts <- results |>
  dplyr::count(status, name = "number_of_MAGs")

total_mags <- nrow(results)
status_summary <- expected_status |>
  dplyr::left_join(observed_counts, by = "status") |>
  dplyr::mutate(
    number_of_MAGs = ifelse(is.na(number_of_MAGs), 0L, number_of_MAGs),
    percent = if (total_mags > 0) number_of_MAGs / total_mags else NA_real_
  )

unexpected_status <- setdiff(unique(results$status), expected_status$status)
if (length(unexpected_status) > 0) {
  unexpected_rows <- data.frame(
    status = unexpected_status,
    interpretation = "Additional status observed in result table; inspect workflow logs before manuscript reporting.",
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(observed_counts, by = "status") |>
    dplyr::mutate(percent = if (total_mags > 0) number_of_MAGs / total_mags else NA_real_)
  status_summary <- dplyr::bind_rows(status_summary, unexpected_rows)
}

qc_summary <- results |>
  dplyr::group_by(status) |>
  dplyr::summarise(
    median_zero_fraction = stats::median(zero_fraction, na.rm = TRUE),
    median_response_variance = stats::median(response_variance, na.rm = TRUE),
    median_subjects = stats::median(n_subjects, na.rm = TRUE),
    median_timepoints = stats::median(n_timepoints, na.rm = TRUE),
    fit_succeeded = sum(fit_succeeded %in% TRUE, na.rm = TRUE),
    primary_analysis = sum(primary_analysis %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  )

table_3 <- status_summary |>
  dplyr::left_join(qc_summary, by = "status") |>
  dplyr::mutate(
    median_zero_fraction = ifelse(is.na(median_zero_fraction), NA_real_, median_zero_fraction),
    median_response_variance = ifelse(is.na(median_response_variance), NA_real_, median_response_variance),
    median_subjects = ifelse(is.na(median_subjects), NA_real_, median_subjects),
    median_timepoints = ifelse(is.na(median_timepoints), NA_real_, median_timepoints),
    fit_succeeded = ifelse(is.na(fit_succeeded), 0L, fit_succeeded),
    primary_analysis = ifelse(is.na(primary_analysis), 0L, primary_analysis)
  ) |>
  dplyr::transmute(
    Status = status,
    `Number of MAGs` = number_of_MAGs,
    Percent = fmt_pct(percent),
    Interpretation = interpretation,
    `Median zero fraction` = ifelse(is.na(median_zero_fraction), "NA", fmt_pct(median_zero_fraction)),
    `Median response variance` = ifelse(is.na(median_response_variance), "NA", formatC(median_response_variance, format = "e", digits = 2)),
    `Median subjects` = ifelse(is.na(median_subjects), "NA", as.character(round(median_subjects, 1))),
    `Median timepoints` = ifelse(is.na(median_timepoints), "NA", as.character(round(median_timepoints, 1))),
    `Fit succeeded` = fit_succeeded,
    `Primary analysis` = primary_analysis
  )

csv_path <- file.path(output_dir, paste0(output_stem, ".csv"))
md_path <- file.path(output_dir, paste0(output_stem, ".md"))
png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

utils::write.csv(table_3, csv_path, row.names = FALSE)

markdown_lines <- c(
  "# Table 3. Model Status Summary",
  "",
  "This table summarizes per-MAG model quality-control status for the full T1-T10 lmerTest LMM workflow. Zero-count expected statuses are retained to document checks that were available but not triggered.",
  "",
  paste(names(table_3), collapse = " | "),
  paste(rep("---", ncol(table_3)), collapse = " | ")
)
for (i in seq_len(nrow(table_3))) {
  row_values <- vapply(table_3[i, ], as.character, character(1))
  row_values <- gsub("\\|", "/", row_values)
  row_values[is.na(row_values)] <- ""
  markdown_lines <- c(markdown_lines, paste(row_values, collapse = " | "))
}
writeLines(markdown_lines, md_path)

display_table <- table_3 |>
  dplyr::select(Status, `Number of MAGs`, Percent, Interpretation) |>
  dplyr::mutate(
    Status = wrap_text(Status, 14),
    `Number of MAGs` = fmt_int(`Number of MAGs`),
    Interpretation = wrap_text(Interpretation, 82)
  )

col_labels <- c(
  "Status", "Number of MAGs", "Percent", "Interpretation"
)
col_widths <- c(0.18, 0.16, 0.12, 0.54)
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
    "Table 3 | Model status summary",
    x = grid::unit(0.022, "npc"),
    y = grid::unit(0.965, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 17, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    paste0(
      "Transparency and quality control for ", fmt_int(total_mags),
      " MAG-level LMM attempts. Full QC medians are retained in the CSV/Markdown outputs."
    ),
    x = grid::unit(0.022, "npc"),
    y = grid::unit(0.925, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9.4, col = pal$muted)
  )

  table_top <- 0.84
  table_bottom <- 0.09
  header_h <- 0.09
  row_h <- (table_top - table_bottom - header_h) / nrow(display_table)

  grid::grid.rect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(table_top - header_h / 2, "npc"),
    width = grid::unit(0.956, "npc"),
    height = grid::unit(header_h, "npc"),
    gp = grid::gpar(fill = pal$header, col = pal$hair, lwd = 0.8)
  )

  for (j in seq_along(col_labels)) {
    x_left <- 0.022 + 0.956 * col_x[j]
    grid::grid.text(
      col_labels[j],
      x = grid::unit(x_left + 0.004, "npc"),
      y = grid::unit(table_top - header_h / 2, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 7.35, fontface = "bold", col = pal$ink, lineheight = 0.86)
    )
  }

  for (i in seq_len(nrow(display_table))) {
    y_top <- table_top - header_h - (i - 1) * row_h
    y_mid <- y_top - row_h / 2
    fill <- if (i %% 2 == 0) "#FFFFFF" else pal$stripe
    grid::grid.rect(
      x = grid::unit(0.5, "npc"),
      y = grid::unit(y_mid, "npc"),
      width = grid::unit(0.956, "npc"),
      height = grid::unit(row_h, "npc"),
      gp = grid::gpar(fill = fill, col = "#E5E7EB", lwd = 0.45)
    )

    row_values <- as.character(display_table[i, ])
    for (j in seq_along(row_values)) {
      x_left <- 0.022 + 0.956 * col_x[j]
      color <- pal$text
      fontface <- "plain"
      if (col_labels[j] == "Status") {
        fontface <- "bold"
        color <- ifelse(row_values[j] == "lmm_ok", pal$blue, pal$orange)
      }
      grid::grid.text(
        row_values[j],
        x = grid::unit(x_left + 0.004, "npc"),
        y = grid::unit(y_mid, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(fontsize = 6.85, col = color, fontface = fontface, lineheight = 0.84)
      )
    }
  }

  for (j in seq_along(col_x)) {
    x_pos <- 0.022 + 0.956 * col_x[j]
    grid::grid.lines(
      x = grid::unit(c(x_pos, x_pos), "npc"),
      y = grid::unit(c(table_bottom, table_top), "npc"),
      gp = grid::gpar(col = "#E5E7EB", lwd = 0.45)
    )
  }

  grid::grid.text(
    "Percent is calculated across all MAG rows in per_mag_lmertest_hpc_results.csv. Expected zero-count status categories are retained.",
    x = grid::unit(0.022, "npc"),
    y = grid::unit(0.036, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.5, col = pal$muted)
  )
}

grDevices::png(filename = png_path, width = 10.8, height = 5.8, units = "in", res = 450, bg = "white", type = "cairo")
draw_table()
grDevices::dev.off()

grDevices::pdf(file = pdf_path, width = 10.8, height = 5.8, bg = "white", useDingbats = FALSE)
draw_table()
grDevices::dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = svg_path, width = 10.8, height = 5.8, bg = "white")
  draw_table()
  grDevices::dev.off()
} else {
  warning("Package 'svglite' is not installed; SVG output was skipped.", call. = FALSE)
}

message("Table 3 saved:")
message("- ", csv_path)
message("- ", md_path)
message("- ", png_path)
message("- ", pdf_path)
if (file.exists(svg_path)) {
  message("- ", svg_path)
}
