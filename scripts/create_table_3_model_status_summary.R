# Create manuscript Table 3: model status summary for per-MAG LMM quality control.
#
# This table reports every expected model status category, including categories
# with zero MAGs in the current run, to make model inclusion/exclusion
# transparent for readers and reviewers.

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
output_dir <- file.path(repo_root, "results", "manuscript_tables")
output_stem <- "table_3_model_status_summary"

if (!file.exists(results_path)) {
  stop("Missing per-MAG LMM result table: ", results_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("MAG", "status", "fit_succeeded", "singular", "too_sparse", "low_variance", "primary_analysis")
missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns) > 0) {
  stop("Missing required result columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_pct <- function(x, digits = 1) paste0(format(round(100 * x, digits), nsmall = digits, trim = TRUE), "%")

wrap_text <- function(x, width) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

status_levels <- c("lmm_ok", "lmm_singular", "too_sparse", "low_variance", "model_failed")
status_interpretation <- c(
  lmm_ok = "Usable model; included in the primary longitudinal testing family.",
  lmm_singular = "Model fit completed but random-effect variance was effectively zero or singular.",
  too_sparse = "Insufficient abundance support for reliable longitudinal modeling.",
  low_variance = "Response was too constant to support a meaningful longitudinal model.",
  model_failed = "Model fitting or output extraction failed."
)

observed_status <- as.character(results$status)
observed_status[is.na(observed_status) | !nzchar(observed_status)] <- "model_failed"
observed_status[!(observed_status %in% status_levels)] <- "model_failed"
results$status_standardized <- observed_status

total_mags <- dplyr::n_distinct(results$MAG)
status_counts <- table(factor(results$status_standardized, levels = status_levels))

table_3 <- data.frame(
  status = status_levels,
  number_of_MAGs = as.integer(status_counts),
  percent = as.integer(status_counts) / total_mags,
  interpretation = unname(status_interpretation[status_levels]),
  stringsAsFactors = FALSE
) %>%
  mutate(
    percent = ifelse(is.na(percent), 0, percent),
    note = case_when(
      status == "lmm_ok" ~ paste0(fmt_int(sum(results$primary_analysis %in% TRUE, na.rm = TRUE)), " MAGs marked primary_analysis == TRUE."),
      status == "lmm_singular" ~ paste0(fmt_int(sum(results$singular %in% TRUE, na.rm = TRUE)), " rows have singular == TRUE."),
      status == "too_sparse" ~ paste0(fmt_int(sum(results$too_sparse %in% TRUE, na.rm = TRUE)), " rows have too_sparse == TRUE."),
      status == "low_variance" ~ paste0(fmt_int(sum(results$low_variance %in% TRUE, na.rm = TRUE)), " rows have low_variance == TRUE."),
      status == "model_failed" ~ paste0(fmt_int(sum(!(results$fit_succeeded %in% TRUE), na.rm = TRUE)), " rows have fit_succeeded != TRUE."),
      TRUE ~ ""
    )
  )

csv_table <- table_3 %>%
  transmute(
    Status = status,
    `Number of MAGs` = number_of_MAGs,
    Percent = fmt_pct(percent, 1),
    Interpretation = interpretation,
    Note = note
  )

csv_path <- file.path(output_dir, paste0(output_stem, ".csv"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))
png_path <- file.path(output_dir, paste0(output_stem, ".png"))
utils::write.csv(csv_table, csv_path, row.names = FALSE)

display_table <- table_3 %>%
  transmute(
    Status = status,
    `Number of MAGs` = fmt_int(number_of_MAGs),
    Percent = fmt_pct(percent, 1),
    Interpretation = wrap_text(interpretation, 54),
    Note = wrap_text(note, 42)
  )

col_labels <- c("Status", "Number of MAGs", "Percent", "Interpretation", "QC note")
col_widths <- c(0.15, 0.14, 0.10, 0.39, 0.22)
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
    "Table 3 | Model status summary",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.965, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 17, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    paste0(
      "Quality-control summary for ", fmt_int(total_mags),
      " MAG-specific longitudinal LMM fits. Percentages use all MAGs in the per-MAG result table as the denominator."
    ),
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.920, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9.4, col = pal$muted)
  )

  table_top <- 0.830
  table_bottom <- 0.170
  header_h <- 0.080
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
      x = grid::unit(x_left + 0.006, "npc"),
      y = grid::unit(table_top - header_h / 2, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 9.0, fontface = "bold", col = pal$ink)
    )
  }

  for (i in seq_len(nrow(display_table))) {
    y_top <- table_top - header_h - (i - 1) * row_h
    y_mid <- y_top - row_h / 2
    fill <- if (i %% 2 == 0) "#FFFFFF" else pal$stripe
    is_usable <- display_table$Status[[i]] == "lmm_ok"

    grid::grid.rect(
      x = grid::unit(0.5, "npc"),
      y = grid::unit(y_mid, "npc"),
      width = grid::unit(0.95, "npc"),
      height = grid::unit(row_h, "npc"),
      gp = grid::gpar(fill = fill, col = "#E5E7EB", lwd = 0.45)
    )

    row_values <- as.character(display_table[i, ])
    for (j in seq_along(row_values)) {
      x_left <- 0.025 + 0.95 * col_x[j]
      color <- if (is_usable && j == 1) pal$blue else pal$text
      fontface <- if (j == 1 || (is_usable && j %in% c(2, 3))) "bold" else "plain"
      grid::grid.text(
        row_values[j],
        x = grid::unit(x_left + 0.006, "npc"),
        y = grid::unit(y_mid, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(fontsize = 8.0, col = color, fontface = fontface, lineheight = 0.90)
      )
    }
  }

  for (j in seq_along(col_x)) {
    x_pos <- 0.025 + 0.95 * col_x[j]
    grid::grid.lines(
      x = grid::unit(c(x_pos, x_pos), "npc"),
      y = grid::unit(c(table_bottom, table_top), "npc"),
      gp = grid::gpar(col = "#E5E7EB", lwd = 0.45)
    )
  }

  grid::grid.text(
    "Notes: Status values are read from per_mag_lmertest_hpc_results.csv. Categories with zero MAGs are retained to document expected QC outcomes.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.100, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.8, col = pal$muted)
  )
  grid::grid.text(
    "Abbreviations: LMM, linear mixed model; MAG, metagenome-assembled genome; QC, quality control.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.065, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.8, col = pal$muted)
  )
}

grDevices::png(filename = png_path, width = 11.8, height = 5.8, units = "in", res = 450, bg = "white", type = "cairo")
draw_table()
grDevices::dev.off()

svglite::svglite(file = svg_path, width = 11.8, height = 5.8, bg = "white")
draw_table()
grDevices::dev.off()

message("Table 3 saved:")
message("- ", csv_path)
message("- ", svg_path)
message("- ", png_path)
