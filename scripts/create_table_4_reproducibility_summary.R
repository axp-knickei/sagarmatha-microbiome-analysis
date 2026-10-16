# Create manuscript Table 4: reproducibility and computational traceability.
#
# The table summarizes the HPC per-MAG LMM execution record using the archived
# hpc_execution_record.csv, session_info.txt, and per_mag_lmertest_hpc_summary.txt.

required_packages <- c("dplyr", "grid", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

library(dplyr)

repo_root <- getwd()
execution_record_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "reproducibility",
  "hpc_execution_record.csv"
)
session_info_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "session",
  "session_info.txt"
)
summary_log_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "logs",
  "per_mag_lmertest_hpc_summary.txt"
)
output_dir <- file.path(repo_root, "results", "manuscript_tables")
output_stem <- "table_4_reproducibility_summary"

required_files <- c(execution_record_path, session_info_path, summary_log_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

execution <- utils::read.csv(execution_record_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(execution) < 1) {
  stop("Execution record is empty: ", execution_record_path, call. = FALSE)
}
record <- execution[1, , drop = FALSE]
session_lines <- readLines(session_info_path, warn = FALSE)
summary_lines <- readLines(summary_log_path, warn = FALSE)

field_value <- function(name, default = "Not recorded") {
  if (!name %in% names(record)) return(default)
  value <- as.character(record[[name]][[1]])
  if (is.na(value) || !nzchar(value)) default else value
}

summary_value <- function(pattern, default = "Not recorded") {
  hit <- grep(pattern, summary_lines, value = TRUE)
  if (length(hit) == 0) return(default)
  trimws(sub("^[^:]+:\\s*", "", hit[[1]]))
}

first_session_line <- function(pattern, default = "Not recorded") {
  hit <- grep(pattern, session_lines, value = TRUE)
  if (length(hit) == 0) return(default)
  trimws(hit[[1]])
}

fmt_runtime <- function(seconds) {
  seconds <- suppressWarnings(as.numeric(seconds))
  if (is.na(seconds)) return("Not recorded")
  if (seconds < 60) return(paste0(round(seconds, 1), " s"))
  paste0(round(seconds / 60, 1), " min")
}

wrap_text <- function(x, width) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

shorten_commit <- function(x) {
  if (is.na(x) || !nzchar(x) || x == "Not recorded") return(x)
  substr(x, 1, 12)
}

main_packages <- field_value("key_R_packages")
main_packages <- gsub(";\\s*", "; ", main_packages)

table_4 <- data.frame(
  field = c(
    "Analysis name",
    "Input file",
    "Script",
    "PBS job",
    "Scheduler / job ID",
    "Run time",
    "R version",
    "Main packages",
    "Conda environment",
    "Git commit",
    "Number of MAGs tested",
    "Primary eligible MAGs",
    "FDR method",
    "Main output table",
    "Session info file",
    "Summary log file"
  ),
  value = c(
    "per-MAG LMM",
    field_value("input_data_file"),
    field_value("script", "scripts/hpc_run_per_mag_lmertest.R"),
    "hpc/job_per_mag_lmertest.pbs",
    paste0(field_value("scheduler"), " / ", field_value("job_id")),
    fmt_runtime(field_value("runtime_seconds")),
    field_value("R_version", first_session_line("^R version")),
    main_packages,
    field_value("conda_environment", "glmm-hpc"),
    shorten_commit(field_value("git_commit")),
    field_value("number_of_tests", summary_value("^MAGs modeled:")),
    summary_value("^lmm_ok:"),
    "Benjamini-Hochberg",
    basename(field_value("main_results_file", "per_mag_lmertest_hpc_results.csv")),
    field_value("session_info_file", "results/hpc_per_mag_lmertest/session/session_info.txt"),
    field_value("summary_log_file", "results/hpc_per_mag_lmertest/logs/per_mag_lmertest_hpc_summary.txt")
  ),
  source = c(
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "Repository workflow",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv; session_info.txt",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv; summary log",
    "per_mag_lmertest_hpc_summary.txt",
    "per_mag_lmertest_hpc_summary.txt",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv"
  ),
  stringsAsFactors = FALSE
)

csv_path <- file.path(output_dir, paste0(output_stem, ".csv"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))
png_path <- file.path(output_dir, paste0(output_stem, ".png"))
utils::write.csv(table_4, csv_path, row.names = FALSE)

display_table <- table_4 %>%
  mutate(
    field = wrap_text(field, 26),
    value = wrap_text(value, 78),
    source = wrap_text(source, 34)
  )

col_labels <- c("Field", "Value", "Source")
col_widths <- c(0.22, 0.56, 0.22)
col_x <- cumsum(c(0, col_widths))

pal <- list(
  ink = "#111827",
  text = "#1F2937",
  muted = "#4B5563",
  hair = "#D1D5DB",
  header = "#E8F1F8",
  stripe = "#F8FAFC",
  blue = "#0072B2"
)

draw_table <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  grid::grid.text(
    "Supplementary Table 4 | Reproducibility summary",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.965, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 17, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    "Computational traceability for the full T1-T10 per-MAG longitudinal LMM workflow.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.925, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9.4, col = pal$muted)
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
      x = grid::unit(x_left + 0.006, "npc"),
      y = grid::unit(table_top - header_h / 2, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 8.5, fontface = "bold", col = pal$ink)
    )
  }

  for (i in seq_len(nrow(display_table))) {
    y_top <- table_top - header_h - (i - 1) * row_h
    y_mid <- y_top - row_h / 2
    fill <- if (i %% 2 == 0) "#FFFFFF" else pal$stripe

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
      color <- if (j == 1) pal$blue else pal$text
      fontface <- if (j == 1) "bold" else "plain"
      grid::grid.text(
        row_values[j],
        x = grid::unit(x_left + 0.006, "npc"),
        y = grid::unit(y_mid, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(fontsize = 6.85, col = color, fontface = fontface, lineheight = 0.86)
      )
    }
  }

  for (j in seq_along(col_x)) {
    x_pos <- 0.025 + 0.95 * col_x[j]
    grid::grid.lines(
      x = grid::unit(c(x_pos, x_pos), "npc"),
      y = grid::unit(c(table_bottom, table_top), "npc"),
      gp = grid::gpar(col = "#E5E7EB", lwd = 0.42)
    )
  }

  grid::grid.text(
    "Notes: The Git commit and job metadata are taken from the archived HPC execution record. A non-clean Git status was recorded at run time and should be considered when interpreting exact reproducibility.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.044, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.1, col = pal$muted)
  )
  grid::grid.text(
    "Abbreviations: FDR, false discovery rate; LMM, linear mixed model; MAG, metagenome-assembled genome; PBS, Portable Batch System.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.024, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.1, col = pal$muted)
  )
}

grDevices::png(filename = png_path, width = 13.2, height = 8.4, units = "in", res = 450, bg = "white", type = "cairo")
draw_table()
grDevices::dev.off()

svglite::svglite(file = svg_path, width = 13.2, height = 8.4, bg = "white")
draw_table()
grDevices::dev.off()

message("Table 4 saved:")
message("- ", csv_path)
message("- ", svg_path)
message("- ", png_path)
