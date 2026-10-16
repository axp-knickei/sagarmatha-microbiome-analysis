# ============================================
# create_supplementary_table_4_class1_analysis.R
# ============================================
# Supplementary Table 4: reproducibility and HPC execution record for the
# predefined Class 1 window LMM analysis.

required_packages <- c("readr", "dplyr", "grid")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
job_script_path <- file.path(repo_root, "hpc", "job_class1_window_lmertest.pbs")
summary_path <- file.path(
  repo_root,
  "results",
  "class1_window_lmertest",
  "logs",
  "class1_window_lmertest_summary.txt"
)
session_info_path <- file.path(
  repo_root,
  "results",
  "class1_window_lmertest",
  "session",
  "session_info.txt"
)
results_table_path <- file.path(
  repo_root,
  "results",
  "class1_window_lmertest",
  "tables",
  "class1_window_lmertest_results.csv"
)
output_dir <- file.path(repo_root, "results", "class1_window_lmertest", "tables")
figure_dir <- file.path(repo_root, "results", "class1_window_lmertest", "figures")
logs_dir <- file.path(repo_root, "logs")

required_files <- c(job_script_path, summary_path, session_info_path, results_table_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "), call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

summary_lines <- readLines(summary_path, warn = FALSE)
job_script_lines <- readLines(job_script_path, warn = FALSE)
session_lines <- readLines(session_info_path, warn = FALSE)
results_table <- readr::read_csv(results_table_path, show_col_types = FALSE)

extract_value <- function(lines, prefix) {
  hit <- grep(paste0("^", prefix), lines, value = TRUE)
  if (length(hit) == 0) {
    return(NA_character_)
  }
  trimws(sub(paste0("^", prefix), "", hit[[1]]))
}

extract_regex <- function(lines, pattern) {
  hit <- grep(pattern, lines, value = TRUE)
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[[1]]
}

job_name_line <- extract_regex(job_script_lines, "^#PBS -N ")
job_name <- trimws(sub("^#PBS -N ", "", job_name_line))

log_pattern <- paste0("^", gsub("\\.", "\\\\.", job_name), "\\..*\\.out$")
out_logs <- list.files(logs_dir, pattern = log_pattern, full.names = TRUE)
job_log_path <- if (length(out_logs) > 0) out_logs[[1]] else NA_character_
job_log_lines <- if (!is.na(job_log_path)) readLines(job_log_path, warn = FALSE) else character()

job_id_line <- extract_regex(job_log_lines, "^PBS job id: ")
job_id <- trimws(sub("^PBS job id: ", "", job_id_line))

completion_line <- extract_regex(job_log_lines, "^=== BlackICE Class 1 window LMM job finished: ")
completion_time <- trimws(
  sub("^=== BlackICE Class 1 window LMM job finished: ", "", sub(" ===$", "", completion_line))
)

project_name <- basename(repo_root)
script_name <- "scripts/hpc_run_class1_window_lmertest.R"
scheduler_name <- "OpenPBS"
input_file <- extract_value(summary_lines, "Input file: ")
windows_file <- extract_value(summary_lines, "Window definition file: ")
conda_environment <- trimws(sub("^Conda environment: ", "", extract_regex(job_log_lines, "^Conda environment: ")))
number_of_tests <- as.character(nrow(results_table))
fdr_family <- "Benjamini-Hochberg FDR across all predefined Class 1 MAG-window tests"
session_info_file <- file.path("results", "class1_window_lmertest", "session", "session_info.txt")

supplementary_table <- tibble::tibble(
  project = project_name,
  script = script_name,
  scheduler = scheduler_name,
  `job ID` = job_id,
  `completion time` = completion_time,
  `input data file` = input_file,
  `window definition file` = windows_file,
  `Conda environment` = conda_environment,
  `number of tests` = number_of_tests,
  `FDR family` = fdr_family,
  `session info file` = session_info_file
)

csv_path <- file.path(output_dir, "supplementary_table_4_class1_reproducibility_hpc_record.csv")
md_path <- file.path(output_dir, "supplementary_table_4_class1_reproducibility_hpc_record.md")
png_path <- file.path(figure_dir, "supplementary_table_4_class1_reproducibility_hpc_record.png")
pdf_path <- file.path(figure_dir, "supplementary_table_4_class1_reproducibility_hpc_record.pdf")
svg_path <- file.path(figure_dir, "supplementary_table_4_class1_reproducibility_hpc_record.svg")

readr::write_csv(supplementary_table, csv_path)

markdown_lines <- c(
  "# Supplementary Table 4. Reproducibility and HPC Execution Record",
  "",
  "This table records the actual HPC execution context for the predefined Class 1 window LMM analysis.",
  "",
  paste(names(supplementary_table), collapse = " | "),
  paste(rep("---", ncol(supplementary_table)), collapse = " | ")
)
row_values <- vapply(supplementary_table[1, ], as.character, character(1))
row_values <- gsub("\\|", "/", row_values)
row_values[is.na(row_values)] <- ""
markdown_lines <- c(markdown_lines, paste(row_values, collapse = " | "))
writeLines(markdown_lines, md_path)

wrap_text <- function(x, width) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

pretty_table <- supplementary_table |>
  dplyr::mutate(
    project = wrap_text(project, 15),
    script = wrap_text(script, 26),
    scheduler = wrap_text(scheduler, 10),
    `job ID` = wrap_text(`job ID`, 16),
    `completion time` = wrap_text(`completion time`, 26),
    `input data file` = wrap_text(`input data file`, 28),
    `window definition file` = wrap_text(`window definition file`, 32),
    `Conda environment` = wrap_text(`Conda environment`, 16),
    `number of tests` = wrap_text(`number of tests`, 10),
    `FDR family` = wrap_text(`FDR family`, 30),
    `session info file` = wrap_text(`session info file`, 30)
  )

display_long <- tibble::tibble(
  field = names(pretty_table),
  value = vapply(pretty_table[1, ], as.character, character(1))
)

draw_table <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text(
    "Supplementary Table 4 | Reproducibility and HPC execution record",
    x = grid::unit(0.03, "npc"),
    y = grid::unit(0.955, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 16, fontface = "bold", col = "#111827")
  )
  grid::grid.text(
    "Execution metadata are pulled from the actual Class 1 PBS wrapper, summary log, session info file, and scheduler output.",
    x = grid::unit(0.03, "npc"),
    y = grid::unit(0.915, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9.1, col = "#4b5563")
  )

  table_top <- 0.84
  table_bottom <- 0.10
  header_h <- 0.07
  row_h <- (table_top - table_bottom - header_h) / nrow(display_long)

  grid::grid.rect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(table_top - header_h / 2, "npc"),
    width = grid::unit(0.94, "npc"),
    height = grid::unit(header_h, "npc"),
    gp = grid::gpar(fill = "#f3f4f6", col = "#d1d5db", lwd = 0.8)
  )

  x_split <- 0.26
  grid::grid.text(
    "Field",
    x = grid::unit(0.03 + 0.94 * (x_split / 2), "npc"),
    y = grid::unit(table_top - header_h / 2, "npc"),
    just = "center",
    gp = grid::gpar(fontsize = 8.2, fontface = "bold", col = "#111827")
  )
  grid::grid.text(
    "Value",
    x = grid::unit(0.03 + 0.94 * (x_split + (1 - x_split) / 2), "npc"),
    y = grid::unit(table_top - header_h / 2, "npc"),
    just = "center",
    gp = grid::gpar(fontsize = 8.2, fontface = "bold", col = "#111827")
  )

  for (i in seq_len(nrow(display_long))) {
    y_top <- table_top - header_h - (i - 1) * row_h
    y_mid <- y_top - row_h / 2
    fill <- if (i %% 2 == 0) "#ffffff" else "#f8fafc"
    grid::grid.rect(
      x = grid::unit(0.5, "npc"),
      y = grid::unit(y_mid, "npc"),
      width = grid::unit(0.94, "npc"),
      height = grid::unit(row_h, "npc"),
      gp = grid::gpar(fill = fill, col = "#e5e7eb", lwd = 0.45)
    )
    grid::grid.text(
      as.character(display_long$field[i]),
      x = grid::unit(0.03 + 0.008, "npc"),
      y = grid::unit(y_mid, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 7.6, fontface = "bold", col = "#1f2937", lineheight = 0.9)
    )
    grid::grid.text(
      as.character(display_long$value[i]),
      x = grid::unit(0.03 + 0.94 * x_split + 0.008, "npc"),
      y = grid::unit(y_mid, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 7.4, col = "#1f2937", lineheight = 0.9)
    )
  }

  for (x_pos in c(0, x_split, 1)) {
    grid::grid.lines(
      x = grid::unit(rep(0.03 + 0.94 * x_pos, 2), "npc"),
      y = grid::unit(c(table_bottom, table_top), "npc"),
      gp = grid::gpar(col = "#e5e7eb", lwd = 0.45)
    )
  }

  grid::grid.text(
    "Run summary: 2646 MAG-window tests across 6 predefined windows; 0 global Class 1 BH-FDR significant results.",
    x = grid::unit(0.03, "npc"),
    y = grid::unit(0.055, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 8.0, col = "#4b5563")
  )
}

png(filename = png_path, width = 15, height = 5.8, units = "in", res = 360, bg = "white")
draw_table()
dev.off()

pdf(file = pdf_path, width = 15, height = 5.8, bg = "white")
draw_table()
dev.off()

svglite::svglite(file = svg_path, width = 15, height = 5.8, bg = "white")
draw_table()
dev.off()

message("Supplementary Table 4 reproducibility outputs created:")
message(" - ", csv_path)
message(" - ", md_path)
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Job ID: ", job_id)
message("Completion time: ", completion_time)
message("Conda environment: ", conda_environment)
message("Number of tests: ", number_of_tests)
