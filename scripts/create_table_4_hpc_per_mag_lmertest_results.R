# ================================================
# create_table_4_hpc_per_mag_lmertest_results.R
# ================================================
# Table 4: reproducibility table for the full T1-T10 per-MAG lmerTest LMM.
#
# This table summarizes computational traceability from the execution record,
# session information, summary log, and primary result table. It does not modify
# or recompute any scientific analysis.

required_packages <- c("grid")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
record_path <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "reproducibility", "hpc_execution_record.csv")
session_path <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "session", "session_info.txt")
summary_path <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "logs", "per_mag_lmertest_hpc_summary.txt")
results_path <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "tables", "per_mag_lmertest_hpc_results.csv")
pbs_path <- "hpc/job_per_mag_lmertest.pbs"
output_dir <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "tables")
output_stem <- "table_4_hpc_per_mag_lmertest_results_20260708_024546"

for (path in c(record_path, session_path, summary_path, results_path)) {
  if (!file.exists(path)) {
    stop("Missing required reproducibility source file: ", path, call. = FALSE)
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

record <- utils::read.csv(record_path, check.names = FALSE, stringsAsFactors = FALSE)
results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
session_lines <- readLines(session_path, warn = FALSE)
summary_lines <- readLines(summary_path, warn = FALSE)

get_record <- function(name, default = NA_character_) {
  if (!name %in% names(record) || nrow(record) == 0) {
    return(default)
  }
  value <- as.character(record[[name]][[1]])
  if (is.na(value) || !nzchar(value)) default else value
}

extract_summary_value <- function(prefix, default = NA_character_) {
  hit <- grep(paste0("^", prefix), summary_lines, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  trimws(sub(paste0("^", prefix), "", hit[[1]]))
}

extract_r_version <- function(default = NA_character_) {
  hit <- grep("^R version ", session_lines, value = TRUE)
  if (length(hit) > 0) {
    return(hit[[1]])
  }
  value <- get_record("R_version", default)
  value
}

shorten_commit <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(NA_character_)
  }
  if (nchar(x) > 12) substr(x, 1, 12) else x
}

wrap_text <- function(x, width) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)

number_of_tests <- get_record("number_of_tests")
if (is.na(number_of_tests)) {
  number_of_tests <- as.character(nrow(results))
}

fdr_method <- extract_summary_value("FDR method:", default = NA_character_)
if (is.na(fdr_method)) {
  fdr_method <- "Benjamini-Hochberg"
}

main_packages <- get_record("key_R_packages")
if (is.na(main_packages)) {
  attached <- grep("^(other attached packages:|\\[1\\])", session_lines, value = TRUE)
  main_packages <- paste(attached, collapse = " ")
}

table_4 <- data.frame(
  Field = c(
    "Analysis name",
    "Input file",
    "Script",
    "PBS job",
    "Scheduler",
    "Job ID",
    "Hostname",
    "R version",
    "Main packages",
    "Conda module",
    "Conda environment",
    "Git branch",
    "Git commit",
    "Git status",
    "Start time",
    "Completion time",
    "Runtime",
    "Number of MAGs tested",
    "FDR method",
    "FDR family",
    "Main output table",
    "Summary log",
    "Session info",
    "Execution record"
  ),
  Value = c(
    get_record("analysis_name", "full_T1_T10_per_MAG_lmertest_LMM"),
    get_record("input_data_file", "data/processed/cleaned_data.rds"),
    get_record("script_path", "scripts/hpc_run_per_mag_lmertest.R"),
    pbs_path,
    get_record("scheduler", "OpenPBS"),
    get_record("job_id"),
    get_record("hostname"),
    extract_r_version(),
    main_packages,
    get_record("conda_module"),
    get_record("conda_environment", "glmm-hpc"),
    get_record("git_branch"),
    get_record("git_commit"),
    get_record("git_status_short"),
    get_record("start_time"),
    get_record("completion_time"),
    paste0(round(as.numeric(get_record("runtime_seconds", NA_character_)), 1), " seconds"),
    number_of_tests,
    fdr_method,
    get_record("fdr_family", get_record("FDR_family")),
    get_record("main_results_file", "results/hpc_per_mag_lmertest/tables/per_mag_lmertest_hpc_results.csv"),
    get_record("summary_log_file", "results/hpc_per_mag_lmertest/logs/per_mag_lmertest_hpc_summary.txt"),
    get_record("session_info_file", "results/hpc_per_mag_lmertest/session/session_info.txt"),
    "results/hpc_per_mag_lmertest/reproducibility/hpc_execution_record.csv"
  ),
  Source = c(
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "repository path",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "session_info.txt",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv / result table fallback",
    "per_mag_lmertest_hpc_summary.txt",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "hpc_execution_record.csv",
    "repository path"
  ),
  stringsAsFactors = FALSE
)

table_4$Value[table_4$Field == "Git commit"] <- shorten_commit(table_4$Value[table_4$Field == "Git commit"])
table_4$Value[table_4$Field == "Git status"] <- ifelse(
  is.na(table_4$Value[table_4$Field == "Git status"]) | !nzchar(table_4$Value[table_4$Field == "Git status"]),
  "clean or not recorded",
  table_4$Value[table_4$Field == "Git status"]
)
table_4$Value[table_4$Field == "Number of MAGs tested"] <- fmt_int(as.numeric(table_4$Value[table_4$Field == "Number of MAGs tested"]))

csv_path <- file.path(output_dir, paste0(output_stem, ".csv"))
md_path <- file.path(output_dir, paste0(output_stem, ".md"))
png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

utils::write.csv(table_4, csv_path, row.names = FALSE)

markdown_lines <- c(
  "# Table 4. Reproducibility Table",
  "",
  "This table summarizes computational traceability for the full T1-T10 per-MAG lmerTest LMM workflow. It is metadata/logging only and does not change the statistical analysis.",
  "",
  paste(names(table_4), collapse = " | "),
  paste(rep("---", ncol(table_4)), collapse = " | ")
)
for (i in seq_len(nrow(table_4))) {
  row_values <- vapply(table_4[i, ], as.character, character(1))
  row_values <- gsub("\\|", "/", row_values)
  row_values[is.na(row_values)] <- ""
  markdown_lines <- c(markdown_lines, paste(row_values, collapse = " | "))
}
writeLines(markdown_lines, md_path)

display_fields <- c(
  "Analysis name", "Input file", "Script", "PBS job", "Scheduler", "Job ID",
  "Hostname", "R version", "Main packages", "Conda environment", "Git commit",
  "Number of MAGs tested", "FDR method", "Main output table"
)

display_table <- table_4[table_4$Field %in% display_fields, c("Field", "Value", "Source")]
display_table$Field <- factor(display_table$Field, levels = display_fields)
display_table <- display_table[order(display_table$Field), ]
display_table$Field <- as.character(display_table$Field)
display_table$Value <- wrap_text(display_table$Value, 76)
display_table$Source <- wrap_text(display_table$Source, 28)

col_labels <- c("Field", "Value", "Source")
col_widths <- c(0.22, 0.57, 0.21)
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
    "Table 4 | Reproducibility table",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.965, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 17, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    "Computational traceability for the full T1-T10 per-MAG lmerTest LMM workflow.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.925, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9.4, col = pal$muted)
  )

  table_top <- 0.86
  table_bottom <- 0.07
  header_h <- 0.055
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
      gp = grid::gpar(fontsize = 8.2, fontface = "bold", col = pal$ink)
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
      gp = grid::gpar(fill = fill, col = "#E5E7EB", lwd = 0.4)
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
        gp = grid::gpar(fontsize = 6.8, col = color, fontface = fontface, lineheight = 0.82)
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
    "Full CSV/Markdown outputs include additional scheduler, runtime, git-status, log, session, and execution-record paths.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.032, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.5, col = pal$muted)
  )
}

grDevices::png(filename = png_path, width = 12.6, height = 8.2, units = "in", res = 450, bg = "white", type = "cairo")
draw_table()
grDevices::dev.off()

grDevices::pdf(file = pdf_path, width = 12.6, height = 8.2, bg = "white", useDingbats = FALSE)
draw_table()
grDevices::dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = svg_path, width = 12.6, height = 8.2, bg = "white")
  draw_table()
  grDevices::dev.off()
} else {
  warning("Package 'svglite' is not installed; SVG output was skipped.", call. = FALSE)
}

message("Table 4 saved:")
message("- ", csv_path)
message("- ", md_path)
message("- ", png_path)
message("- ", pdf_path)
if (file.exists(svg_path)) {
  message("- ", svg_path)
}
