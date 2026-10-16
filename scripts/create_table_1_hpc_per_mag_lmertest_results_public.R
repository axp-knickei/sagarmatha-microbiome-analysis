# Public renderer: values come from external inputs; see docs/public-rendering.md.
# ================================================
# create_table_1_hpc_per_mag_lmertest_results.R
# ================================================
# Table 1: sample and metadata summary for the full T1-T10 per-MAG lmerTest
# LMM workflow.
#
# This table describes the repeated-measures study design and the model
# inclusion outcome. It does not modify, filter, or recompute the primary LMM
# result table.

required_packages <- c("dplyr", "grid")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
data_path <- file.path(repo_root, "data", "processed", "cleaned_data.rds")
results_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "tables",
  "per_mag_lmertest_hpc_results.csv"
)
output_dir <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "tables")
output_stem <- "table_1_hpc_per_mag_lmertest_results"

# External inputs and isolated output destination (run from repository root).
data_path <- Sys.getenv("PUBLIC_CLEANED_DATA", unset = data_path)
results_path <- Sys.getenv("PUBLIC_LMM_RESULTS", unset = results_path)
output_dir <- Sys.getenv("PUBLIC_RENDER_OUTPUT_DIR", unset = output_dir)

if (!file.exists(data_path)) {
  stop("Missing cleaned data file: ", data_path, call. = FALSE)
}
if (!file.exists(results_path)) {
  stop("Missing per-MAG LMM result table: ", results_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cleaned <- readRDS(data_path)
results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)

required_data_columns <- c("sampleID", "MAG", "Subject", "Timepoint", "age", "gender", "altitude")
required_result_columns <- c("MAG", "status", "primary_analysis", "fit_succeeded", "primary_significant", "time_p_value")
missing_data_columns <- setdiff(required_data_columns, names(cleaned))
missing_result_columns <- setdiff(required_result_columns, names(results))
if (length(missing_data_columns) > 0) {
  stop("Missing required cleaned-data columns: ", paste(missing_data_columns, collapse = ", "), call. = FALSE)
}
if (length(missing_result_columns) > 0) {
  stop("Missing required result columns: ", paste(missing_result_columns, collapse = ", "), call. = FALSE)
}

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_num <- function(x, digits = 1) format(round(x, digits), nsmall = digits, trim = TRUE)

collapse_named_counts <- function(counts) {
  paste(paste0(names(counts), ": ", as.integer(counts)), collapse = "; ")
}

wrap_text <- function(x, width) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

sample_metadata <- cleaned |>
  dplyr::distinct(sampleID, Subject, Timepoint, age, gender, altitude)

subject_metadata <- sample_metadata |>
  dplyr::distinct(Subject, age, gender)

timepoint_levels <- paste0("T", seq_len(10))
observed_timepoints <- sort(unique(as.character(sample_metadata$Timepoint)))
if (all(timepoint_levels %in% observed_timepoints)) {
  ordered_timepoints <- timepoint_levels
} else {
  ordered_timepoints <- observed_timepoints
}

samples_per_timepoint <- table(factor(sample_metadata$Timepoint, levels = ordered_timepoints))
n_subjects <- dplyr::n_distinct(subject_metadata$Subject)
missing_per_timepoint <- pmax(n_subjects - as.integer(samples_per_timepoint), 0L)
names(missing_per_timepoint) <- names(samples_per_timepoint)

age_values <- unique(subject_metadata[, c("Subject", "age")])$age
gender_counts <- table(subject_metadata$gender, useNA = "ifany")
status_counts <- table(results$status, useNA = "ifany")

total_mags <- dplyr::n_distinct(results$MAG)
primary_mags <- sum(results$primary_analysis %in% TRUE, na.rm = TRUE)
lmm_ok_mags <- sum(results$status == "lmm_ok", na.rm = TRUE)
excluded_or_failed <- total_mags - primary_mags
model_failed_or_non_ok <- sum(!(results$status %in% "lmm_ok"), na.rm = TRUE)
nominal_mags <- sum(results$primary_analysis %in% TRUE & results$time_p_value < 0.05, na.rm = TRUE)
fdr_mags <- sum(results$primary_significant %in% TRUE, na.rm = TRUE)

table_1 <- data.frame(
  section = c(
    rep("Study design", 5),
    rep("Participant metadata", 3),
    rep("MAG analysis set", 4)
  ),
  metric = c(
    "Number of subjects",
    "Number of samples",
    "Number of timepoints",
    "Samples per timepoint",
    "Missing samples per timepoint",
    "Age summary",
    "Gender distribution",
    "Altitude range",
    "Total MAGs",
    "MAGs modeled in primary family",
    "MAGs excluded or failed primary modeling",
    "Exploratory/FDR result counts"
  ),
  value = c(
    fmt_int(n_subjects),
    fmt_int(dplyr::n_distinct(sample_metadata$sampleID)),
    fmt_int(length(unique(sample_metadata$Timepoint))),
    collapse_named_counts(samples_per_timepoint),
    collapse_named_counts(missing_per_timepoint),
    paste0(
      "median ", fmt_num(stats::median(age_values, na.rm = TRUE), 1),
      " years; range ", fmt_num(min(age_values, na.rm = TRUE), 0),
      "-", fmt_num(max(age_values, na.rm = TRUE), 0)
    ),
    collapse_named_counts(gender_counts),
    paste0(fmt_num(min(sample_metadata$altitude, na.rm = TRUE), 0), "-", fmt_num(max(sample_metadata$altitude, na.rm = TRUE), 0), " m"),
    fmt_int(total_mags),
    fmt_int(primary_mags),
    paste0(fmt_int(excluded_or_failed), " excluded from primary family; ", fmt_int(model_failed_or_non_ok), " non-ok status"),
    paste0(fmt_int(nominal_mags), " nominal raw p < 0.05; ", fmt_int(fdr_mags), " BH-FDR q < 0.05")
  ),
  notes = c(
    "Unique Subject identifiers in cleaned_data.rds.",
    "Unique sampleID values after removing MAG-level duplication.",
    "Observed expedition stages represented in sample metadata.",
    "Counts are unique samples, not MAG rows.",
    paste0("Calculated as ", n_subjects, " observed subjects minus samples available at each timepoint."),
    "Calculated once per subject to avoid repeated-measure duplication.",
    "Calculated once per subject to avoid repeated-measure duplication.",
    "Altitude values in sample metadata; includes baseline/return altitude.",
    "Unique MAG identifiers in per-MAG result table.",
    "Rows with primary_analysis == TRUE; this is the BH-FDR testing family.",
    paste0("Status distribution: ", collapse_named_counts(status_counts), "."),
    paste0("Nominal candidates are exploratory; ", fdr_mags, " MAGs pass primary BH-FDR in the supplied results.")
  ),
  stringsAsFactors = FALSE
)

csv_path <- file.path(output_dir, paste0(output_stem, ".csv"))
md_path <- file.path(output_dir, paste0(output_stem, ".md"))
png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

utils::write.csv(table_1, csv_path, row.names = FALSE)

markdown_lines <- c(
  "# Table 1. Sample and Metadata Summary",
  "",
  "This table summarizes the repeated-measures study design and the MAG modeling set used for the full T1-T10 per-MAG lmerTest LMM workflow.",
  "",
  paste(names(table_1), collapse = " | "),
  paste(rep("---", ncol(table_1)), collapse = " | ")
)
for (i in seq_len(nrow(table_1))) {
  row_values <- vapply(table_1[i, ], as.character, character(1))
  row_values <- gsub("\\|", "/", row_values)
  row_values[is.na(row_values)] <- ""
  markdown_lines <- c(markdown_lines, paste(row_values, collapse = " | "))
}
writeLines(markdown_lines, md_path)

display_table <- table_1 |>
  dplyr::mutate(
    section = wrap_text(section, 17),
    metric = wrap_text(metric, 28),
    value = wrap_text(value, 36),
    notes = wrap_text(notes, 58)
  )

col_labels <- c("Section", "Metric", "Value", "Definition / note")
col_widths <- c(0.16, 0.25, 0.24, 0.35)
col_x <- cumsum(c(0, col_widths))

pal <- list(
  ink = "#111827",
  text = "#1F2937",
  muted = "#4B5563",
  hair = "#D1D5DB",
  faint = "#F8FAFC",
  header = "#E8F1F8",
  stripe = "#F8FAFC",
  blue = "#0072B2",
  orange = "#D55E00"
)

draw_table <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  grid::grid.text(
    "Table 1 | Sample and metadata summary",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.965, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 17, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    paste0(
      "Repeated-measures gut microbiome study across T1-T10 expedition stages. ",
      "Counts are derived from cleaned_data.rds and the full per-MAG lmerTest result table."
    ),
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.925, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9.4, col = pal$muted)
  )

  table_top <- 0.86
  table_bottom <- 0.07
  header_h <- 0.065
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
    x_right <- 0.025 + 0.95 * col_x[j + 1]
    grid::grid.text(
      col_labels[j],
      x = grid::unit(x_left + 0.006, "npc"),
      y = grid::unit(table_top - header_h / 2, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 9.0, fontface = "bold", col = pal$ink)
    )
  }

  section_seen <- character()
  for (i in seq_len(nrow(display_table))) {
    y_top <- table_top - header_h - (i - 1) * row_h
    y_mid <- y_top - row_h / 2
    fill <- if (i %% 2 == 0) "#FFFFFF" else pal$stripe
    grid::grid.rect(
      x = grid::unit(0.5, "npc"),
      y = grid::unit(y_mid, "npc"),
      width = grid::unit(0.95, "npc"),
      height = grid::unit(row_h, "npc"),
      gp = grid::gpar(fill = fill, col = "#E5E7EB", lwd = 0.45)
    )

    section_value <- as.character(display_table$section[[i]])
    original_section <- table_1$section[[i]]
    if (original_section %in% section_seen) {
      section_value <- ""
    } else {
      section_seen <- c(section_seen, original_section)
    }

    row_values <- c(section_value, display_table$metric[[i]], display_table$value[[i]], display_table$notes[[i]])
    for (j in seq_along(row_values)) {
      x_left <- 0.025 + 0.95 * col_x[j]
      fontface <- if (j == 1 && nzchar(row_values[j])) "bold" else "plain"
      color <- if (j == 1 && nzchar(row_values[j])) pal$blue else pal$text
      grid::grid.text(
        row_values[j],
        x = grid::unit(x_left + 0.006, "npc"),
        y = grid::unit(y_mid, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(fontsize = 7.65, col = color, fontface = fontface, lineheight = 0.88)
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
    "Abbreviations: MAG, metagenome-assembled genome; LMM, linear mixed model; BH-FDR, Benjamini-Hochberg false discovery rate.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.034, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.7, col = pal$muted)
  )
}

grDevices::png(filename = png_path, width = 12.8, height = 7.8, units = "in", res = 450, bg = "white", type = "cairo")
draw_table()
grDevices::dev.off()

grDevices::pdf(file = pdf_path, width = 12.8, height = 7.8, bg = "white", useDingbats = FALSE)
draw_table()
grDevices::dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = svg_path, width = 12.8, height = 7.8, bg = "white")
  draw_table()
  grDevices::dev.off()
} else {
  warning("Package 'svglite' is not installed; SVG output was skipped.", call. = FALSE)
}

message("Table 1 saved:")
message("- ", csv_path)
message("- ", md_path)
message("- ", png_path)
message("- ", pdf_path)
if (file.exists(svg_path)) {
  message("- ", svg_path)
}
