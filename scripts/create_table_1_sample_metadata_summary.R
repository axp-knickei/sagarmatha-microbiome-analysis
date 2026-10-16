# Create manuscript Table 1: sample and metadata summary for the gut MAG dataset.
#
# The table describes the gut abundance dataset used for alpha diversity and
# the available per-MAG LMM result table used to summarize MAG model inclusion.

required_packages <- c("dplyr", "grid", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

library(dplyr)

repo_root <- getwd()
metadata_path <- file.path(repo_root, "data", "raw", "sample_info.csv")
abundance_path <- file.path(repo_root, "data", "raw", "gut.csv")
taxonomy_path <- file.path(repo_root, "data", "raw", "GDTB.csv")
model_results_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "tables",
  "per_mag_lmertest_hpc_results.csv"
)

output_dir <- file.path(repo_root, "results", "manuscript_tables")
output_stem <- "table_1_sample_metadata_summary"

required_files <- c(metadata_path, abundance_path, taxonomy_path, model_results_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- utils::read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
abundance <- utils::read.csv(abundance_path, check.names = FALSE, stringsAsFactors = FALSE)
taxonomy <- utils::read.csv(taxonomy_path, check.names = FALSE, stringsAsFactors = FALSE)
model_results <- utils::read.csv(model_results_path, check.names = FALSE, stringsAsFactors = FALSE)

required_metadata_columns <- c("sampleID", "group", "Timepoint", "date", "gender", "age", "altitude")
required_model_columns <- c("MAG", "status", "primary_analysis")
missing_metadata_columns <- setdiff(required_metadata_columns, names(metadata))
missing_model_columns <- setdiff(required_model_columns, names(model_results))
if (!"sampleID" %in% names(abundance)) {
  stop("Abundance table must contain a sampleID column.", call. = FALSE)
}
if (length(missing_metadata_columns) > 0) {
  stop("Metadata table is missing required columns: ", paste(missing_metadata_columns, collapse = ", "), call. = FALSE)
}
if (length(missing_model_columns) > 0) {
  stop("Model result table is missing required columns: ", paste(missing_model_columns, collapse = ", "), call. = FALSE)
}

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_num <- function(x, digits = 1) format(round(x, digits), nsmall = digits, trim = TRUE)
fmt_pct <- function(x, digits = 0) paste0(fmt_num(100 * x, digits), "%")

natural_timepoint_order <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  numeric_tp <- x[grepl("^T[0-9]+$", x)]
  other_tp <- setdiff(x, numeric_tp)
  c(numeric_tp[order(as.integer(sub("^T", "", numeric_tp)))], sort(other_tp))
}

collapse_named_counts <- function(counts) {
  paste(paste0(names(counts), ": ", as.integer(counts)), collapse = "; ")
}

wrap_text <- function(x, width) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

metadata$sampleID <- as.character(metadata$sampleID)
abundance$sampleID <- as.character(abundance$sampleID)
matched_metadata <- metadata %>%
  filter(sampleID %in% abundance$sampleID)

if (nrow(matched_metadata) == 0) {
  stop("No sampleID values match between metadata and gut abundance table.", call. = FALSE)
}

sample_metadata <- matched_metadata %>%
  distinct(sampleID, group, Timepoint, date, gender, age, altitude)
subject_metadata <- sample_metadata %>%
  arrange(group, sampleID) %>%
  distinct(group, .keep_all = TRUE)

timepoint_levels <- natural_timepoint_order(sample_metadata$Timepoint)
samples_per_timepoint <- table(factor(sample_metadata$Timepoint, levels = timepoint_levels))
n_subjects <- n_distinct(sample_metadata$group)
missing_per_timepoint <- pmax(n_subjects - as.integer(samples_per_timepoint), 0L)
names(missing_per_timepoint) <- names(samples_per_timepoint)

age_values <- suppressWarnings(as.numeric(subject_metadata$age))
gender_counts <- table(subject_metadata$gender, useNA = "ifany")
altitude_values <- suppressWarnings(as.numeric(sample_metadata$altitude))

total_mags <- length(setdiff(names(abundance), "sampleID"))
taxonomy_matches <- if ("Genome" %in% names(taxonomy)) sum(setdiff(names(abundance), "sampleID") %in% taxonomy$Genome) else NA_integer_
modeled_mags <- sum(model_results$primary_analysis %in% TRUE, na.rm = TRUE)
model_result_mags <- n_distinct(model_results$MAG)
excluded_or_failed <- total_mags - modeled_mags
status_counts <- table(model_results$status, useNA = "ifany")

table_1 <- data.frame(
  section = c(
    rep("Study design", 5),
    rep("Participant metadata", 3),
    rep("MAG analysis set", 3)
  ),
  metric = c(
    "Number of subjects",
    "Number of gut samples",
    "Number of timepoints",
    "Samples per timepoint",
    "Missing samples per timepoint",
    "Age summary",
    "Gender distribution",
    "Altitude range",
    "Total MAGs",
    "MAGs modeled",
    "MAGs excluded or failed"
  ),
  value = c(
    fmt_int(n_subjects),
    fmt_int(n_distinct(sample_metadata$sampleID)),
    fmt_int(length(timepoint_levels)),
    collapse_named_counts(samples_per_timepoint),
    collapse_named_counts(missing_per_timepoint),
    paste0(
      "median ", fmt_num(stats::median(age_values, na.rm = TRUE), 1),
      " years; IQR ", fmt_num(stats::quantile(age_values, 0.25, na.rm = TRUE), 1),
      "-", fmt_num(stats::quantile(age_values, 0.75, na.rm = TRUE), 1),
      "; range ", fmt_num(min(age_values, na.rm = TRUE), 0),
      "-", fmt_num(max(age_values, na.rm = TRUE), 0)
    ),
    paste0(
      paste(
        paste0(names(gender_counts), ": ", as.integer(gender_counts), " (", fmt_pct(as.integer(gender_counts) / sum(gender_counts), 0), ")"),
        collapse = "; "
      )
    ),
    paste0(fmt_num(min(altitude_values, na.rm = TRUE), 0), "-", fmt_num(max(altitude_values, na.rm = TRUE), 0), " m"),
    paste0(fmt_int(total_mags), if (!is.na(taxonomy_matches)) paste0(" (", fmt_int(taxonomy_matches), " matched to GTDB taxonomy)") else ""),
    fmt_int(modeled_mags),
    paste0(fmt_int(excluded_or_failed), " total; status: ", collapse_named_counts(status_counts))
  ),
  note = c(
    "Unique subject IDs in metadata records matched to gut abundance samples.",
    "Unique sampleID values present in both sample_info.csv and gut.csv.",
    "Observed timepoint labels in matched gut samples, including auxiliary labels when present.",
    "Counts are unique gut samples, not MAG-level long-format rows.",
    paste0("Calculated as ", n_subjects, " possible subjects minus observed samples at each timepoint."),
    "Age summarized once per subject to avoid repeated-measure duplication.",
    "Gender summarized once per subject to avoid repeated-measure duplication.",
    "Range across matched gut sample metadata.",
    "MAG feature columns in gut.csv; taxonomy match uses GDTB.csv Genome identifiers.",
    "Primary modeled MAGs are rows with primary_analysis == TRUE in per_mag_lmertest_hpc_results.csv.",
    "Excluded/failed count is total gut MAGs minus primary modeled MAGs; status distribution is from the per-MAG result table."
  ),
  stringsAsFactors = FALSE
)

csv_path <- file.path(output_dir, paste0(output_stem, ".csv"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))
png_path <- file.path(output_dir, paste0(output_stem, ".png"))

utils::write.csv(table_1, csv_path, row.names = FALSE)

display_table <- table_1 %>%
  mutate(
    section = wrap_text(section, 18),
    metric = wrap_text(metric, 28),
    value = wrap_text(value, 44),
    note = wrap_text(note, 64)
  )

col_labels <- c("Section", "Metric", "Value", "Definition / note")
col_widths <- c(0.15, 0.23, 0.27, 0.35)
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
    "Table 1 | Sample and metadata summary",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.965, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 17, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    "Matched gut MAG abundance samples, participant metadata, and per-MAG model inclusion summary.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.925, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 9.4, col = pal$muted)
  )

  table_top <- 0.865
  table_bottom <- 0.078
  header_h <- 0.062
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
      gp = grid::gpar(fontsize = 9.1, fontface = "bold", col = pal$ink)
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

    section_value <- display_table$section[[i]]
    if (table_1$section[[i]] %in% section_seen) {
      section_value <- ""
    } else {
      section_seen <- c(section_seen, table_1$section[[i]])
    }

    row_values <- c(section_value, display_table$metric[[i]], display_table$value[[i]], display_table$note[[i]])
    for (j in seq_along(row_values)) {
      x_left <- 0.025 + 0.95 * col_x[j]
      fontface <- if (j == 1 && nzchar(row_values[j])) "bold" else "plain"
      color <- if (j == 1 && nzchar(row_values[j])) pal$blue else pal$text
      grid::grid.text(
        row_values[j],
        x = grid::unit(x_left + 0.006, "npc"),
        y = grid::unit(y_mid, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(fontsize = 7.8, col = color, fontface = fontface, lineheight = 0.88)
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
    "Abbreviations: GTDB, Genome Taxonomy Database; IQR, interquartile range; MAG, metagenome-assembled genome.",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.038, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.7, col = pal$muted)
  )
}

grDevices::png(filename = png_path, width = 12.8, height = 7.4, units = "in", res = 450, bg = "white", type = "cairo")
draw_table()
grDevices::dev.off()

svglite::svglite(file = svg_path, width = 12.8, height = 7.4, bg = "white")
draw_table()
grDevices::dev.off()

message("Table 1 saved:")
message("- ", csv_path)
message("- ", svg_path)
message("- ", png_path)
