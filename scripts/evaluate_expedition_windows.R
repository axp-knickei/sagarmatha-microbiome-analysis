# ============================================
# evaluate_expedition_windows.R
# ============================================
# Metadata-only evaluation of biologically predefined expedition windows for
# secondary longitudinal MAG abundance analysis. This script does not modify
# raw data and does not run MAG abundance models.

has_yaml <- requireNamespace("yaml", quietly = TRUE)
execution_start_time <- Sys.time()

execution_record_helper <- file.path("scripts", "utils", "write_hpc_execution_record.R")
if (file.exists(execution_record_helper)) {
  tryCatch(
    source(execution_record_helper),
    error = function(e) warning("Could not source HPC execution recorder; continuing metadata evaluation: ", conditionMessage(e), call. = FALSE)
  )
}

is_absolute_path <- function(path) {
  grepl("^/", path) || grepl("^[A-Za-z]:[/\\\\]", path)
}

translate_windows_path <- function(path) {
  if (.Platform$OS.type != "unix" || !grepl("^[A-Za-z]:[/\\\\]", path)) {
    return(path)
  }

  drive <- tolower(substr(path, 1, 1))
  rest <- substring(path, 3)
  rest <- gsub("\\\\", "/", rest)
  file.path("/mnt", drive, rest)
}

resolve_path <- function(path) {
  path <- translate_windows_path(path)
  if (is_absolute_path(path)) path else file.path(getwd(), path)
}

env_value <- function(...) {
  keys <- c(...)
  for (key in keys) {
    value <- Sys.getenv(key, unset = NA_character_)
    if (!is.na(value) && nzchar(value)) {
      return(value)
    }
  }
  NULL
}

parse_date_value <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }
  x_chr <- as.character(x)
  parsed <- as.Date(x_chr)
  need_compact <- is.na(parsed) & grepl("^[0-9]{8}$", x_chr)
  parsed[need_compact] <- as.Date(x_chr[need_compact], format = "%Y%m%d")
  parsed
}

median_date <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(as.Date(NA))
  }
  as.Date(stats::median(as.numeric(x)), origin = "1970-01-01")
}

timepoint_sort_key <- function(x) {
  ifelse(grepl("^T[0-9]+$", x), as.integer(sub("^T", "", x)), NA_integer_)
}

format_date <- function(x) {
  ifelse(is.na(x), NA_character_, format(as.Date(x), "%Y-%m-%d"))
}

collapse_timepoints <- function(x) {
  paste(x, collapse = " -> ")
}

altitude_note <- function(n_unique, min_alt, max_alt) {
  if (is.na(n_unique) || n_unique == 0) {
    return("altitude missing")
  }
  if (n_unique == 1) {
    return("consistent altitude within timepoint")
  }
  sprintf("variable altitude within timepoint (%s unique values; range %s m)", n_unique, max_alt - min_alt)
}

transition_direction <- function(delta_altitude) {
  if (is.na(delta_altitude)) {
    return("unknown")
  }
  if (delta_altitude >= 500) {
    return("ascent")
  }
  if (delta_altitude <= -500) {
    return("descent")
  }
  "plateau/slight change"
}

transition_interpretation <- function(from_tp, to_tp, day_gap, alt_delta, direction) {
  gap_text <- ifelse(is.na(day_gap), "unknown date gap", paste0(day_gap, " day gap"))
  alt_text <- ifelse(is.na(alt_delta), "unknown altitude change", paste0(ifelse(alt_delta > 0, "+", ""), alt_delta, " m"))
  if (direction == "ascent") {
    return(sprintf("%s from %s to %s with altitude increase (%s).", gap_text, from_tp, to_tp, alt_text))
  }
  if (direction == "descent") {
    return(sprintf("%s from %s to %s with altitude decrease (%s).", gap_text, from_tp, to_tp, alt_text))
  }
  sprintf("%s from %s to %s with stable or modest altitude change (%s).", gap_text, from_tp, to_tp, alt_text)
}

timepoint_stage_context <- function() {
  data.frame(
    Timepoint = paste0("T", 1:10),
    corrected_stage_label = c(
      "baseline_low_altitude_pre_expedition",
      "early_ascent_first_mountain_exposure",
      "higher_ascent_near_base_camp_level",
      "intermediate_or_lower_altitude_after_T3",
      "high_climb_exposure_around_6500m_not_summit",
      "short_interval_post_T5_descent_sample",
      "pre_summit_high_altitude_staging",
      "post_summit_post_extreme_altitude_sample",
      "return_to_low_altitude_post_expedition_descent",
      "later_low_altitude_follow_up_recovery_like"
    ),
    figure_informed_exposure_context = c(
      "Baseline or low-altitude pre-expedition sampling.",
      "Early ascent or first mountain exposure.",
      "Higher ascent phase, near or approaching Base Camp-level altitude.",
      "Sampling altitude is lower/intermediate relative to surrounding expedition movement; interpret with date and trajectory context.",
      "High-altitude sampling around 6500 m; important high-climb exposure but not the summit.",
      "Short-interval post-T5 sampling after descent from the T5 high-altitude point.",
      "High-altitude staging before summit push, likely around Base Camp-level altitude.",
      "Post-summit or post-extreme-altitude sampling; the major summit-push exposure occurs before/around this sample.",
      "Return to low altitude or post-expedition descent.",
      "Later low-altitude follow-up or recovery-like phase."
    ),
    stringsAsFactors = FALSE
  )
}

transition_exposure_context <- function() {
  data.frame(
    transition = paste0("T", 1:9, " -> T", 2:10),
    figure_informed_interpretation = c(
      "Transition from baseline/low-altitude sampling into early ascent or first mountain exposure.",
      "Continued ascent toward higher mountain or Base Camp-level exposure.",
      "Figure-informed context suggests expedition movement after higher ascent; do not infer biology from sampling altitude alone.",
      "High-climb exposure interval leading to T5 around 6500 m; T5 is important but is not the summit.",
      "Short-interval descent after the T5 high-altitude sampling point.",
      "Pre-summit high-altitude staging interval before the major summit-push period.",
      "Major summit-push interval: figure trajectory rises sharply toward 8848.86 m and descends before/around T8.",
      "Post-summit/post-extreme-altitude period followed by return to low altitude.",
      "Later low-altitude follow-up or recovery-like interval."
    ),
    figure_suggests_between_sample_exposure = c(
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE
    ),
    between_sample_exposure_note = c(
      "Exposure trajectory begins moving away from low-altitude baseline.",
      "Exposure trajectory continues upward beyond the prior sampling altitude.",
      "Figure context indicates expedition movement between samples; sampling altitude alone may understate exposure dynamics.",
      "High-climb exposure around 6500 m is represented by T5 sampling altitude, but this is not summit exposure.",
      "Descent between sampling points is biologically meaningful despite short date gap.",
      "Staging before summit push may matter even when sampling altitude remains around high-altitude/Base Camp level.",
      "Summit-push exposure to 8848.86 m is inferred from the figure, not encoded as a metadata sampling altitude.",
      "Sampling altitude drops to low altitude after prior extreme-altitude exposure.",
      "Sampling and exposure are both low-altitude/recovery-like in the available evidence."
    ),
    stringsAsFactors = FALSE
  )
}

count_complete_subjects <- function(sample_metadata, timepoints) {
  subset <- sample_metadata[sample_metadata$Timepoint %in% timepoints, , drop = FALSE]
  if (nrow(subset) == 0) {
    return(0L)
  }
  tab <- table(subset$Subject, subset$Timepoint)
  if (!all(timepoints %in% colnames(tab))) {
    return(0L)
  }
  sum(rowSums(tab[, timepoints, drop = FALSE] > 0) == length(timepoints))
}

window_caveat <- function(complete_subjects, min_subjects_per_timepoint, altitude_span, date_span_days) {
  caveats <- character(0)
  if (complete_subjects < 5) {
    caveats <- c(caveats, sprintf("only %s subjects have all three timepoints", complete_subjects))
  }
  if (min_subjects_per_timepoint < 5) {
    caveats <- c(caveats, sprintf("minimum per-timepoint subject count is %s", min_subjects_per_timepoint))
  }
  if (altitude_span >= 5000) {
    caveats <- c(caveats, "altitude and expedition date are strongly coupled across this window")
  }
  if (date_span_days >= 30) {
    caveats <- c(caveats, sprintf("wide date span of %s days", date_span_days))
  }
  caveats <- c(caveats, "sampling altitude is not the full exposure altitude trajectory")
  if (length(caveats) == 0) {
    return("standard longitudinal caveats: repeated measures, missing samples, and date-altitude coupling")
  }
  paste(caveats, collapse = "; ")
}

write_markdown_table <- function(df, columns) {
  df <- df[, columns, drop = FALSE]
  out <- character(0)
  out <- c(out, paste(names(df), collapse = " | "))
  out <- c(out, paste(rep("---", ncol(df)), collapse = " | "))
  for (i in seq_len(nrow(df))) {
    values <- vapply(df[i, , drop = FALSE], as.character, character(1))
    values[is.na(values)] <- ""
    out <- c(out, paste(values, collapse = " | "))
  }
  out
}

config_path <- file.path(getwd(), "config", "config.yml")
config <- list()
if (has_yaml && file.exists(config_path)) {
  config <- yaml::read_yaml(config_path)
}

configured_raw_dir <- config$data$raw_dir
if (is.null(configured_raw_dir) || !nzchar(configured_raw_dir)) {
  configured_raw_dir <- "data/raw"
}
metadata_filename <- config$data$metadata
if (is.null(metadata_filename) || !nzchar(metadata_filename)) {
  metadata_filename <- "sample_info.csv"
}
processed_path <- config$data$processed
if (is.null(processed_path) || !nzchar(processed_path)) {
  processed_path <- "data/processed/cleaned_data.rds"
}

raw_dir <- env_value("BLACKICE_RAW_DIR", "RAW_DATA_DIR")
if (is.null(raw_dir)) {
  raw_dir <- configured_raw_dir
}
metadata_path <- file.path(resolve_path(raw_dir), metadata_filename)
processed_path_resolved <- resolve_path(processed_path)

metadata_source <- NA_character_
metadata_source_note <- NA_character_

if (file.exists(metadata_path)) {
  sample_metadata <- utils::read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
  metadata_source <- metadata_path
  metadata_source_note <- "raw sample_info.csv"
} else if (file.exists(processed_path_resolved)) {
  processed <- readRDS(processed_path_resolved)
  required_processed_cols <- c("sampleID", "Timepoint", "date", "altitude")
  missing_processed_cols <- setdiff(required_processed_cols, names(processed))
  if (length(missing_processed_cols) > 0) {
    stop("Processed data is missing required metadata columns: ", paste(missing_processed_cols, collapse = ", "), call. = FALSE)
  }
  subject_col <- if ("Subject" %in% names(processed)) "Subject" else if ("group" %in% names(processed)) "group" else NA_character_
  if (is.na(subject_col)) {
    stop("Processed data does not contain Subject or group column.", call. = FALSE)
  }
  keep_cols <- unique(c("sampleID", subject_col, "group", "Subject", "Timepoint", "date", "date_converted", "gender", "age", "altitude"))
  keep_cols <- intersect(keep_cols, names(processed))
  sample_metadata <- unique(as.data.frame(processed[, keep_cols, drop = FALSE]))
  metadata_source <- processed_path_resolved
  metadata_source_note <- "unique sample-level metadata extracted from processed cleaned_data.rds because raw sample_info.csv was not available locally"
} else {
  stop("No raw metadata file or processed data file found.", call. = FALSE)
}

if (!"Subject" %in% names(sample_metadata)) {
  if ("group" %in% names(sample_metadata)) {
    sample_metadata$Subject <- sample_metadata$group
  } else {
    stop("Metadata does not contain group or Subject column.", call. = FALSE)
  }
}
if (!"group" %in% names(sample_metadata)) {
  sample_metadata$group <- sample_metadata$Subject
}

sample_metadata$sampleID <- as.character(sample_metadata$sampleID)
sample_metadata$Subject <- as.character(sample_metadata$Subject)
sample_metadata$group <- as.character(sample_metadata$group)
sample_metadata$Timepoint <- as.character(sample_metadata$Timepoint)
sample_metadata$altitude <- suppressWarnings(as.numeric(sample_metadata$altitude))
if ("date_converted" %in% names(sample_metadata)) {
  sample_metadata$date_parsed <- parse_date_value(sample_metadata$date_converted)
} else {
  sample_metadata$date_parsed <- parse_date_value(sample_metadata$date)
}

sample_metadata <- sample_metadata[!duplicated(sample_metadata$sampleID), , drop = FALSE]

output_dir <- file.path("results", "metadata_expedition_windows")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

timepoint_values <- sort(unique(sample_metadata$Timepoint), na.last = TRUE)
observed_timepoints <- timepoint_values[!is.na(timepoint_values)]
numeric_order <- timepoint_sort_key(observed_timepoints)

timepoint_summary <- do.call(
  rbind,
  lapply(observed_timepoints, function(tp) {
    rows <- sample_metadata[sample_metadata$Timepoint == tp, , drop = FALSE]
    min_alt <- min(rows$altitude, na.rm = TRUE)
    med_alt <- stats::median(rows$altitude, na.rm = TRUE)
    max_alt <- max(rows$altitude, na.rm = TRUE)
    if (!is.finite(min_alt)) min_alt <- NA_real_
    if (!is.finite(med_alt)) med_alt <- NA_real_
    if (!is.finite(max_alt)) max_alt <- NA_real_
    n_alt <- length(unique(rows$altitude[!is.na(rows$altitude)]))
    data.frame(
      Timepoint = tp,
      timepoint_numeric_order = timepoint_sort_key(tp),
      n_samples = length(unique(rows$sampleID)),
      n_subjects = length(unique(rows$Subject)),
      earliest_date = format_date(min(rows$date_parsed, na.rm = TRUE)),
      latest_date = format_date(max(rows$date_parsed, na.rm = TRUE)),
      median_date = format_date(median_date(rows$date_parsed)),
      min_altitude = min_alt,
      median_altitude = med_alt,
      max_altitude = max_alt,
      altitude_range = max_alt - min_alt,
      n_unique_altitude_values = n_alt,
      altitude_consistency_note = altitude_note(n_alt, min_alt, max_alt),
      stringsAsFactors = FALSE
    )
  })
)

timepoint_summary$median_date_parsed <- as.Date(timepoint_summary$median_date)
timepoint_summary <- timepoint_summary[order(timepoint_summary$median_date_parsed, timepoint_summary$median_altitude), , drop = FALSE]
timepoint_summary$chronological_order <- seq_len(nrow(timepoint_summary))
stage_context <- timepoint_stage_context()
timepoint_summary <- merge(timepoint_summary, stage_context, by = "Timepoint", all.x = TRUE, sort = FALSE)
timepoint_summary <- timepoint_summary[order(timepoint_summary$chronological_order), , drop = FALSE]

date_order_matches_timepoint_order <- TRUE
if (all(!is.na(timepoint_summary$timepoint_numeric_order))) {
  date_order_matches_timepoint_order <- identical(
    timepoint_summary$Timepoint,
    timepoint_summary$Timepoint[order(timepoint_summary$timepoint_numeric_order)]
  )
}

timepoint_summary_out <- timepoint_summary[, c(
  "chronological_order",
  "Timepoint",
  "timepoint_numeric_order",
  "n_samples",
  "n_subjects",
  "earliest_date",
  "latest_date",
  "median_date",
  "min_altitude",
  "median_altitude",
  "max_altitude",
  "altitude_range",
  "n_unique_altitude_values",
  "altitude_consistency_note",
  "corrected_stage_label",
  "figure_informed_exposure_context"
)]
observed_timepoints_chronological <- timepoint_summary_out$Timepoint

transitions <- data.frame()
if (nrow(timepoint_summary) >= 2) {
  transition_context <- transition_exposure_context()
  transition_rows <- vector("list", nrow(timepoint_summary) - 1)
  for (i in seq_len(nrow(timepoint_summary) - 1)) {
    from <- timepoint_summary[i, ]
    to <- timepoint_summary[i + 1, ]
    day_gap <- as.integer(as.Date(to$median_date) - as.Date(from$median_date))
    alt_delta <- to$median_altitude - from$median_altitude
    direction <- transition_direction(alt_delta)
    transition_rows[[i]] <- data.frame(
      from_timepoint = from$Timepoint,
      to_timepoint = to$Timepoint,
      transition = paste(from$Timepoint, "->", to$Timepoint),
      from_median_date = from$median_date,
      to_median_date = to$median_date,
      date_difference_days = day_gap,
      from_median_altitude = from$median_altitude,
      to_median_altitude = to$median_altitude,
      altitude_difference_m = alt_delta,
      sampling_altitude_direction = direction,
      sampling_altitude_interpretation = transition_interpretation(from$Timepoint, to$Timepoint, day_gap, alt_delta, direction),
      stringsAsFactors = FALSE
    )
  }
  transitions <- do.call(rbind, transition_rows)
  transitions <- merge(transitions, transition_context, by = "transition", all.x = TRUE, sort = FALSE)
  transitions <- transitions[match(paste(timepoint_summary$Timepoint[-nrow(timepoint_summary)], "->", timepoint_summary$Timepoint[-1]), transitions$transition), , drop = FALSE]
}

summary_by_tp <- timepoint_summary
rownames(summary_by_tp) <- summary_by_tp$Timepoint

window_specs <- list(
  list(
    label = "initial_ascent_to_high_altitude",
    timepoints = c("T1", "T2", "T3"),
    biological_interpretation = "baseline to early ascent/adaptation and first higher mountain exposure",
    exposure_interpretation = "Sampling altitude rises from 0 m to 3650 m and then 5200 m; figure context supports this as early expedition ascent rather than simple linear calendar time.",
    justification = "Dates and altitude move from 0 m at T1 to 3650 m at T2 and 5200 m at T3, with a large pre-expedition-to-expedition date gap.",
    response_detected = "MAGs that respond during initial expedition entry, ascent, and first high-altitude exposure.",
    evidence_type = "metadata_plus_figure",
    evidence_type_explanation = "Early ascent is supported by metadata sampling altitude increase from 0 m to 3650 m to 5200 m and by the figure trajectory."
  ),
  list(
    label = "ascent_and_intermediate_adjustment",
    timepoints = c("T2", "T3", "T4"),
    biological_interpretation = "early ascent, higher ascent, and intermediate adjustment after T3",
    exposure_interpretation = "Sampling altitude increases from 3650 m to 5200 m and remains 5200 m at T4, but figure context cautions that expedition movement between sampling points may not be a true plateau.",
    justification = "The window captures the first high-altitude transition and the following timepoint without using MAG results; T4 should be interpreted as intermediate adjustment rather than a simple plateau endpoint.",
    response_detected = "MAGs showing early high-altitude adaptation, delayed response to ascent, or adjustment after initial high-altitude exposure.",
    evidence_type = "metadata_plus_figure",
    evidence_type_explanation = "Ascent and intermediate adjustment are supported by sampling altitude/date metadata plus figure-informed movement context between samples."
  ),
  list(
    label = "high_climb_6500m_exposure",
    timepoints = c("T4", "T5", "T6"),
    biological_interpretation = "high-climb exposure around 6500 m and short-interval descent afterward",
    exposure_interpretation = "T5 samples around 6500 m, then T6 samples after descent to 5200 m two days later. T5 is a high-climb exposure but not the summit.",
    justification = "T5 is the only metadata sampling point at 6500 m and is flanked by 5200 m samples, defining a biologically defensible high-climb exposure window distinct from the later summit push.",
    response_detected = "MAGs sensitive to high-climb stress around 6500 m, immediate descent, or short-lag post-exposure response.",
    evidence_type = "metadata_plus_figure",
    evidence_type_explanation = "The T5 high-climb exposure around 6500 m is directly visible in sampling metadata and supported by figure context; it is not the summit."
  ),
  list(
    label = "pre_summit_staging_to_post_summit",
    timepoints = c("T6", "T7", "T8"),
    biological_interpretation = "post-T5 descent, pre-summit staging, and post-summit/post-extreme-altitude sampling",
    exposure_interpretation = "Sampling altitude is 5200 m at T6-T8, but figure context indicates the major summit-push exposure occurs between T7 and T8 and reaches 8848.86 m before descent around T8.",
    justification = "This replaces the oversimplified plateau label: T6-T8 are not chosen as flat movement, but because T7 to T8 brackets the figure-informed summit-push interval.",
    response_detected = "MAGs responding to summit-push preparation, extreme-altitude exposure, or early post-summit effects.",
    evidence_type = "figure_informed_critical_interval",
    evidence_type_explanation = "Sampling altitude appears unchanged at 5200 m across T6-T8, but the key biological exposure is the T7 -> T8 summit-push interval inferred from the figure trajectory."
  ),
  list(
    label = "summit_push_and_descent_recovery_transition",
    timepoints = c("T7", "T8", "T9"),
    biological_interpretation = "pre-summit staging, post-summit sampling, and return to low altitude",
    exposure_interpretation = "The biologically extreme interval is T7 -> T8 in the figure, followed by low-altitude sampling at T9. Sampling altitude alone shows 5200 m, 5200 m, 0 m.",
    justification = "This window directly captures the figure-informed summit-push interval and the subsequent descent/recovery transition without treating T7-T8 as a plateau.",
    response_detected = "MAGs that respond to extreme-altitude summit exposure, post-summit descent, or early low-altitude recovery.",
    evidence_type = "figure_informed_plus_sampling_descent",
    evidence_type_explanation = "The summit-push exposure is figure-informed, while the drop from T8 to T9 low altitude is directly supported by sampling metadata."
  ),
  list(
    label = "post_expedition_low_altitude_recovery",
    timepoints = c("T8", "T9", "T10"),
    biological_interpretation = "post-summit/post-extreme-altitude sampling followed by low-altitude recovery-like follow-up",
    exposure_interpretation = "T8 is interpreted as post-summit/post-extreme-altitude sampling, then T9 and T10 are low-altitude samples at 0 m.",
    justification = "The window tests whether post-extreme-altitude changes persist, rebound, or normalize after return to low altitude.",
    response_detected = "MAGs that rebound, persist, or normalize after return to low altitude.",
    evidence_type = "metadata_plus_figure",
    evidence_type_explanation = "Low-altitude post-expedition and recovery-like sampling is supported by metadata altitude at T9-T10 and figure context that T8 follows extreme exposure."
  )
)

windows <- data.frame()
for (spec in window_specs) {
  tps <- spec$timepoints
  if (!all(tps %in% timepoint_summary$Timepoint)) {
    next
  }
  rows <- summary_by_tp[tps, , drop = FALSE]
  date_start <- min(as.Date(rows$median_date), na.rm = TRUE)
  date_end <- max(as.Date(rows$median_date), na.rm = TRUE)
  min_alt <- min(rows$min_altitude, na.rm = TRUE)
  max_alt <- max(rows$max_altitude, na.rm = TRUE)
  complete_subjects <- count_complete_subjects(sample_metadata, tps)
  min_subjects <- min(rows$n_subjects, na.rm = TRUE)
  window_row <- data.frame(
    window_label = spec$label,
    timepoints = paste(tps, collapse = ","),
    forward_timepoint_path = collapse_timepoints(tps),
    date_span = paste(format_date(date_start), "to", format_date(date_end)),
    date_span_days = as.integer(date_end - date_start),
    altitude_span = paste0(min_alt, " to ", max_alt, " m"),
    sampling_altitude_span = paste0(min_alt, " to ", max_alt, " m"),
    sampling_altitude_span_m = max_alt - min_alt,
    n_complete_subjects_all_three_timepoints = complete_subjects,
    min_subjects_per_timepoint = min_subjects,
    lmm_sample_subject_support = ifelse(complete_subjects >= 5 && min_subjects >= 5, "adequate_for_secondary_LMM_screening", "limited_for_LMM_review_before_modeling"),
    evidence_type = spec$evidence_type,
    evidence_type_explanation = spec$evidence_type_explanation,
    biological_interpretation = spec$biological_interpretation,
    figure_informed_altitude_exposure_interpretation = spec$exposure_interpretation,
    why_justified = spec$justification,
    possible_mag_response = spec$response_detected,
    caveats = window_caveat(complete_subjects, min_subjects, max_alt - min_alt, as.integer(date_end - date_start)),
    stringsAsFactors = FALSE
  )
  windows <- rbind(windows, window_row)
}

missing_counts <- data.frame(
  field = c("Timepoint", "date", "altitude", "Subject"),
  missing_values = c(
    sum(is.na(sample_metadata$Timepoint) | !nzchar(sample_metadata$Timepoint)),
    sum(is.na(sample_metadata$date_parsed)),
    sum(is.na(sample_metadata$altitude)),
    sum(is.na(sample_metadata$Subject) | !nzchar(sample_metadata$Subject))
  ),
  stringsAsFactors = FALSE
)

subject_counts <- aggregate(sampleID ~ Subject, sample_metadata, function(x) length(unique(x)))
names(subject_counts) <- c("Subject", "n_samples")
subjects_with_few_samples <- subject_counts[subject_counts$n_samples < 3, , drop = FALSE]

raw_row_loss_report <- file.path("data", "processed", "row_loss_investigation.txt")
row_loss_text <- if (file.exists(raw_row_loss_report)) paste(readLines(raw_row_loss_report, warn = FALSE), collapse = "\n") else ""
ta_tb_evidence <- grepl("\\bTa\\b", row_loss_text) || grepl("\\bTb\\b", row_loss_text)
ta_tb_observed_in_inspected_metadata <- any(sample_metadata$Timepoint %in% c("Ta", "Tb"), na.rm = TRUE)

write.csv(timepoint_summary_out, file.path(output_dir, "timepoint_altitude_date_summary.csv"), row.names = FALSE)
write.csv(transitions, file.path(output_dir, "consecutive_timepoint_transitions.csv"), row.names = FALSE)
write.csv(windows, file.path(output_dir, "proposed_biological_windows.csv"), row.names = FALSE)

metadata_summary_lines <- c(
  sprintf("- Configured metadata path from config/data-prep: `%s`", metadata_path),
  sprintf("- Metadata source inspected: `%s`", metadata_source),
  sprintf("- Source note: %s", metadata_source_note),
  sprintf("- Rows inspected: %s", nrow(sample_metadata)),
  sprintf("- Unique sampleID: %s", length(unique(sample_metadata$sampleID))),
  sprintf("- Unique subjects from group/Subject: %s", length(unique(sample_metadata$Subject))),
  sprintf("- Observed Timepoint values in chronological order: %s", paste(observed_timepoints_chronological, collapse = ", ")),
  sprintf("- Date range: %s to %s", format_date(min(sample_metadata$date_parsed, na.rm = TRUE)), format_date(max(sample_metadata$date_parsed, na.rm = TRUE))),
  sprintf("- Altitude range: %s to %s m", min(sample_metadata$altitude, na.rm = TRUE), max(sample_metadata$altitude, na.rm = TRUE)),
  sprintf("- Ta/Tb observed in inspected metadata: %s", ifelse(ta_tb_observed_in_inspected_metadata, "yes", "no")),
  sprintf("- Ta/Tb evidence from row-loss report: %s", ifelse(ta_tb_evidence, "yes; raw metadata included Ta/Tb but current processed data excluded them", "not found")),
  sprintf("- Date order matches numeric T1-T10 order: %s", ifelse(date_order_matches_timepoint_order, "yes", "no"))
)

missing_lines <- paste0("- ", missing_counts$field, ": ", missing_counts$missing_values)
few_sample_lines <- if (nrow(subjects_with_few_samples) == 0) {
  "- No subjects have fewer than three inspected samples."
} else {
  paste0("- Subjects with fewer than three inspected samples: ", paste(paste0(subjects_with_few_samples$Subject, " (", subjects_with_few_samples$n_samples, ")"), collapse = ", "))
}

large_gaps <- transitions[!is.na(transitions$date_difference_days) & transitions$date_difference_days >= 14, , drop = FALSE]
large_gap_line <- if (nrow(large_gaps) == 0) {
  "- No consecutive median-date gaps of 14 days or more."
} else {
  paste0("- Large date gaps: ", paste(paste0(large_gaps$transition, " (", large_gaps$date_difference_days, " days)"), collapse = "; "))
}

variable_altitude_tps <- timepoint_summary_out[timepoint_summary_out$n_unique_altitude_values > 1, , drop = FALSE]
variable_alt_line <- if (nrow(variable_altitude_tps) == 0) {
  "- Altitude is consistent within every inspected timepoint."
} else {
  paste0("- Timepoints with within-timepoint altitude variation: ", paste(variable_altitude_tps$Timepoint, collapse = ", "))
}
summit_push_line <- "- Figure-informed prior: T7 -> T8 contains the major summit-push exposure toward 8848.86 m and descent before/around T8; this exposure is not directly encoded as metadata sampling altitude."
stage_lines <- paste0(
  "- ",
  timepoint_summary_out$Timepoint,
  ": ",
  timepoint_summary_out$corrected_stage_label,
  " — ",
  timepoint_summary_out$figure_informed_exposure_context
)

md <- c(
  "# Expedition Window Metadata Evaluation",
  "",
  sprintf("Report written: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "Evaluate project metadata and propose a small number of biologically predefined forward three-timepoint windows for secondary longitudinal MAG abundance analysis. Window choices are based on date and altitude structure, not MAG p-values.",
  "",
  "## Data Source Inspected",
  "",
  metadata_summary_lines,
  "",
  "Raw data were not modified. When raw `sample_info.csv` was unavailable locally, the script used unique sample-level metadata from `data/processed/cleaned_data.rds`.",
  "",
  "## Metadata Summary",
  "",
  "Missing values in key metadata fields:",
  "",
  missing_lines,
  "",
  few_sample_lines,
  "",
  "## Sampling Altitude Versus Exposure Trajectory",
  "",
  "This report distinguishes two altitude concepts:",
  "",
  "- Sampling altitude: the altitude recorded in metadata at the sampling timepoint.",
  "- Exposure altitude trajectory: the altitude path between sampling points shown by the supporting expedition figure.",
  "",
  "A sampling timepoint can therefore have a moderate recorded altitude even if the interval immediately before it included more extreme altitude exposure. The key correction is that T6-T8 should not be interpreted as a simple 5200 m plateau: the figure-informed prior places the major summit-push interval between T7 and T8, with exposure toward 8848.86 m before descent around T8.",
  "",
  "## Timepoint Altitude/Date Summary",
  "",
  write_markdown_table(
    timepoint_summary_out,
    c("Timepoint", "n_samples", "n_subjects", "earliest_date", "latest_date", "median_date", "min_altitude", "median_altitude", "max_altitude", "n_unique_altitude_values", "altitude_consistency_note", "corrected_stage_label")
  ),
  "",
  "## Consecutive Transition Interpretation",
  "",
  large_gap_line,
  variable_alt_line,
  summit_push_line,
  "",
  write_markdown_table(
    transitions,
    c("transition", "date_difference_days", "altitude_difference_m", "sampling_altitude_direction", "sampling_altitude_interpretation", "figure_informed_interpretation", "figure_suggests_between_sample_exposure", "between_sample_exposure_note")
  ),
  "",
  "## Corrected Expedition-Stage Interpretation",
  "",
  stage_lines,
  "",
  "Corrections relative to the first draft:",
  "",
  "- T6-T8 is not labeled as a simple high-altitude plateau or flat movement.",
  "- T5 is labeled as high-climb exposure around 6500 m, not the summit.",
  "- T7 -> T8 is marked as the summit-push interval containing figure-informed 8848.86 m exposure.",
  "- T8 is labeled post-summit/post-extreme-altitude sampling.",
  "",
  "## Proposed Class 1 Biologically Predefined Windows",
  "",
  "These windows are predefined from expedition chronology, metadata sampling altitude, and figure-informed altitude-exposure trajectory. They should be treated as secondary biological tests, not data-mined p-value windows.",
  "Evidence type labels are included to make clear which windows are directly supported by sampling metadata and which depend partly on the figure-informed altitude exposure trajectory.",
  "",
  "## Evidence Type Definitions",
  "",
  "- `metadata_only`: window definition is supported only by sampling metadata, without needing figure trajectory context.",
  "- `figure_informed`: window definition depends mainly on the expedition trajectory figure.",
  "- `metadata_plus_figure`: window definition is supported by both metadata sampling altitude/date and the figure trajectory.",
  "- `figure_informed_critical_interval`: window contains a biologically critical between-sample exposure that is not directly encoded as sampling altitude.",
  "- `figure_informed_plus_sampling_descent`: window combines a figure-informed extreme exposure interval with a metadata-supported descent or recovery transition.",
  "",
  write_markdown_table(
    windows,
    c("window_label", "forward_timepoint_path", "evidence_type", "evidence_type_explanation", "date_span", "sampling_altitude_span", "n_complete_subjects_all_three_timepoints", "lmm_sample_subject_support", "biological_interpretation", "figure_informed_altitude_exposure_interpretation", "why_justified", "possible_mag_response", "caveats")
  ),
  "",
  "## Optional Exploratory-Window Note",
  "",
  "All possible forward three-timepoint windows may be useful only as exploratory discovery scanning. If used, inference must correct globally across all MAG x window tests, not separately within each window, to avoid p-value hunting.",
  "",
  "## Limitations",
  "",
  "- `group` is interpreted as `Subject`; no treatment group is assumed.",
  "- Ta and Tb were excluded from the current processed dataset because `scripts/01_data_prep.R` coerces `Timepoint` to T1-T10 before filtering. The row-loss report shows Ta at 3800 m on 2024-04-16 and Tb at 5800 m on 2024-04-23 in raw metadata.",
  "- Some subjects have missing timepoints, and some proposed windows have few complete subjects across all three timepoints.",
  "- Date and altitude are strongly confounded in ascent, peak-exposure, descent, and recovery-like windows.",
  "- Sampling altitude and exposure altitude trajectory differ conceptually; T7 -> T8 summit-push exposure is inferred from the supporting figure and is not directly encoded as a metadata altitude value.",
  "- T8 is interpreted as post-summit/post-extreme-altitude sampling from figure context, even though the metadata sampling altitude is 5200 m.",
  "- Per-timepoint altitude appears consistent in the inspected processed metadata, but this should be rechecked if raw metadata are restored or Ta/Tb are included.",
  "- Windows with fewer than five complete subjects should be considered unstable for LMM-based inference.",
  "",
  "## Recommended Next Modeling Step",
  "",
  "Use these Class 1 predefined windows for secondary MAG abundance analysis only after finalizing whether Ta/Tb should be included and confirming the figure-informed summit-push interpretation. For each window, subset to subjects with enough repeated observations, fit the same response-scale decision framework, and control FDR across all MAGs and all predefined windows tested."
)

writeLines(md, file.path(output_dir, "proposed_biological_windows.md"))

cat("Metadata expedition window evaluation complete\n")
cat("============================================\n")
cat("Timepoints found:", length(observed_timepoints), "\n")
cat("Observed timepoints:", paste(observed_timepoints_chronological, collapse = ", "), "\n")
cat("Corrected stage labels:\n")
for (i in seq_len(nrow(timepoint_summary_out))) {
  cat(" - ", timepoint_summary_out$Timepoint[[i]], ": ", timepoint_summary_out$corrected_stage_label[[i]], "\n", sep = "")
}
cat("Proposed biological windows:", paste(windows$window_label, collapse = ", "), "\n")
cat("Evidence types assigned:\n")
for (i in seq_len(nrow(windows))) {
  cat(" - ", windows$forward_timepoint_path[[i]], ": ", windows$evidence_type[[i]], "\n", sep = "")
}
cat("Output files:\n")
cat(" - ", file.path(output_dir, "timepoint_altitude_date_summary.csv"), "\n", sep = "")
cat(" - ", file.path(output_dir, "consecutive_timepoint_transitions.csv"), "\n", sep = "")
cat(" - ", file.path(output_dir, "proposed_biological_windows.csv"), "\n", sep = "")
cat(" - ", file.path(output_dir, "proposed_biological_windows.md"), "\n", sep = "")
cat("Warnings for review:\n")
if (!file.exists(metadata_path)) {
  cat(" - Raw metadata not available locally; used processed cleaned_data.rds metadata.\n")
}
if (ta_tb_evidence && !ta_tb_observed_in_inspected_metadata) {
  cat(" - Ta/Tb exist in raw metadata evidence but are excluded from processed data.\n")
}
if (nrow(subjects_with_few_samples) > 0) {
  cat(" - Subjects with fewer than three samples:", paste(subjects_with_few_samples$Subject, collapse = ", "), "\n")
}
if (any(windows$n_complete_subjects_all_three_timepoints < 5)) {
  cat(" - Some proposed windows have fewer than five complete subjects.\n")
}
if (any(windows$sampling_altitude_span_m >= 5000)) {
  cat(" - Some windows strongly confound expedition date and altitude.\n")
}
cat(" - Sampling altitude is not the full exposure altitude trajectory; T7 -> T8 summit-push exposure is figure-informed and not directly encoded in metadata.\n")

tryCatch(
  {
    if (exists("write_hpc_execution_record", mode = "function")) {
      write_hpc_execution_record(
        analysis_name = "metadata_expedition_window_definition",
        script_path = "scripts/evaluate_expedition_windows.R",
        result_dir = output_dir,
        input_data_file = metadata_source,
        window_definition_file = NA_character_,
        config_files = c("config/config.yml")[file.exists("config/config.yml")],
        number_of_tests = NA,
        fdr_family = "NA; metadata-only window definition, no MAG hypothesis testing",
        main_results_file = file.path(output_dir, "proposed_biological_windows.csv"),
        summary_log_file = file.path(output_dir, "proposed_biological_windows.md"),
        session_info_file = NA_character_,
        notes = c(
          "Metadata-only expedition window definition; no MAG model fitting.",
          paste0("Metadata source note: ", metadata_source_note),
          paste0("Proposed Class 1 windows: ", nrow(windows))
        ),
        start_time = execution_start_time
      )
    }
  },
  error = function(e) warning("Failed to write HPC execution record; continuing metadata evaluation: ", conditionMessage(e), call. = FALSE)
)
