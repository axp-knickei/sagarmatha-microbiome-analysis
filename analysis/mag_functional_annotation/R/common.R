suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

required_packages <- function(extra = character()) {
  unique(c("data.table", "yaml", "ggplot2", "scales", "svglite", "Matrix", extra))
}

assert_packages <- function(extra = character()) {
  packages <- required_packages(extra)
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(packages)
}

parse_stage_args <- function(default_config = "config/mag_functional_annotation.yml") {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(config = default_config, force = FALSE)
  i <- 1L
  while (i <= length(args)) {
    if (args[[i]] == "--config" && i < length(args)) {
      out$config <- args[[i + 1L]]
      i <- i + 2L
    } else if (args[[i]] == "--force") {
      out$force <- TRUE
      i <- i + 1L
    } else {
      stop("Unknown argument: ", args[[i]], call. = FALSE)
    }
  }
  out
}

read_magfunc_config <- function(path) {
  if (!file.exists(path)) stop("Configuration not found: ", path, call. = FALSE)
  cfg <- yaml::read_yaml(path)
  required <- c("inputs", "columns", "identifiers", "quality_control", "analysis", "ontologies", "outputs")
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0) stop("Configuration is missing sections: ", paste(missing, collapse = ", "), call. = FALSE)
  cfg$.config_path <- normalizePath(path, mustWork = TRUE)
  cfg
}

resolve_input_paths <- function(cfg) {
  paths <- cfg$inputs
  overrides <- cfg$input_environment_overrides %||% list()
  for (key in names(paths)) {
    env_name <- overrides[[key]] %||% ""
    if (nzchar(env_name)) {
      env_value <- Sys.getenv(env_name, unset = "")
      if (nzchar(env_value)) paths[[key]] <- env_value
    }
  }
  paths
}

resolve_output_dir <- function(cfg) {
  Sys.getenv("MAGFUNC_OUTPUT_DIR", unset = cfg$outputs$output_dir)
}

output_paths <- function(cfg, create = TRUE) {
  root <- resolve_output_dir(cfg)
  paths <- list(
    root = root,
    qc = file.path(root, "00_quality_control"),
    convergence = file.path(root, "01_functional_convergence"),
    unknown = file.path(root, "02_unknown_by_cluster_size"),
    taxonomy = file.path(root, "03_taxonomic_distribution"),
    longitudinal = file.path(root, "04_longitudinal_functional_potential"),
    hysteresis = file.path(root, "05_ascent_descent_hysteresis"),
    targeted_validation = file.path(root, "06_targeted_candidate_validation"),
    exploratory_followup = file.path(root, "07_pathway_unknown_exploration"),
    carrier_decomposition = file.path(root, "08_shared_unique_carrier_decomposition"),
    shared_component_followup = file.path(root, "09_shared_component_followup"),
    candidate_curation = file.path(root, "10_candidate_annotation_audit"),
    figures = file.path(root, "figures"),
    tables = file.path(root, "tables"),
    intermediate = file.path(root, "intermediate"),
    logs = file.path(root, "logs"),
    reproducibility = file.path(root, "reproducibility")
  )
  if (create) invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
  paths
}

guard_stage <- function(marker, args, cfg) {
  overwrite <- isTRUE(args$force) || isTRUE(cfg$outputs$overwrite)
  if (file.exists(marker) && !overwrite) {
    stop("Stage completion marker already exists: ", marker, ". Use --force only when intentional.", call. = FALSE)
  }
  if (file.exists(marker) && overwrite) {
    removal_status <- unlink(marker)
    if (removal_status != 0L || file.exists(marker)) {
      stop("Could not remove the previous stage completion marker: ", marker, call. = FALSE)
    }
  }
  invisible(overwrite)
}

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data.table::as.data.table(x), path, na = "")
  invisible(path)
}

read_table <- function(path, ...) {
  if (!file.exists(path)) stop("Required table not found: ", path, call. = FALSE)
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) return(data.table::fread(path, ...))

  gzip <- unname(Sys.which("gzip"))
  if (!nzchar(gzip)) {
    stop("The gzip command is required to read compressed result tables: ", path, call. = FALSE)
  }
  command <- paste(shQuote(gzip), "-cd", shQuote(normalizePath(path, mustWork = TRUE)))
  data.table::fread(cmd = command, ...)
}

save_object <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path, compress = "gzip")
  invisible(path)
}

read_object <- function(path) {
  if (!file.exists(path)) stop("Required intermediate file not found: ", path, call. = FALSE)
  readRDS(path)
}

parse_mag_id <- function(x, pattern) {
  x <- as.character(x)
  matched <- grepl(pattern, x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  out[matched] <- sub(pattern, "\\1", x[matched], perl = TRUE)
  out
}

count_file_lines <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  output <- tryCatch(system2("wc", c("-l", path), stdout = TRUE, stderr = FALSE), error = function(e) character())
  if (length(output) == 0) return(NA_real_)
  suppressWarnings(as.numeric(strsplit(trimws(output[[1]]), "[[:space:]]+")[[1]][[1]]))
}

annotation_header <- function(path, max_lines = 500L) {
  lines <- readLines(path, n = max_lines, warn = FALSE)
  index <- grep("^#?query\\t", lines, perl = TRUE)
  if (length(index) == 0) stop("Could not find the eggNOG annotation header beginning with query or #query.", call. = FALSE)
  index <- index[[1]]
  names <- strsplit(sub("^#", "", lines[[index]]), "\\t", fixed = FALSE)[[1]]
  list(line = index, names = names)
}

read_annotation <- function(path, cfg, select_columns = NULL, nrows = Inf) {
  header <- annotation_header(path)
  if (is.null(select_columns)) select_columns <- header$names
  absent <- setdiff(select_columns, header$names)
  if (length(absent) > 0) stop("Annotation file is missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  select_index <- match(select_columns, header$names)
  dt <- withCallingHandlers(
    data.table::fread(
      path,
      sep = "\t",
      quote = "",
      skip = header$line,
      header = FALSE,
      col.names = header$names[select_index],
      select = select_index,
      nrows = nrows,
      fill = TRUE,
      na.strings = unlist(cfg$quality_control$missing_tokens),
      showProgress = interactive()
    ),
    warning = function(w) {
      stop("Annotation parsing warning: ", conditionMessage(w), call. = FALSE)
    }
  )
  query_col <- cfg$columns$annotation_query
  dt <- dt[!is.na(get(query_col)) & !startsWith(get(query_col), "#")]
  dt
}

read_cluster <- function(path, nrows = Inf) {
  dt <- data.table::fread(
    path,
    sep = "\t",
    quote = "",
    header = FALSE,
    nrows = nrows,
    fill = FALSE,
    showProgress = interactive()
  )
  if (ncol(dt) != 2L) stop("Cluster file must have exactly two tab-separated columns; observed ", ncol(dt), ".", call. = FALSE)
  data.table::setnames(dt, c("representative_id", "member_id"))
  dt
}

ontology_specs <- function(cfg) {
  lapply(cfg$ontologies, function(x) {
    list(
      id = as.character(x$id),
      column = as.character(x$column),
      split_mode = as.character(x$split_mode),
      strip_prefix = as.character(x$strip_prefix %||% "")
    )
  })
}

split_annotation_terms <- function(annotation, query_col, spec) {
  values <- annotation[, .(query = get(query_col), value = get(spec$column))]
  values <- values[!is.na(value) & nzchar(trimws(value))]
  if (nrow(values) == 0) return(data.table(query = character(), term = character()))
  if (identical(spec$split_mode, "characters")) {
    out <- values[, .(term = unlist(strsplit(value, "", fixed = TRUE))), by = query]
  } else if (identical(spec$split_mode, "comma")) {
    out <- values[, .(term = unlist(strsplit(value, ",", fixed = TRUE))), by = query]
  } else {
    stop("Unsupported split_mode for ", spec$id, ": ", spec$split_mode, call. = FALSE)
  }
  out[, term := trimws(term)]
  if (nzchar(spec$strip_prefix)) out[, term := sub(spec$strip_prefix, "", term, perl = TRUE)]
  out <- out[!is.na(term) & nzchar(term) & term != "-"]
  unique(out, by = c("query", "term"))
}

register_results <- function(cfg, records) {
  paths <- output_paths(cfg)
  registry_path <- file.path(paths$intermediate, "results_registry.csv")
  records <- data.table::as.data.table(records)
  required <- c("result_number", "category", "question", "artifact_type", "path", "description", "stage")
  absent <- setdiff(required, names(records))
  if (length(absent) > 0) stop("Result registry records are missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  existing <- if (file.exists(registry_path)) data.table::fread(registry_path) else data.table()
  combined <- data.table::rbindlist(list(existing, records), fill = TRUE, use.names = TRUE)
  combined <- unique(combined, by = "path", fromLast = TRUE)
  data.table::setorder(combined, result_number, path)
  write_table(combined, registry_path)
  invisible(registry_path)
}

relative_to_root <- function(path) {
  root <- normalizePath(getwd(), mustWork = TRUE)
  normalized <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(root, .Platform$file.sep)
  if (startsWith(normalized, prefix)) substring(normalized, nchar(prefix) + 1L) else path
}

theme_nature <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.2, colour = "#4D4D4D", margin = ggplot2::margin(b = 7)),
      axis.title = ggplot2::element_text(face = "plain", colour = "#222222"),
      axis.text = ggplot2::element_text(colour = "#333333"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

magfunc_palette <- c(
  blue = "#3B6FB6",
  gold = "#D39B2A",
  orange = "#D96C3F",
  olive = "#788C43",
  pink = "#B05A7A",
  ink = "#222222",
  muted = "#7A7A7A",
  pale = "#D9E2EF"
)

save_publication_plot <- function(plot, stem, cfg, width = NULL, height = NULL) {
  paths <- output_paths(cfg)
  width <- width %||% cfg$outputs$figure_width_in %||% 7.2
  height <- height %||% cfg$outputs$figure_height_in %||% 5.2
  dpi <- cfg$outputs$png_dpi %||% 600
  pdf_path <- file.path(paths$figures, paste0(stem, ".pdf"))
  png_path <- file.path(paths$figures, paste0(stem, ".png"))
  svg_path <- file.path(paths$figures, paste0(stem, ".svg"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(svg_path, plot, width = width, height = height, device = svglite::svglite, bg = "white")
  c(pdf_path, png_path, svg_path)
}

parse_metadata_date <- function(x) {
  x <- as.character(x)
  out <- as.Date(x, format = "%Y%m%d")
  missing <- is.na(out)
  if (any(missing)) out[missing] <- as.Date(x[missing])
  out
}

prepare_metadata <- function(path, cfg) {
  dt <- data.table::fread(path)
  cols <- cfg$columns
  required <- unlist(cols[c("sample_id", "subject_id", "timepoint", "date", "altitude")], use.names = FALSE)
  absent <- setdiff(required, names(dt))
  if (length(absent) > 0) stop("Metadata is missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  dt[, analysis_date := parse_metadata_date(get(cols$date))]
  dt[, Subject := as.character(get(cols$subject_id))]
  dt[, Timepoint := as.character(get(cols$timepoint))]
  dt[, altitude_numeric := as.numeric(get(cols$altitude))]
  dt <- dt[!Timepoint %in% unlist(cfg$analysis$excluded_timepoints)]
  dt[, Timepoint := factor(Timepoint, levels = unlist(cfg$analysis$ordered_timepoints), ordered = TRUE)]
  if (any(is.na(dt$Timepoint))) stop("Metadata contains timepoints outside the configured T1-T10 analysis after exclusions.", call. = FALSE)
  dt[, days_from_baseline_subject := as.numeric(analysis_date - min(analysis_date)), by = Subject]
  dt
}

write_session_info <- function(cfg, script_name) {
  paths <- output_paths(cfg)
  path <- file.path(paths$reproducibility, paste0(script_name, "_session_info.txt"))
  capture.output(sessionInfo(), file = path)
  invisible(path)
}

safe_bh <- function(p) {
  out <- rep(NA_real_, length(p))
  keep <- is.finite(p)
  out[keep] <- stats::p.adjust(p[keep], method = "BH")
  out
}
