`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

parse_beta_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(config = "config/beta_diversity.yml", force = FALSE, metric = NULL)
  i <- 1
  while (i <= length(args)) {
    if (args[[i]] == "--config") { i <- i + 1; out$config <- args[[i]]
    } else if (args[[i]] == "--force") { out$force <- TRUE
    } else if (args[[i]] == "--metric") { i <- i + 1; out$metric <- args[[i]]
    } else stop("Unknown argument: ", args[[i]], call. = FALSE)
    i <- i + 1
  }
  out
}

require_beta_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

load_beta_config <- function(path) {
  require_beta_packages("yaml")
  if (!file.exists(path)) stop("Configuration not found: ", path, call. = FALSE)
  config <- yaml::read_yaml(path)
  config$.config_path <- normalizePath(path, mustWork = TRUE)
  config
}

beta_output_dir <- function(config, ...) file.path(config$outputs$output_dir, ...)

ensure_beta_dirs <- function(config) {
  dirs <- c("audit", "qc", "prepared_data", "distances", "ordination", "statistics",
            "figures", "tables", "reports", "logs", "provenance")
  invisible(vapply(file.path(config$outputs$output_dir, dirs), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

assert_writable_output <- function(path, force = FALSE) {
  if (file.exists(path) && !force) {
    stop("Refusing to overwrite existing output: ", path, ". Re-run with --force or use a new output directory.", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

read_delimited_table <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  suffix <- tolower(tools::file_ext(path))
  if (suffix == "csv") utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  else if (suffix %in% c("tsv", "txt")) utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  else stop("Unsupported tabular format: ", suffix, call. = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}

write_matrix_tsv <- function(x, path) {
  write_tsv(data.frame(sample_id = rownames(x), x, check.names = FALSE), path)
}

canonicalize_metadata <- function(metadata, config) {
  cols <- config$columns
  metadata$sample_id <- as.character(metadata[[cols$sample_id]])
  metadata$subject_id <- as.character(metadata[[cols$subject_id]])
  metadata$timepoint <- factor(metadata[[cols$timepoint]], levels = unlist(cols$ordered_time_levels), ordered = TRUE)
  metadata$sample_date <- as.Date(as.character(metadata[[cols$time]]), format = "%Y%m%d")
  metadata$days_from_subject_baseline <- ave(as.numeric(metadata$sample_date), metadata$subject_id, FUN = function(x) x - min(x, na.rm = TRUE))
  metadata$altitude <- suppressWarnings(as.numeric(metadata$altitude))
  metadata$age <- suppressWarnings(as.numeric(metadata$age))
  metadata$gender <- factor(metadata$gender)
  metadata
}
