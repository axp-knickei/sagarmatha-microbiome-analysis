suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

required_dynamics_packages <- function(extra = character()) {
  unique(c("data.table", "yaml", "vegan", "permute", "lme4", "ggplot2", "scales", "svglite", extra))
}

assert_dynamics_packages <- function(extra = character()) {
  packages <- required_dynamics_packages(extra)
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(packages)
}

parse_dynamics_args <- function(default_config = "config/functional_community_dynamics.yml") {
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

read_dynamics_config <- function(path) {
  if (!file.exists(path)) stop("Configuration not found: ", path, call. = FALSE)
  cfg <- yaml::read_yaml(path)
  required <- c("inputs", "columns", "analysis", "outputs")
  absent <- setdiff(required, names(cfg))
  if (length(absent)) stop("Configuration is missing sections: ", paste(absent, collapse = ", "), call. = FALSE)
  cfg$.config_path <- normalizePath(path, mustWork = TRUE)
  cfg
}

resolve_dynamics_inputs <- function(cfg) {
  paths <- cfg$inputs
  overrides <- cfg$input_environment_overrides %||% list()
  for (key in names(paths)) {
    env_name <- overrides[[key]] %||% ""
    if (nzchar(env_name)) {
      value <- Sys.getenv(env_name, unset = "")
      if (nzchar(value)) paths[[key]] <- value
    }
  }
  paths
}

dynamics_paths <- function(cfg, create = TRUE) {
  root <- Sys.getenv("FUNCDYN_OUTPUT_DIR", unset = cfg$outputs$output_dir)
  out <- list(
    root = root,
    qc = file.path(root, "00_quality_control"),
    prepared = file.path(root, "01_prepared_profiles"),
    multivariate = file.path(root, "02_multivariate_analysis"),
    longitudinal = file.path(root, "03_longitudinal_stability"),
    figures = file.path(root, "figures"),
    reports = file.path(root, "reports"),
    reproducibility = file.path(root, "reproducibility"),
    intermediate = file.path(root, "intermediate"),
    logs = file.path(root, "logs")
  )
  if (create) invisible(lapply(out, dir.create, recursive = TRUE, showWarnings = FALSE))
  out
}

guard_dynamics_stage <- function(marker, args, cfg) {
  overwrite <- isTRUE(args$force) || isTRUE(cfg$outputs$overwrite)
  if (file.exists(marker) && !overwrite) {
    stop("Stage completion marker already exists: ", marker, ". Use --force only when intentional.", call. = FALSE)
  }
  if (file.exists(marker) && overwrite && unlink(marker) != 0L) {
    stop("Could not remove previous completion marker: ", marker, call. = FALSE)
  }
  invisible(overwrite)
}

read_dynamics_table <- function(path, ...) {
  if (!file.exists(path)) stop("Required input not found: ", path, call. = FALSE)
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) return(data.table::fread(path, ...))
  gzip <- unname(Sys.which("gzip"))
  if (!nzchar(gzip)) stop("gzip is required to read: ", path, call. = FALSE)
  data.table::fread(cmd = paste(shQuote(gzip), "-cd", shQuote(normalizePath(path, mustWork = TRUE))), ...)
}

write_dynamics_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data.table::as.data.table(x), path, na = "")
  invisible(path)
}

write_dynamics_session <- function(cfg, stage) {
  path <- file.path(dynamics_paths(cfg)$reproducibility, paste0(stage, "_session_info.txt"))
  capture.output(sessionInfo(), file = path)
  invisible(path)
}

safe_bh <- function(x) {
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  out[keep] <- stats::p.adjust(x[keep], method = "BH")
  out
}

participant_transition_sensitivity <- function(matched, layer, role) {
  required <- c(
    "subject_id", "taxonomic_standardized_distance",
    "functional_standardized_distance", "log_functional_to_taxonomic_ratio"
  )
  absent <- setdiff(required, names(matched))
  if (length(absent)) {
    stop("Transition sensitivity input is missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  }

  participant <- matched[, .(
    transitions = .N,
    spearman_rho = if (
      .N >= 3L &&
      stats::sd(taxonomic_standardized_distance) > 0 &&
      stats::sd(functional_standardized_distance) > 0
    ) {
      suppressWarnings(stats::cor(
        taxonomic_standardized_distance,
        functional_standardized_distance,
        method = "spearman"
      ))
    } else {
      NA_real_
    },
    median_log_functional_to_taxonomic_ratio = stats::median(
      log_functional_to_taxonomic_ratio,
      na.rm = TRUE
    )
  ), by = subject_id]
  participant[, `:=`(layer = layer, role = role)]

  summarize_signed_rank <- function(values, family, interpretation) {
    values <- values[is.finite(values)]
    if (length(values) < 3L) {
      return(data.table(
        layer = layer,
        role = role,
        sensitivity_family = family,
        participant_summary = NA_real_,
        ratio = NA_real_,
        positive_direction_participants = sum(values > 0),
        participants = length(values),
        p_value = NA_real_,
        status = "insufficient_participants",
        interpretation = interpretation
      ))
    }
    test <- suppressWarnings(stats::wilcox.test(
      values,
      mu = 0,
      alternative = "two.sided",
      exact = FALSE
    ))
    estimate <- stats::median(values)
    data.table(
      layer = layer,
      role = role,
      sensitivity_family = family,
      participant_summary = estimate,
      ratio = if (family == "participant_median_relative_movement") exp(estimate) else NA_real_,
      positive_direction_participants = sum(values > 0),
      participants = length(values),
      p_value = unname(test$p.value),
      status = "fitted",
      interpretation = interpretation
    )
  }

  summary <- rbindlist(list(
    summarize_signed_rank(
      participant$spearman_rho,
      "participant_spearman_transition_coupling",
      "Median participant-level Spearman correlation; positive values mean taxonomic and functional transitions move together within participants."
    ),
    summarize_signed_rank(
      participant$median_log_functional_to_taxonomic_ratio,
      "participant_median_relative_movement",
      "Median of participant-specific median log ratios; ratios below one support functional smoothing and ratios above one do not."
    )
  ), fill = TRUE)

  list(participant = participant, summary = summary)
}

close_rows <- function(x) {
  totals <- rowSums(x)
  if (any(!is.finite(totals) | totals <= 0)) stop("Every profile row must have a finite positive total.", call. = FALSE)
  x / totals
}

hellinger_transform <- function(x) sqrt(close_rows(x))

multiplicative_replace <- function(x, fraction = 0.65) {
  x <- close_rows(x)
  replaced <- t(apply(x, 1L, function(row) {
    zeros <- row == 0
    k <- sum(zeros)
    if (!k) return(row)
    positive <- row[!zeros]
    if (!length(positive)) stop("Aitchison zero replacement encountered an all-zero row.", call. = FALSE)
    delta <- min(fraction * min(positive), 0.5 / k)
    row[zeros] <- delta
    row[!zeros] <- positive * ((1 - k * delta) / sum(positive))
    row
  }))
  dimnames(replaced) <- dimnames(x)
  replaced
}

clr_transform <- function(x) {
  if (any(!is.finite(x) | x <= 0)) stop("CLR input must be finite and strictly positive.", call. = FALSE)
  logged <- log(x)
  logged - rowMeans(logged)
}

profile_distance <- function(x, metric) {
  if (metric == "hellinger") return(stats::dist(hellinger_transform(x), method = "euclidean"))
  if (metric == "bray_curtis") return(vegan::vegdist(close_rows(x), method = "bray"))
  if (metric == "aitchison") return(stats::dist(clr_transform(multiplicative_replace(x)), method = "euclidean"))
  stop("Unsupported distance metric: ", metric, call. = FALSE)
}

distance_summary <- function(distance) {
  values <- as.numeric(distance)
  data.table(
    minimum = min(values), median = stats::median(values), mean = mean(values),
    maximum = max(values), missing_values = sum(!is.finite(values))
  )
}

restricted_permutations <- function(metadata, nperm) {
  control <- permute::how(nperm = nperm)
  permute::setBlocks(control) <- factor(metadata$subject_id)
  control
}

tidy_adonis <- function(fit, layer, metric, nperm, n) {
  tab <- as.data.frame(fit)
  tab$hypothesis <- rownames(tab)
  rownames(tab) <- NULL
  names(tab) <- sub("Pr\\(>F\\)", "p_value", names(tab))
  data.table(
    layer = layer,
    metric = metric,
    hypothesis = tab$hypothesis,
    permutations = nperm,
    df = tab$Df,
    sum_squares = tab$SumOfSqs,
    r_squared = tab$R2,
    statistic = tab$F,
    p_value = tab$p_value,
    sample_count = n
  )
}

make_pcoa <- function(distance, metadata, layer, metric) {
  k <- min(2L, attr(distance, "Size") - 1L)
  fit <- vegan::wcmdscale(distance, k = k, eig = TRUE, add = "lingoes")
  points <- as.data.table(fit$points, keep.rownames = "sample_id")
  setnames(points, old = names(points)[2:ncol(points)], new = paste0("PCoA", seq_len(ncol(points) - 1L)))
  if (!"PCoA2" %in% names(points)) points[, PCoA2 := 0]
  points <- merge(points, metadata, by = "sample_id", all.x = TRUE, sort = FALSE)
  total_positive <- sum(fit$eig[fit$eig > 0])
  variance <- if (total_positive > 0) 100 * fit$eig[seq_len(k)] / total_positive else rep(NA_real_, k)
  points[, `:=`(layer = layer, metric = metric, axis1_percent = variance[[1]], axis2_percent = variance[[min(2L, length(variance))]])]
  points
}

timepoint_levels <- function(cfg) as.character(unlist(cfg$analysis$ordered_timepoints))

prepare_dynamics_metadata <- function(path, cfg) {
  metadata <- data.table::fread(path)
  cols <- cfg$columns
  required <- unlist(cols[c("sample_id", "subject_id", "timepoint", "date")], use.names = FALSE)
  absent <- setdiff(required, names(metadata))
  if (length(absent)) stop("Metadata is missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  metadata[, sample_id := as.character(get(cols$sample_id))]
  metadata[, subject_id := as.character(get(cols$subject_id))]
  metadata[, timepoint := as.character(get(cols$timepoint))]
  date_text <- as.character(metadata[[cols$date]])
  parsed_date <- as.IDate(date_text, format = "%Y%m%d")
  parsed_date[is.na(parsed_date)] <- as.IDate(date_text[is.na(parsed_date)])
  metadata[, sample_date := parsed_date]
  if (anyNA(metadata$sample_date)) stop("Metadata contains unparseable sample dates.", call. = FALSE)
  metadata <- metadata[!timepoint %in% unlist(cfg$analysis$excluded_timepoints)]
  metadata[, timepoint := factor(timepoint, levels = timepoint_levels(cfg), ordered = TRUE)]
  if (anyNA(metadata$timepoint)) stop("Unexpected timepoint remained after exclusions.", call. = FALSE)
  if (anyDuplicated(metadata$sample_id)) stop("Metadata sample identifiers are not unique.", call. = FALSE)
  metadata[, days_from_earliest_observation := as.numeric(sample_date - min(sample_date)), by = subject_id]
  metadata
}

configured_layers <- function(cfg) {
  rbindlist(lapply(cfg$analysis$functional_layers, function(x) {
    data.table(layer = x$id, role = x$role, file = x$file)
  }))
}

configured_metrics <- function(cfg) {
  rbindlist(lapply(cfg$analysis$metrics, function(x) {
    data.table(metric = x$id, metric_role = x$role)
  }))
}

fit_dynamics_lmer <- function(formula, data) {
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      lme4::lmer(formula, data = data, REML = FALSE),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = identity
  )
  if (inherits(fit, "error")) {
    return(list(status = paste0("model_error: ", conditionMessage(fit)), fit = NULL,
                singular = NA, converged = FALSE, messages = paste(unique(warnings), collapse = "; ")))
  }
  messages <- c(warnings, unlist(fit@optinfo$conv$lme4$messages %||% character()))
  list(
    status = "fitted", fit = fit, singular = lme4::isSingular(fit),
    converged = !length(fit@optinfo$conv$lme4$messages %||% character()),
    messages = paste(unique(messages), collapse = "; ")
  )
}

checksum_file <- function(path) {
  unname(tools::md5sum(path))
}

checksum_object <- function(x) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 3)
  checksum_file(path)
}

extract_longitudinal_records <- function(distance, metadata, layer, metric) {
  matrix <- as.matrix(distance)
  metadata <- metadata[match(rownames(matrix), sample_id)]
  baseline <- rbindlist(lapply(split(metadata, metadata$subject_id), function(subject) {
    setorder(subject, sample_date, sample_id)
    reference <- subject$sample_id[[1]]
    data.table(
      record_type = "distance_to_earliest_observation",
      layer = layer,
      metric = metric,
      subject_id = subject$subject_id,
      sample_id = subject$sample_id,
      from_sample = reference,
      to_sample = subject$sample_id,
      from_timepoint = as.character(subject$timepoint[[1]]),
      to_timepoint = as.character(subject$timepoint),
      transition_days = subject$days_from_earliest_observation,
      distance = matrix[subject$sample_id, reference]
    )
  }))
  adjacent <- rbindlist(lapply(split(metadata, metadata$subject_id), function(subject) {
    setorder(subject, sample_date, sample_id)
    if (nrow(subject) < 2L) return(NULL)
    from <- head(subject$sample_id, -1L)
    to <- tail(subject$sample_id, -1L)
    data.table(
      record_type = "adjacent_observed_transition",
      layer = layer,
      metric = metric,
      subject_id = subject$subject_id[-1L],
      sample_id = to,
      from_sample = from,
      to_sample = to,
      from_timepoint = as.character(head(subject$timepoint, -1L)),
      to_timepoint = as.character(tail(subject$timepoint, -1L)),
      transition_days = as.numeric(diff(subject$sample_date)),
      distance = mapply(function(a, b) matrix[a, b], from, to)
    )
  }), fill = TRUE)
  list(baseline = baseline, adjacent = adjacent)
}

recovery_records <- function(distance, metadata, layer, metric) {
  matrix <- as.matrix(distance)
  eligible <- metadata[, .(has_t1 = any(timepoint == "T1"), has_t10 = any(timepoint == "T10")), by = subject_id][has_t1 & has_t10, subject_id]
  if (!length(eligible)) return(data.table())
  rbindlist(lapply(eligible, function(id) {
    subject <- metadata[subject_id == id]
    setorder(subject, sample_date, sample_id)
    t1 <- subject[timepoint == "T1", sample_id][[1]]
    t10 <- subject[timepoint == "T10", sample_id][[1]]
    distances <- matrix[subject$sample_id, t1]
    peak <- max(distances)
    t10_distance <- matrix[t10, t1]
    data.table(
      record_type = "exploratory_t1_recovery",
      layer = layer,
      metric = metric,
      subject_id = id,
      sample_id = t10,
      from_sample = t1,
      to_sample = t10,
      from_timepoint = "T1",
      to_timepoint = "T10",
      transition_days = as.numeric(subject[sample_id == t10, sample_date] - subject[sample_id == t1, sample_date]),
      distance = t10_distance,
      peak_distance_from_t1 = peak,
      recovery_fraction_from_peak = if (peak > 0) 1 - t10_distance / peak else NA_real_
    )
  }))
}

theme_dynamics <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(colour = "#4D4D4D", margin = ggplot2::margin(b = 7)),
      axis.text = ggplot2::element_text(colour = "#333333"),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

dynamics_palette <- c(
  T1 = "#332288", T2 = "#117733", T3 = "#44AA99", T4 = "#88CCEE", T5 = "#DDCC77",
  T6 = "#CC6677", T7 = "#AA4499", T8 = "#882255", T9 = "#999933", T10 = "#661100"
)

save_dynamics_plot <- function(plot, stem, cfg, width = NULL, height = NULL) {
  paths <- dynamics_paths(cfg)
  width <- width %||% cfg$outputs$figure_width_in %||% 7.4
  height <- height %||% cfg$outputs$figure_height_in %||% 5.4
  dpi <- cfg$outputs$png_dpi %||% 400
  pdf <- file.path(paths$figures, paste0(stem, ".pdf"))
  png <- file.path(paths$figures, paste0(stem, ".png"))
  svg <- file.path(paths$figures, paste0(stem, ".svg"))
  ggplot2::ggsave(pdf, plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(png, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(svg, plot, width = width, height = height, device = svglite::svglite, bg = "white")
  c(pdf, png, svg)
}
