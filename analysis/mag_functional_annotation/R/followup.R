# Shared helpers for the prespecified MAG functional-annotation follow-up.

configured_followup_candidates <- function(cfg) {
  candidates <- cfg$followup$candidates %||% list()
  if (length(candidates) == 0L) {
    stop("Configuration must define at least one followup candidate.", call. = FALSE)
  }
  out <- data.table::rbindlist(lapply(candidates, function(x) {
    data.table::data.table(
      ontology = as.character(x$ontology),
      term = as.character(x$term),
      selection_basis = as.character(x$selection_basis %||% "prespecified follow-up")
    )
  }), fill = TRUE)
  if (anyNA(out$ontology) || anyNA(out$term) || any(!nzchar(out$ontology)) || any(!nzchar(out$term))) {
    stop("Every followup candidate requires non-empty ontology and term values.", call. = FALSE)
  }
  if (anyDuplicated(out, by = c("ontology", "term"))) {
    stop("Followup candidates must be unique by ontology and term.", call. = FALSE)
  }
  out[, candidate_id := paste(ontology, term, sep = ":")]
  data.table::setcolorder(out, c("candidate_id", "ontology", "term", "selection_basis"))
  out
}

configured_candidate_curation <- function(cfg) {
  candidates <- cfg$followup$candidates %||% list()
  review <- cfg$followup$candidate_curation %||% list()
  review_date <- as.character(review$review_date %||% "")
  relationship_statement <- as.character(review$relationship_statement %||% "")
  relationship_claim_boundary <- as.character(review$relationship_claim_boundary %||% "")

  out <- data.table::rbindlist(lapply(candidates, function(x) {
    curation <- x$curation %||% list()
    sources <- curation$authoritative_sources %||% list()
    source_names <- vapply(sources, function(source) {
      as.character(source$name %||% "")
    }, character(1))
    source_urls <- vapply(sources, function(source) {
      as.character(source$url %||% "")
    }, character(1))
    data.table::data.table(
      candidate_id = paste(as.character(x$ontology), as.character(x$term), sep = ":"),
      ontology = as.character(x$ontology),
      term = as.character(x$term),
      accepted_name = as.character(curation$accepted_name %||% ""),
      gene_symbol = as.character(curation$gene_symbol %||% ""),
      entity_type = as.character(curation$entity_type %||% ""),
      reaction = as.character(curation$reaction %||% ""),
      authoritative_source_names = paste(source_names, collapse = "; "),
      authoritative_source_urls = paste(source_urls, collapse = "; "),
      review_date = review_date,
      candidate_claim_boundary = as.character(curation$claim_boundary %||% ""),
      relationship_statement = relationship_statement,
      relationship_claim_boundary = relationship_claim_boundary
    )
  }), fill = TRUE)

  required_nonempty <- c(
    "candidate_id", "ontology", "term", "accepted_name", "gene_symbol",
    "entity_type", "authoritative_source_names", "authoritative_source_urls",
    "review_date", "candidate_claim_boundary", "relationship_statement",
    "relationship_claim_boundary"
  )
  empty <- required_nonempty[vapply(required_nonempty, function(field) {
    any(is.na(out[[field]]) | !nzchar(trimws(out[[field]])))
  }, logical(1))]
  if (length(empty) > 0L) {
    stop(
      "Candidate curation metadata contains empty required fields: ",
      paste(empty, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(out$candidate_id)) {
    stop("Candidate curation metadata must be unique by candidate ID.", call. = FALSE)
  }
  out
}

align_gut_metadata <- function(cfg) {
  inputs <- resolve_input_paths(cfg)
  gut <- data.table::fread(inputs$abundance_file)
  sample_col <- cfg$columns$sample_id
  gut_mags <- setdiff(names(gut), sample_col)
  metadata_all <- data.table::fread(inputs$metadata_file)
  missing_metadata <- setdiff(as.character(gut[[sample_col]]), as.character(metadata_all[[sample_col]]))
  if (length(missing_metadata) > 0L) {
    stop(length(missing_metadata), " gut sample IDs have no metadata row.", call. = FALSE)
  }
  metadata <- prepare_metadata(inputs$metadata_file, cfg)
  eligible_ids <- intersect(as.character(gut[[sample_col]]), as.character(metadata[[sample_col]]))
  gut <- gut[get(sample_col) %in% eligible_ids]
  metadata <- metadata[get(sample_col) %in% eligible_ids]
  data.table::setorderv(gut, sample_col)
  data.table::setorderv(metadata, sample_col)
  if (!identical(as.character(gut[[sample_col]]), as.character(metadata[[sample_col]]))) {
    stop("Gut and metadata sample order could not be aligned.", call. = FALSE)
  }
  abundance <- as.matrix(gut[, ..gut_mags])
  storage.mode(abundance) <- "double"
  if (!identical(cfg$analysis$abundance_row_normalization, "none")) {
    stop("The approved primary analysis requires abundance_row_normalization: none.", call. = FALSE)
  }
  list(
    gut = gut,
    gut_mags = gut_mags,
    metadata = metadata,
    abundance = abundance,
    row_sums = rowSums(abundance),
    sample_ids = as.character(gut[[sample_col]]),
    sample_col = sample_col
  )
}

fit_candidate_timepoint_model <- function(data, cfg, pseudocount, capture_adjusted = FALSE) {
  d <- data.table::copy(data)
  if (!is.finite(pseudocount) || pseudocount <= 0) {
    return(list(status = "invalid_pseudocount", model = NULL, null_model = NULL))
  }
  d[, response := log10(potential + pseudocount)]
  d[, Timepoint_unordered := factor(
    as.character(Timepoint),
    levels = unlist(cfg$analysis$ordered_timepoints)
  )]
  if (capture_adjusted) {
    d[, captured_log10_centered := log10(abundance_sum_as_provided) - mean(log10(abundance_sum_as_provided))]
  }
  captured_warnings <- character()
  tryCatch(withCallingHandlers({
    control <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
    null_formula <- if (capture_adjusted) {
      stats::as.formula("response ~ captured_log10_centered + (1 | Subject)")
    } else {
      stats::as.formula("response ~ 1 + (1 | Subject)")
    }
    full_formula <- if (capture_adjusted) {
      stats::as.formula("response ~ captured_log10_centered + Timepoint_unordered + (1 | Subject)")
    } else {
      stats::as.formula("response ~ Timepoint_unordered + (1 | Subject)")
    }
    null_model <- lme4::lmer(null_formula, data = d, REML = FALSE, control = control)
    full_model <- lme4::lmer(full_formula, data = d, REML = FALSE, control = control)
    comparison <- stats::anova(null_model, full_model)
    convergence_messages <- unique(c(
      captured_warnings,
      as.character(full_model@optinfo$conv$lme4$messages %||% character()),
      if (!is.null(full_model@optinfo$conv$opt) && full_model@optinfo$conv$opt != 0L) {
        paste0("optimizer_code_", full_model@optinfo$conv$opt)
      } else character()
    ))
    list(
      status = "fitted",
      model = full_model,
      null_model = null_model,
      chisq = as.numeric(comparison$Chisq[[2L]]),
      df = as.numeric(comparison$Df[[2L]]),
      p_value = as.numeric(comparison$`Pr(>Chisq)`[[2L]]),
      singular = lme4::isSingular(full_model, tol = 1e-4),
      converged = length(convergence_messages) == 0L,
      convergence_messages = paste(convergence_messages, collapse = " | "),
      data = d
    )
  }, warning = function(w) {
    captured_warnings <<- unique(c(captured_warnings, conditionMessage(w)))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    list(status = paste0("model_error: ", conditionMessage(e)), model = NULL, null_model = NULL)
  })
}

fit_timepoint_response_model <- function(data, cfg, response_column) {
  d <- data.table::copy(data)
  if (!response_column %in% names(d)) {
    return(list(status = paste0("missing_response_column: ", response_column), model = NULL, null_model = NULL))
  }
  d[, response := as.numeric(get(response_column))]
  if (any(!is.finite(d$response))) {
    return(list(status = "nonfinite_response", model = NULL, null_model = NULL))
  }
  d[, Timepoint_unordered := factor(
    as.character(Timepoint),
    levels = unlist(cfg$analysis$ordered_timepoints)
  )]
  captured_warnings <- character()
  tryCatch(withCallingHandlers({
    control <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
    null_model <- lme4::lmer(response ~ 1 + (1 | Subject), data = d, REML = FALSE, control = control)
    full_model <- lme4::lmer(response ~ Timepoint_unordered + (1 | Subject), data = d, REML = FALSE, control = control)
    comparison <- stats::anova(null_model, full_model)
    convergence_messages <- unique(c(
      captured_warnings,
      as.character(full_model@optinfo$conv$lme4$messages %||% character()),
      if (!is.null(full_model@optinfo$conv$opt) && full_model@optinfo$conv$opt != 0L) {
        paste0("optimizer_code_", full_model@optinfo$conv$opt)
      } else character()
    ))
    list(
      status = "fitted",
      model = full_model,
      null_model = null_model,
      chisq = as.numeric(comparison$Chisq[[2L]]),
      df = as.numeric(comparison$Df[[2L]]),
      p_value = as.numeric(comparison$`Pr(>Chisq)`[[2L]]),
      singular = lme4::isSingular(full_model, tol = 1e-4),
      converged = length(convergence_messages) == 0L,
      convergence_messages = paste(convergence_messages, collapse = " | "),
      data = d
    )
  }, warning = function(w) {
    captured_warnings <<- unique(c(captured_warnings, conditionMessage(w)))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    list(status = paste0("model_error: ", conditionMessage(e)), model = NULL, null_model = NULL)
  })
}

fixed_timepoint_estimates <- function(model, timepoints) {
  timepoints <- as.character(timepoints)
  new_data <- data.frame(Timepoint_unordered = factor(timepoints, levels = timepoints))
  fixed_names <- names(lme4::fixef(model))
  if ("captured_log10_centered" %in% fixed_names) new_data$captured_log10_centered <- 0
  fixed_formula <- if ("captured_log10_centered" %in% fixed_names) {
    stats::as.formula("~ captured_log10_centered + Timepoint_unordered")
  } else {
    stats::as.formula("~ Timepoint_unordered")
  }
  matrix <- stats::model.matrix(fixed_formula, new_data)
  absent <- setdiff(fixed_names, colnames(matrix))
  if (length(absent) > 0L) {
    matrix <- cbind(matrix, matrix(0, nrow = nrow(matrix), ncol = length(absent), dimnames = list(NULL, absent)))
  }
  matrix <- matrix[, fixed_names, drop = FALSE]
  beta <- lme4::fixef(model)
  covariance <- as.matrix(stats::vcov(model))[fixed_names, fixed_names, drop = FALSE]
  estimate <- as.vector(matrix %*% beta)
  standard_error <- sqrt(pmax(0, diag(matrix %*% covariance %*% t(matrix))))
  data.table::data.table(
    Timepoint = timepoints,
    estimate_log10 = estimate,
    standard_error = standard_error,
    lower_95 = estimate - 1.96 * standard_error,
    upper_95 = estimate + 1.96 * standard_error,
    potential_scale_estimate = 10^estimate
  )
}

candidate_pairwise_contrasts <- function(model, timepoints) {
  timepoints <- as.character(timepoints)
  new_data <- data.frame(Timepoint_unordered = factor(timepoints, levels = timepoints))
  fixed_names <- names(lme4::fixef(model))
  if ("captured_log10_centered" %in% fixed_names) new_data$captured_log10_centered <- 0
  fixed_formula <- if ("captured_log10_centered" %in% fixed_names) {
    stats::as.formula("~ captured_log10_centered + Timepoint_unordered")
  } else {
    stats::as.formula("~ Timepoint_unordered")
  }
  matrix <- stats::model.matrix(fixed_formula, new_data)
  absent <- setdiff(fixed_names, colnames(matrix))
  if (length(absent) > 0L) {
    matrix <- cbind(matrix, matrix(0, nrow = nrow(matrix), ncol = length(absent), dimnames = list(NULL, absent)))
  }
  matrix <- matrix[, fixed_names, drop = FALSE]
  beta <- lme4::fixef(model)
  covariance <- as.matrix(stats::vcov(model))[fixed_names, fixed_names, drop = FALSE]
  pairs <- utils::combn(seq_along(timepoints), 2L)
  rows <- lapply(seq_len(ncol(pairs)), function(i) {
    before_index <- pairs[1L, i]
    after_index <- pairs[2L, i]
    contrast <- matrix[after_index, ] - matrix[before_index, ]
    estimate <- as.numeric(sum(contrast * beta))
    standard_error <- sqrt(max(0, as.numeric(t(contrast) %*% covariance %*% contrast)))
    z_value <- if (standard_error > 0) estimate / standard_error else NA_real_
    p_value <- if (is.finite(z_value)) 2 * stats::pnorm(-abs(z_value)) else NA_real_
    data.table::data.table(
      before_timepoint = timepoints[[before_index]],
      after_timepoint = timepoints[[after_index]],
      estimate_log10_after_minus_before = estimate,
      standard_error = standard_error,
      lower_95 = estimate - 1.96 * standard_error,
      upper_95 = estimate + 1.96 * standard_error,
      ratio_after_vs_before = 10^estimate,
      z_value = z_value,
      p_value = p_value
    )
  })
  out <- data.table::rbindlist(rows)
  out[, p_adj_bh_within_candidate := safe_bh(p_value)]
  out
}

paired_change_drivers <- function(data, before_timepoint, after_timepoint,
                                  id_columns, value_column = "mag_value") {
  required <- c("Subject", "Timepoint", id_columns, value_column)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Paired driver data are missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  before_timepoint <- as.character(before_timepoint)
  after_timepoint <- as.character(after_timepoint)
  d <- data.table::copy(data)[as.character(Timepoint) %in% c(before_timepoint, after_timepoint)]
  if (nrow(d) == 0L) return(data.table::data.table())
  d[, Timepoint := as.character(Timepoint)]
  grouped <- d[, .(value = sum(as.numeric(get(value_column)))), by = c("Subject", "Timepoint", id_columns)]
  formula <- stats::as.formula(paste(
    paste(c("Subject", id_columns), collapse = " + "), "~ Timepoint"
  ))
  wide <- data.table::dcast(grouped, formula, value.var = "value")
  if (!all(c(before_timepoint, after_timepoint) %in% names(wide))) return(data.table::data.table())
  wide <- wide[is.finite(get(before_timepoint)) & is.finite(get(after_timepoint))]
  if (nrow(wide) == 0L) return(data.table::data.table())
  wide[, `:=`(
    before_value = as.numeric(get(before_timepoint)),
    after_value = as.numeric(get(after_timepoint))
  )]
  wide[, delta_after_minus_before := after_value - before_value]
  out <- wide[, .(
    n_pairs = .N,
    mean_before = mean(before_value),
    mean_after = mean(after_value),
    mean_delta = mean(delta_after_minus_before),
    median_delta = stats::median(delta_after_minus_before),
    standard_error_delta = if (.N > 1L) stats::sd(delta_after_minus_before) / sqrt(.N) else NA_real_,
    positive_subjects = sum(delta_after_minus_before > 0),
    negative_subjects = sum(delta_after_minus_before < 0),
    unchanged_subjects = sum(delta_after_minus_before == 0)
  ), by = id_columns]
  out[, `:=`(
    before_timepoint = before_timepoint,
    after_timepoint = after_timepoint
  )]
  total_absolute_change <- sum(abs(out$mean_delta))
  net_change <- sum(out$mean_delta)
  out[, `:=`(
    absolute_change_share = if (total_absolute_change > 0) abs(mean_delta) / total_absolute_change else NA_real_,
    signed_fraction_of_net_change = if (abs(net_change) > .Machine$double.eps) mean_delta / net_change else NA_real_,
    net_change = net_change,
    total_absolute_change = total_absolute_change
  )]
  out
}

within_subject_correlation <- function(x, y, subject) {
  centered_x <- x - ave(x, subject, FUN = function(z) mean(z, na.rm = TRUE))
  centered_y <- y - ave(y, subject, FUN = function(z) mean(z, na.rm = TRUE))
  keep <- is.finite(centered_x) & is.finite(centered_y)
  if (sum(keep) < 3L || stats::sd(centered_x[keep]) == 0 || stats::sd(centered_y[keep]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(centered_x[keep], centered_y[keep], method = "spearman"))
}

parse_mag_construction_provenance <- function(mag) {
  mag <- as.character(mag)
  specimen <- sub("^.*-[0-9]{2}([BCD])_bin\\..*$", "\\1", mag, perl = TRUE)
  specimen[specimen == mag] <- NA_character_
  participant <- sub("^.*-([0-9]{2})[BCD]_bin\\..*$", "\\1", mag, perl = TRUE)
  participant[participant == mag] <- NA_character_
  project <- sub("^([A-Za-z]+).*$", "\\1", mag, perl = TRUE)
  project[project == mag & !grepl("^[A-Za-z]+", mag)] <- NA_character_
  tier <- sub("^.*\\.(strict|permissive|orig)$", "\\1", mag, perl = TRUE)
  tier[tier == mag] <- NA_character_
  data.table::data.table(
    MAG = mag,
    construction_project = project,
    construction_participant = participant,
    construction_specimen_code = specimen,
    construction_specimen = unname(c(B = "fecal", C = "saliva", D = "skin")[specimen]),
    construction_tier = tier
  )
}

resolve_optional_followup_input <- function(cfg, key, environment_variable) {
  configured <- cfg$followup$optional_inputs[[key]] %||% ""
  value <- Sys.getenv(environment_variable, unset = as.character(configured))
  if (!nzchar(value)) return(NA_character_)
  value
}
