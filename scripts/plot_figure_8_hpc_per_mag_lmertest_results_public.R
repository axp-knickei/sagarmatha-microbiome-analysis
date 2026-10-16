# Public renderer: values come from external inputs; see docs/public-rendering.md.
# ============================================
# plot_figure_8_hpc_per_mag_lmertest_results.R
# ============================================
# Figure 8: diagnostic panel for selected top per-MAG LMMs.
#
# This script refits the same LMM structure used by the HPC workflow for a
# small set of top exploratory MAGs, then saves residual, Q-Q, random-effect,
# and observed-vs-fitted diagnostics. It does not alter the primary model
# results table.

required_packages <- c("ggplot2", "dplyr", "lme4")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
input_path <- file.path(repo_root, "data", "processed", "cleaned_data.rds")
results_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "tables",
  "per_mag_lmertest_hpc_results.csv"
)
output_dir <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "figures")
output_stem <- "figure_8_hpc_per_mag_lmertest_results"

# External inputs and isolated output destination (run from repository root).
input_path <- Sys.getenv("PUBLIC_CLEANED_DATA", unset = input_path)
results_path <- Sys.getenv("PUBLIC_LMM_RESULTS", unset = results_path)
output_dir <- Sys.getenv("PUBLIC_RENDER_OUTPUT_DIR", unset = output_dir)

if (!file.exists(input_path)) {
  stop("Missing processed input data: ", input_path, call. = FALSE)
}
if (!file.exists(results_path)) {
  stop("Missing per-MAG LMM result table: ", results_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- readRDS(input_path)
results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)

required_data_columns <- c("MAG", "Y", "Subject", "sampleID", "days_from_baseline_subject")
missing_data_columns <- setdiff(required_data_columns, names(df))
if (length(missing_data_columns) > 0) {
  stop("Missing required columns in cleaned_data.rds: ", paste(missing_data_columns, collapse = ", "), call. = FALSE)
}

required_result_columns <- c(
  "MAG", "status", "primary_analysis", "time_estimate", "time_p_value",
  "time_p_adj_bh_primary", "zero_fraction", "singular"
)
missing_result_columns <- setdiff(required_result_columns, names(results))
if (length(missing_result_columns) > 0) {
  stop("Missing required result columns: ", paste(missing_result_columns, collapse = ", "), call. = FALSE)
}

smallest_nonzero_y <- min(df$Y[df$Y > 0], na.rm = TRUE)
if (!is.finite(smallest_nonzero_y)) {
  stop("Cannot compute pseudocount because no positive Y values were found.", call. = FALSE)
}
pseudocount <- smallest_nonzero_y / 2

selected <- results |>
  dplyr::filter(
    primary_analysis %in% TRUE,
    status == "lmm_ok",
    is.finite(time_p_value),
    time_p_value < 0.05
  ) |>
  dplyr::arrange(time_p_value) |>
  dplyr::slice_head(n = 4) |>
  dplyr::mutate(
    mag_label = ifelse(nchar(MAG) > 32, paste0(substr(MAG, 1, 29), "..."), MAG),
    panel_label = paste0(
      dplyr::row_number(), ". ", mag_label,
      "\nbeta = ", sprintf("%+.3f", time_estimate),
      "; p = ", ifelse(time_p_value < 0.001, formatC(time_p_value, format = "e", digits = 1), sprintf("%.3f", time_p_value)),
      "; q = ", sprintf("%.3f", time_p_adj_bh_primary)
    )
  )

if (nrow(selected) == 0) {
  stop("No nominal primary lmm_ok MAGs available for Figure 8 diagnostics.", call. = FALSE)
}

fit_one_mag <- function(mag_id) {
  mag_df <- df[df$MAG == mag_id, , drop = FALSE]
  mag_df <- mag_df[stats::complete.cases(mag_df[, c("Y", "days_from_baseline_subject", "Subject")]), , drop = FALSE]
  mag_df$Subject <- droplevels(factor(mag_df$Subject))
  mag_df$response_value <- log10(mag_df$Y + pseudocount)

  warning_messages <- character(0)
  error_messages <- character(0)
  fit <- tryCatch(
    withCallingHandlers(
      lme4::lmer(
        response_value ~ days_from_baseline_subject + (1 | Subject),
        data = mag_df,
        control = lme4::lmerControl(optimizer = "bobyqa")
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_messages <<- c(error_messages, conditionMessage(e))
      NULL
    }
  )

  list(fit = fit, data = mag_df, warnings = warning_messages, errors = error_messages)
}

fit_list <- setNames(lapply(selected$MAG, fit_one_mag), selected$MAG)

diagnostics <- dplyr::bind_rows(lapply(seq_len(nrow(selected)), function(i) {
  mag_id <- selected$MAG[[i]]
  fit_obj <- fit_list[[mag_id]]$fit
  if (is.null(fit_obj)) return(data.frame())
  fit_df <- fit_list[[mag_id]]$data
  residuals <- stats::residuals(fit_obj, type = "pearson")
  fitted_values <- stats::fitted(fit_obj)
  qq <- stats::qqnorm(residuals, plot.it = FALSE)
  data.frame(
    MAG = mag_id,
    panel_label = selected$panel_label[[i]],
    Subject = fit_df$Subject,
    observed = fit_df$response_value,
    fitted = fitted_values,
    residual = residuals,
    qq_theoretical = qq$x,
    qq_sample = qq$y,
    days_from_baseline_subject = fit_df$days_from_baseline_subject,
    stringsAsFactors = FALSE
  )
}))

random_effects <- dplyr::bind_rows(lapply(seq_len(nrow(selected)), function(i) {
  mag_id <- selected$MAG[[i]]
  fit_obj <- fit_list[[mag_id]]$fit
  if (is.null(fit_obj)) return(data.frame())
  ran <- lme4::ranef(fit_obj, condVar = TRUE)$Subject
  post_var <- attr(lme4::ranef(fit_obj, condVar = TRUE)$Subject, "postVar")
  se <- sqrt(as.numeric(post_var[1, 1, ]))
  data.frame(
    MAG = mag_id,
    panel_label = selected$panel_label[[i]],
    Subject = rownames(ran),
    random_intercept = ran[, "(Intercept)"],
    conf_low = ran[, "(Intercept)"] - 1.96 * se,
    conf_high = ran[, "(Intercept)"] + 1.96 * se,
    stringsAsFactors = FALSE
  )
}))

diagnostics$panel_label <- factor(diagnostics$panel_label, levels = selected$panel_label)
random_effects$panel_label <- factor(random_effects$panel_label, levels = selected$panel_label)
random_effects <- random_effects |>
  dplyr::group_by(panel_label) |>
  dplyr::arrange(random_intercept, .by_group = TRUE) |>
  dplyr::mutate(Subject_ordered = factor(Subject, levels = Subject)) |>
  dplyr::ungroup()

fit_summary <- data.frame(
  MAG = selected$MAG,
  panel_label = selected$panel_label,
  time_estimate = selected$time_estimate,
  time_p_value = selected$time_p_value,
  time_p_adj_bh_primary = selected$time_p_adj_bh_primary,
  zero_fraction = selected$zero_fraction,
  singular_primary_table = selected$singular,
  refit_singular = vapply(fit_list, function(x) if (is.null(x$fit)) NA else lme4::isSingular(x$fit), logical(1)),
  warning_count = vapply(fit_list, function(x) length(x$warnings), integer(1)),
  error_count = vapply(fit_list, function(x) length(x$errors), integer(1)),
  warnings = vapply(fit_list, function(x) paste(unique(x$warnings), collapse = " | "), character(1)),
  errors = vapply(fit_list, function(x) paste(unique(x$errors), collapse = " | "), character(1)),
  stringsAsFactors = FALSE
)

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_sci <- function(x) formatC(x, format = "e", digits = 2)

pal <- list(
  ink = "#111827",
  text = "#374151",
  muted = "#6B7280",
  hair = "#D1D5DB",
  faint = "#F8FAFC",
  blue = "#0072B2",
  sky = "#DCEEF8",
  orange = "#D55E00",
  red = "#CC3311",
  green = "#009E73",
  grey = "#9CA3AF"
)

base_theme <- ggplot2::theme_minimal(base_size = 8.8, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 10.2, color = pal$ink),
    plot.subtitle = ggplot2::element_text(size = 7.5, color = pal$muted, margin = ggplot2::margin(b = 4)),
    axis.title = ggplot2::element_text(size = 7.8, color = pal$text),
    axis.text = ggplot2::element_text(size = 6.8, color = pal$text),
    strip.text = ggplot2::element_text(face = "bold", size = 6.7, color = pal$ink, lineheight = 0.88),
    strip.background = ggplot2::element_rect(fill = pal$faint, color = pal$hair, linewidth = 0.30),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.20),
    legend.position = "none",
    plot.margin = ggplot2::margin(4, 5, 4, 5)
  )

panel_a <- ggplot2::ggplot(diagnostics, ggplot2::aes(x = fitted, y = residual)) +
  ggplot2::geom_hline(yintercept = 0, color = pal$hair, linewidth = 0.45) +
  ggplot2::geom_point(color = pal$blue, alpha = 0.70, size = 1.15) +
  ggplot2::facet_wrap(~panel_label, nrow = 1, scales = "free_x") +
  ggplot2::labs(
    title = "A  Residuals vs fitted values",
    subtitle = "Checks residual spread and fitted-value trends on the modeled response scale",
    x = "Fitted value",
    y = "Pearson residual"
  ) +
  base_theme

panel_b <- ggplot2::ggplot(diagnostics, ggplot2::aes(x = qq_theoretical, y = qq_sample)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, color = pal$hair, linewidth = 0.45) +
  ggplot2::geom_point(color = pal$green, alpha = 0.72, size = 1.05) +
  ggplot2::facet_wrap(~panel_label, nrow = 1, scales = "free") +
  ggplot2::labs(
    title = "B  Residual Q-Q diagnostics",
    subtitle = "Compares Pearson residuals with the expected Gaussian quantiles",
    x = "Theoretical quantile",
    y = "Sample residual quantile"
  ) +
  base_theme

panel_c <- ggplot2::ggplot(random_effects, ggplot2::aes(x = random_intercept, y = Subject_ordered)) +
  ggplot2::geom_vline(xintercept = 0, color = pal$hair, linewidth = 0.45) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = conf_low, xmax = conf_high),
    orientation = "y",
    width = 0,
    color = pal$grey,
    linewidth = 0.35
  ) +
  ggplot2::geom_point(color = pal$orange, size = 1.35) +
  ggplot2::facet_wrap(~panel_label, nrow = 1, scales = "free_x") +
  ggplot2::labs(
    title = "C  Subject random-intercept estimates",
    subtitle = "Caterpillar plot shows subject-specific baseline shifts with approximate 95% intervals",
    x = "Random intercept",
    y = "Subject"
  ) +
  base_theme +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5.8))

panel_d <- ggplot2::ggplot(diagnostics, ggplot2::aes(x = fitted, y = observed)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, color = pal$hair, linewidth = 0.45) +
  ggplot2::geom_point(ggplot2::aes(color = days_from_baseline_subject), alpha = 0.78, size = 1.15) +
  ggplot2::scale_color_gradient(low = pal$sky, high = pal$blue) +
  ggplot2::facet_wrap(~panel_label, nrow = 1, scales = "free") +
  ggplot2::labs(
    title = "D  Observed vs fitted response",
    subtitle = "Color intensity follows days from subject baseline",
    x = "Fitted value",
    y = "Observed log10(Y + pseudocount)"
  ) +
  base_theme

summary_text <- paste0(
  "Diagnostic selection: top ", fmt_int(nrow(selected)),
  paste0(" nominal lmm_ok MAGs by raw p-value; ", sum(results$primary_analysis %in% TRUE & results$time_p_adj_bh_primary < 0.05, na.rm = TRUE), " primary MAGs passed BH-FDR q < 0.05. "),
  "Models were refit with the same fixed/random structure and bobyqa optimizer used by the HPC workflow. ",
  "Pseudocount = ", fmt_sci(pseudocount), "; refit warnings = ",
  fmt_int(sum(fit_summary$warning_count)), "; refit errors = ", fmt_int(sum(fit_summary$error_count)), "."
)

draw_figure <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  layout <- grid::grid.layout(
    nrow = 6,
    ncol = 1,
    heights = grid::unit(c(0.50, 1, 1, 1, 1, 0.42), c("in", "null", "null", "null", "null", "in"))
  )
  grid::pushViewport(grid::viewport(layout = layout))

  grid::grid.text(
    "Figure 8 | Diagnostic panel for selected MAG models",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.88, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 15.5, fontface = "bold", col = pal$ink),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
  )
  grid::grid.text(
    "Top exploratory per-MAG LMMs were checked for residual behavior, fitted-value patterns, and subject-level random effects.",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.20, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 9.0, col = pal$text),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
  )

  print(panel_a, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panel_b, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
  print(panel_c, vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1))
  print(panel_d, vp = grid::viewport(layout.pos.row = 5, layout.pos.col = 1))

  grid::grid.text(
    paste(strwrap(summary_text, width = 175), collapse = "\n"),
    x = grid::unit(0, "npc"),
    y = grid::unit(0.58, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 7.2, col = pal$muted, lineheight = 0.95),
    vp = grid::viewport(layout.pos.row = 6, layout.pos.col = 1)
  )

  grid::popViewport()
}

png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

png(filename = png_path, width = 14.0, height = 12.4, units = "in", res = 420, bg = "white")
draw_figure()
dev.off()

pdf(file = pdf_path, width = 14.0, height = 12.4, bg = "white", useDingbats = FALSE)
draw_figure()
dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = svg_path, width = 14.0, height = 12.4, bg = "white")
  draw_figure()
  dev.off()
}

diagnostics_path <- file.path(output_dir, paste0(output_stem, "_model_diagnostics.csv"))
random_effects_path <- file.path(output_dir, paste0(output_stem, "_random_intercepts.csv"))
summary_path <- file.path(output_dir, paste0(output_stem, "_diagnostic_model_summary.csv"))

utils::write.csv(diagnostics, diagnostics_path, row.names = FALSE)
utils::write.csv(random_effects, random_effects_path, row.names = FALSE)
utils::write.csv(fit_summary, summary_path, row.names = FALSE)

message("Figure 8 diagnostic outputs created:")
message(" - ", png_path)
message(" - ", pdf_path)
if (file.exists(svg_path)) message(" - ", svg_path)
message(" - ", diagnostics_path)
message(" - ", random_effects_path)
message(" - ", summary_path)
