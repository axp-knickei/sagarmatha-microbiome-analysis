# Public renderer: values come from external inputs; see docs/public-rendering.md.
# ============================================
# plot_figure_5_hpc_per_mag_lmertest_results.R
# ============================================
# Figure 5: top longitudinal MAG trajectories from the full T1-T10 per-MAG
# lmerTest LMM workflow.
#
# Selection rules are preserved; significance counts are calculated from supplied results.

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
output_stem <- "figure_5_hpc_per_mag_lmertest_results"

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

required_data_columns <- c("sampleID", "MAG", "Y", "Subject", "Timepoint", "days_from_baseline_subject", "altitude")
missing_data_columns <- setdiff(required_data_columns, names(df))
if (length(missing_data_columns) > 0) {
  stop("Missing required columns in cleaned_data.rds: ", paste(missing_data_columns, collapse = ", "), call. = FALSE)
}

required_result_columns <- c(
  "MAG", "status", "primary_analysis", "time_estimate", "time_std_error",
  "time_p_value", "time_p_adj_bh_primary", "zero_fraction", "n_subjects", "n_timepoints"
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

# Select balanced exploratory candidates. The thresholds avoid showing panels
# dominated by extreme missingness while preserving the top longitudinal signals.
candidate_pool <- results |>
  dplyr::filter(
    primary_analysis %in% TRUE,
    status == "lmm_ok",
    is.finite(time_estimate),
    is.finite(time_p_value),
    time_p_value < 0.05,
    n_subjects >= 8,
    n_timepoints >= 8,
    zero_fraction <= 0.75
  ) |>
  dplyr::mutate(
    direction = dplyr::if_else(time_estimate > 0, "Increasing", "Decreasing"),
    abs_time_estimate = abs(time_estimate)
  )

if (nrow(candidate_pool) == 0) {
  stop("No exploratory candidate MAGs passed the Figure 5 selection rule.", call. = FALSE)
}

top_per_direction <- 4L
selected <- candidate_pool |>
  dplyr::group_by(direction) |>
  dplyr::arrange(time_p_adj_bh_primary, dplyr::desc(abs_time_estimate), time_p_value, .by_group = TRUE) |>
  dplyr::slice_head(n = top_per_direction) |>
  dplyr::ungroup()

if (nrow(selected) < 2 * top_per_direction) {
  selected <- candidate_pool |>
    dplyr::arrange(time_p_adj_bh_primary, dplyr::desc(abs_time_estimate), time_p_value) |>
    dplyr::slice_head(n = min(8L, nrow(candidate_pool)))
}

selected <- selected |>
  dplyr::mutate(
    direction = factor(direction, levels = c("Increasing", "Decreasing")),
    direction_rank = dplyr::row_number(),
    mag_short = ifelse(nchar(MAG) > 32, paste0(substr(MAG, 1, 29), "..."), MAG),
    panel_label = paste0(
      mag_short,
      "\n",
      ifelse(time_estimate > 0, "beta = +", "beta = "),
      sprintf("%.3f", time_estimate),
      "; q = ",
      sprintf("%.3f", time_p_adj_bh_primary)
    )
  ) |>
  dplyr::arrange(direction, dplyr::desc(abs_time_estimate))

selected$panel_label <- factor(selected$panel_label, levels = selected$panel_label)

plot_df <- df |>
  dplyr::filter(MAG %in% selected$MAG) |>
  dplyr::inner_join(
    selected |>
      dplyr::select(MAG, direction, panel_label, time_estimate, time_p_value, time_p_adj_bh_primary, zero_fraction),
    by = "MAG"
  ) |>
  dplyr::mutate(
    Subject = factor(Subject, levels = sort(unique(Subject))),
    Timepoint = factor(Timepoint, levels = paste0("T", 1:10)),
    response_value = log10(Y + pseudocount)
  )

make_prediction <- function(mag_id) {
  mag_df <- plot_df[plot_df$MAG == mag_id, , drop = FALSE]
  mag_df$Subject <- droplevels(factor(mag_df$Subject))
  fit <- tryCatch(
    lme4::lmer(
      response_value ~ days_from_baseline_subject + (1 | Subject),
      data = mag_df,
      control = lme4::lmerControl(optimizer = "bobyqa")
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame())
  }

  day_seq <- seq(
    min(mag_df$days_from_baseline_subject, na.rm = TRUE),
    max(mag_df$days_from_baseline_subject, na.rm = TRUE),
    length.out = 80
  )
  pred_df <- data.frame(days_from_baseline_subject = day_seq)
  xmat <- stats::model.matrix(~ days_from_baseline_subject, data = pred_df)
  beta <- lme4::fixef(fit)
  vc <- as.matrix(stats::vcov(fit))
  fit_value <- as.numeric(xmat %*% beta)
  se_value <- sqrt(rowSums((xmat %*% vc) * xmat))

  meta <- selected[selected$MAG == mag_id, , drop = FALSE]
  data.frame(
    MAG = mag_id,
    days_from_baseline_subject = day_seq,
    fitted = fit_value,
    conf_low = fit_value - 1.96 * se_value,
    conf_high = fit_value + 1.96 * se_value,
    direction = meta$direction[[1]],
    panel_label = meta$panel_label[[1]],
    stringsAsFactors = FALSE
  )
}

prediction_df <- dplyr::bind_rows(lapply(selected$MAG, make_prediction))
prediction_df$panel_label <- factor(prediction_df$panel_label, levels = levels(selected$panel_label))
prediction_df$direction <- factor(prediction_df$direction, levels = c("Increasing", "Decreasing"))

timepoint_summary <- plot_df |>
  dplyr::distinct(Timepoint, altitude) |>
  dplyr::group_by(Timepoint) |>
  dplyr::summarise(altitude = stats::median(altitude, na.rm = TRUE), .groups = "drop")

total_candidates <- sum(results$primary_analysis %in% TRUE & results$time_p_value < 0.05, na.rm = TRUE)
fdr_significant <- sum(results$primary_significant %in% TRUE, na.rm = TRUE)

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_sci <- function(x) formatC(x, format = "e", digits = 2)

pal <- list(
  ink = "#111827",
  text = "#374151",
  muted = "#6B7280",
  hair = "#D1D5DB",
  faint = "#F3F4F6",
  blue = "#0072B2",
  sky = "#DCEEF8",
  orange = "#D55E00",
  peach = "#FFF2E8",
  red = "#CC3311",
  rose = "#FFF0EE",
  green = "#009E73"
)

direction_palette <- c("Increasing" = pal$blue, "Decreasing" = pal$red)
direction_fill <- c("Increasing" = pal$sky, "Decreasing" = pal$rose)

base_theme <- ggplot2::theme_minimal(base_size = 9.4, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11.0, color = pal$ink),
    plot.subtitle = ggplot2::element_text(size = 8.0, color = pal$muted, margin = ggplot2::margin(b = 5)),
    axis.title = ggplot2::element_text(size = 8.4, color = pal$text),
    axis.text = ggplot2::element_text(size = 7.3, color = pal$text),
    strip.text = ggplot2::element_text(face = "bold", size = 7.1, color = pal$ink, lineheight = 0.9),
    strip.background = ggplot2::element_rect(fill = "#F8FAFC", color = pal$hair, linewidth = 0.35),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.22),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 8),
    legend.text = ggplot2::element_text(size = 7.5),
    plot.margin = ggplot2::margin(5, 6, 5, 6)
  )

trajectory_plot <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = days_from_baseline_subject, y = response_value)
) +
  ggplot2::geom_line(
    ggplot2::aes(group = Subject),
    color = "#9CA3AF",
    alpha = 0.34,
    linewidth = 0.35
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = direction),
    shape = 21,
    size = 1.7,
    color = "white",
    stroke = 0.25,
    alpha = 0.78
  ) +
  ggplot2::geom_ribbon(
    data = prediction_df,
    ggplot2::aes(
      x = days_from_baseline_subject,
      y = fitted,
      ymin = conf_low,
      ymax = conf_high,
      fill = direction
    ),
    inherit.aes = FALSE,
    alpha = 0.22,
    color = NA
  ) +
  ggplot2::geom_line(
    data = prediction_df,
    ggplot2::aes(x = days_from_baseline_subject, y = fitted, color = direction),
    inherit.aes = FALSE,
    linewidth = 0.85
  ) +
  ggplot2::facet_wrap(~panel_label, nrow = 2, scales = "free_y") +
  ggplot2::scale_color_manual(values = direction_palette, name = "Selected trend") +
  ggplot2::scale_fill_manual(values = direction_fill, name = "Selected trend") +
  ggplot2::scale_x_continuous(breaks = seq(0, 90, by = 15)) +
  ggplot2::labs(
    title = "A  Exploratory top longitudinal MAG trajectories",
    subtitle = "Thin grey lines are subjects; colored line and ribbon are fixed-effect LMM trend and 95% CI",
    x = "Days from subject baseline",
    y = "log10(Y + pseudocount)"
  ) +
  base_theme +
  ggplot2::theme(
    panel.spacing = grid::unit(0.62, "lines"),
    legend.position = "bottom"
  )

raw_summary <- plot_df |>
  dplyr::group_by(direction, panel_label, Timepoint) |>
  dplyr::summarise(
    mean_raw = mean(Y, na.rm = TRUE),
    se_raw = stats::sd(Y, na.rm = TRUE) / sqrt(sum(!is.na(Y))),
    days = stats::median(days_from_baseline_subject, na.rm = TRUE),
    .groups = "drop"
  )

raw_plot <- ggplot2::ggplot(raw_summary, ggplot2::aes(x = days, y = mean_raw, color = direction, fill = direction)) +
  ggplot2::geom_line(linewidth = 0.58) +
  ggplot2::geom_point(shape = 21, size = 1.65, color = "white", stroke = 0.25) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = pmax(0, mean_raw - se_raw), ymax = mean_raw + se_raw),
    width = 1.0,
    linewidth = 0.28,
    alpha = 0.78
  ) +
  ggplot2::facet_wrap(~panel_label, nrow = 2, scales = "free_y") +
  ggplot2::scale_color_manual(values = direction_palette, guide = "none") +
  ggplot2::scale_fill_manual(values = direction_fill, guide = "none") +
  ggplot2::scale_x_continuous(breaks = seq(0, 90, by = 30)) +
  ggplot2::labs(
    title = "B  Raw abundance summaries for the same MAGs",
    subtitle = "Points show timepoint mean raw relative abundance; error bars show standard error",
    x = "Days from subject baseline",
    y = "Mean raw abundance"
  ) +
  base_theme +
  ggplot2::theme(
    panel.spacing = grid::unit(0.62, "lines"),
    strip.text = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank()
  )

selection_text <- paste(
  strwrap(
    paste0(
      "Exploratory nominal-candidate selection ",
      "(lmm_ok, raw p < 0.05, n_subjects >= 8, n_timepoints >= 8, zero_fraction <= 0.75), ",
      "balanced across increasing and decreasing estimates. Nominal candidates = ",
      fmt_int(total_candidates), "; FDR-significant MAGs = ", fmt_int(fdr_significant),
      "; pseudocount = ", fmt_sci(pseudocount), "."
    ),
    width = 170
  ),
  collapse = "\n"
)

draw_figure <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  layout <- grid::grid.layout(
    nrow = 4,
    ncol = 1,
    heights = grid::unit(c(0.50, 1.42, 0.92, 0.46), c("in", "null", "null", "in"))
  )
  grid::pushViewport(grid::viewport(layout = layout))

  grid::grid.text(
    "Figure 5 | Top longitudinal MAG trajectories",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.88, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 15.5, fontface = "bold", col = pal$ink),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
  )
  grid::grid.text(
    "Selected exploratory MAGs show interpretable subject-level trajectories across expedition time; these are not FDR-confirmed discoveries.",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.20, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 9.0, col = pal$text),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
  )

  print(trajectory_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(raw_plot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))

  grid::grid.text(
    selection_text,
    x = grid::unit(0, "npc"),
    y = grid::unit(0.58, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 7.1, col = pal$muted, lineheight = 0.95),
    vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1)
  )

  grid::popViewport()
}

png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

png(filename = png_path, width = 13.2, height = 11.0, units = "in", res = 420, bg = "white")
draw_figure()
dev.off()

pdf(file = pdf_path, width = 13.2, height = 11.0, bg = "white", useDingbats = FALSE)
draw_figure()
dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = svg_path, width = 13.2, height = 11.0, bg = "white")
  draw_figure()
  dev.off()
}

selected_path <- file.path(output_dir, paste0(output_stem, "_selected_mags.csv"))
summary_path <- file.path(output_dir, paste0(output_stem, "_trajectory_summary.csv"))

utils::write.csv(
  selected |>
    dplyr::select(
      MAG, direction, time_estimate, time_std_error, time_p_value,
      time_p_adj_bh_primary, zero_fraction, n_subjects, n_timepoints
    ),
  selected_path,
  row.names = FALSE
)

utils::write.csv(
  data.frame(
    metric = c(
      "selected_mags",
      "nominal_candidate_mags",
      "fdr_significant_mags",
      "pseudocount",
      "selection_status",
      "minimum_n_subjects",
      "minimum_n_timepoints",
      "maximum_zero_fraction"
    ),
    value = c(
      nrow(selected),
      total_candidates,
      fdr_significant,
      pseudocount,
      "exploratory_nominal_not_fdr_confirmed",
      8,
      8,
      0.75
    )
  ),
  summary_path,
  row.names = FALSE
)

message("Figure 5 trajectory outputs created:")
message(" - ", png_path)
message(" - ", pdf_path)
if (file.exists(svg_path)) message(" - ", svg_path)
message(" - ", selected_path)
message(" - ", summary_path)
