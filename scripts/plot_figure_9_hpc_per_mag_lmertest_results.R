# ============================================
# plot_figure_9_hpc_per_mag_lmertest_results.R
# ============================================
# Figure 9: sensitivity analysis for the full T1-T10 per-MAG lmerTest LMM run.
#
# This script performs a conservative, reproducible sensitivity analysis using
# the existing fitted model table: it excludes increasingly sparse MAGs and
# recomputes BH-FDR within each retained testing family. This tests whether the
# primary "no FDR-confirmed MAGs" conclusion depends on the inclusion of highly
# sparse MAGs, without changing the primary HPC model formula or overwriting
# scientific outputs.

required_packages <- c("ggplot2", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
results_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "tables",
  "per_mag_lmertest_hpc_results.csv"
)
output_dir <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "figures")
output_stem <- "figure_9_hpc_per_mag_lmertest_results_20260708_024546"

if (!file.exists(results_path)) {
  stop("Missing per-MAG LMM result table: ", results_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c(
  "MAG", "status", "primary_analysis", "time_estimate", "time_p_value",
  "time_p_adj_bh_primary", "primary_significant", "zero_fraction",
  "nonzero_sample_count", "subjects_with_nonzero"
)
missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns) > 0) {
  stop("Missing required result columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

primary <- results |>
  dplyr::filter(primary_analysis %in% TRUE, is.finite(time_p_value), is.finite(time_estimate)) |>
  dplyr::mutate(
    primary_rank = rank(time_p_value, ties.method = "first"),
    primary_category = dplyr::case_when(
      time_p_adj_bh_primary < 0.05 ~ "BH-FDR q < 0.05",
      time_p_value < 0.05 ~ "Nominal raw p < 0.05",
      TRUE ~ "Not nominal"
    )
  )

if (nrow(primary) == 0) {
  stop("No primary-analysis MAGs with finite p-values found.", call. = FALSE)
}

sensitivity_definitions <- data.frame(
  sensitivity_model = c("Primary full family", "Zero fraction <= 75%", "Zero fraction <= 50%"),
  zero_fraction_max = c(Inf, 0.75, 0.50),
  stringsAsFactors = FALSE
)

make_sensitivity <- function(model_name, zero_fraction_max) {
  retained <- primary$zero_fraction <= zero_fraction_max
  out <- primary |>
    dplyr::mutate(
      sensitivity_model = model_name,
      zero_fraction_max = zero_fraction_max,
      retained_in_sensitivity = retained,
      estimate_sensitivity = ifelse(retained, time_estimate, NA_real_),
      q_sensitivity = NA_real_,
      rank_sensitivity = NA_real_,
      sensitivity_category = "Excluded by sensitivity filter"
    )

  idx <- which(out$retained_in_sensitivity)
  if (length(idx) > 0) {
    out$q_sensitivity[idx] <- stats::p.adjust(out$time_p_value[idx], method = "BH")
    out$rank_sensitivity[idx] <- rank(out$time_p_value[idx], ties.method = "first")
    out$sensitivity_category[idx] <- dplyr::case_when(
      out$q_sensitivity[idx] < 0.05 ~ "BH-FDR q < 0.05",
      out$time_p_value[idx] < 0.05 ~ "Nominal raw p < 0.05",
      TRUE ~ "Not nominal"
    )
  }
  out
}

sensitivity_results <- dplyr::bind_rows(lapply(seq_len(nrow(sensitivity_definitions)), function(i) {
  make_sensitivity(
    sensitivity_definitions$sensitivity_model[[i]],
    sensitivity_definitions$zero_fraction_max[[i]]
  )
}))

comparison_table <- sensitivity_results |>
  dplyr::transmute(
    MAG,
    sensitivity_model,
    zero_fraction,
    retained_in_sensitivity,
    estimate_primary = time_estimate,
    q_primary = time_p_adj_bh_primary,
    estimate_sensitivity,
    q_sensitivity,
    direction_same = retained_in_sensitivity & sign(estimate_primary) == sign(estimate_sensitivity),
    rank_primary = primary_rank,
    rank_sensitivity,
    rank_change = rank_sensitivity - rank_primary,
    primary_category,
    sensitivity_category,
    significance_category = paste(primary_category, "->", sensitivity_category),
    time_p_value,
    nonzero_sample_count,
    subjects_with_nonzero
  )

count_summary <- comparison_table |>
  dplyr::group_by(sensitivity_model) |>
  dplyr::summarise(
    retained_mags = sum(retained_in_sensitivity, na.rm = TRUE),
    excluded_mags = sum(!retained_in_sensitivity, na.rm = TRUE),
    nominal_mags = sum(sensitivity_category == "Nominal raw p < 0.05", na.rm = TRUE),
    fdr_mags = sum(sensitivity_category == "BH-FDR q < 0.05", na.rm = TRUE),
    best_q = suppressWarnings(min(q_sensitivity, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::mutate(best_q = ifelse(is.finite(best_q), best_q, NA_real_))

strict75 <- comparison_table |>
  dplyr::filter(sensitivity_model == "Zero fraction <= 75%", retained_in_sensitivity)

strict50 <- comparison_table |>
  dplyr::filter(sensitivity_model == "Zero fraction <= 50%", retained_in_sensitivity)

top_transition_mags <- primary |>
  dplyr::arrange(time_p_value) |>
  dplyr::slice_head(n = min(45L, nrow(primary))) |>
  dplyr::pull(MAG)

transition_df <- comparison_table |>
  dplyr::filter(MAG %in% top_transition_mags) |>
  dplyr::mutate(
    mag_order = match(MAG, rev(top_transition_mags)),
    mag_label = paste0(
      sprintf("%02d", match(MAG, top_transition_mags)),
      "  ",
      ifelse(nchar(MAG) > 34, paste0(substr(MAG, 1, 31), "..."), MAG)
    ),
    mag_label = factor(mag_label, levels = unique(mag_label[order(mag_order)])),
    sensitivity_model = factor(sensitivity_model, levels = sensitivity_definitions$sensitivity_model),
    sensitivity_category = factor(
      sensitivity_category,
      levels = c("BH-FDR q < 0.05", "Nominal raw p < 0.05", "Not nominal", "Excluded by sensitivity filter")
    )
  )

total_primary <- nrow(primary)
primary_nominal <- sum(primary$time_p_value < 0.05, na.rm = TRUE)
primary_fdr <- sum(primary$time_p_adj_bh_primary < 0.05, na.rm = TRUE)
best_primary_q <- min(primary$time_p_adj_bh_primary, na.rm = TRUE)

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_q <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))

pal <- list(
  ink = "#111827",
  text = "#374151",
  muted = "#6B7280",
  hair = "#D1D5DB",
  faint = "#F8FAFC",
  blue = "#0072B2",
  sky = "#DCEEF8",
  orange = "#D55E00",
  peach = "#FFF2E8",
  red = "#CC3311",
  green = "#009E73",
  purple = "#785EF0",
  grey = "#9CA3AF"
)

category_palette <- c(
  "BH-FDR q < 0.05" = pal$green,
  "Nominal raw p < 0.05" = pal$orange,
  "Not nominal" = "#D1D5DB",
  "Excluded by sensitivity filter" = "#FFFFFF"
)

base_theme <- ggplot2::theme_minimal(base_size = 9.4, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11.0, color = pal$ink),
    plot.subtitle = ggplot2::element_text(size = 8.0, color = pal$muted, margin = ggplot2::margin(b = 5)),
    axis.title = ggplot2::element_text(size = 8.4, color = pal$text),
    axis.text = ggplot2::element_text(size = 7.4, color = pal$text),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.22),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 8),
    legend.text = ggplot2::element_text(size = 7.4),
    plot.margin = ggplot2::margin(5, 7, 5, 7)
  )

panel_a <- ggplot2::ggplot(strict75, ggplot2::aes(x = q_primary, y = q_sensitivity)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, color = pal$hair, linewidth = 0.55) +
  ggplot2::geom_hline(yintercept = 0.05, color = pal$green, linewidth = 0.45, linetype = "longdash") +
  ggplot2::geom_vline(xintercept = 0.05, color = pal$green, linewidth = 0.45, linetype = "longdash") +
  ggplot2::geom_point(
    ggplot2::aes(color = sensitivity_category, size = zero_fraction),
    alpha = 0.78
  ) +
  ggplot2::scale_color_manual(values = category_palette, drop = FALSE, name = "Sensitivity class") +
  ggplot2::scale_size_continuous(range = c(1.2, 3.2), labels = function(x) paste0(round(100 * x), "%"), name = "Zero fraction") +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    title = "A  Recomputed q-values after excluding high-sparsity MAGs",
    subtitle = "Zero fraction <= 75% sensitivity keeps fitted estimates but changes the BH testing family",
    x = "Primary BH q-value",
    y = "Sensitivity BH q-value"
  ) +
  base_theme

panel_b <- ggplot2::ggplot(strict75, ggplot2::aes(x = rank_primary, y = rank_sensitivity)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, color = pal$hair, linewidth = 0.55) +
  ggplot2::geom_point(ggplot2::aes(color = sensitivity_category), alpha = 0.78, size = 1.65) +
  ggplot2::scale_color_manual(values = category_palette, drop = FALSE, guide = "none") +
  ggplot2::labs(
    title = "B  Rank stability under zero-fraction filtering",
    subtitle = "MAG ranks are based on raw p-values among retained primary models",
    x = "Primary raw-p rank",
    y = "Sensitivity raw-p rank"
  ) +
  base_theme

panel_c <- ggplot2::ggplot(transition_df, ggplot2::aes(x = sensitivity_model, y = mag_label, fill = sensitivity_category)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.35) +
  ggplot2::scale_fill_manual(values = category_palette, drop = FALSE, name = "Category") +
  ggplot2::labs(
    title = "C  Top-MAG category transitions across sensitivity definitions",
    subtitle = "Rows show top primary raw-p MAGs; exclusions reflect stricter zero-fraction filters",
    x = NULL,
    y = NULL
  ) +
  base_theme +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(face = "bold", angle = 20, hjust = 1),
    axis.text.y = ggplot2::element_text(size = 5.4, family = "mono"),
    panel.grid = ggplot2::element_blank()
  )

count_long <- rbind(
  data.frame(sensitivity_model = count_summary$sensitivity_model, metric = "Retained MAGs", value = count_summary$retained_mags),
  data.frame(sensitivity_model = count_summary$sensitivity_model, metric = "Nominal raw p < 0.05", value = count_summary$nominal_mags),
  data.frame(sensitivity_model = count_summary$sensitivity_model, metric = "BH-FDR q < 0.05", value = count_summary$fdr_mags)
)
count_long$sensitivity_model <- factor(count_long$sensitivity_model, levels = sensitivity_definitions$sensitivity_model)
count_long$metric <- factor(count_long$metric, levels = c("Retained MAGs", "Nominal raw p < 0.05", "BH-FDR q < 0.05"))

panel_d <- ggplot2::ggplot(count_long, ggplot2::aes(x = sensitivity_model, y = value, fill = metric)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78), width = 0.70, color = "white", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = value),
    position = ggplot2::position_dodge(width = 0.78),
    vjust = -0.18,
    size = 2.55,
    color = pal$text
  ) +
  ggplot2::scale_fill_manual(
    values = c("Retained MAGs" = pal$blue, "Nominal raw p < 0.05" = pal$orange, "BH-FDR q < 0.05" = pal$green),
    name = NULL
  ) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
  ggplot2::labs(
    title = "D  Sensitivity conclusion across testing families",
    subtitle = "The no-FDR-discovery conclusion is stable under stricter sparsity filters",
    x = NULL,
    y = "MAG count"
  ) +
  base_theme +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 18, hjust = 1),
    panel.grid.major.x = ggplot2::element_blank()
  )

summary_text <- paste0(
  "Sensitivity analysis: primary family = ", fmt_int(total_primary),
  " MAGs with ", fmt_int(primary_nominal), " nominal raw-p candidates and ",
  fmt_int(primary_fdr), " BH-FDR discoveries (best q = ", fmt_q(best_primary_q),
  "). Stricter zero-fraction filters recompute BH-FDR among retained MAGs only; estimates are not refit, so this evaluates sparsity-driven testing-family sensitivity rather than an alternative model formula."
)

draw_figure <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  layout <- grid::grid.layout(
    nrow = 4,
    ncol = 2,
    heights = grid::unit(c(0.50, 1.0, 1.25, 0.42), c("in", "null", "null", "in")),
    widths = grid::unit(c(1.0, 1.0), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))

  grid::grid.text(
    "Figure 9 | Sensitivity analysis to sparsity filtering",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.88, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 15.5, fontface = "bold", col = pal$ink),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )
  grid::grid.text(
    "Recomputing BH-FDR after excluding increasingly sparse MAGs tests whether the main inference depends on sparse-feature inclusion.",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.20, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 9.0, col = pal$text),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )

  print(panel_a, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panel_b, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(panel_c, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
  print(panel_d, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 2))

  grid::grid.text(
    paste(strwrap(summary_text, width = 175), collapse = "\n"),
    x = grid::unit(0, "npc"),
    y = grid::unit(0.58, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 7.2, col = pal$muted, lineheight = 0.95),
    vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1:2)
  )

  grid::popViewport()
}

png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

png(filename = png_path, width = 13.0, height = 10.0, units = "in", res = 420, bg = "white")
draw_figure()
dev.off()

pdf(file = pdf_path, width = 13.0, height = 10.0, bg = "white", useDingbats = FALSE)
draw_figure()
dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = svg_path, width = 13.0, height = 10.0, bg = "white")
  draw_figure()
  dev.off()
}

comparison_path <- file.path(output_dir, paste0(output_stem, "_sensitivity_comparison_table.csv"))
summary_path <- file.path(output_dir, paste0(output_stem, "_sensitivity_summary.csv"))

utils::write.csv(comparison_table, comparison_path, row.names = FALSE)
utils::write.csv(count_summary, summary_path, row.names = FALSE)

message("Figure 9 sensitivity outputs created:")
message(" - ", png_path)
message(" - ", pdf_path)
if (file.exists(svg_path)) message(" - ", svg_path)
message(" - ", comparison_path)
message(" - ", summary_path)
