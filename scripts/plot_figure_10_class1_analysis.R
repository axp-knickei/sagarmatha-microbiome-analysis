# ============================================
# plot_figure_10_class1_analysis.R
# ============================================
# Supplementary/Main Figure 10: effect-size heatmap for nominal Class 1 MAGs.
# Rows are MAGs with at least one nominal raw-p candidate. Columns are the six
# predefined Class 1 windows. Cell color is the window_order estimate from the
# Class 1 LMM. Dots mark raw p < 0.05; stars would mark global q < 0.05.

required_packages <- c("readr", "dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- getwd()
results_path <- file.path(
  repo_root,
  "results",
  "class1_window_lmertest",
  "tables",
  "class1_window_lmertest_results.csv"
)
windows_path <- file.path(
  repo_root,
  "results",
  "metadata_expedition_windows",
  "proposed_biological_windows.csv"
)
output_dir <- file.path(repo_root, "results", "class1_window_lmertest", "figures")

if (!file.exists(results_path)) {
  stop("Missing Class 1 results table: ", results_path, call. = FALSE)
}
if (!file.exists(windows_path)) {
  stop("Missing Class 1 window definition table: ", windows_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

class1_results <- readr::read_csv(results_path, show_col_types = FALSE)
class1_windows <- readr::read_csv(windows_path, show_col_types = FALSE)

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_q <- function(x) ifelse(is.na(x), "NA", sprintf("%.4f", x))

window_order <- class1_windows$window_label
compact_window_names <- c(
  initial_ascent_to_high_altitude = "Initial ascent",
  ascent_and_intermediate_adjustment = "Ascent adjustment",
  high_climb_6500m_exposure = "6500 m climb",
  pre_summit_staging_to_post_summit = "Pre/post summit",
  summit_push_and_descent_recovery_transition = "Summit + descent",
  post_expedition_low_altitude_recovery = "Recovery"
)
window_labels <- setNames(
  paste0(
    "W", seq_along(window_order), "\n",
    gsub(" -> ", "-", class1_windows$forward_timepoint_path),
    "\n",
    unname(compact_window_names[window_order])
  ),
  window_order
)

# No taxonomy columns are currently present in the Class 1 output or processed
# data, so MAG identifiers are used as row labels. This block is deliberately
# simple and will continue to work if taxonomy columns are added later.
taxonomy_columns <- grep(
  "taxonomy|taxon|phylum|family|genus|species|^class$|^order$",
  names(class1_results),
  ignore.case = TRUE,
  value = TRUE
)
has_taxonomy <- length(taxonomy_columns) > 0

nominal_results <- class1_results |>
  dplyr::filter(primary_analysis %in% TRUE, nominal_p_lt_0_05 %in% TRUE)

if (nrow(nominal_results) == 0) {
  stop("No nominal Class 1 candidates found for effect-size heatmap.", call. = FALSE)
}

nominal_mags <- unique(nominal_results$MAG)

row_order <- class1_results |>
  dplyr::filter(MAG %in% nominal_mags) |>
  dplyr::group_by(MAG) |>
  dplyr::summarise(
    best_q = suppressWarnings(min(q_global_class1_bh[primary_analysis %in% TRUE], na.rm = TRUE)),
    best_raw_p = suppressWarnings(min(raw_p_value[primary_analysis %in% TRUE], na.rm = TRUE)),
    max_abs_effect = suppressWarnings(max(abs(estimate_window_order[primary_analysis %in% TRUE]), na.rm = TRUE)),
    nominal_count = sum(primary_analysis %in% TRUE & nominal_p_lt_0_05 %in% TRUE, na.rm = TRUE),
    best_window = window_label[which.min(ifelse(primary_analysis %in% TRUE, raw_p_value, Inf))][1],
    .groups = "drop"
  ) |>
  dplyr::mutate(
    best_q = ifelse(is.infinite(best_q), NA_real_, best_q),
    best_raw_p = ifelse(is.infinite(best_raw_p), NA_real_, best_raw_p),
    max_abs_effect = ifelse(is.infinite(max_abs_effect), NA_real_, max_abs_effect),
    best_window = factor(best_window, levels = window_order)
  ) |>
  dplyr::arrange(best_window, best_q, best_raw_p, dplyr::desc(nominal_count), dplyr::desc(max_abs_effect), MAG)

mag_levels <- row_order$MAG

plot_data <- class1_results |>
  dplyr::filter(MAG %in% nominal_mags) |>
  dplyr::mutate(
    MAG = factor(MAG, levels = rev(mag_levels)),
    window_label = factor(window_label, levels = window_order),
    display_estimate = ifelse(primary_analysis %in% TRUE, estimate_window_order, NA_real_),
    nominal = primary_analysis %in% TRUE & nominal_p_lt_0_05 %in% TRUE,
    fdr_significant = primary_analysis %in% TRUE & q_global_class1_bh < 0.05
  )

n_nominal_mags <- length(nominal_mags)
n_nominal_tests <- sum(plot_data$nominal, na.rm = TRUE)
n_fdr <- sum(plot_data$fdr_significant, na.rm = TRUE)
best_q <- min(plot_data$q_global_class1_bh[plot_data$primary_analysis %in% TRUE], na.rm = TRUE)
if (!is.finite(best_q)) {
  best_q <- NA_real_
}

# Symmetric color limits keep effect-size sign interpretable. Extreme values are
# squished by the scale limits rather than changing the palette range.
effect_limit <- stats::quantile(abs(plot_data$display_estimate), probs = 0.98, na.rm = TRUE)
effect_limit <- max(0.25, as.numeric(effect_limit))
effect_limit <- ceiling(effect_limit * 10) / 10

row_label_df <- row_order |>
  dplyr::mutate(
    MAG = factor(MAG, levels = rev(mag_levels)),
    row_label = as.character(MAG)
  )

audit_text <- paste0(
  fmt_int(n_nominal_tests),
  " nominal MAG-window cells across ",
  fmt_int(n_nominal_mags),
  " unique MAGs; ",
  fmt_int(n_fdr),
  " global FDR-significant cells; best q = ",
  fmt_q(best_q)
)

taxonomy_note <- if (has_taxonomy) {
  paste0("Taxonomy-like columns detected: ", paste(taxonomy_columns, collapse = ", "))
} else {
  "No taxonomy annotation columns detected in current Class 1 outputs; rows are labeled by MAG identifier."
}

figure <- ggplot2::ggplot(plot_data, ggplot2::aes(x = window_label, y = MAG)) +
  ggplot2::geom_tile(
    ggplot2::aes(fill = display_estimate),
    color = "white",
    linewidth = 0.2,
    width = 0.96,
    height = 0.96
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(plot_data, nominal),
    shape = 21,
    size = 1.65,
    stroke = 0.55,
    fill = "white",
    color = "#111827"
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(plot_data, fdr_significant),
    shape = 8,
    size = 2.4,
    stroke = 0.85,
    color = "#111827"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#D55E00",
    mid = "white",
    high = "#0072B2",
    midpoint = 0,
    limits = c(-effect_limit, effect_limit),
    na.value = "#e5e7eb",
    name = "window_order\nestimate"
  ) +
  ggplot2::scale_x_discrete(labels = unname(window_labels[window_order]), position = "top") +
  ggplot2::scale_y_discrete(labels = function(x) x) +
  ggplot2::labs(
    title = "Main Figure 10 | Effect-size heatmap for nominal Class 1 MAGs",
    subtitle = paste0(
      "Cell color shows the Class 1 LMM window_order estimate; dots mark raw p < 0.05; stars would mark global q < 0.05.\n",
      audit_text
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      taxonomy_note,
      "\nGrey cells indicate no primary lmm_ok estimate for that MAG-window. ",
      "Absence of stars is expected here because no Class 1 test passed global BH-FDR correction."
    )
  ) +
  ggplot2::theme_minimal(base_size = 10.2) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = 15.5, color = "#111827"),
    plot.subtitle = ggplot2::element_text(size = 9.4, color = "#4b5563", lineheight = 0.95, margin = ggplot2::margin(b = 8)),
    axis.text.x = ggplot2::element_text(size = 8.4, face = "bold", color = "#111827", lineheight = 0.86, margin = ggplot2::margin(b = 8)),
    axis.text.y = ggplot2::element_text(size = 4.2, color = "#374151"),
    legend.position = "right",
    legend.title = ggplot2::element_text(face = "bold", size = 8.8),
    legend.text = ggplot2::element_text(size = 8),
    plot.caption = ggplot2::element_text(size = 8.0, color = "#4b5563", hjust = 0, margin = ggplot2::margin(t = 10)),
    plot.margin = ggplot2::margin(14, 18, 16, 18)
  )

png_path <- file.path(output_dir, "figure_10_class1_analysis_effect_size_heatmap_nominal_mags.png")
pdf_path <- file.path(output_dir, "figure_10_class1_analysis_effect_size_heatmap_nominal_mags.pdf")
svg_path <- file.path(output_dir, "figure_10_class1_analysis_effect_size_heatmap_nominal_mags.svg")

ggplot2::ggsave(png_path, figure, width = 10.8, height = 15.2, dpi = 360, bg = "white", limitsize = FALSE)
ggplot2::ggsave(pdf_path, figure, width = 10.8, height = 15.2, bg = "white", limitsize = FALSE)
ggplot2::ggsave(svg_path, figure, width = 10.8, height = 15.2, bg = "white", limitsize = FALSE, device = svglite::svglite)

message("Main Figure 10 Class 1 nominal-MAG effect-size heatmap created:")
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Inputs:")
message(" - ", results_path)
message(" - ", windows_path)
message("Summary:")
message(" - Nominal MAG-window cells: ", fmt_int(n_nominal_tests))
message(" - Unique nominal MAGs: ", fmt_int(n_nominal_mags))
message(" - Global Class 1 BH-FDR significant cells: ", fmt_int(n_fdr))
message(" - Best q_global_class1_bh: ", fmt_q(best_q))
message(" - Effect color limit: +/-", effect_limit)
message(" - ", taxonomy_note)
