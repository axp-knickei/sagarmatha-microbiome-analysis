# ============================================
# plot_figure_3_class1_analysis.R
# ============================================
# Main Figure 3: MAG-window testing matrix for the Class 1 analysis.
# Each cell is one MAG x predefined-window test, colored by model status.
# The goal is model auditability: where were LMMs usable, singular, or too sparse?

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

fmt_int <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)

window_order <- class1_windows$window_label
compact_window_names <- c(
  initial_ascent_to_high_altitude = "Initial ascent",
  ascent_and_intermediate_adjustment = "Ascent adjustment",
  high_climb_6500m_exposure = "6500 m climb",
  pre_summit_staging_to_post_summit = "Pre/post summit",
  summit_push_and_descent_recovery_transition = "Summit + descent",
  post_expedition_low_altitude_recovery = "Recovery"
)
window_label_map <- setNames(
  paste0(
    "W", seq_along(class1_windows$window_label),
    "\n",
    gsub(" -> ", "-", class1_windows$forward_timepoint_path),
    "\n",
    unname(compact_window_names[class1_windows$window_label])
  ),
  class1_windows$window_label
)

status_levels <- c("lmm_ok", "lmm_singular", "too_sparse", "low_variance", "model_failed")
status_palette <- c(
  lmm_ok = "#0072B2",       # Okabe-Ito blue
  lmm_singular = "#E69F00", # Okabe-Ito orange
  too_sparse = "#D55E00",   # vermillion
  low_variance = "#CC79A7", # reddish purple
  model_failed = "#4B5563"  # neutral dark grey
)

plot_data <- class1_results |>
  dplyr::mutate(
    status = factor(status, levels = status_levels),
    window_label = factor(window_label, levels = window_order)
  )

# Order MAG rows by modelability so the matrix reveals sparse-test structure.
mag_order <- plot_data |>
  dplyr::group_by(MAG) |>
  dplyr::summarise(
    n_lmm_ok = sum(status == "lmm_ok", na.rm = TRUE),
    n_singular = sum(status == "lmm_singular", na.rm = TRUE),
    n_sparse = sum(status == "too_sparse", na.rm = TRUE),
    n_failed_or_other = sum(status %in% c("low_variance", "model_failed"), na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(n_lmm_ok),
    n_sparse,
    n_singular,
    n_failed_or_other,
    MAG
  ) |>
  dplyr::mutate(mag_row = dplyr::row_number())

plot_data <- plot_data |>
  dplyr::left_join(mag_order |> dplyr::select(MAG, mag_row), by = "MAG") |>
  dplyr::mutate(
    window_index = as.integer(window_label)
  )

status_counts <- plot_data |>
  dplyr::count(status, name = "n") |>
  dplyr::mutate(status = as.character(status))

status_lookup <- setNames(status_counts$n, status_counts$status)
count_status <- function(status) {
  idx <- match(status, names(status_lookup))
  if (is.na(idx)) {
    return(0L)
  }
  value <- status_lookup[[idx]]
  if (is.null(value) || is.na(value)) 0L else as.integer(value)
}

n_mags <- dplyr::n_distinct(plot_data$MAG)
n_windows <- dplyr::n_distinct(plot_data$window_label)
n_tests <- nrow(plot_data)
n_ok <- count_status("lmm_ok")
n_singular <- count_status("lmm_singular")
n_sparse <- count_status("too_sparse")
n_low_variance <- count_status("low_variance")
n_failed <- count_status("model_failed")

window_summary <- plot_data |>
  dplyr::group_by(window_label, window_index) |>
  dplyr::summarise(
    lmm_ok = sum(status == "lmm_ok", na.rm = TRUE),
    too_sparse = sum(status == "too_sparse", na.rm = TRUE),
    lmm_singular = sum(status == "lmm_singular", na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    label = paste0("OK ", lmm_ok, "\nSparse ", too_sparse),
    y = n_mags + 22
  )

legend_df <- data.frame(
  status = factor(status_levels, levels = status_levels),
  x = seq_along(status_levels),
  y = 1,
  label = c(
    paste0("lmm_ok\n", fmt_int(n_ok)),
    paste0("singular\n", fmt_int(n_singular)),
    paste0("too_sparse\n", fmt_int(n_sparse)),
    paste0("low_variance\n", fmt_int(n_low_variance)),
    paste0("model_failed\n", fmt_int(n_failed))
  )
)

audit_text <- paste0(
  "Total tests attempted: ", fmt_int(n_tests),
  "   |   lmm_ok: ", fmt_int(n_ok),
  "   |   lmm_singular: ", fmt_int(n_singular),
  "   |   too_sparse: ", fmt_int(n_sparse),
  "   |   low_variance: ", fmt_int(n_low_variance),
  "   |   model_failed: ", fmt_int(n_failed)
)

figure <- ggplot2::ggplot() +
  ggplot2::annotate(
    "rect",
    xmin = 0.45,
    xmax = n_windows + 0.55,
    ymin = 0.5,
    ymax = n_mags + 0.5,
    fill = "#f8fafc",
    color = "#d1d5db",
    linewidth = 0.35
  ) +
  ggplot2::geom_tile(
    data = plot_data,
    ggplot2::aes(x = window_index, y = mag_row, fill = status),
    width = 0.92,
    height = 0.98,
    color = NA
  ) +
  ggplot2::geom_vline(
    xintercept = seq(1.5, n_windows - 0.5, by = 1),
    color = "white",
    linewidth = 0.65
  ) +
  ggplot2::geom_text(
    data = window_summary,
    ggplot2::aes(x = window_index, y = y, label = label),
    size = 3.1,
    lineheight = 0.9,
    color = "#374151",
    fontface = "bold"
  ) +
  ggplot2::annotate(
    "segment",
    x = 0.45,
    xend = n_windows + 0.55,
    y = n_mags + 10,
    yend = n_mags + 10,
    color = "#d1d5db",
    linewidth = 0.35
  ) +
  ggplot2::annotate(
    "text",
    x = 0.45,
    y = n_mags + 48,
    label = "Main Figure 3 | Class 1 MAG-window testing matrix",
    hjust = 0,
    size = 5.2,
    fontface = "bold",
    color = "#111827"
  ) +
  ggplot2::annotate(
    "text",
    x = 0.45,
    y = n_mags + 35,
    label = "Rows are MAGs ordered by modelability; columns are biologically predefined Class 1 windows.",
    hjust = 0,
    size = 3.35,
    color = "#4b5563"
  ) +
  ggplot2::annotate(
    "text",
    x = 0.45,
    y = -28,
    label = audit_text,
    hjust = 0,
    size = 3.15,
    color = "#374151"
  ) +
  ggplot2::geom_tile(
    data = legend_df,
    ggplot2::aes(x = x * 0.75 + 0.1, y = -52, fill = status),
    width = 0.18,
    height = 8.5,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_text(
    data = legend_df,
    ggplot2::aes(x = x * 0.75 + 0.26, y = -52, label = label),
    hjust = 0,
    size = 2.8,
    lineheight = 0.86,
    color = "#374151",
    inherit.aes = FALSE
  ) +
  ggplot2::scale_fill_manual(values = status_palette, drop = FALSE, name = "Model status") +
  ggplot2::scale_x_continuous(
    breaks = seq_along(window_order),
    labels = unname(window_label_map[window_order]),
    limits = c(0.35, n_windows + 0.65),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = c(1, 100, 200, 300, 400, n_mags),
    labels = c("1", "100", "200", "300", "400", fmt_int(n_mags)),
    limits = c(-72, n_mags + 58),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    x = NULL,
    y = "MAG index after modelability ordering",
    caption = paste0(
      "Primary-analysis cells are lmm_ok. Flagged cells are retained for auditability but excluded from primary significance claims. ",
      "\nThe matrix shows the complete predefined Class 1 testing universe: ",
      fmt_int(n_mags), " MAGs x ", fmt_int(n_windows), " windows = ", fmt_int(n_tests), " tests."
    )
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 8.6, color = "#111827", lineheight = 0.86, margin = ggplot2::margin(t = 8)),
    axis.text.y = ggplot2::element_text(size = 8.8, color = "#4b5563"),
    axis.title.y = ggplot2::element_text(size = 10.2, color = "#374151", margin = ggplot2::margin(r = 10)),
    legend.position = "none",
    plot.caption = ggplot2::element_text(size = 8.5, color = "#4b5563", hjust = 0, margin = ggplot2::margin(t = 12)),
    plot.margin = ggplot2::margin(18, 20, 16, 18)
  )

png_path <- file.path(output_dir, "figure_3_class1_analysis_testing_matrix.png")
pdf_path <- file.path(output_dir, "figure_3_class1_analysis_testing_matrix.pdf")
svg_path <- file.path(output_dir, "figure_3_class1_analysis_testing_matrix.svg")

ggplot2::ggsave(png_path, figure, width = 10.8, height = 8.0, dpi = 360, bg = "white")
ggplot2::ggsave(pdf_path, figure, width = 10.8, height = 8.0, bg = "white")
ggplot2::ggsave(svg_path, figure, width = 10.8, height = 8.0, bg = "white", device = svglite::svglite)

message("Main Figure 3 Class 1 testing matrix created:")
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Inputs:")
message(" - ", results_path)
message(" - ", windows_path)
message("Summary:")
message(" - MAGs: ", fmt_int(n_mags))
message(" - Windows: ", fmt_int(n_windows))
message(" - Tests attempted: ", fmt_int(n_tests))
message(" - lmm_ok: ", fmt_int(n_ok))
message(" - lmm_singular: ", fmt_int(n_singular))
message(" - too_sparse: ", fmt_int(n_sparse))
message(" - low_variance: ", fmt_int(n_low_variance))
message(" - model_failed: ", fmt_int(n_failed))
