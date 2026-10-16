# ============================================
# plot_figure_8_class1_analysis.R
# ============================================
# Supplementary/Main Figure 8: raw p-value and global q-value distributions.
# The purpose is to show why Class 1 nominal structure did not survive global
# BH-FDR correction.

required_packages <- c("readr", "dplyr", "ggplot2", "grid")
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
output_dir <- file.path(repo_root, "results", "class1_window_lmertest", "figures")

if (!file.exists(results_path)) {
  stop("Missing Class 1 results table: ", results_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

class1_results <- readr::read_csv(results_path, show_col_types = FALSE)

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_q <- function(x) ifelse(is.na(x), "NA", sprintf("%.4f", x))

primary_results <- class1_results |>
  dplyr::filter(primary_analysis %in% TRUE) |>
  dplyr::filter(!is.na(raw_p_value), !is.na(q_global_class1_bh))

if (nrow(primary_results) == 0) {
  stop("No primary Class 1 lmm_ok tests found for p/q distribution plotting.", call. = FALSE)
}

n_primary <- nrow(primary_results)
n_nominal <- sum(primary_results$raw_p_value < 0.05, na.rm = TRUE)
n_fdr <- sum(primary_results$q_global_class1_bh < 0.05, na.rm = TRUE)
best_q <- min(primary_results$q_global_class1_bh, na.rm = TRUE)
best_raw_p <- min(primary_results$raw_p_value, na.rm = TRUE)

ranked_q <- primary_results |>
  dplyr::arrange(q_global_class1_bh, raw_p_value) |>
  dplyr::mutate(rank = dplyr::row_number())

ranked_q_zoom <- ranked_q |>
  dplyr::filter(q_global_class1_bh <= 0.12)
if (nrow(ranked_q_zoom) == 0) {
  ranked_q_zoom <- ranked_q |>
    dplyr::slice_head(n = min(50L, nrow(ranked_q)))
}

palette <- list(
  ink = "#111827",
  muted = "#4b5563",
  grid = "#e5e7eb",
  blue = "#0072B2",
  vermillion = "#D55E00",
  amber = "#E69F00",
  red = "#dc2626"
)

panel_a <- ggplot2::ggplot(primary_results, ggplot2::aes(x = raw_p_value)) +
  ggplot2::geom_histogram(
    bins = 40,
    boundary = 0,
    closed = "left",
    fill = "#bfdbfe",
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::annotate(
    "rect",
    xmin = 0,
    xmax = 0.05,
    ymin = -Inf,
    ymax = Inf,
    fill = "#fee2e2",
    alpha = 0.45
  ) +
  ggplot2::geom_histogram(
    bins = 40,
    boundary = 0,
    closed = "left",
    fill = "#0072B2",
    color = "white",
    linewidth = 0.25,
    alpha = 0.78
  ) +
  ggplot2::geom_vline(xintercept = 0.05, color = palette$red, linetype = "22", linewidth = 0.65) +
  ggplot2::annotate(
    "text",
    x = 0.055,
    y = Inf,
    label = "raw p = 0.05",
    hjust = 0,
    vjust = 1.7,
    size = 3.0,
    color = "#991b1b"
  ) +
  ggplot2::annotate(
    "label",
    x = 0.98,
    y = Inf,
    label = paste0(fmt_int(n_nominal), " nominal raw-p candidates\namong ", fmt_int(n_primary), " lmm_ok tests"),
    hjust = 1,
    vjust = 1.25,
    size = 3.15,
    lineheight = 0.92,
    label.r = grid::unit(0.08, "inches"),
    label.size = 0.35,
    fill = "white",
    color = palette$ink
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), expand = c(0, 0)) +
  ggplot2::labs(
    title = "A  Raw p-value distribution",
    x = "Raw p-value from window_order term",
    y = "Number of MAG-window tests"
  ) +
  ggplot2::theme_minimal(base_size = 10.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12.2, color = palette$ink),
    axis.title = ggplot2::element_text(size = 10.2, color = "#374151"),
    axis.text = ggplot2::element_text(size = 8.8, color = palette$muted),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = palette$grid, linewidth = 0.35)
  )

panel_b <- ggplot2::ggplot(ranked_q_zoom, ggplot2::aes(x = rank, y = q_global_class1_bh)) +
  ggplot2::annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 0,
    ymax = 0.05,
    fill = "#fee2e2",
    alpha = 0.45
  ) +
  ggplot2::geom_point(color = palette$blue, alpha = 0.42, size = 1.0) +
  ggplot2::geom_line(color = palette$blue, alpha = 0.42, linewidth = 0.35) +
  ggplot2::geom_hline(yintercept = 0.05, color = palette$red, linetype = "22", linewidth = 0.65) +
  ggplot2::geom_point(
    data = dplyr::slice_min(ranked_q, q_global_class1_bh, n = 1, with_ties = FALSE),
    color = palette$red,
    fill = "white",
    shape = 21,
    size = 3.2,
    stroke = 1.1
  ) +
  ggplot2::annotate(
    "label",
    x = max(ranked_q_zoom$rank) * 0.96,
    y = 0.092,
    label = paste0("best q = ", fmt_q(best_q), "\nq < 0.05 tests = ", fmt_int(n_fdr)),
    hjust = 1,
    vjust = 0,
    size = 3.15,
    lineheight = 0.92,
    label.r = grid::unit(0.08, "inches"),
    label.size = 0.35,
    fill = "white",
    color = palette$ink
  ) +
  ggplot2::annotate(
    "text",
    x = max(ranked_q_zoom$rank) * 0.04,
    y = 0.052,
    label = "global BH-FDR threshold q = 0.05",
    hjust = 0,
    vjust = 0,
    size = 3.0,
    color = "#991b1b"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, max(ranked_q_zoom$rank), by = 50),
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 0.12),
    breaks = seq(0, 0.12, by = 0.02),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = "B  Lowest ranked global Class 1 BH-FDR q-values",
    x = paste0("Top ", nrow(ranked_q_zoom), " MAG-window tests ranked by q_global_class1_bh"),
    y = "q_global_class1_bh"
  ) +
  ggplot2::theme_minimal(base_size = 10.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12.2, color = palette$ink),
    axis.title = ggplot2::element_text(size = 10.2, color = "#374151"),
    axis.text = ggplot2::element_text(size = 8.8, color = palette$muted),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = palette$grid, linewidth = 0.35)
  )

draw_combined <- function() {
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 4,
    ncol = 2,
    heights = grid::unit(c(0.62, 4.6, 0.16, 0.42), c("in", "null", "in", "in")),
    widths = grid::unit(c(1, 1), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))

  grid::grid.text(
    "Main Figure 8 | P-value and global q-value distribution",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.72, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 17, fontface = "bold", col = palette$ink),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )
  grid::grid.text(
    "Nominal Class 1 structure was present, but evidence did not survive global BH-FDR correction.",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.2, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 10.5, col = palette$muted),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )

  print(panel_a, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panel_b, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))

  grid::grid.text(
    paste0(
      "Primary lmm_ok tests: ", fmt_int(n_primary),
      " | nominal raw p < 0.05: ", fmt_int(n_nominal),
      " | best raw p: ", formatC(best_raw_p, format = "e", digits = 2),
      " | best global q: ", fmt_q(best_q),
      " | FDR-significant tests: ", fmt_int(n_fdr)
    ),
    x = grid::unit(0, "npc"),
    y = grid::unit(0.55, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 9.5, col = "#374151"),
    vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1:2)
  )
  grid::grid.text(
    "BH-FDR controls the expected false discovery proportion across the predefined Class 1 test family; raw p-value candidates remain hypothesis-generating unless they pass global correction.",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.55, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 8.5, col = palette$muted),
    vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1:2)
  )
  grid::popViewport()
}

png_path <- file.path(output_dir, "figure_8_class1_analysis_pvalue_qvalue_distribution.png")
pdf_path <- file.path(output_dir, "figure_8_class1_analysis_pvalue_qvalue_distribution.pdf")
svg_path <- file.path(output_dir, "figure_8_class1_analysis_pvalue_qvalue_distribution.svg")

png(filename = png_path, width = 12.2, height = 6.5, units = "in", res = 360, bg = "white")
draw_combined()
dev.off()

pdf(file = pdf_path, width = 12.2, height = 6.5, bg = "white")
draw_combined()
dev.off()

svglite::svglite(file = svg_path, width = 12.2, height = 6.5, bg = "white")
draw_combined()
dev.off()

message("Main Figure 8 Class 1 p/q distribution diagnostic created:")
message(" - ", png_path)
message(" - ", pdf_path)
message(" - ", svg_path)
message("Input:")
message(" - ", results_path)
message("Summary:")
message(" - Primary lmm_ok tests: ", fmt_int(n_primary))
message(" - Nominal raw p < 0.05 candidates: ", fmt_int(n_nominal))
message(" - FDR-significant tests: ", fmt_int(n_fdr))
message(" - Best raw p-value: ", formatC(best_raw_p, format = "e", digits = 2))
message(" - Best q_global_class1_bh: ", fmt_q(best_q))
