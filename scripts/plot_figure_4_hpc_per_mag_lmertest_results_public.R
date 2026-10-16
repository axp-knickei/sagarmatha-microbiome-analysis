# Public renderer: values come from external inputs; see docs/public-rendering.md.
# ============================================
# plot_figure_4_hpc_per_mag_lmertest_results.R
# ============================================
# Figure 4: global per-MAG LMM result landscape for the full T1-T10 HPC run.
#
# This figure visualizes the statistical landscape across all primary modeled
# MAGs. It distinguishes nominal raw-p signals from BH-FDR-supported findings
# and avoids treating volcano-plot associations as causal evidence.

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
output_stem <- "figure_4_hpc_per_mag_lmertest_results"

# External inputs and isolated output destination (run from repository root).
results_path <- Sys.getenv("PUBLIC_LMM_RESULTS", unset = results_path)
output_dir <- Sys.getenv("PUBLIC_RENDER_OUTPUT_DIR", unset = output_dir)

if (!file.exists(results_path)) {
  stop("Missing per-MAG LMM result table: ", results_path, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c(
  "MAG", "status", "primary_analysis", "time_estimate", "time_std_error",
  "time_p_value", "time_p_adj_bh_primary", "primary_significant",
  "zero_fraction", "n_subjects", "n_timepoints"
)
missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns) > 0) {
  stop("Missing required result columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

primary <- results |>
  dplyr::filter(primary_analysis %in% TRUE) |>
  dplyr::mutate(
    time_p_value = pmax(time_p_value, .Machine$double.xmin),
    time_p_adj_bh_primary = pmax(time_p_adj_bh_primary, .Machine$double.xmin),
    neg_log10_p = -log10(time_p_value),
    neg_log10_q = -log10(time_p_adj_bh_primary),
    effect_strength = abs(time_estimate),
    effect_direction = dplyr::case_when(
      time_estimate > 0 ~ "Increasing",
      time_estimate < 0 ~ "Decreasing",
      TRUE ~ "No directional estimate"
    ),
    result_class = dplyr::case_when(
      primary_significant %in% TRUE ~ "BH-FDR q < 0.05",
      time_p_value < 0.05 ~ "Nominal raw p < 0.05",
      TRUE ~ "Not nominal"
    ),
    result_class = factor(
      result_class,
      levels = c("BH-FDR q < 0.05", "Nominal raw p < 0.05", "Not nominal")
    )
  ) |>
  dplyr::arrange(time_p_value) |>
  dplyr::mutate(
    p_rank = dplyr::row_number(),
    effect_rank = dplyr::row_number(dplyr::desc(effect_strength))
  )

if (nrow(primary) == 0) {
  stop("No primary-analysis MAGs found in result table.", call. = FALSE)
}

top_n <- min(12L, nrow(primary))
top_mag_table <- primary |>
  dplyr::arrange(time_p_value) |>
  dplyr::slice_head(n = top_n) |>
  dplyr::mutate(
    volcano_label = as.character(dplyr::row_number()),
    short_mag = sub("^(QML|EVR)", "\\1", MAG),
    short_mag = ifelse(nchar(short_mag) > 33, paste0(substr(short_mag, 1, 30), "..."), short_mag)
  )

ranked_effect <- primary |>
  dplyr::arrange(time_estimate) |>
  dplyr::mutate(effect_order = dplyr::row_number())

status_counts <- as.data.frame(table(results$status), stringsAsFactors = FALSE)
names(status_counts) <- c("status", "n")

total_mags <- nrow(results)
primary_mags <- nrow(primary)
nominal_mags <- sum(primary$time_p_value < 0.05, na.rm = TRUE)
fdr_mags <- sum(primary$primary_significant %in% TRUE, na.rm = TRUE)
increasing_nominal <- sum(primary$time_p_value < 0.05 & primary$time_estimate > 0, na.rm = TRUE)
decreasing_nominal <- sum(primary$time_p_value < 0.05 & primary$time_estimate < 0, na.rm = TRUE)
best_p <- min(primary$time_p_value, na.rm = TRUE)
best_q <- min(primary$time_p_adj_bh_primary, na.rm = TRUE)
effect_min <- min(primary$time_estimate, na.rm = TRUE)
effect_max <- max(primary$time_estimate, na.rm = TRUE)

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_p <- function(x) {
  if (is.na(x) || !is.finite(x)) return("NA")
  if (x < 0.001) formatC(x, format = "e", digits = 2) else sprintf("%.3f", x)
}

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
  green = "#009E73",
  mint = "#EAF7F2",
  red = "#CC3311",
  purple = "#785EF0",
  grey = "#9CA3AF"
)

class_palette <- c(
  "BH-FDR q < 0.05" = pal$green,
  "Nominal raw p < 0.05" = pal$orange,
  "Not nominal" = pal$grey
)

direction_palette <- c(
  "Increasing" = pal$blue,
  "Decreasing" = pal$red,
  "No directional estimate" = pal$grey
)

base_theme <- ggplot2::theme_minimal(base_size = 9.6, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11.1, color = pal$ink),
    plot.subtitle = ggplot2::element_text(size = 8.2, color = pal$muted, margin = ggplot2::margin(b = 6)),
    plot.caption = ggplot2::element_text(size = 7.0, color = pal$muted, hjust = 0),
    axis.title = ggplot2::element_text(size = 8.4, color = pal$text),
    axis.text = ggplot2::element_text(size = 7.6, color = pal$text),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.25),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 8),
    legend.text = ggplot2::element_text(size = 7.5),
    plot.margin = ggplot2::margin(6, 7, 6, 7)
  )

panel_a <- ggplot2::ggplot(
  primary,
  ggplot2::aes(x = time_estimate, y = neg_log10_p)
) +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = pal$orange, linewidth = 0.45, linetype = "longdash") +
  ggplot2::geom_vline(xintercept = 0, color = pal$hair, linewidth = 0.45) +
  ggplot2::geom_point(
    ggplot2::aes(color = result_class, size = zero_fraction),
    alpha = 0.82
  ) +
  ggplot2::geom_label(
    data = top_mag_table,
    ggplot2::aes(label = volcano_label),
    size = 2.45,
    label.padding = grid::unit(0.08, "lines"),
    label.r = grid::unit(0.06, "lines"),
    fill = "white",
    color = pal$ink,
    linewidth = 0.25,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = class_palette, drop = FALSE, name = "Result class") +
  ggplot2::scale_size_continuous(
    range = c(1.2, 3.4),
    labels = function(x) paste0(round(100 * x), "%"),
    name = "Zero fraction"
  ) +
  ggplot2::labs(
    title = "A  Volcano landscape of per-MAG longitudinal effects",
    subtitle = paste0(
      nominal_mags, " nominal raw-p candidates; ",
      fdr_mags, " MAGs pass primary BH-FDR q < 0.05"
    ),
    x = "Time effect estimate, log10(Y + pseudocount) per day",
    y = expression(-log[10]("raw p-value"))
  ) +
  base_theme +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

panel_b <- ggplot2::ggplot(
  ranked_effect,
  ggplot2::aes(x = effect_order, y = time_estimate)
) +
  ggplot2::geom_hline(yintercept = 0, color = pal$hair, linewidth = 0.5) +
  ggplot2::geom_point(
    ggplot2::aes(color = result_class),
    size = 1.75,
    alpha = 0.82
  ) +
  ggplot2::scale_color_manual(values = class_palette, drop = FALSE, name = "Result class") +
  ggplot2::labs(
    title = "B  Ranked longitudinal effect sizes",
    subtitle = paste0(
      increasing_nominal, " nominal increasing MAGs; ",
      decreasing_nominal, " nominal decreasing MAGs"
    ),
    x = "MAG rank ordered by time-effect estimate",
    y = "Time effect estimate"
  ) +
  base_theme +
  ggplot2::theme(legend.position = "none")

panel_c <- ggplot2::ggplot(primary, ggplot2::aes(x = time_p_value)) +
  ggplot2::geom_histogram(binwidth = 0.05, boundary = 0, fill = pal$sky, color = "white", linewidth = 0.25) +
  ggplot2::geom_vline(xintercept = 0.05, color = pal$orange, linewidth = 0.55, linetype = "longdash") +
  ggplot2::annotate(
    "text",
    x = 0.05,
    y = Inf,
    label = "raw p = 0.05",
    hjust = -0.08,
    vjust = 1.45,
    size = 2.7,
    color = pal$orange
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ggplot2::labs(
    title = "C  Raw p-value distribution",
    subtitle = paste0(nominal_mags, " nominal candidates; ", fdr_mags, " primary BH-FDR discoveries"),
    x = "Raw p-value for days_from_baseline_subject",
    y = "Primary MAG count"
  ) +
  base_theme

top_table_text <- paste(
  sprintf(
    "%2s  %-33s  beta=% .3f  p=%s  q=%s",
    top_mag_table$volcano_label,
    top_mag_table$short_mag,
    top_mag_table$time_estimate,
    vapply(top_mag_table$time_p_value, fmt_p, character(1)),
    vapply(top_mag_table$time_p_adj_bh_primary, fmt_p, character(1))
  ),
  collapse = "\n"
)

summary_text <- paste0(
  "Primary LMM result landscape: ", fmt_int(primary_mags), " primary modeled MAGs from ",
  fmt_int(total_mags), " total result rows; best raw p = ", fmt_p(best_p),
  "; best primary BH q = ", fmt_p(best_q), ". Effect range = ",
  sprintf("%.3f", effect_min), " to ", sprintf("%.3f", effect_max),
  " log10 units per day. Associations are longitudinal model signals, not causal claims."
)

draw_top_table <- function() {
  grid::grid.roundrect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.50, "npc"),
    width = grid::unit(0.98, "npc"),
    height = grid::unit(0.88, "npc"),
    r = grid::unit(0.035, "snpc"),
    gp = grid::gpar(fill = "white", col = pal$hair, lwd = 0.8)
  )
  grid::grid.text(
    "D  Top raw-p ranked MAGs mapped to volcano labels",
    x = grid::unit(0.04, "npc"),
    y = grid::unit(0.90, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 10.1, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    paste0(sum(top_mag_table$time_p_adj_bh_primary < 0.05, na.rm = TRUE), " displayed MAGs have primary BH-FDR q < 0.05; nominal-only rows are exploratory."),
    x = grid::unit(0.04, "npc"),
    y = grid::unit(0.81, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 7.7, col = pal$muted)
  )
  grid::grid.text(
    top_table_text,
    x = grid::unit(0.04, "npc"),
    y = grid::unit(0.72, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 6.4, col = pal$text, fontfamily = "mono", lineheight = 0.95)
  )
}

draw_figure <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  layout <- grid::grid.layout(
    nrow = 4,
    ncol = 2,
    heights = grid::unit(c(0.50, 1.25, 1.0, 0.35), c("in", "null", "null", "in")),
    widths = grid::unit(c(1.15, 1.0), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))

  grid::grid.text(
    "Figure 4 | Global per-MAG LMM result landscape",
    x = grid::unit(0, "npc"),
    y = grid::unit(0.88, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 15.5, fontface = "bold", col = pal$ink),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )
  grid::grid.text(
    paste0("Full T1-T10 screening: ", nominal_mags, " nominal candidates and ", fdr_mags, " primary BH-FDR discoveries."),
    x = grid::unit(0, "npc"),
    y = grid::unit(0.20, "npc"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 9.0, col = pal$text),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )

  print(panel_a, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panel_b, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(panel_c, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
  grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
  draw_top_table()
  grid::popViewport()

  grid::grid.text(
    summary_text,
    x = grid::unit(0, "npc"),
    y = grid::unit(0.57, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 7.5, col = pal$muted),
    vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1:2)
  )

  grid::popViewport()
}

png_path <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_path <- file.path(output_dir, paste0(output_stem, ".svg"))

png(filename = png_path, width = 12.8, height = 9.2, units = "in", res = 420, bg = "white")
draw_figure()
dev.off()

pdf(file = pdf_path, width = 12.8, height = 9.2, bg = "white", useDingbats = FALSE)
draw_figure()
dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = svg_path, width = 12.8, height = 9.2, bg = "white")
  draw_figure()
  dev.off()
}

summary_path <- file.path(output_dir, paste0(output_stem, "_result_landscape_summary.csv"))
top_table_path <- file.path(output_dir, paste0(output_stem, "_top_raw_p_mags.csv"))

utils::write.csv(
  data.frame(
    metric = c(
      "total_result_rows",
      "primary_modeled_mags",
      "nominal_raw_p_lt_0_05",
      "bh_fdr_significant",
      "nominal_increasing",
      "nominal_decreasing",
      "best_raw_p",
      "best_primary_bh_q",
      "minimum_time_estimate",
      "maximum_time_estimate"
    ),
    value = c(
      total_mags,
      primary_mags,
      nominal_mags,
      fdr_mags,
      increasing_nominal,
      decreasing_nominal,
      best_p,
      best_q,
      effect_min,
      effect_max
    )
  ),
  summary_path,
  row.names = FALSE
)

utils::write.csv(
  top_mag_table |>
    dplyr::select(
      volcano_label, MAG, time_estimate, time_std_error, time_p_value,
      time_p_adj_bh_primary, zero_fraction, n_subjects, n_timepoints
    ),
  top_table_path,
  row.names = FALSE
)

message("Figure 4 result-landscape outputs created:")
message(" - ", png_path)
message(" - ", pdf_path)
if (file.exists(svg_path)) message(" - ", svg_path)
message(" - ", summary_path)
message(" - ", top_table_path)
