beta_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
                  "#F0E442", "#000000", "#8C6D31", "#6A3D9A", "#1B9E77", "#7570B3")
beta_shapes <- c(16, 17, 15, 18, 3, 7, 8, 0, 1, 2, 4, 5)

theme_beta <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.35),
                   axis.ticks = ggplot2::element_line(linewidth = 0.35),
                   legend.key = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title.position = "plot")
}

save_beta_plot <- function(plot, stem, width = 180, height = 140) {
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "mm", device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, units = "mm", dpi = 400, bg = "white")
}
