suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

parse_synthesis_args <- function(default_config = "config/manuscript_functional_synthesis.yml") {
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

read_synthesis_config <- function(path) {
  if (!file.exists(path)) stop("Configuration not found: ", path, call. = FALSE)
  cfg <- yaml::read_yaml(path)
  required <- c("inputs", "analysis", "outputs")
  absent <- setdiff(required, names(cfg))
  if (length(absent)) {
    stop("Configuration is missing sections: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  cfg$.config_path <- normalizePath(path, mustWork = TRUE)
  cfg
}

resolve_synthesis_inputs <- function(cfg) {
  paths <- cfg$inputs
  overrides <- cfg$input_environment_overrides %||% list()
  for (key in names(paths)) {
    env_name <- overrides[[key]] %||% ""
    env_value <- if (nzchar(env_name)) Sys.getenv(env_name, unset = "") else ""
    if (nzchar(env_value)) paths[[key]] <- env_value
  }
  paths
}

synthesis_paths <- function(cfg, create = TRUE) {
  root <- Sys.getenv("MANSYN_OUTPUT_DIR", unset = cfg$outputs$output_dir)
  paths <- list(
    root = root,
    qc = file.path(root, "00_quality_control"),
    main = file.path(root, "main"),
    reproducibility = file.path(root, "reproducibility")
  )
  if (create) invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
  paths
}

guard_synthesis <- function(marker, args, cfg) {
  overwrite <- isTRUE(args$force) || isTRUE(cfg$outputs$overwrite)
  if (file.exists(marker) && !overwrite) {
    stop("Completion marker already exists: ", marker, ". Use --force only when intentional.", call. = FALSE)
  }
  if (file.exists(marker) && overwrite) unlink(marker)
  invisible(overwrite)
}

require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Required synthesis inputs are missing:\n", paste("-", missing, collapse = "\n"), call. = FALSE)
  }
  invisible(paths)
}

read_result_table <- function(path) {
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) return(data.table::fread(path))
  gzip <- unname(Sys.which("gzip"))
  if (!nzchar(gzip)) stop("gzip is required to read: ", path, call. = FALSE)
  data.table::fread(cmd = paste(shQuote(gzip), "-cd", shQuote(normalizePath(path, mustWork = TRUE))))
}

write_result_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data.table::as.data.table(x), path, na = "")
  invisible(path)
}

theme_manuscript <- function(base_size = 7.6) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "#202124"),
      axis.text = ggplot2::element_text(colour = "#35373A"),
      axis.title = ggplot2::element_text(colour = "#202124"),
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "#202124"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, colour = "#202124"),
      legend.position = "top",
      legend.justification = "left",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.height = grid::unit(3.2, "mm"),
      legend.key.width = grid::unit(5.0, "mm"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 0.8, margin = ggplot2::margin(b = 3, l = 10)),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.4, colour = "#55585C", margin = ggplot2::margin(b = 4, l = 10)),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 3),
      plot.tag.position = c(0, 1),
      plot.margin = ggplot2::margin(5, 6, 5, 6)
    )
}

save_composite <- function(grob, stem, cfg, paths) {
  width <- as.numeric(cfg$outputs$figure_width_mm) / 25.4
  height <- as.numeric(cfg$outputs$figure_height_mm) / 25.4
  png_path <- file.path(paths$main, paste0(stem, ".png"))
  pdf_path <- file.path(paths$main, paste0(stem, ".pdf"))
  svg_path <- file.path(paths$main, paste0(stem, ".svg"))
  ggplot2::ggsave(png_path, grob, width = width, height = height, units = "in", dpi = as.integer(cfg$outputs$png_dpi), bg = "white")
  ggplot2::ggsave(pdf_path, grob, width = width, height = height, units = "in", device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(svg_path, grob, width = width, height = height, units = "in", device = svglite::svglite, bg = "white")
  c(pdf_path, png_path, svg_path)
}

wrap_table_text <- function(x, width) {
  vapply(as.character(x), function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

save_table_preview <- function(x, stem, cfg, paths, widths, height) {
  preview <- data.frame(lapply(seq_along(x), function(i) wrap_table_text(x[[i]], widths[[i]])), check.names = FALSE)
  names(preview) <- names(x)
  theme <- gridExtra::ttheme_minimal(
    base_size = 7.4,
    base_family = "sans",
    core = list(
      fg_params = list(hjust = 0, x = 0.02, col = "#202124"),
      bg_params = list(fill = rep(c("#FFFFFF", "#F4F6F8"), length.out = nrow(preview)), col = "#D9DDE2", lwd = 0.35)
    ),
    colhead = list(
      fg_params = list(col = "white", fontface = "bold", hjust = 0, x = 0.02),
      bg_params = list(fill = "#443F90", col = "#443F90")
    )
  )
  grob <- gridExtra::tableGrob(preview, rows = NULL, theme = theme)
  grob$widths <- grid::unit(widths / sum(widths), "npc")
  width <- as.numeric(cfg$outputs$figure_width_mm) / 25.4
  png_path <- file.path(paths$main, paste0(stem, ".png"))
  pdf_path <- file.path(paths$main, paste0(stem, ".pdf"))
  ggplot2::ggsave(png_path, grob, width = width, height = height, dpi = as.integer(cfg$outputs$table_preview_dpi), bg = "white")
  ggplot2::ggsave(pdf_path, grob, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  c(pdf_path, png_path)
}
