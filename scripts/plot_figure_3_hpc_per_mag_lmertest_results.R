# ============================================
# plot_figure_3_hpc_per_mag_lmertest_results.R
# ============================================
# Figure 3: model workflow schematic from MAG abundance table to per-MAG LMM.
#
# This figure is a conceptual workflow diagram, but key labels are populated
# from the actual HPC per-MAG LMM outputs and reproducibility record.

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
execution_record_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "reproducibility",
  "hpc_execution_record.csv"
)
session_info_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "session",
  "session_info.txt"
)
summary_log_path <- file.path(
  repo_root,
  "results",
  "hpc_per_mag_lmertest",
  "logs",
  "per_mag_lmertest_hpc_summary.txt"
)
output_dir <- file.path(repo_root, "results", "hpc_per_mag_lmertest", "figures")
output_stem <- "figure_3_hpc_per_mag_lmertest_results_20260708_024546"

required_files <- c(results_path, execution_record_path, session_info_path, summary_log_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required workflow output file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- utils::read.csv(results_path, check.names = FALSE, stringsAsFactors = FALSE)
execution_record <- utils::read.csv(execution_record_path, check.names = FALSE, stringsAsFactors = FALSE)

needed_cols <- c("MAG", "status", "primary_analysis", "time_p_value", "time_p_adj_bh_primary", "primary_significant")
missing_cols <- setdiff(needed_cols, names(results))
if (length(missing_cols) > 0) {
  stop("Missing required columns in per-MAG result table: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

record <- execution_record[1, , drop = FALSE]

total_mags <- nrow(results)
primary_mags <- sum(results$primary_analysis %in% TRUE, na.rm = TRUE)
nominal_mags <- sum(results$primary_analysis %in% TRUE & results$time_p_value < 0.05, na.rm = TRUE)
fdr_mags <- sum(results$primary_significant %in% TRUE, na.rm = TRUE)
singular_mags <- sum(results$status == "lmm_singular", na.rm = TRUE)
too_sparse_mags <- sum(results$status == "too_sparse", na.rm = TRUE)
model_failed_mags <- sum(results$status == "model_failed", na.rm = TRUE)
best_q <- suppressWarnings(min(results$time_p_adj_bh_primary[results$primary_analysis %in% TRUE], na.rm = TRUE))
if (!is.finite(best_q)) best_q <- NA_real_

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_q <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
coalesce_text <- function(x, fallback = "not recorded") {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) fallback else x
}

conda_env <- coalesce_text(record$conda_environment)
conda_module <- coalesce_text(record$conda_module)
r_version <- sub(" \\(.*$", "", coalesce_text(record$R_version))
job_id <- coalesce_text(record$job_id)
runtime <- suppressWarnings(as.numeric(record$runtime_seconds))
runtime_label <- if (is.finite(runtime)) paste0(sprintf("%.1f", runtime), " s") else "not recorded"
git_commit <- coalesce_text(record$git_commit)
git_short <- if (nchar(git_commit) >= 8) substr(git_commit, 1, 8) else git_commit

pal <- list(
  ink = "#111827",
  text = "#374151",
  muted = "#6B7280",
  hair = "#D1D5DB",
  faint = "#F8FAFC",
  blue = "#0072B2",
  sky = "#DCEEF8",
  green = "#009E73",
  mint = "#EAF7F2",
  orange = "#D55E00",
  peach = "#FFF2E8",
  purple = "#785EF0",
  lavender = "#F1EEFF",
  red = "#CC3311",
  rose = "#FFF0EE",
  gold = "#E69F00",
  cream = "#FFF8E6"
)

workflow_nodes <- data.frame(
  id = 1:8,
  x = rep(3.2, 8),
  y = seq(7.75, 0.95, length.out = 8),
  fill = c(pal$faint, pal$sky, pal$mint, pal$lavender, pal$peach, pal$cream, pal$mint, pal$rose),
  accent = c(pal$muted, pal$blue, pal$green, pal$purple, pal$orange, pal$gold, pal$green, pal$red),
  kicker = c("INPUT", "EXTRACT", "TRANSFORM", "MODEL", "QC STATUS", "MULTIPLE TESTING", "RESULT TABLE", "INTERPRETATION"),
  title = c(
    "Processed MAG abundance table",
    "One MAG-specific abundance vector at a time",
    "Primary response scale",
    "Per-MAG longitudinal LMM",
    "Model status classification",
    "Benjamini-Hochberg FDR",
    "Ranked scientific output",
    "Diagnostics and biological review"
  ),
  body = c(
    paste0("data/processed/cleaned_data.rds\n", fmt_int(total_mags), " MAGs with repeated subject sampling"),
    "Loop over MAG IDs; retain Subject, sampleID, Timepoint,\nand days_from_baseline_subject",
    "log10(Y + pseudocount)\npseudocount = half of smallest non-zero abundance",
    "log10(Y + c) ~ days_from_baseline_subject + (1 | Subject)\nlmerTest Satterthwaite p-values",
    paste0(primary_mags, " lmm_ok primary MAGs; ", singular_mags, " singular; ", too_sparse_mags, " too sparse\nflagged rows retained but excluded from primary claims"),
    paste0("One correction family across modeled MAGs\n", nominal_mags, " raw p < 0.05; ", fdr_mags, " q < 0.05; best q = ", fmt_q(best_q)),
    "results/hpc_per_mag_lmertest/tables/\nper_mag_lmertest_hpc_results.csv",
    "Candidate summaries, response-scale caveats,\nmanual plot review, biological plausibility"
  ),
  stringsAsFactors = FALSE
)

hpc_nodes <- data.frame(
  x = 8.25,
  y = c(7.65, 6.35, 5.05, 3.75, 2.45),
  fill = c(pal$faint, pal$sky, pal$mint, pal$cream, pal$lavender),
  accent = c(pal$muted, pal$blue, pal$green, pal$gold, pal$purple),
  kicker = c("SCHEDULER", "ENVIRONMENT", "RUNTIME", "LOGGING", "REPRODUCIBILITY"),
  title = c(
    "OpenPBS wrapper",
    "Conda + R environment",
    "Batch execution",
    "Operational logs",
    "Execution record"
  ),
  body = c(
    "hpc/job_per_mag_lmertest.pbs\nnodes=1:ppn=4; mem=16gb; walltime=24h",
    paste0("module: ", conda_module, "\nenv: ", conda_env, "; ", r_version),
    paste0("job_id: ", job_id, "\nruntime recorded: ", runtime_label),
    "logs/<job>.out and logs/<job>.err\nscheduler logs are operational, not scientific tables",
    paste0("hpc_execution_record.csv + session_info.txt\ngit commit: ", git_short)
  ),
  stringsAsFactors = FALSE
)

draw_box <- function(x, y, width, height, fill, accent, kicker, title, body, title_size = 8.8, body_size = 7.2) {
  grid::grid.roundrect(
    x = grid::unit(x, "native"),
    y = grid::unit(y, "native"),
    width = grid::unit(width, "native"),
    height = grid::unit(height, "native"),
    r = grid::unit(0.05, "snpc"),
    gp = grid::gpar(fill = fill, col = pal$hair, lwd = 0.8)
  )
  grid::grid.rect(
    x = grid::unit(x - width / 2 + 0.055, "native"),
    y = grid::unit(y, "native"),
    width = grid::unit(0.11, "native"),
    height = grid::unit(height * 0.82, "native"),
    gp = grid::gpar(fill = accent, col = NA)
  )
  grid::grid.text(
    kicker,
    x = grid::unit(x - width / 2 + 0.22, "native"),
    y = grid::unit(y + height * 0.30, "native"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 6.2, fontface = "bold", col = accent)
  )
  grid::grid.text(
    title,
    x = grid::unit(x - width / 2 + 0.22, "native"),
    y = grid::unit(y + height * 0.08, "native"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = title_size, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    body,
    x = grid::unit(x - width / 2 + 0.22, "native"),
    y = grid::unit(y - height * 0.24, "native"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = body_size, col = pal$text, lineheight = 0.92)
  )
}

draw_arrow <- function(x0, y0, x1, y1, col = pal$muted, lwd = 1.0) {
  grid::grid.segments(
    x0 = grid::unit(x0, "native"),
    y0 = grid::unit(y0, "native"),
    x1 = grid::unit(x1, "native"),
    y1 = grid::unit(y1, "native"),
    arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed"),
    gp = grid::gpar(col = col, lwd = lwd, lineend = "round")
  )
}

draw_figure <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  grid::pushViewport(grid::viewport(
    xscale = c(0, 10.5),
    yscale = c(0, 9.25),
    width = grid::unit(1, "npc") - grid::unit(0.18, "in"),
    height = grid::unit(1, "npc") - grid::unit(0.18, "in")
  ))

  grid::grid.text(
    "Figure 3 | HPC per-MAG mixed-model workflow",
    x = grid::unit(0.0, "native"),
    y = grid::unit(9.05, "native"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 16, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    "One longitudinal LMM is fitted per MAG, then status filters, BH-FDR correction, and reproducibility records separate scientific results from scheduler logs.",
    x = grid::unit(0.0, "native"),
    y = grid::unit(8.68, "native"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 8.8, col = pal$text)
  )

  grid::grid.text(
    "Scientific analysis path",
    x = grid::unit(0.72, "native"),
    y = grid::unit(8.24, "native"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 9.0, fontface = "bold", col = pal$blue)
  )
  grid::grid.text(
    "HPC execution and reproducibility layer",
    x = grid::unit(6.42, "native"),
    y = grid::unit(8.24, "native"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 9.0, fontface = "bold", col = pal$purple)
  )

  for (i in seq_len(nrow(workflow_nodes))) {
    draw_box(
      workflow_nodes$x[[i]],
      workflow_nodes$y[[i]],
      width = 4.95,
      height = 0.72,
      fill = workflow_nodes$fill[[i]],
      accent = workflow_nodes$accent[[i]],
      kicker = workflow_nodes$kicker[[i]],
      title = workflow_nodes$title[[i]],
      body = workflow_nodes$body[[i]]
    )
    if (i < nrow(workflow_nodes)) {
      draw_arrow(workflow_nodes$x[[i]], workflow_nodes$y[[i]] - 0.43, workflow_nodes$x[[i + 1]], workflow_nodes$y[[i + 1]] + 0.43)
    }
  }

  for (i in seq_len(nrow(hpc_nodes))) {
    draw_box(
      hpc_nodes$x[[i]],
      hpc_nodes$y[[i]],
      width = 3.75,
      height = 0.78,
      fill = hpc_nodes$fill[[i]],
      accent = hpc_nodes$accent[[i]],
      kicker = hpc_nodes$kicker[[i]],
      title = hpc_nodes$title[[i]],
      body = hpc_nodes$body[[i]],
      title_size = 8.3,
      body_size = 6.9
    )
    if (i < nrow(hpc_nodes)) {
      draw_arrow(hpc_nodes$x[[i]], hpc_nodes$y[[i]] - 0.46, hpc_nodes$x[[i + 1]], hpc_nodes$y[[i + 1]] + 0.46, col = pal$purple)
    }
  }

  # Cross-links show how the scheduler layer launches the scientific workflow
  # and where the workflow writes reproducibility artifacts.
  draw_arrow(6.38, 7.65, 5.78, 7.75, col = pal$purple, lwd = 0.9)
  draw_arrow(6.38, 5.05, 5.78, 3.86, col = pal$purple, lwd = 0.9)
  draw_arrow(6.38, 2.45, 5.78, 1.92, col = pal$purple, lwd = 0.9)

  grid::grid.roundrect(
    x = grid::unit(8.25, "native"),
    y = grid::unit(1.15, "native"),
    width = grid::unit(3.75, "native"),
    height = grid::unit(1.08, "native"),
    r = grid::unit(0.04, "snpc"),
    gp = grid::gpar(fill = "#FFFFFF", col = pal$hair, lwd = 0.8)
  )
  grid::grid.text(
    "Key distinction",
    x = grid::unit(6.55, "native"),
    y = grid::unit(1.50, "native"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 8.2, fontface = "bold", col = pal$ink)
  )
  grid::grid.text(
    "Scientific output: per_mag_lmertest_hpc_results.csv\nOperational output: PBS stdout/stderr and scheduler logs\nReproducibility output: session_info + hpc_execution_record",
    x = grid::unit(6.55, "native"),
    y = grid::unit(1.05, "native"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 6.9, col = pal$text, lineheight = 0.95)
  )

  grid::grid.text(
    paste0(
      "Actual run summarized: ", fmt_int(total_mags), " MAG tests; ",
      fmt_int(primary_mags), " primary lmm_ok; ", fmt_int(nominal_mags),
      " nominal raw p < 0.05; ", fmt_int(fdr_mags),
      " BH-FDR significant. Input, result, session, and execution-record paths are read from repository outputs."
    ),
    x = grid::unit(0.0, "native"),
    y = grid::unit(0.22, "native"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontsize = 7.3, col = pal$muted)
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

summary_path <- file.path(output_dir, paste0(output_stem, "_workflow_summary.csv"))
utils::write.csv(
  data.frame(
    metric = c(
      "total_mag_tests",
      "primary_lmm_ok_tests",
      "nominal_raw_p_lt_0_05",
      "bh_fdr_significant",
      "singular_lmm",
      "too_sparse",
      "model_failed",
      "best_primary_bh_q",
      "job_id",
      "conda_environment",
      "main_results_file",
      "session_info_file",
      "execution_record_file"
    ),
    value = c(
      total_mags,
      primary_mags,
      nominal_mags,
      fdr_mags,
      singular_mags,
      too_sparse_mags,
      model_failed_mags,
      best_q,
      job_id,
      conda_env,
      "results/hpc_per_mag_lmertest/tables/per_mag_lmertest_hpc_results.csv",
      "results/hpc_per_mag_lmertest/session/session_info.txt",
      "results/hpc_per_mag_lmertest/reproducibility/hpc_execution_record.csv"
    )
  ),
  summary_path,
  row.names = FALSE
)

message("Figure 3 workflow outputs created:")
message(" - ", png_path)
message(" - ", pdf_path)
if (file.exists(svg_path)) message(" - ", svg_path)
message(" - ", summary_path)
