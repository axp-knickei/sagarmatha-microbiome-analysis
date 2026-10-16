# HPC Per-MAG LMM Plots

This page records the plotting outputs for the per-MAG longitudinal LMM abundance workflow and the scripts that generate them.

## Integrated HPC Entrypoint

The per-MAG analysis and plotting workflow is now wired through `hpc/job_per_mag_lmertest.pbs`.

- `RUN_STAGE=auto` or `RUN_STAGE=full` runs the per-MAG analysis and then launches the plotting orchestrator when the job is non-chunked.
- `RUN_STAGE=analysis` runs only the model-fitting stage.
- `RUN_STAGE=postprocess` combines chunk outputs, restores the canonical result files, and runs the plotting scripts.
- The plotting orchestrator is [`scripts/hpc_plot_per_mag_lmertest_results.R`](../scripts/hpc_plot_per_mag_lmertest_results.R).
- Chunk aggregation is handled by [`scripts/combine_hpc_per_mag_chunks.R`](../scripts/combine_hpc_per_mag_chunks.R).

## Main Per-MAG Result Figures

### Figure 1: Study design and sampling structure
- Script: [`scripts/plot_figure_1_hpc_per_mag_lmertest_results.R`](../scripts/plot_figure_1_hpc_per_mag_lmertest_results.R)
- Purpose: metadata-only sampling design schematic for the longitudinal MAG study.

### Figure 2: Abundance transformation rationale
- Script: [`scripts/plot_figure_2_hpc_per_mag_lmertest_results.R`](../scripts/plot_figure_2_hpc_per_mag_lmertest_results.R)
- Purpose: shows sparsity, zero inflation, and the `log10(Y + pseudocount)` modeling scale.

### Figure 3: HPC workflow schematic
- Script: [`scripts/plot_figure_3_hpc_per_mag_lmertest_results.R`](../scripts/plot_figure_3_hpc_per_mag_lmertest_results.R)
- Purpose: workflow diagram from processed MAG table through HPC execution and reproducibility logging.

### Figure 4: Result landscape
- Script: [`scripts/plot_figure_4_hpc_per_mag_lmertest_results_public.R`](../scripts/plot_figure_4_hpc_per_mag_lmertest_results_public.R)
- Purpose: summarizes top raw p-value MAGs and the result landscape.

### Figure 5: Trajectory summary
- Script: [`scripts/plot_figure_5_hpc_per_mag_lmertest_results_public.R`](../scripts/plot_figure_5_hpc_per_mag_lmertest_results_public.R)
- Purpose: selected MAG abundance trajectories and trend summaries.

### Figure 6: Taxonomy summary
- Script: [`scripts/plot_figure_6_hpc_per_mag_lmertest_results_public.R`](../scripts/plot_figure_6_hpc_per_mag_lmertest_results_public.R)
- Purpose: taxonomy context for significant and nominal MAGs.

### Figure 7: Model status audit
- Script: [`scripts/plot_figure_7_hpc_per_mag_lmertest_results_public.R`](../scripts/plot_figure_7_hpc_per_mag_lmertest_results_public.R)
- Purpose: model status summary and flagged MAG audit table.

### Figure 8: Diagnostics
- Script: [`scripts/plot_figure_8_hpc_per_mag_lmertest_results_public.R`](../scripts/plot_figure_8_hpc_per_mag_lmertest_results_public.R)
- Purpose: residual and random-intercept diagnostics for selected MAGs.

### Figure 9: Sensitivity analysis
- Script: [`scripts/plot_figure_9_hpc_per_mag_lmertest_results.R`](../scripts/plot_figure_9_hpc_per_mag_lmertest_results.R)
- Purpose: compares sensitivity checks for the per-MAG result set.

### Figure 10: Selected-MAG heatmap
- Script: [`scripts/plot_figure_10_hpc_per_mag_lmertest_results_public.R`](../scripts/plot_figure_10_hpc_per_mag_lmertest_results_public.R)
- Purpose: heatmap-style view of selected MAGs and their values across timepoints.

## Figure Outputs

Each figure is written as:

- `.png`
- `.pdf`
- `.svg`

Most figures also emit one or more companion CSV summary tables in `results/hpc_per_mag_lmertest/figures/`.

## Supporting Plotting Scripts

- [`scripts/summarize_hpc_candidate_signals.R`](../scripts/summarize_hpc_candidate_signals.R) - candidate signal heatmaps and trajectories
- [`scripts/diagnose_response_scale_models.R`](../scripts/diagnose_response_scale_models.R) - response-scale diagnostics and trajectory checks

## Notes

- These plots belong to the per-MAG abundance LMM workflow, not the alpha-diversity workflow.
- The result figures are descriptive and diagnostic. They support interpretation of the per-MAG model table rather than replacing it.
- The primary scientific table remains `results/hpc_per_mag_lmertest/tables/per_mag_lmertest_hpc_results.csv`.
