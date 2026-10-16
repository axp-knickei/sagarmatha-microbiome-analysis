# Alpha-Diversity Plots

This page records the alpha-diversity figures that have been generated so far and the scripts that produce them.

## Publication-Style Alpha-Diversity Figures

Script:
- [`scripts/alpha_diversity/05_plot_alpha_diversity_publication.R`](../scripts/alpha_diversity/05_plot_alpha_diversity_publication.R)

Figures:
- `alpha_diversity_boxplot_jitter_by_timepoint`
- `alpha_diversity_violin_boxplot_jitter_by_timepoint`
- `alpha_diversity_longitudinal_spaghetti_by_subject`
- `alpha_diversity_faceted_multimetric_panel`
- `alpha_diversity_model_trend_days_from_baseline`
- `alpha_diversity_subject_timepoint_heatmap`

## Heatmap Figures

Script:
- [`scripts/alpha_diversity/06_plot_alpha_diversity_heatmap.R`](../scripts/alpha_diversity/06_plot_alpha_diversity_heatmap.R)

Figures:
- `alpha_diversity_heatmap_shannon_raw`
- `alpha_diversity_heatmap_shannon_zscore`
- `alpha_diversity_heatmap_simpson_raw`
- `alpha_diversity_heatmap_simpson_zscore`
- `alpha_diversity_heatmap_observed_richness_raw`
- `alpha_diversity_heatmap_observed_richness_zscore`
- `alpha_diversity_heatmap_multimetric_raw`
- `alpha_diversity_heatmap_multimetric_zscore`

## Alpha-Diversity LMM Trend Figures

Script:
- [`scripts/alpha_diversity_lmm/03_plot_alpha_diversity_lmm_trends.R`](../scripts/alpha_diversity_lmm/03_plot_alpha_diversity_lmm_trends.R)

Figures:
- `alpha_lmm_observed_model_trend_shannon`
- `alpha_lmm_observed_model_trend_simpson`
- `alpha_lmm_observed_model_trend_observed_richness`
- `alpha_lmm_observed_model_trend_chao1`
- `alpha_lmm_observed_model_trend_ace`
- `alpha_lmm_observed_model_trend_all_metrics`

## Alpha-Diversity LMM Diagnostic Figures

Script:
- [`scripts/alpha_diversity_lmm/04_alpha_lmm_diagnostics.R`](../scripts/alpha_diversity_lmm/04_alpha_lmm_diagnostics.R)

Figures:
- `alpha_lmm_residuals_vs_fitted_shannon`
- `alpha_lmm_qq_shannon`
- `alpha_lmm_random_intercepts_shannon`
- `alpha_lmm_residuals_vs_fitted_simpson`
- `alpha_lmm_qq_simpson`
- `alpha_lmm_random_intercepts_simpson`
- `alpha_lmm_residuals_vs_fitted_observed_richness`
- `alpha_lmm_qq_observed_richness`
- `alpha_lmm_random_intercepts_observed_richness`

## Notes

- The publication figures summarize the sample-level alpha-diversity workflow.
- The heatmaps are descriptive visualizations of subject-timepoint patterns.
- The LMM trend and diagnostic figures summarize the fitted sample-level longitudinal models.
