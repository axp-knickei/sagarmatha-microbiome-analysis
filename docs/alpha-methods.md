# Alpha-diversity methods and execution

Use the [alpha schematic](alpha-diversity-schematic.md) and
[plot guide](alpha-diversity-plots.md) with the [workflow commands](public-workflow.md#2-diversity-branches).
Preparation joins authorized metadata to the raw sample-by-MAG abundance table.
Alpha uses the raw-input branch, not the primary cleaned cohort automatically.

Shannon, Simpson and observed richness summarize sample-level diversity.
Count-based Chao1/ACE require suitable count-like input; preserve unavailable
metric flags rather than manufacturing counts from relative-abundance-like values.
R vegan and Python scikit-bio calculations are compared by the documented stages.
Heatmaps include raw values and within-subject z-scores; missing combinations
remain missing. Descriptive plots are not independent statistical tests.

Alpha LMMs use `metric_value ~ days_from_baseline + (1 | subjectID)`; the optional
adjusted model adds altitude, age and gender only when metadata completeness and
collinearity checks permit it. Inspect status tables and residual/QQ diagnostics.
Inference depends on the available fitting packages as described in the
[environment guide](public-environment.md). Window alpha fitting also requires
the locally regenerated window table and prepared alpha metadata.

```bash
Rscript scripts/alpha_diversity/06_plot_alpha_diversity_heatmap.R
Rscript scripts/alpha_diversity/08_create_alpha_diversity_publication_summary_table.R
qsub hpc/job_alpha_diversity_lmm.pbs
qsub hpc/job_alpha_diversity_window_lmm.pbs
```

Generated tables, figures and Markdown reports under `results/` and `docs/reports/`
are private analytical outputs and excluded from the public archive. Run rendering
in an isolated reproduction checkout; do not overwrite selected private evidence.
