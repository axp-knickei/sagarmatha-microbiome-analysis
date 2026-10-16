# Alpha-Diversity Schematic

This page captures the confirmed high-level alpha-diversity workflow for the repository.

## Confirmed Flow

```text
gut.csv / abundance table
    -> calculate alpha-diversity per sample
    -> alpha_diversity_with_metadata.csv
    -> fit LMM for each alpha metric
    -> plot observed trajectories + model-estimated trend
```

## Implementation Map

- Input abundance table:
  - `data/raw/gut.csv`
- Alpha-diversity preparation and metric calculation:
  - `scripts/alpha_diversity/01_prepare_alpha_diversity_inputs.R`
  - `scripts/alpha_diversity/02_calculate_alpha_diversity_vegan.R`
  - `scripts/alpha_diversity/03_calculate_alpha_diversity_scikitbio.py`
- Metadata join and comparison output:
  - `scripts/alpha_diversity/04_compare_alpha_diversity_tools.R`
  - output: `results/alpha_diversity/tables/alpha_diversity_with_metadata.csv`
- Alpha-diversity LMM fitting:
  - `scripts/alpha_diversity_lmm/01_prepare_alpha_lmm_input.R`
  - `scripts/alpha_diversity_lmm/02_fit_alpha_diversity_lmm.R`
- Trajectory plotting and diagnostics:
  - `scripts/alpha_diversity_lmm/03_plot_alpha_diversity_lmm_trends.R`
  - `scripts/alpha_diversity_lmm/04_alpha_lmm_diagnostics.R`

## What This Is

This is a sample-level alpha-diversity workflow. It is not the per-MAG longitudinal LMM workflow used in `hpc/job_per_mag_lmertest.pbs`.

The outputs are intended to support:

- observed-by-timepoint summaries,
- subject-level longitudinal trajectories,
- model-estimated trends by alpha metric,
- QC and diagnostic review.

## Related Documents

- [Release Workflow](public-workflow.md)
- [Environment Guide](public-environment.md)
- [Alpha-Diversity Plots](alpha-diversity-plots.md)
- [Alpha-Diversity LMM Report](public-workflow.md)
