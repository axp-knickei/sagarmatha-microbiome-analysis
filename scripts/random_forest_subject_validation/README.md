# Random Forest subject-validation manuscript outputs

## Purpose

This post-processing workflow converts one completed subject-aware Random
Forest validation run into a reader-first manuscript package. It reads only the
saved audit, prediction, comparison, and validation files produced by the
canonical workflow. It does not read raw data, alter processed data, refit the
Random Forest, select a better strategy after seeing held-out performance, or
reinterpret feature importance as biological association.

The main figure has three panels:

1. the pooled leave-one-subject-out confusion matrix;
2. recall within each expedition-linked category; and
3. comparison of overall accuracy, balanced accuracy, and macro F1 with the
   majority and stratified-random baselines.

Two main tables provide exact overall/baseline metrics and per-category
precision, recall, and F1. Feature-importance rankings remain supplementary
because the validation model did not show reliable out-of-subject
discrimination and native importance is not a biological effect estimate.

## Entrypoint

Run from the repository root and name the completed validation run explicitly:

```bash
Rscript scripts/random_forest_subject_validation/01_build_manuscript_outputs.R \
  --run-dir results/random_forest_subject_validation/development_20260724_150733
```

On Nebula, submit the tracked OpenPBS wrapper:

```bash
qsub -v RFMS_RUN_DIR=results/random_forest_subject_validation/development_20260724_150733 \
  hpc/job_random_forest_subject_validation_manuscript.pbs
```

The script verifies that the selected run manifest is complete, requires every
canonical source-output check to have passed, and independently recomputes the
confusion-derived overall and class-specific metrics before rendering.
Production and development runs are accepted. Smoke runs are rejected unless
an explicit run directory and `--allow-smoke` are supplied for layout testing.

## Software environment

All required packages are already declared in `environment.yml`. On Nebula,
update the existing `glmm-hpc` environment after pulling this commit:

```bash
qsub hpc/update_glmm_env.pbs
```

Wait for that job to finish successfully before submitting the manuscript job.
The publication render uses `ggplot2`, `gridExtra`, `scales`, `svglite`, and
`ragg`. PDF and SVG are vector outputs; PNG and LZW-compressed TIFF are exported
at 600 dpi. No interactive package installation on a compute node is needed.

## Outputs

Outputs are written below the selected run:

```text
results/random_forest_subject_validation/<run-id>/manuscript/
  figures/
    main_figure_random_forest_subject_validation.{pdf,svg,png,tiff}
  tables/
    main_table_1_model_performance.{csv,tsv,md,pdf,svg,png,tiff}
    main_table_2_class_performance.{csv,tsv,md,pdf,svg,png,tiff}
  panel_data/
    figure_panel_a_confusion_matrix.csv
    figure_panel_b_class_recall.csv
    figure_panel_c_baseline_comparison.csv
  metadata/
    manuscript_output_validation.csv
    manuscript_render_manifest.yml
    manuscript_run_summary.csv
    session_info.txt
```

The machine-readable CSV/TSV tables retain numeric values. Presentation files
use reader-facing rounding. The palette values come from the Ordnance Survey
GeoDataViz Toolkit colour palettes v0.9.1; labels, values, category ordering,
and model-specific shapes keep the figure interpretable without color.

## Safe reruns

Existing outputs are protected by default. After reviewing a previous render,
replace only the selected run's `manuscript/` directory with:

```bash
qsub -v RFMS_RUN_DIR=results/random_forest_subject_validation/development_20260724_150733,RFMS_FORCE=1 \
  hpc/job_random_forest_subject_validation_manuscript.pbs
```

Never use a broad cleanup command against `results/`. Review `git status`
before and after the job because source result directories may contain selected
tracked evidence.
