# Random Forest subject-aware validation

## Scientific question and scope

This workflow asks whether gut MAG abundance profiles predict the four
expedition-linked altitude categories (`0`, `3650`, `5200`, `6500`) for subjects
excluded from model training. The four-class target is retained for continuity
with the legacy Random Forest analysis.

Altitude is perfectly determined by timepoint and sampling date in the current
dataset. The workflow therefore evaluates prediction among confounded
altitude/time expedition categories. It cannot isolate altitude from chronology,
expedition phase, diet, treatment, or other time-varying conditions.

## Authoritative and historical files

The current entry point is
`analysis/random_forest_subject_validation/run_random_forest_subject_validation.R`;
its authoritative configuration is
`config/random_forest_subject_validation.yml`.

The older `scripts/run_random_forest_model.R`,
`config/random_forest_config.yml`, `hpc/job_random_forest.pbs`, and
`results/random_forest/rf_*` directories are historical. They are audited but
never modified by this workflow.

## Workflow map

```text
cleaned_data.rds (long: sampleID × MAG)
  │
  ├─ 01 input, legacy-split, class, coverage, and confounding audit
  ├─ 03 LOSO majority and repeated stratified-random baselines
  ├─ 04 primary unweighted LOSO Random Forest
  ├─ 05 secondary repeated grouped five-fold CV
  ├─ 06 LOSO class-imbalance strategy comparison
  ├─ 07 LOSO training-variance feature-count comparison
  ├─ 11 LOSO tree-count stability
  ├─ 08 implemented/disabled outcome status table
  ├─ 09 reports, subject-macro summaries, and importance-rank stability
  ├─ 10 PNG/PDF/SVG figures generated from run CSV files
  └─ 12 programmatic output-integrity validation
       └─ manifest status changes from running to complete
```

Every stage is a fresh non-interactive R process. Run from the repository root;
repository-relative paths are intentional.

## Input requirements

Required columns are configured, currently `sampleID`, `Subject`, `MAG`, `Y`,
and `altitude`; `Timepoint`, `date_converted`, and
`days_from_baseline_subject` support the audit.

The input must contain one unique row for every sample–MAG combination. All
samples must contain the same MAG identifiers, abundance must be numeric,
finite, and non-negative, and metadata must be constant within a sample.
Missing combinations are rejected; they are not assumed to be biological
zeros. Current sample abundance sums are close to, but not always, one, so the
scale is described as relative-abundance-like. Exact upstream normalization is
not encoded in the processed RDS.

## Leakage-safe preprocessing

Inside each outer training fold, the workflow:

1. calculates missingness, sample prevalence, subject prevalence, mean
   abundance, and variance;
2. retains features meeting configured thresholds;
3. optionally ranks eligible MAGs by training-fold variance;
4. estimates missing-value medians from training samples;
5. applies the configured `log10(x + pseudocount)` transformation;
6. derives class weights or sampling targets from training outcomes; and
7. fits the Random Forest.

The held-out subject cannot affect these operations. Identifiers, outcomes,
timepoints, dates, and other metadata are not predictors. The pseudocount is a
fixed, prespecified configuration value inherited from the earlier data audit;
it is not re-estimated during resampling.

## Resampling, baselines, and metrics

LOSO-CV is primary: each included subject is assessed once. Pooled metrics
weight samples equally. Subject-macro accuracy weights subjects equally.

Subject folds often lack altitude classes. Their `balanced_accuracy` is
explicitly a present-class recall average and is descriptive only.
`balanced_accuracy_all_classes` is `NA` unless all four classes occur. A fold with one observation can have perfect descriptive accuracy without
providing strong evidence of generalization.

Metrics include accuracy, balanced accuracy, macro precision/recall/F1,
multiclass log loss, ordinal-rank MAE and RMSE, within-one-category accuracy,
and quadratic weighted kappa. Probability rows must be finite, bounded,
class-ordered, and sum to one.

The majority baseline is learned in each outer training fold. The
stratified-random baseline draws from each training fold’s class proportions
using deterministic repeated seeds; its mean and Monte Carlo standard deviation
are reported. Both baselines use the training prevalence distribution for log
loss; their hard class predictions differ.

Inverse-frequency weights and balanced within-tree sampling are derived
separately in each outer training fold. Their combination intentionally applies
both corrections and can over-emphasize minority classes; poor combined
performance is not by itself evidence of a programming error.

## Configuration modes

- `smoke`: 50-tree primary fits, 25/50-tree stability, one grouped-CV repeat,
  and 100 stratified-random repetitions.
- `development`: 500-tree primary fits, 500–5,000-tree stability, ten grouped
  repeats, and 500 baseline repetitions.
- `production`: 2,000-tree primary fits, 500–10,000-tree stability, ten grouped
  repeats, and 1,000 baseline repetitions.

All modes execute the same stages. The resolved mode is propagated to every
subprocess and the resolved configuration is copied into the run directory.

Hyperparameters are fixed, not nested-tuned. Strategy comparisons reuse outer
LOSO assessment predictions; choosing a best strategy afterward would make its
reported performance optimistic. Confirmatory tuning would require grouped
inner resampling inside every LOSO training fold and an independently locked
final workflow.

## Public synthetic validation

The public test at `tests/testthat/test-random-forest-subject-validation.R` is
assembled from its approved public replacement; the private original is excluded.
Run all public suites from the repository root:

```bash
Rscript tests/public_release/run-public-tests.R /path/to/local/test-logs
```

The RF fixture is generated inside a temporary directory from an explicit formula:
six `SYNTHETIC_SUBJECT_*` identifiers, eight `SYNTHETIC_MAG_*` identifiers, four
fictional numeric labels (11, 22, 33, 44), two visits per label and subject, and
artificial nonnegative abundances. The resulting complete sample-by-feature grid
has 48 samples and 384 rows; these are fixture choices, not study cohort counts.
Numeric labels meet the existing ordered-class contract and are not real altitudes.
One feature is deliberately constant to exercise filtering. No study input or
scientific output is read, and no fixture values come from those files.

Only test-local paths, the fixed fixture seed, smoke mode and fictional labels
are set for these checks. Production configuration and all model, preprocessing,
class-balancing and resampling implementations remain unchanged. Checks cover
input rejection, subject separation, training-only feature fitting, assessment
perturbation, baseline and RF reproducibility, LOSO and grouped-CV prediction
coverage/probabilities, class-balancing modes, metrics and output collision safety.
Hand-calculated metric expectations describe artificial prediction examples.
No accuracy target or feature ranking is prescribed for the fitted synthetic RF.
The existing input validator rejects missing abundances with its nonfinite-value
error; the test records this behavior without changing the validator.

These software checks do not certify study performance, real-data class coverage,
biological feature stability or end-to-end production audit/report completeness.
Validate the full orchestrator and production outputs with separately authorized
inputs and record diagnostics/environment provenance before release approval.
The public runner can generate synthetic results/plot artifacts; those artifacts
must remain outside the final source inventory. Public staging tests are a local
checkpoint, not frozen-snapshot or scientific-validation approval.

## Local execution and outputs

From the repository root:

```bash
Rscript -e 'testthat::test_file(
  "tests/testthat/test-random-forest-subject-validation.R",
  reporter = "summary"
)'

Rscript analysis/random_forest_subject_validation/run_random_forest_subject_validation.R \
  smoke smoke_YYYYMMDD_review
```

Run IDs must be unique. Outputs are stored under:

```text
results/random_forest_subject_validation/<run-id>/
├── audit/
├── predictions/
├── comparisons/
├── features/
├── figures/{png,pdf,svg}/
└── reports/
```

`reports/run_manifest.yml` records status, timestamps, hashes, Git state,
package versions, seed, and completed/failed stages. Only `status: complete`
means every stage and output-integrity check succeeded.

## Manuscript post-processing

A completed run can be converted into one reader-first three-panel main figure,
two exact main tables, panel-data exports, and rendering validation metadata
without refitting the Random Forest:

```bash
Rscript scripts/random_forest_subject_validation/01_build_manuscript_outputs.R \
  --run-dir results/random_forest_subject_validation/<run-id>
```

The post-processing script reads only saved subject-validation outputs,
independently reconstructs confusion-derived metrics, rejects incomplete runs,
and protects existing manuscript outputs unless `--force` is supplied. Nebula
execution and the complete output inventory are documented in
[`scripts/random_forest_subject_validation/README.md`](../scripts/random_forest_subject_validation/README.md).

## HPC execution

The prepared job `hpc/job_random_forest_subject_validation.pbs` requests one
CPU because the implementation is sequential, 8 GB memory, and 12 hours. These
are conservative limits; smoke runtime alone is insufficient for precise
production resource prediction.

Submit manually from the Nebula repository root:

```bash
qsub -v RF_VALIDATION_MODE=production \
  hpc/job_random_forest_subject_validation.pbs
```

Do not submit until a development-mode tree-stability run has been reviewed.

## Implemented, disabled, and unsupported

- Implemented: four-class classification and ordinal metrics from the same
  predictions.
- Exploratory comparisons: imbalance strategies, training-variance feature
  counts, tree counts, and repeated grouped CV.
- Disabled and not implemented: numeric-altitude regression and within-subject
  baseline change.
- Disabled: phase prediction because no authoritative mapping exists.
- Unsupported: taxonomic aggregation because no reproducible mapping exists.

## Troubleshooting

- “Refusing to reuse output directory”: choose a new run ID.
- “Incomplete sample-by-feature grid”: regenerate valid processed data; do not
  fill unknown measurements with zero.
- “Observed outcome classes do not match”: check configured levels and metadata.
- A `timedatectl` warning during `sessionInfo()` does not affect model fitting.
- Failed runs retain `status: failed` and `failed_stage` in their manifest.


## Evaluation boundaries


The four-class analysis is retained as the primary comparison with the legacy
model. Exact accuracy must be read alongside balanced accuracy, macro F1,
per-class recall, ordinal rank error, and majority/stratified-random baselines.

LOSO folds are subjects, not equally informative experimental replicates:
subjects contribute different numbers of samples and altitude classes, and some participants
may have only one observation. Pooled out-of-fold metrics answer sample-weighted
generalization; fold distributions show heterogeneity but should not be treated
as equally precise independent estimates. Individual samples must not be bootstrapped
as independent observations.

Feature importance is exploratory predictive importance. A high or stable MAG
rank does not establish univariate association, longitudinal change,
statistical significance, biological mechanism, or causality. Candidate MAGs
require separate, multiplicity-controlled longitudinal testing and taxonomy
review.

Native importance magnitudes are not directly comparable across forests with
different subjects or predictor sets. The workflow also summarizes within-fold
importance ranks, but rank recurrence remains sensitive to prevalence,
variance, and correlated predictors.

The recommended scientific framing is an exploratory comparison of exact
four-class and ordinal evaluation. Phase, numeric-altitude, and
within-subject-change analyses are currently disabled. Within-subject change
would require baseline data from a new subject and answers a different
prediction question from prediction without subject-specific history.


## Publication figure execution



From the repository root, select a completed run explicitly:

```bash
conda run -n glmm-hpc Rscript scripts/create_random_forest_publication_figure.R \
  --run-dir results/random_forest/rf_YYYYMMDD_HHMMSS
```

The default output is `RUN_DIR/publication/`. Existing output is protected; pass `--force` only when replacement is intentional.

For a faster formatting check:

```bash
conda run -n glmm-hpc Rscript scripts/create_random_forest_publication_figure.R \
  --run-dir results/random_forest/rf_YYYYMMDD_HHMMSS \
  --bootstrap-replicates 100 \
  --output-dir /tmp/blackice-rf-figure-check
```

## Nebula OpenPBS execution

Submit from the repository root after the Random Forest production job completes:

```bash
mkdir -p logs
qsub -v RF_RUN_DIR=results/random_forest/rf_YYYYMMDD_HHMMSS \
  hpc/job_random_forest_publication.pbs
```

The default uses 1,000 subject-bootstrap replicates. Override it when required:

```bash
qsub -v RF_RUN_DIR=results/random_forest/rf_YYYYMMDD_HHMMSS,RF_BOOTSTRAPS=2000 \
  hpc/job_random_forest_publication.pbs
```



## Feature stability



Eligibility filtering was estimated independently from each training set
before its held-out participant was assessed. Because `preprocessing.top_n` is
disabled, every eligible MAG entered its fold-specific model. Inclusion
frequency is consequently a filtering quality-control measure, not repeated
identification as an important feature.

Taxonomy is joined by the complete, unchanged MAG identifier to `GDTB.csv`
column `Genome`. The project intentionally retains the source filename
`GDTB.csv`, although the database is normally abbreviated GTDB. The plotting
workflow validates duplicate taxonomy keys, exact match coverage, and unchanged
row counts before producing figures.

Native `MeanDecreaseAccuracy` is summarized by its median and interquartile
range across eligible folds. Within-fold ranks are independently recomputed and
summarized by their median and IQR. These are algorithm-specific predictive
rankings, not effect sizes or significance tests. Predictive importance does not establish association, significance, mechanism or causality.


## Public feature-stability reporting

Use the canonical feature-stability renderer with completed validation outputs;
its public-safe replacement preserves the private original, calculations, ranking
rules and input consistency checks. It does not train a model. Performance text
reports supplied quantitative values without preset discrimination judgments.
Absent or incomplete summaries remain unavailable. The existing input/layout
limitations and synthetic reporting checks are described in the
[rendering contract](public-rendering.md#random-forest-feature-stability).
