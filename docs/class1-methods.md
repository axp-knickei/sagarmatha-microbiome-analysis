# Class 1 window methods and validation

The secondary Class 1 screen fits predefined three-timepoint windows. Regenerate
the excluded window CSV using the [authorized-metadata instructions](public-workflow.md#window-csv-regeneration).
The six definitions and rationale remain unchanged in the generator.

## Core Design

Each Class 1 window contains exactly three forward timepoints, such as:

```text
T1 -> T2 -> T3
T2 -> T3 -> T4
T4 -> T5 -> T6
T6 -> T7 -> T8
T7 -> T8 -> T9
T8 -> T9 -> T10
```

For each MAG-window test, the model is:

```r
log10(Y + pseudocount) ~ window_order + (1 | Subject)
```

where `window_order` is coded inside each window:

```text
first timepoint = 0
second timepoint = 1
third timepoint = 2
```

The fixed effect of `window_order` tests whether the MAG tends to increase or decrease across that predefined biological window.

## Run Full Class 1 Screen

Use the repository-declared `glmm-hpc` environment:

```bash
conda run -n glmm-hpc Rscript scripts/hpc_run_class1_window_lmertest.R
```

On Nebula HPC, load the Conda module first:

```bash
module load intelpython/2022.2.0
conda run -n glmm-hpc Rscript scripts/hpc_run_class1_window_lmertest.R
```

Or submit the PBS job:

```bash
qsub hpc/job_class1_window_lmertest.pbs
```

## First Outputs To Inspect

After the full run, inspect model behavior before biological interpretation:

```text
results/class1_window_lmertest/logs/class1_window_lmertest_summary.txt
results/class1_window_lmertest/tables/class1_window_direction_summary.csv
results/class1_window_lmertest/tables/class1_window_lmertest_results.csv
results/class1_window_lmertest/tables/class1_window_lmertest_results_primary_only.csv
```

Check these items in order:

1. Status counts per window.
2. Best `q_global_class1_bh`.
3. Number of `q_global_class1_bh < 0.05`.
4. Number of nominal raw `p < 0.05` rows.
5. Direction balance among nominal candidates.
6. Complete-subject support per window.

Do not start with biological stories or candidate interpretation. First check whether the model behavior is stable.

## Required Status Table

Before interpreting candidates, create or inspect a status table with:

```text
window_label
total_MAGs
lmm_ok
too_sparse
lmm_singular
low_variance
model_failed
percent_lmm_ok
```

This table tells which predefined windows are statistically usable.

The most important question is:

```text
Are the summit-related windows losing too many MAGs or subjects because of sparsity?
```

Pay special attention to:

```text
T6 -> T7 -> T8
T7 -> T8 -> T9
```

These windows are biologically important because the T7 -> T8 interval is figure-informed as the summit-push interval toward 8848.86 m. They may also be statistically fragile if many MAGs are too sparse.

## Interpretation Rules

Use conservative language:

```text
stage-local association
predefined biological window
secondary analysis
exploratory candidate
MAG-window signal
hypothesis-generating result
```

Avoid unsupported claims:

```text
discovery
mediation
altitude causes MAG change
summit MAG signature
```

Only consider stronger language if a MAG-window passes global Class 1 BH-FDR and later diagnostics, taxonomy review, and trajectory review support the interpretation.

## Overlapping Window Warning

Class 1 windows partially overlap. For example, `T6 -> T7 -> T8` and `T7 -> T8 -> T9` share two timepoints.

Therefore:

```text
Because predefined windows partially overlap, window-level tests are interpreted as related secondary stage-local analyses rather than independent biological discoveries. Benjamini-Hochberg FDR correction is applied globally across all predefined Class 1 MAG-window tests.
```

## Correction and review

Use `q_global_class1_bh` across all eligible MAG-window tests, not raw p-values
or correction within one selected window. Retain `lmm_ok`, singular, sparse,
low-variance and failed statuses. Review sample/subject support, diagnostics,
taxonomy and trajectories before interpreting candidates. Overlapping windows
are related secondary analyses. Recorded sampling altitude is distinct from
between-sample exposure; preserve the predefined biological rationale in code.
A later detection-model companion is a separate scientific decision, not a
release-packaging change.

## Accepted production reproducibility limitation

On 2026-10-04, the owner accepted the Class 1 production validation as
**passed with numerically negligible differences**. The workflow reproduced all
model classifications and scientific conclusions. Three MAG-window fits showed
small numerical differences from the historical reference that exceeded the
predefined strict comparison tolerance: absolute `1e-12` plus relative `1e-8`.
That strict comparison remains failed; its tolerance is unchanged.

Targeted reruns were deterministic and showed no convergence or singularity issue
for those three fits. Effect directions, eligibility, nominal and BH significance
decisions, multiple-testing conclusions and interpretation were unchanged.
Historical input bytes were not fully fingerprinted, fitted model objects were
unavailable, exact compiled numerical-library builds were unavailable, and CPU/
thread settings were incompletely recorded. Historical Rocky Linux 8.10 and
OpenBLAS 0.3.33 differ from validation Ubuntu 24.04.4 and OpenBLAS 0.3.30.
These are plausible computational contributors; the exact cause is unproven.

Detailed input-level evidence, model identities, comparison tables and fit objects
remain private. This acceptance applies to Class 1 only and does not certify other
production workflows or the final frozen snapshot.
