# Analysis Pipelines

This directory contains self-contained staged workflows whose implementation,
helpers, and tests are kept together. User-facing launchers may live in
`scripts/` or `hpc/`.

## Beta diversity

Location: `analysis/beta_diversity/`

Stages:

1. `00_inspect_repository.py` inventories candidate inputs without loading
   large files fully.
2. `01_validate_inputs.R` validates configuration and input structure.
3. `02_prepare_data.R` builds aligned abundance, metadata, and taxonomy data.
4. `03_compute_distances.R` computes configured distance matrices.
5. `04_run_ordination.R` runs ordination.
6. `05_run_permanova.R` runs restricted-permutation tests.
7. `06_run_longitudinal_analysis.R` evaluates longitudinal distance summaries.
8. `07_make_figures.R` creates figures.
9. `08_build_report.R` assembles the workflow report.
10. `09_create_smoke_test_data.R` creates deterministic synthetic test fixtures.

Run the workflow with:

```bash
bash scripts/run_beta_diversity_local.sh --config config/beta_diversity.yml
```

Use `config/beta_diversity_test.yml` and `make beta-test` for the smoke-test
path. Workflow-local helpers are under `utils/`; tests are under `tests/`.
Detailed methods and output descriptions are in
[`docs/beta_diversity/workflow.md`](../docs/beta_diversity/workflow.md).

## Subject-aware Random Forest validation

Location: `analysis/random_forest_subject_validation/`

`run_random_forest_subject_validation.R` orchestrates audit, grouped
resampling, baseline comparisons, leave-one-subject-out Random Forest,
sensitivity comparisons, summaries, figures, and output validation. Shared
functions are in `R/common.R`.

Run from the repository root:

```bash
Rscript analysis/random_forest_subject_validation/run_random_forest_subject_validation.R smoke
```

The workflow reads `config/random_forest_subject_validation.yml` and writes a
named run directory under `results/random_forest_subject_validation/`. See
[`docs/random-forest-methods.md`](../docs/random-forest-methods.md)
for interpretation and limitations.

After a run is complete, build its reader-first manuscript figure and tables
without refitting the model with
`scripts/random_forest_subject_validation/01_build_manuscript_outputs.R` or
`hpc/job_random_forest_subject_validation_manuscript.pbs`.

## MAG functional annotation

Location: `analysis/mag_functional_annotation/`

This staged, QC-gated OpenPBS workflow audits global MMseqs protein clusters,
propagates eggNOG representative annotations to member MAGs, profiles catalogue
convergence and unknown proteins, joins gut MAG taxonomy, constructs
unnormalized longitudinal predicted functional potential, and tests
prespecified equal-altitude paired contrasts. See
[`docs/functional-methods.md`](../docs/functional-methods.md).

## Functional community dynamics

Location: `analysis/functional_community_dynamics/`

This independent workflow compares predicted KO, KEGG-module, EC, and COG
profile restructuring with the 441-MAG taxonomic profile on the exact same
T1–T10 sample set. It runs participant-restricted PERMANOVA with dispersion
checks, matched transition models, participant-level transition sensitivity,
and an explicitly exploratory T1-to-T10 recovery summary. See
[`docs/functional-community-methods.md`](../docs/functional-community-methods.md).

## Manuscript functional synthesis

Location: `analysis/manuscript_functional_synthesis/`

This post-processing workflow reads the two completed functional-analysis
result trees without changing them. It creates a five-panel, 183-mm-wide main
figure, two editable main tables with PDF/PNG previews, exact panel data and a
completion marker. See the
[workflow README](manuscript_functional_synthesis/README.md) for required
source results, invocation, overrides and interpretation limits. The related
draft manuscript remains private and is excluded from v1.0.0; it is not a
workflow input.

## Why these workflows are not under `scripts/`

These analyses were developed as cohesive staged pipelines with their own
helpers and tests. This directory boundary is historical but functional. Do
not move either pipeline without updating shell runners, scheduler wrappers,
`source()` calls, tests, and documentation together.

## Public synthesis replacements

Functional and NetMoss synthesis runners use `_public.R` renderers and validators.
Supply completed outputs externally; see the [public rendering guide](../docs/public-rendering.md).
