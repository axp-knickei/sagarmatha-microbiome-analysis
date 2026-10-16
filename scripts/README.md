# Script Guide

This directory contains executable analysis scripts, post-processing scripts,
publication-output builders, diagnostics, and shared utilities. Run scripts
from the repository root unless their documentation states otherwise.

## Canonical workflow entrypoints

| Workflow | Entrypoint | Configuration | Main outputs | Status |
|---|---|---|---|---|
| Data preparation | `01_data_prep.R` | `config/config.yml` | `data/processed/cleaned_data.rds` | Shared active prerequisite |
| Per-MAG longitudinal LMM | `hpc_run_per_mag_lmertest.R` | `config/model_decision_rules.yml`, `config/response_scale_decision.yml` | `results/hpc_per_mag_lmertest/` | Primary inference workflow |
| Class 1 predefined-window LMM | `hpc_run_class1_window_lmertest.R` | Window definitions under `results/metadata_expedition_windows/` | `results/class1_window_lmertest/` | Secondary validation workflow |
| Alpha diversity | `alpha_diversity/01_prepare_alpha_diversity_inputs.R` | Environment variables and raw-data defaults documented in the scripts | `results/alpha_diversity/` | Supporting analysis |
| Alpha-diversity LMM | `alpha_diversity_lmm/01_prepare_alpha_lmm_input.R` | Prepared alpha-diversity outputs | `results/alpha_diversity_lmm/` | Supporting analysis |
| Alpha-diversity window LMM | `alpha_diversity_window_lmm/fit_plot_alpha_window_lmm.R` | Prepared alpha-diversity and window definitions | `results/alpha_diversity_window_lmm/` | Supporting analysis |
| Beta diversity | `run_beta_diversity_local.sh` | `config/beta_diversity.yml` | `results/beta_diversity/` | Active supporting pipeline; implementation is under `analysis/` |
| Random Forest | `run_random_forest_model.R` | `config/random_forest_config.yml` | `results/random_forest/` | Exploratory prediction workflow |
| Subject-aware Random Forest validation | `analysis/random_forest_subject_validation/run_random_forest_subject_validation.R` | `config/random_forest_subject_validation.yml` | `results/random_forest_subject_validation/` | Validation workflow |
| Subject-validation manuscript outputs | `random_forest_subject_validation/01_build_manuscript_outputs.R` | Completed subject-validation result files | `<run-id>/manuscript/` | Reader-first post-processing; does not refit the model |
| Random Forest MAG feature stability | `plot_random_forest_feature_stability.R` | Completed subject-validation feature tables plus `data/raw/GDTB.csv` | `<run-id>/figures/feature_stability/` | Publication post-processing; eligibility QC and exploratory native-importance stability |
| NetMoss | `netmoss/run_netmoss_pipeline.R` | `config/netmoss_analysis.yml` | `results/netmoss/` | Exploratory supporting analysis |
| NetMoss altitude contrast | `netmoss/run_netmoss_pipeline.R` | `config/netmoss_altitude_analysis.yml` | `results/netmoss_altitude/` | Separate exploratory T1/T2/T9/T10 versus T3–T8 analysis |
| NetMoss altitude stability | `netmoss/run_altitude_stability.R` | `config/netmoss_altitude_stability.yml` | `results/netmoss_altitude_stability/` | Matched-participant bootstrap and leave-one-out internal reproducibility analysis |
| Matched-timepoint NetMoss companion | `netmoss/run_matched_timepoint_companion.R` | `config/netmoss_matched_timepoint_companion.yml` | `results/netmoss_matched_timepoint_companion/` | Exploratory paired abundance analysis; not NetMoss network inference |
| MAG functional annotation | `run_mag_functional_annotation.sh` | `config/mag_functional_annotation.yml` | `results/mag_functional_annotation/` | Active staged HPC analysis |
| Functional community dynamics | `run_functional_community_dynamics.sh` | `config/functional_community_dynamics.yml` | `results/functional_community_dynamics/` | Active whole-profile taxonomic–functional comparison |
| Manuscript functional synthesis | `run_manuscript_functional_synthesis.sh` | `config/manuscript_functional_synthesis.yml` | `results/manuscript_functional_synthesis/` | Reader-first post-processing of the two completed functional workflows |
| Diversity/per-MAG manuscript synthesis | `manuscript_diversity_per_mag/01_build_main_figures_tables.R` | Completed alpha-, beta-, and per-MAG result tables | `results/manuscript_diversity_per_mag/` | Reader-first post-processing; does not rerun the source analyses |

HPC and scheduler entrypoints are indexed in [hpc/README.md](../hpc/README.md).

## Directory layout

- `alpha_diversity/` — prepare inputs, calculate metrics with vegan and
  scikit-bio, compare implementations, and build plots and summary tables.
  No stage `07` is present in the current tracked sequence; preserve the
  existing numbering and review history before filling or renumbering the gap.
- `alpha_diversity_lmm/` — prepare, fit, plot, and diagnose longitudinal
  alpha-diversity models.
- `alpha_diversity_window_lmm/` — fit and plot predefined biological-window
  models.
- `random_forest_subject_validation/` — validate and render the main manuscript
  figure, exact tables, panel data, and reproducibility metadata from one
  completed subject-aware Random Forest run without refitting it.
- `inspection/` — read-only or report-writing input inspection utilities.
  `inspection/report_abundance_totals.R` audits sample totals, summarizes the
  unrepresented fraction by subject and timepoint, and creates a total-sum-
  normalized sensitivity table without modifying raw input.
  `inspection/harmonize_pbmc_metadata.py` parses PBMC `rank` identifiers,
  reconciles source and microbiome dates, audits participant-timepoint and gut
  coverage, and writes a harmonized metadata crosswalk without modifying raw
  input.
- `netmoss/` — installation, environment validation, analysis, publication
  output, and workflow-local tests.
- `../analysis/mag_functional_annotation/` — full-file QC, catalogue analyses,
  longitudinal functional potential, paired contrasts, publication figures,
  and result indexing; launched through `run_mag_functional_annotation.sh`.
- `../analysis/functional_community_dynamics/` — identical-cohort functional
  and taxonomic profile distances, restricted-permutation tests, longitudinal
  movement, and exploratory recovery; launched through
  `run_functional_community_dynamics.sh`.
- `../analysis/manuscript_functional_synthesis/` — five-panel main figure,
  two compact main tables and exact panel-data exports; launched through
  `run_manuscript_functional_synthesis.sh` after both source workflows finish.
- `manuscript_diversity_per_mag/` — three manuscript main figures, three main
  tables, and an exhaustive supplementary evidence manifest from the completed
  alpha-diversity, beta-diversity, and primary per-MAG workflows.
- `utils/` — shared helpers sourced by multiple scripts.

The beta-diversity and subject-validation implementations are under
[`analysis/`](../analysis/README.md) because each is a self-contained staged
pipeline with local helpers and tests.

## Post-processing and publication scripts

Root-level `plot_figure_*`, `create_table_*`, and `create_supplementary_table_*`
scripts consume existing result tables and write publication-oriented files.
Their number is meaningful only within the analysis name embedded in the
filename.

The following Figure 2 scripts are retained as design history and must not be
treated as interchangeable without checking their output names:

- `plot_figure_2_class1_analysis.R`
- `plot_figure_2_class1_analysis_nature_style.R`
- `plot_figure_2_class1_analysis_nature_style_v2.R`

No variant is deleted or renamed by repository maintenance alone.

## Legacy and exploratory scripts

- `02_run_glmm.R` and `03_results.R` implement the original configurable
  single-model pipeline. They are maintained for compatibility but are not the
  primary per-MAG inference path.
- `run_per_mag_glmm.R`, `run_per_mag_lmertest.R`, `test_per_mag_models.R`, and
  related diagnostic scripts document earlier exploratory model development.
- `hpc_run_per_mag_lmertest.R` remains the canonical per-MAG implementation.

## Conventions

- Run from the repository root so relative paths and `source()` calls resolve.
- Keep scientific decisions in the appropriate YAML file rather than changing
  constants casually in scripts.
- Do not overwrite raw inputs.
- Treat `results/` as generated or retained evidence according to
  [results/README.md](../results/README.md).
- Add purpose, inputs, outputs, dependencies, usage, and workflow relationship
  to new user-facing script headers.

## v1.0.0 release boundary

Generated scientific figures and result tables are excluded from v1.0.0; their
producer scripts remain in the proposed source archive. The original
`results/metadata_expedition_windows/proposed_biological_windows.csv` is also
excluded because it includes aggregate study metadata. Its six predefined
windows and rationale remain in `evaluate_expedition_windows.R`. Before using
the Class 1/alpha-window models or their plotting/table scripts, follow
[window CSV regeneration and the complete dependency list](../docs/public-workflow.md#window-csv-regeneration)
with separately obtained authorized metadata matching the original 86-sample
processed selection. Default `results/` paths identify locally generated inputs,
not files supplied in the public archive.

The [release workflow](../docs/public-workflow.md) identifies current-publication
cohorts and exact entrypoints. The beta fixture generator now creates wholly
synthetic data without reading raw study inputs. Use
`config/netmoss_historical_publication.yml` only for the preserved historical
Ta-inclusive NetMoss scenario, with an isolated output override.

The untracked PBMC harmonization and abundance-total/normalization utilities
listed in this working guide remain outside v1.0.0; preserve them locally and
omit their navigation entries from the final public snapshot.

## Public rendering entrypoints

Per-MAG release renderers use `_public.R` suffixes. Synthesis runners call public
renderers and validators. Private originals are excluded from the proposed archive.
See [external input and validation instructions](../docs/public-rendering.md) and
the [replacement map](../docs/releases/v1.0.0-public-replacements.csv).


The combined diversity/per-MAG renderer is assembled from its approved public
replacement at the canonical entrypoint shown above. The private original is
excluded. See its [external-bundle and optional-evidence contract](manuscript_diversity_per_mag/README.md#public-replacement-and-external-bundles).

The RF feature-stability public replacement also occupies its canonical plotting
path. See [performance reporting and required inputs](../docs/public-rendering.md#random-forest-feature-stability).
