# Public analytical workflow

Run from the repository root in a separate reproduction checkout with authorized
inputs. This guide maps existing behavior; no scientific results were regenerated
during release preparation. The [manifest](release-workflow-manifest.csv) records
entrypoints and output roots. All generated research outputs remain outside this code archive.

## 1. Inputs, cleaning and primary inference

```bash
Rscript scripts/01_data_prep.R
Rscript scripts/validate_processed_data.R
Rscript scripts/hpc_run_per_mag_lmertest.R
Rscript scripts/summarize_hpc_candidate_signals.R
```

`01_data_prep.R` uses `config/config.yml` for inputs. Its legacy model section
is not the primary LMM authority: consult `model_decision_rules.yml` and
`response_scale_decision.yml`. The primary per-MAG result table is
`results/hpc_per_mag_lmertest/tables/per_mag_lmertest_hpc_results.csv`.
Diagnostics, logs and session records are sibling directories. Both all-model
and primary-eligible BH columns are retained; significance uses the latter.
Chunked runs require `scripts/combine_hpc_per_mag_chunks.R` to recompute global
FDR. Never interpret chunk-local p-values as final corrected results.

`investigate_row_loss.R`, `check_model_structure.R`, response-scale diagnostic
scripts and manual-review preparation retain decision provenance. They are not
replacement primary workflows. Earlier `test_per_mag_models.R` and
`summarize_per_mag_model_test.R` describe exploratory development. They need
not be rerun to package the existing publication evidence.

Class 1 window validation depends on a locally regenerated window CSV; the
original is excluded from v1.0.0. Complete [window CSV regeneration](#window-csv-regeneration)
with authorized metadata using `scripts/evaluate_expedition_windows.R`, followed by
`scripts/hpc_run_class1_window_lmertest.R`. Follow
[class1-methods.md](class1-methods.md) for eligibility,
global FDR and overlapping-window caveats. Preserve existing figure variants as
design history; select manuscript figures explicitly rather than choosing by name.

### Window CSV regeneration

**Release decision:** exclude the original
`results/metadata_expedition_windows/proposed_biological_windows.csv` from the
public GitHub release and Zenodo archive. It combines predefined analysis
settings with aggregate study metadata: date/altitude spans, participant counts,
model-support assessments and related caveats. It has no participant identifiers
or MAG measurements, but aggregation does not establish redistribution permission.
The public release does not supply the research metadata or processed dataset.
Generated scientific figures and result tables are also excluded; their producer
scripts remain part of the proposed source archive.

The six ordered windows and scientific rationale are retained in `window_specs`
in `scripts/evaluate_expedition_windows.R`: T1/T2/T3, T2/T3/T4, T4/T5/T6,
T6/T7/T8, T7/T8/T9 and T8/T9/T10. The generator adds summaries from metadata;
the complete CSV cannot be recreated from those definitions alone.

#### Required authorized inputs and cohort

The generator reads `config/config.yml`. It first looks for `sample_info.csv`
under `BLACKICE_RAW_DIR` (or `RAW_DATA_DIR`, then `data.raw_dir`, normally
`data/raw/`). If that file is absent, it reads `data.processed`, normally
`data/processed/cleaned_data.rds`. It does **not** honor the preprocessing
script's `BLACKICE_PROCESSED_PATH` override; use the configured processed path.

The preserved generation report records the fallback route: unique sample-level
metadata extracted from the original 86-sample processed dataset. That selection
has 13 subjects and T1–T10 observations, excluding Ta/Tb. Reproduce that selection,
including the same sample IDs, subject assignments, dates and altitudes. A count
of 86 alone is insufficient. The unfiltered raw metadata has 119 rows; allowing
the generator to select it can change the summaries even though the six windows
stay predefined. Do not substitute the 103-sample diversity cohort.

For the fallback route, the processed RDS must contain `sampleID`, `Timepoint`,
`date`, `altitude`, and `Subject` or `group`. Retain `date_converted` from
preprocessing when present. The script deduplicates by sample ID and calculates
summaries; it does not fit MAG models. Obtain the original authorized processed
dataset separately, or generate it with `scripts/01_data_prep.R` from the original
authorized `gut.csv` and `sample_info.csv`, following the [input contract](release-data-requirements.md).
Verify the selection privately against the original study records before using
it as historical reproduction evidence.

#### Commands and generated files

Use an isolated reproduction checkout with the documented R environment and
`yaml` installed so the configured input paths are read. Place the authorized
processed RDS at `data/processed/cleaned_data.rds` with the current configuration.
Check its cohort before running the generator; do not run these generation
commands in a checkout containing selected evidence you need to preserve.

```bash
# Optional preparation when starting with the original authorized raw inputs:
# place gut.csv and sample_info.csv in data/raw/, then run:
Rscript scripts/01_data_prep.R

# Check the processed cohort; this does not certify sample identity or input hashes.
Rscript -e 'd <- readRDS("data/processed/cleaned_data.rds"); stopifnot(length(unique(d$sampleID)) == 86L, length(unique(d$Subject)) == 13L, all(as.character(d$Timepoint) %in% paste0("T", 1:10)))'

# An empty temporary directory prevents the generator from preferring raw metadata.
WINDOW_RAW_GUARD_DIR="$(mktemp -d)"
BLACKICE_RAW_DIR="$WINDOW_RAW_GUARD_DIR" Rscript scripts/evaluate_expedition_windows.R
rmdir "$WINDOW_RAW_GUARD_DIR"
```

Skip the optional preparation command if an authorized matching processed RDS
has already been supplied. With a custom `data.processed` value, use that same
path in the cohort check. Do not overwrite or modify the authorized raw files.

The generator creates `results/metadata_expedition_windows/`, writes the window
CSV, `timepoint_altitude_date_summary.csv`, `consecutive_timepoint_transitions.csv`
and a Markdown report, and may write private execution records through the shared
helper. Its output directory is fixed; it has no output-directory override.
These generated research materials remain outside the public archive. Verify six
rows in the preserved order and compare metadata selection, input hashes and CSV
values with the original private evidence. No new full-table reproduction was
performed for this documentation change; the checklist's reproduction gate remains
open. Deterministic calculations require the same inputs, code and compatible
runtime, and do not alone certify byte-for-byte agreement with a historical file.

#### Dependent workflows and launchers

The following defaults refer to a file generated locally, **not a bundled CSV**.
Regenerate it before running these scripts; model and alpha input datasets are
also obtained or generated separately with authorization.

| Consumer | Required preparation or existing override |
|---|---|
| `scripts/hpc_run_class1_window_lmertest.R` | Regenerated CSV at the default path, or `--windows-csv PATH` |
| `scripts/alpha_diversity_window_lmm/fit_plot_alpha_window_lmm.R` | Regenerated CSV, or `BLACKICE_WINDOW_DEFINITION_PATH`; also prepared alpha metadata via `BLACKICE_ALPHA_METADATA_PATH` |
| `scripts/plot_main_figure_1_class1_analysis.R` | CSV and `timepoint_altitude_date_summary.csv` at their generated default paths |
| `scripts/plot_figure_2_class1_analysis.R`, `scripts/plot_figure_2_class1_analysis_nature_style.R`, `scripts/plot_figure_2_class1_analysis_nature_style_v2.R` | CSV at the default path plus completed Class 1 results |
| `scripts/plot_figure_3_class1_analysis.R`, `plot_figure_4_class1_analysis.R`, `plot_figure_5_class1_analysis.R`, `plot_figure_6_class1_analysis.R`, `plot_figure_9_class1_analysis.R`, `plot_figure_10_class1_analysis.R` | CSV at the default path plus completed Class 1 results |
| `scripts/create_supplementary_table_1_class1_analysis.R` | CSV at the default path; reports its metadata summaries and rationale |
| `scripts/create_supplementary_table_2_class1_analysis.R` | CSV at the default path plus completed Class 1 results |
| `hpc/job_class1_window_lmertest.pbs` | Generate the CSV before submission; `BLACKICE_CLASS1_WINDOWS` overrides the CSV path and `BLACKICE_INPUT` the processed model input |
| `hpc/job_alpha_diversity_window_lmm.pbs` | Generate the CSV before submission; `BLACKICE_WINDOW_DEFINITION_PATH` and `BLACKICE_ALPHA_METADATA_PATH` override its inputs |

The PBS jobs do not invoke the generator. Plot/table readers use fixed default
paths; supplying a model override does not redirect them. The alpha and reporting
scripts need additional CSV columns, so a minimal table of window labels and
timepoints is not a replacement for the regenerated full table. Retain script
behavior and scientific definitions; no synthetic or sanitized replacement is
provided by this release decision.

## 2. Diversity branches

```bash
# Raw-input alpha metrics, tool comparison and publication plots
make alpha-diversity
Rscript scripts/alpha_diversity/06_plot_alpha_diversity_heatmap.R
Rscript scripts/alpha_diversity_lmm/01_prepare_alpha_lmm_input.R
Rscript scripts/alpha_diversity_lmm/02_fit_alpha_diversity_lmm.R
Rscript scripts/alpha_diversity_lmm/03_plot_alpha_diversity_lmm_trends.R
Rscript scripts/alpha_diversity_lmm/04_alpha_lmm_diagnostics.R
Rscript scripts/alpha_diversity/08_create_alpha_diversity_publication_summary_table.R
# Requires the expedition-window definitions as well:
Rscript scripts/alpha_diversity_window_lmm/fit_plot_alpha_window_lmm.R

# Separate raw-input beta branch: stages 01 through 08
bash scripts/run_beta_diversity_local.sh --config config/beta_diversity.yml
```

These are generation commands, not permission to replace release evidence.
Alpha and the production beta input branch retain raw-input timepoints, including
Ta/Tb. The owner confirmed both as supporting cohorts for the current manuscript.
The beta run and figure/table mapping still require sign-off. Primary per-MAG
inference uses the configured cleaned-data selection.

**Alpha reproduction blocker:** the selected committed outputs use the
normal-approximation fallback, while an environment with lmerTest chooses
Satterthwaite inference. Recover the historical environment and approve the
comparison before running this sequence to replace manuscript artifacts.
Do not edit p-values manually or silently change the method.

Beta uses subject-restricted permutations, primary Bray–Curtis, binary Jaccard
and Aitchison distances plus sensitivities. With no tree, UniFrac is explicitly
skipped. Inspect dispersion and longitudinal diagnostics alongside PERMANOVA.

Once approved source results exist:

```bash
Rscript scripts/manuscript_diversity_per_mag/01_build_main_figures_tables.R
```

This builds main figures/tables and a supplementary evidence manifest under
`results/manuscript_diversity_per_mag/`. It also refits up to three qualifying MAGs per direction
for display trends; it does not replace the primary screening table. Because
current local alpha outputs differ from the selected committed evidence, do not
assume the existing synthesis bundle matches the release candidate.

## 3. Random Forest and subject-aware validation

```bash
Rscript scripts/run_random_forest_model.R
Rscript analysis/random_forest_subject_validation/run_random_forest_subject_validation.R development
```

The first command is exploratory sample-level classification. The second runs
subject-aware baselines, LOSO, repeated grouped CV, feature/balancing/outcome
sensitivities and validation. Production uses more trees and is not independent
replication. Preprocessing is fitted within training folds; evaluation subjects
must remain disjoint. No metadata or subject identifier should become a MAG
predictor.

Select an authorized completed run explicitly; no saved research run is supplied.

```bash
Rscript scripts/random_forest_subject_validation/01_build_manuscript_outputs.R \
  --run-dir results/random_forest_subject_validation/RUN_ID
Rscript scripts/plot_random_forest_feature_stability.R \
  --feature-results-dir results/random_forest_subject_validation/RUN_ID/features \
  --taxonomy-file data/raw/GDTB.csv \
  --output-dir results/random_forest_subject_validation/RUN_ID/figures/feature_stability
```

Existing outputs can trigger overwrite protection. Review prior bundles before
using any force option. Native importance is not causal or independent biomarker
validation. The model's altitude classes do not identify altitude effects apart
from expedition time and other confounding.

## 4. NetMoss and companions

Install pinned NetMoss dependencies as described in [environment evidence](public-environment.md).
The separate historical config is recovered from `fa1f2f3`, matching the
preserved manifest's T1/T2/Ta versus T8/T9/T10 scenario. It retains original
settings, including the original output default. For a future reviewed rerun,
use an isolated output path:

```bash
Rscript scripts/netmoss/run_netmoss_pipeline.R \
  --config config/netmoss_historical_publication.yml \
  --output-dir results/netmoss_historical_reproduction \
  --acknowledge-exploratory
```

The historical configuration is distinct from the current Ta/Tb-excluded
configuration. Supply authorized completed outputs separately before rendering:

```bash
bash scripts/run_netmoss_manuscript_synthesis.sh
```

The current `netmoss_analysis.yml` is a DIFFERENT Ta/Tb-excluded scenario and
must not regenerate the historical publication path accidentally. The altitude
contrast, altitude stability, and matched-timepoint companion have separate
configs/output roots. Their respective entrypoints are
`run_netmoss_pipeline.R`, `run_altitude_stability.R`, and
`run_matched_timepoint_companion.R` under `scripts/netmoss/`; use their documented
config arguments. The companion evaluates paired abundance, not network scores.
All are exploratory. Do not infer direct microbial interactions or validated
rewiring from thresholded correlations.

## 5. Functional workflows and synthesis

These stages require the upstream cluster/annotation files and substantial
memory/time. Full processing is normally run through the optional HPC wrappers.

```bash
bash scripts/run_mag_functional_annotation.sh --config config/mag_functional_annotation.yml
bash scripts/run_functional_community_dynamics.sh --config config/functional_community_dynamics.yml
bash scripts/run_manuscript_functional_synthesis.sh --config config/manuscript_functional_synthesis.yml
```

The first runs stages 0–10, including QC, catalogue and longitudinal analyses,
same-data candidate checks and annotation audit. The second consumes completed
functional profiles and matches taxonomic and functional cohorts. The third
renders existing result tables; it is not a new discovery analysis.

The original gut abundances remain unnormalized in primary functional potential;
existing normalized sensitivity outputs are labeled separately. BH families,
exact MAG joins, candidate-selection caveats and the distinction between genomic
potential and biological activity remain unchanged. Read the workflow-specific
README files for detailed methods and stage restart controls.


## Validation and external outputs

Input files and prior-stage outputs named here are obtained separately or generated
locally; they are not bundled. Run in an isolated checkout with authorized inputs,
record input/code/configuration hashes and package versions, inspect status and
completion markers, and retain failed/singular model flags. Syntax checks and
synthetic tests do not establish a complete scientific workflow reproduction.

See [public release gates](public-release-checklist.md). No full scientific rerun
was performed during this documentation review. The legacy `make report` target
is omitted from the public Makefile by owner approval. No report template is needed by the documented analytical workflows.

Public rendering entrypoints, external input overrides and synthetic checks are
documented in the [rendering guide](public-rendering.md).
