# HPC Workflow Guide

This directory contains environment setup, OpenPBS jobs, Slurm launchers, and
scheduler-specific validation helpers. Submit jobs from the repository root.
Canonical checkout-relative input and result paths are indexed in
[`docs/data-locations.md`](../docs/data-locations.md).

## OpenPBS jobs

| Job | Workflow | Status |
|---|---|---|
| `job_per_mag_lmertest.pbs` | Primary per-MAG LMM, including dry-run and chunked modes | Canonical primary job |
| `job_class1_window_lmertest.pbs` | Class 1 predefined-window LMM | Active secondary job |
| `job_alpha_diversity.pbs` | Alpha-diversity calculation and plots | Active supporting job |
| `job_alpha_diversity_heatmap.pbs` | Alpha-diversity heatmaps | Active supporting job |
| `job_alpha_diversity_lmm.pbs` | Alpha-diversity LMM workflow | Active supporting job |
| `job_alpha_diversity_publication_summary_table.pbs` | Alpha-diversity publication summary | Active supporting job |
| `job_alpha_diversity_window_lmm.pbs` | Alpha-diversity biological-window LMM | Active supporting job |
| `beta_diversity_full.pbs` | Full or smoke-test beta-diversity pipeline | Active supporting job |
| `job_random_forest.pbs` | Exploratory Random Forest | Active exploratory job |
| `job_random_forest_publication.pbs` | Random Forest publication outputs | Post-processing job |
| `job_random_forest_subject_validation.pbs` | Subject-aware Random Forest validation | Validation job |
| `job_random_forest_subject_validation_manuscript.pbs` | Reader-first subject-validation main figure and tables | Manuscript post-processing job; does not refit the model |
| `job_plot_random_forest_feature_stability.pbs` | Four-panel MAG eligibility, native-importance and rank-stability figure | Publication post-processing job; does not refit the model |
| `job_netmoss_readiness.pbs` | NetMoss readiness assessment | Diagnostic job |
| `job_netmoss_smoke_test.pbs` | NetMoss wrapper smoke tests | Test job |
| `job_netmoss_analysis.pbs` | Exploratory NetMoss analysis | Explicitly exploratory |
| `job_netmoss_altitude_analysis.pbs` | Low/initial versus high-altitude NetMoss comparison | Separate explicitly exploratory scenario |
| `job_netmoss_altitude_stability.pbs` | Matched-participant genus/family NetMoss stability analysis | Internal reproducibility only; explicitly exploratory |
| `job_netmoss_altitude_stability_smoke_test.pbs` | NetMoss altitude stability smoke test | Test-data-only job |
| `job_netmoss_publication.pbs` | NetMoss publication outputs | Post-processing job |
| `job_netmoss_manuscript_synthesis.pbs` | Reader-first NetMoss main figure and tables | Post-processing only; does not rerun NetMoss2 |
| `job_netmoss_matched_timepoint_companion.pbs` | Advanced all-pairs matched-participant genus-abundance companion | Exploratory; not NetMoss network inference |
| `job_netmoss_matched_timepoint_smoke_test.pbs` | Matched-timepoint companion smoke test | Test-data-only job |
| `job_mag_functional_annotation.pbs` | Full MAG protein-cluster and functional annotation workflow | Active staged OpenPBS workflow |
| `job_functional_community_dynamics.pbs` | Functional profile restructuring, taxonomic comparison, and exploratory recovery | Active staged OpenPBS workflow |
| `job_manuscript_functional_synthesis.pbs` | Reader-first main figure and tables from the two completed functional workflows | Manuscript post-processing job |
| `job_manuscript_diversity_per_mag.pbs` | Reader-first alpha-, beta-, and primary per-MAG main figures and tables | Manuscript post-processing job |
| `job_template.pbs` | Original three-stage single-model pipeline | Legacy only |

`update_glmm_env.pbs` and `update_netmoss_env.pbs` update or validate the
configured Conda environment on a compute node. `setup_hpc.sh` performs initial
checkout setup. Review site modules and environment names before submission.

## Window-job input prerequisite for v1.0.0

`job_class1_window_lmertest.pbs` and `job_alpha_diversity_window_lmm.pbs` require
a locally regenerated `proposed_biological_windows.csv`. The original CSV and
research metadata are excluded from the public code archive. Neither job creates
the table. Follow [window CSV regeneration](../docs/public-workflow.md#window-csv-regeneration)
with separately obtained authorized metadata matching the original 86-sample
processed selection before submission. The Class 1 job accepts
`BLACKICE_CLASS1_WINDOWS`; the alpha-window job accepts
`BLACKICE_WINDOW_DEFINITION_PATH`. Default result paths do not imply bundled data.

## Slurm launchers

The beta-diversity workflow also includes `beta_diversity_full.slurm`,
`beta_diversity_metric_array.slurm`, and `beta_diversity_report.slurm`. They are
alternatives for a Slurm site, not interchangeable directives for OpenPBS.

## Per-MAG helpers

- `run_hpc_lmertest_dryrun.sh` runs the per-MAG implementation with a small MAG
  limit.
- `run_hpc_lmertest_full.sh` runs the full per-MAG implementation.

The OpenPBS job remains the documented Nebula submission path because it also
records scheduler context and supports chunk arrays.

## Local notes and logs

Files matching `hpc/*.local.*` are intentionally ignored and may contain
machine-specific paths. They must not be copied into tracked documentation.
Scheduler stdout, stderr, downloaded logs, and result archives are operational
artifacts rather than source code.

The matched-timepoint companion uses 8 CPU workers and 16 GB RAM on `workq`
for 10,000-replicate participant bootstraps across the 45 approved timepoint
pairs after excluding Ta and Tb. It
does not request GPUs because the observed Nebula PBS inventory reports
`ngpus 0/0`, and it does not reserve all 96 cores because this dataset would
not use an entire node efficiently. See
[`docs/netmoss-methods.md`](../docs/netmoss-methods.md).

The MAG functional-annotation workflow includes follow-up Stages 6–10 for the
same-data validation of `K04028` and EC `1.2.1.87`, captured-fraction
sensitivity, functional co-carriage, unknown-singleton evidence readiness, and
shared-versus-candidate-specific carrier decomposition. Stage 9 identifies the
BH-supported shared-component timepoint contrasts, checks participant
influence, and attributes paired changes exactly to MAGs and genera.
Stage 10 then audits the representative-protein evidence, exact carrier-MAG
concordance, reviewed candidate identities, and the boundary between MAG-level
co-carriage and unsupported pathway claims. It reads tracked curation metadata
offline and does not require internet access on a compute node.
After completed Stages 0–5 have been synchronized, submit only the follow-up
with `MAGFUNC_FROM_STAGE=6,MAGFUNC_TO_STAGE=10`. If Stage 9 is already complete,
submit only Stage 10 with `MAGFUNC_FROM_STAGE=10,MAGFUNC_TO_STAGE=10`. The detailed contract and
success checks are in
[`docs/functional-methods.md`](../docs/functional-methods.md).

For safe repository updates and result synchronization, follow
[`docs/nebula-repository-update.md`](../docs/nebula-repository-update.md) and
[`docs/nebula-results-download.md`](../docs/nebula-results-download.md).

The functional-community-dynamics workflow is separate from the completed
two-candidate pipeline. It reads the existing predicted-potential files, but
recomputes taxonomic and functional distances on an identical T1–T10 subset.
Submit it with `qsub hpc/job_functional_community_dynamics.pbs`; bounded reruns
use `FUNCDYN_FROM_STAGE` and `FUNCDYN_TO_STAGE` from 0 through 4.

After both functional workflows are complete, the manuscript synthesis can be
submitted with `qsub hpc/job_manuscript_functional_synthesis.pbs`. Optional
source-result overrides are `MANSYN_MAGFUNC_RESULTS_DIR` and
`MANSYN_FUNCDYN_RESULTS_DIR`. The job uses one CPU and 16 GB because it reads
completed compact result tables and renders figures; it does not refit the
scientific models.
