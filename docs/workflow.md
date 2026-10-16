# Workflow Guide

## Operating model

GitHub is the source of truth for tracked code, configuration, documentation,
and selected analysis evidence. Raw data remain outside version control. Local
or HPC runs create model-ready data and workflow-specific outputs under
`data/processed/` and `results/`.

```text
local development and smoke tests
              |
              | reviewed commits
              v
           GitHub
              |
              | intentional pull/checkout
              v
        scheduler-managed HPC
              |
              | dry run -> full run -> recorded outputs
              v
        reviewed result evidence
```

Do not run `git pull` from inside a scheduled job. Select and record the code
revision before submission.

## 1. Prepare the environment

From a clean clone:

```bash
conda env create -f environment.yml
conda activate glmm-hpc
```

On Nebula, follow [Nebula Repository Update](nebula-repository-update.md) for
the site module and environment-update procedure. The Conda environment is the
primary cross-language runtime; see [config/README.md](../config/README.md) for
the relationship to `renv.lock`.

## 2. Provide and validate inputs

Place `sample_info.csv`, `gut.csv`, and any workflow-specific taxonomy input in
`data/raw/`, or use the documented external-path override. Never commit raw
data or private storage paths.

For the shared model-ready table:

```bash
Rscript scripts/01_data_prep.R
Rscript scripts/validate_processed_data.R
Rscript scripts/validate_formula.R
```

Review row counts, sample-ID matching, missingness, model variables, and the
generated validation reports before submitting model jobs.

## 3. Choose the workflow

| Workflow | Local entrypoint | HPC entrypoint | Status |
|---|---|---|---|
| Per-MAG longitudinal LMM | `scripts/hpc_run_per_mag_lmertest.R` | `hpc/job_per_mag_lmertest.pbs` | Primary inference workflow |
| Class 1 predefined-window LMM | `scripts/hpc_run_class1_window_lmertest.R` | `hpc/job_class1_window_lmertest.pbs` | Secondary validation |
| Alpha diversity | `make alpha-diversity` or individual staged scripts | `hpc/job_alpha_diversity.pbs` and related jobs | Supporting analysis |
| Beta diversity | `scripts/run_beta_diversity_local.sh` | `hpc/beta_diversity_full.pbs` | Supporting pipeline |
| Random Forest | `scripts/run_random_forest_model.R` | `hpc/job_random_forest.pbs` | Exploratory prediction |
| Subject-aware Random Forest validation | `analysis/random_forest_subject_validation/run_random_forest_subject_validation.R` | `hpc/job_random_forest_subject_validation.pbs` | Validation workflow |
| NetMoss | `scripts/netmoss/run_netmoss_pipeline.R` | `hpc/job_netmoss_analysis.pbs` | Exploratory supporting analysis |

See [scripts/README.md](../scripts/README.md),
[analysis/README.md](../analysis/README.md), and
[hpc/README.md](../hpc/README.md) for stage and wrapper details.

The original `make all` and `hpc/job_template.pbs` path remains available for
legacy single-model checks; it is not the current per-MAG inference path.

## 4. Test before scaling

Use the smallest workflow-specific test mode before a production run:

```bash
# Primary per-MAG dry run
qsub -v RUN_MODE=dryrun hpc/job_per_mag_lmertest.pbs

# Beta-diversity smoke test
make beta-test

# Subject-aware Random Forest smoke run
Rscript analysis/random_forest_subject_validation/run_random_forest_subject_validation.R smoke

# NetMoss wrapper tests
bash hpc/tests/test_netmoss_pbs.sh
```

Inspect logs, validation tables, session information, and scheduler resource
usage before increasing scope.

## 5. Run and monitor on HPC

Example primary submissions:

```bash
# Full per-MAG run
qsub hpc/job_per_mag_lmertest.pbs

# Chunked per-MAG run
qsub -J 1-10 -v CHUNK_COUNT=10 hpc/job_per_mag_lmertest.pbs
```

Use `qstat` to monitor OpenPBS jobs. Do not execute heavy workflows directly on
the login node. Site-specific access details belong in ignored local notes,
not in tracked documentation.

## 6. Review and synchronize results

Confirm that expected tables, logs, diagnostics, session information, and
execution records exist. For chunked per-MAG analysis, combine chunks and
recompute BH-FDR across the complete testing family:

```bash
Rscript scripts/combine_hpc_per_mag_chunks.R
```

Follow [Nebula Results Download](nebula-results-download.md) for the canonical
dry-run-first `rsync` procedure. Do not use `--delete` unless deletions are
explicitly approved.

## 7. Preserve reproducibility

- Record the Git branch and commit used for the run.
- Keep configuration, input checksums, scheduler context, session information,
  and primary outputs together when supported.
- Review `git status` before staging anything; some selected result evidence is
  tracked and may be regenerated locally.
- Separate executable changes from documentation and generated-result changes.
- Treat candidate signals as exploratory unless they satisfy the documented
  multiple-testing rule and scientific review requirements.

## Troubleshooting map

| Situation | First reference |
|---|---|
| Environment or package failure | [HPC Guide](../hpc/README.md) and [Environment Guide](public-environment.md) |
| Input or row-loss concern | [Local Data Preflight](local-data-preflight.md) and [Data Adaptation Checklist](release-data-requirements.md) |
| Scheduler or run-record concern | [HPC Execution Inspection](hpc-execution-inspection-guide-2026-07-14.md) |
| Unsure which script is canonical | [Script Guide](../scripts/README.md) |
| Unsure whether an output should be committed | [Results Policy](../results/README.md) |
