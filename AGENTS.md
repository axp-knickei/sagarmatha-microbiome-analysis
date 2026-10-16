# Repository Guidelines

## Project Structure and Workflow Status

This repository is a reproducible longitudinal microbiome analysis project for
local machines and scheduler-managed HPC systems. It began as a configurable
three-stage GLMM template, but the active primary analysis is now a per-MAG
longitudinal LMM. Do not treat the legacy template as the only or current
workflow.

Use the directory indexes before editing:

- `scripts/README.md` maps canonical, supporting, exploratory, and legacy
  scripts.
- `analysis/README.md` describes the staged beta-diversity and subject-aware
  Random Forest validation pipelines.
- `config/README.md` identifies the owner and authority of each YAML file.
- `hpc/README.md` maps OpenPBS and Slurm wrappers to workflows.
- `docs/README.md` indexes operational and scientific documentation.
- `results/README.md` explains generated outputs versus selected tracked
  evidence.

Core inputs belong in `data/raw/` or an explicitly documented external raw-data
directory. Generated model-ready data belongs in `data/processed/`. Never
modify raw inputs as part of code maintenance.

## Canonical Entry Points

- Data preparation: `scripts/01_data_prep.R`
- Primary per-MAG LMM: `scripts/hpc_run_per_mag_lmertest.R`
- Primary OpenPBS wrapper: `hpc/job_per_mag_lmertest.pbs`
- Class 1 window LMM: `scripts/hpc_run_class1_window_lmertest.R`
- Alpha diversity: `scripts/alpha_diversity/`
- Beta diversity: `scripts/run_beta_diversity_local.sh`
- Exploratory Random Forest: `scripts/run_random_forest_model.R`
- Subject-aware Random Forest validation:
  `analysis/random_forest_subject_validation/run_random_forest_subject_validation.R`
- NetMoss: `scripts/netmoss/`

`scripts/02_run_glmm.R`, `scripts/03_results.R`, and `hpc/job_template.pbs`
belong to the legacy single-model path. Preserve them unless a task explicitly
authorizes compatibility changes or removal.

## Build, Test, and Development Commands

- `make setup`: create the `glmm-hpc` Conda environment.
- `make data`: prepare `data/processed/cleaned_data.rds`.
- `make analyze`: run the legacy configured single-model fit.
- `make results`: create legacy tables and diagnostics.
- `make all`: run the complete legacy three-stage pipeline.
- `make alpha-diversity`: run the alpha-diversity stages listed in the Makefile.
- `make beta-env`: validate beta-diversity runtime dependencies.
- `make beta-test`: create and run the beta-diversity smoke-test workflow.
- `make beta-diversity`: run the configured production beta-diversity workflow.
- `qsub -v RUN_MODE=dryrun hpc/job_per_mag_lmertest.pbs`: submit a primary
  per-MAG dry run on OpenPBS.

The public Makefile omits the obsolete `make report` target by owner approval.
No report template is supplied; this omission does not change analytical recipes.
The private original Makefile remains unchanged; `Makefile.public` is the approved
assembly source for the public root `Makefile`.

## Coding and Naming Conventions

Run scripts from the repository root. Preserve sequential `NN_description.R`
names within staged pipelines. Use `here()` or validated repository-relative
paths where the workflow already follows that convention. Use two-space YAML
indentation and conventional R spacing around assignments and operators.

New user-facing scripts should document purpose, expected inputs, generated
outputs, dependencies, invocation, and relationship to adjacent stages.
Comments should explain analysis intent and workflow constraints rather than
restate individual expressions.

Do not rename or move files merely to normalize naming. Exact paths are used by
shell runners, scheduler wrappers, `source()` calls, documentation, and
reproducibility records.

## Testing Guidelines

There are workflow-local tests under `analysis/beta_diversity/tests/` and
`scripts/netmoss/tests/`, plus repository-level Random Forest tests under
`tests/testthat/`. Validate the smallest relevant workflow first.

For documentation-only or comment-only work, run syntax parsing, shell syntax
checks, Markdown link checks, `git diff --check`, and Make dry runs. Do not run
scientific pipelines unnecessarily because they may regenerate tracked result
evidence. For executable changes, run the applicable smoke test and then the
full workflow required for analysis handoff.

## Result and Data Safety

The repository contains both ignored generated outputs and selected tracked
result evidence. Always inspect `git status` before and after running a script.
Do not overwrite, delete, untrack, or reformat result artifacts during
repository maintenance. Do not commit raw data, processed RDS files, model
objects, scheduler output, credentials, or private machine paths unless an
explicit evidence policy says otherwise.

Some existing scripts may regenerate tracked evidence. Preserve unrelated
working-tree changes and stage only files that belong to the approved task.

## Commit and Pull Request Guidelines

Use concise imperative summaries such as `Document reproducible environment
setup` and `Harden HPC GLMM workflow`. Keep navigation documentation, script
comment changes, and executable behavior changes in separate commits where
practical. Pull requests should state workflow impact, list validation commands,
mention configuration changes, and include figure or table excerpts when
outputs change.

## Security and Machine-Specific Configuration

Do not commit credentials, `.env` files, SSH keys, raw datasets, scheduler
output, private HPC paths, or ignored `hpc/*.local.*` notes. Keep reusable
defaults in tracked templates and machine-specific overrides in ignored files.

## Nebula Result Synchronization

When asked how to download, copy, or synchronize data produced by HPC jobs,
follow `docs/nebula-results-download.md` as the canonical response. Obtain the
SSH alias and current HPC checkout path from the ignored
`hpc/hpc_access.local.md` note, but never reproduce credentials, private keys,
passphrases, or other secrets. Always present the dry run before the real
`rsync` commands, preserve replaced local files with `--backup`, and do not use
`--delete` unless the user explicitly requests and approves deletions.

## Raw MAG Data Interpretation

Before interpreting `data/raw/gut.csv` MAG identifiers or joining taxonomy,
read `docs/data/input-interpretation.md`. Treat each full MAG column name as the stable key
and join it exactly to `GDTB.csv` column `Genome`. The `B`, `C`, and `D` letters
record the specimen type associated with the retained representative bin; they
do not override the later gut-specific MAG selection represented by membership
in `gut.csv`. Embedded dates and altitude-like values are construction
provenance and do not consistently equal `sample_info.csv`; observation
metadata must come from the row-level `sampleID` join. Do not infer undocumented
meanings for `EVR`, `QML`, `strict`, `permissive`, or `orig`.

## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

When the user types `/graphify`, use the installed graphify skill or instructions before doing anything else.

Rules:
- For codebase questions, first run `graphify query "<question>"` when graphify-out/graph.json exists. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- Dirty graphify-out/ files are expected after hooks or incremental updates; dirty graph files are not a reason to skip graphify. Only skip graphify if the task is about stale or incorrect graph output, or the user explicitly says not to use it.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).
