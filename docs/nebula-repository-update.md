# Updating the Repository on Nebula

Use this runbook to safely update an existing public analytical checkout on the
Nebula OpenPBS cluster. Run the commands from an SSH session on Nebula. Keep
credentials, private paths, raw data, results, and scheduler output out of Git.

## 1. Enter the repository and inspect it

```bash
cd /path/to/sagarmatha-microbiome-analysis
git status
git branch --show-current
```

Do not discard a modified working tree. HPC-only scripts, environment changes,
and generated exports may need to be preserved before updating.

## 2. Fetch and review the incoming commits

```bash
git fetch origin
git status
git log --oneline HEAD..origin/$(git branch --show-current)
```

An empty final command means the current branch has no incoming commits. Confirm
that the current branch is the branch that should be updated before continuing.

## 3. Protect local HPC files

`git stash -u` includes untracked files. If `exports/` exists, move it outside
the checkout first so large archives are not written into the Git stash:

```bash
mv exports ../sagarmatha-microbiome-analysis-exports-backup
git stash push -u -m "HPC local files before repository update"
```

Skip the `mv` command when `exports/` does not exist. If the backup destination
already exists, choose another name instead of overwriting it.

## 4. Fast-forward the current branch

```bash
git pull --ff-only
git stash pop
```

Restore exports if they were moved:

```bash
mv ../sagarmatha-microbiome-analysis-exports-backup exports
```

`git stash pop` may automatically merge `environment.yml`. Verify that it did
not leave conflict markers or whitespace errors:

```bash
git diff -- environment.yml
git diff --check
grep -nE '^(<<<<<<<|=======|>>>>>>>)' environment.yml
```

The last two commands should print nothing. If conflict markers are reported,
resolve the file before updating the environment or submitting analyses.

## 5. Verify the repository update

```bash
git status
git log -1 --oneline
```

Local HPC files may still appear as modified or untracked. This is expected when
they are intentionally machine-specific; do not commit credentials, raw data,
results, exports, or scheduler logs.

## 6. Update the existing Conda environment with OpenPBS

The tracked update job loads Nebula's `intelpython/2022.2.0` module, updates the
existing `glmm-hpc` environment from `environment.yml`, and runs the beta-
diversity environment validation. Submit it from the repository root:

```bash
qsub hpc/update_glmm_env.pbs
```

To override the defaults when necessary:

```bash
qsub -v CONDA_MODULE_NAME=intelpython/2022.2.0,CONDA_ENV_NAME=glmm-hpc hpc/update_glmm_env.pbs
```

Save the returned job ID and monitor it:

```bash
qstat <JOB_ID>
qstat -u "$USER"
```

When the job leaves the queue, inspect its combined log:

```bash
less update_glmm_env.log
```

A successful job ends with both:

```text
COMPLETED_ENVIRONMENT_VALIDATION
=== Environment update completed: ... ===
```

Do not use `hpc/create_conda_env.pbs` to refresh an existing environment. An
environment-creation job is only appropriate when `glmm-hpc` does not exist.

## 7. Submit analysis jobs

Only submit analysis after the environment update succeeds. Nebula uses
OpenPBS and `qsub`. Submit the beta-diversity workflow with
`qsub hpc/beta_diversity_full.pbs`. The tracked `*.slurm` beta-diversity files
remain available for Slurm clusters and must not be submitted on Nebula.

## Recovery notes

- If `git pull --ff-only` refuses to update, do not use `git reset --hard`.
  Recheck the current branch, `git status`, and the stash before proceeding.
- If `git stash pop` reports conflicts, keep both the incoming dependency
  changes and required HPC-specific settings while resolving them.
- Git updates are intentionally performed interactively before job submission;
  PBS jobs should not pull changing source code while an analysis is starting.
