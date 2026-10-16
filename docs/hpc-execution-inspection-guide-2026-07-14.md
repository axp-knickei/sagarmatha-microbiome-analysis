# HPC Execution Inspection Guide

**Recorded:** 2026-07-14 20:43:53 WIB (UTC+07:00)
**Purpose:** Inspect which repository version, scripts, and scheduler jobs were executed on an HPC system after cloning `sagarmatha-microbiome-analysis`.

## 1. Check the Current Repository Version

From the cloned repository:

```bash
cd /path/to/sagarmatha-microbiome-analysis
git log -10 --oneline --decorate
git status
```

This identifies the checked-out commit and shows whether the working tree contains local changes.

## 2. Inspect Shell Command History

Search the current account's shell history for likely analysis and scheduler commands:

```bash
history | grep -E 'qsub|sbatch|Rscript|run_.*\.sh|make '
```

To display the 50 most recent history entries:

```bash
history 50
```

Shell history is supporting evidence only. It may be incomplete if commands were run from another shell, account, session, or non-interactive job.

## 3. Check Currently Running PBS Jobs

```bash
qstat -u "$USER"
```

For detailed information about one job:

```bash
qstat -f <JOB_ID>
```

Replace `<JOB_ID>` with the scheduler job identifier.

## 4. Check Completed Scheduler Jobs

Available commands depend on the cluster and scheduler configuration. For OpenPBS/PBS systems, try:

```bash
qstat -x -u "$USER"
```

or:

```bash
tracejob <JOB_ID>
```

Some clusters provide a compatible accounting command such as:

```bash
qacct -j <JOB_ID>
```

`qacct` is commonly associated with Grid Engine rather than OpenPBS, so consult the cluster documentation if it is unavailable.

## 5. Find Scheduler and Workflow Logs

From the repository root, list common log files:

```bash
find . -type f \( -name '*.log' -o -name '*.out' -o -name '*.err' \) -print
```

Search those files for script execution, scheduler metadata, completion messages, and errors:

```bash
rg -n 'qsub|PBS|JOB_ID|Git commit|started|completed|finished|ERROR|FAILED' \
  --glob '*.log' --glob '*.out' --glob '*.err' .
```

## 6. Inspect Main Workflow Logs

When these files exist, inspect their most recent messages:

```bash
tail -100 update_glmm_env.log
tail -100 blackice_beta_diversity.pbs.log
find logs -maxdepth 2 -type f -print
```

Check the environment update log for successful validation or failure messages:

```bash
grep -nE 'COMPLETED_ENVIRONMENT_VALIDATION|Environment update completed|ERROR|failed' \
  update_glmm_env.log
```

## 7. Check Generated Outputs and Timestamps

List the 50 most recently modified processed-data and result files:

```bash
find results data/processed -type f -printf '%TY-%Tm-%Td %TH:%TM %p\n' \
  2>/dev/null | sort -r | head -50
```

For beta-diversity outputs specifically:

```bash
find results/beta_diversity results/beta_diversity_test \
  -type f -printf '%TY-%Tm-%Td %TH:%TM %p\n' \
  2>/dev/null | sort -r | head -50
```

File modification times can show when outputs were written, but they do not independently prove which code version produced them.

## 8. Check for Recorded Git and PBS Metadata

Search text-based records for commit and scheduler metadata:

```bash
rg -n 'git rev-parse|Git commit|PBS_JOBID|PBS_JOBNAME' \
  logs results --glob '*.log' --glob '*.txt' --glob '*.md' 2>/dev/null
```

The repository's automatic execution-recording approach is documented separately in [`hpc_execution_recording.md`](hpc_execution_recording.md).

## Combined Inspection

The following block provides a compact overview:

```bash
echo '=== CURRENT COMMIT ==='
git log -1 --oneline

echo '=== PBS JOBS ==='
qstat -u "$USER"

echo '=== RECENT HISTORY ==='
history | grep -E 'qsub|sbatch|Rscript|run_.*\.sh|make ' | tail -30

echo '=== RECENT LOGS ==='
find . -type f \( -name '*.log' -o -name '*.out' -o -name '*.err' \) \
  -printf '%TY-%Tm-%Td %TH:%TM %p\n' 2>/dev/null | sort -r | head -30
```

## Evidence Priority

For the most reliable reconstruction, compare:

1. Scheduler state and accounting information from `qstat` or `tracejob`.
2. Scheduler and workflow logs.
3. Generated output timestamps and recorded reproducibility metadata.
4. Shell history as supplementary evidence.

No single source is guaranteed to be complete. Confidence is highest when the scheduler records, logs, Git commit, and output timestamps agree.
