# HPC Execution Recording

This repository records standardized reproducibility metadata for HPC-oriented analysis runs. The record is intended to support thesis or manuscript supplementary tables by preserving the runtime environment, Git state, input files, result files, and FDR/testing scope used for each analysis.

## Workflows

Execution records are written by:

- `scripts/hpc_run_per_mag_lmertest.R`
- `scripts/hpc_run_class1_window_lmertest.R`
- `scripts/evaluate_expedition_windows.R`

This is metadata/logging only. It does not change model formulas, pseudocount rules, filtering thresholds, FDR correction logic, candidate definitions, or output table schemas.

## Output Location

Each workflow writes:

- `results/<analysis_dir>/reproducibility/hpc_execution_record.csv`
- `results/<analysis_dir>/reproducibility/hpc_execution_record.md`

The recorder captures scheduler metadata when available, including OpenPBS job IDs, hostname, user, working directory, Git branch/commit/status, Conda environment, loaded modules, R version, key R package versions, input/config files, result files, test counts, FDR family, and notes.

## Combining Records

After one or more workflows have run, combine available records with:

```bash
Rscript scripts/summarize_hpc_execution_records.R
```

This writes:

- `results/hpc_execution_records_combined.csv`
- `results/hpc_execution_records_combined.md`

The combined table is formatted with stable columns suitable for supplementary material review.

## Reproducibility Value

The execution record links each results directory to the exact script, inputs, configuration files, Git commit, runtime environment, and testing/FDR scope. This helps distinguish scientific results from machine-specific execution details and provides an audit trail for rerunning or reviewing analyses.

## Safety

The recorder is fail-safe. If execution recording fails, the analysis prints a warning and continues. Missing environment variables, unavailable Git metadata, or absent R packages are recorded as blank/`NA` rather than causing workflow failure.
