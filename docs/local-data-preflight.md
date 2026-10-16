# Local Data Preflight Before HPC Analysis

Purpose: make sure the raw dataset, data wrangling, model formula, and analysis scripts work on a local computer before running the larger job on HPC.

Do not start with `scripts/02_run_glmm.R`. Start with the raw data structure, then make `scripts/01_data_prep.R` convert the raw file into the exact model-ready table needed by the statistical model.

## 1. Inspect the Raw Dataset

First identify whether the data is CSV, TSV, or another text format.

For CSV:

```bash
head -5 data/raw/your_file.csv
```

For TSV:

```bash
head -5 data/raw/your_file.tsv
```

Then inspect the file in R.

For CSV:

```r
readr::read_csv("data/raw/your_file.csv", n_max = 5)
```

For TSV:

```r
readr::read_tsv("data/raw/your_file.tsv", n_max = 5)
```

Answer these questions before changing the model:

- What is the response variable?
- What is the subject, sample, or replicate ID column?
- What is the time or sampling point column?
- What is the group, condition, treatment, or comparison column?
- Are values numeric, categorical, dates, or IDs?
- How are missing values encoded: blank, `NA`, `N/A`, `-`, `NULL`, or something else?
- Is the data already in long format, or does it need reshaping from wide to long?

## 2. Define the Model-Ready Schema

The default model in `config/config.yml` is:

```yaml
model:
  formula: "Y ~ Time * Group + (1 + Time | Subject)"
```

This means the processed data must contain at least these columns:

```text
Y
Time
Group
Subject
```

Those names are template names. The real raw data may use different names, for example:

```text
relative_abundance
sampling_day
treatment
sample_id
```

Use `scripts/01_data_prep.R` to rename and transform raw columns into stable model columns. Keep `scripts/02_run_glmm.R` focused on fitting the model.

## 3. Make a Data Dictionary

Create a small mapping before editing scripts.

```text
Raw column name       Meaning                  Model column
sample_id             Subject/sample ID         Subject
day                   Time point                Time
treatment             Experimental group        Group
relative_abundance    Response variable         Y
```

This mapping is the contract between the dataset and the analysis scripts.

## 4. Put Wrangling in `01_data_prep.R`

All data cleanup should happen in `scripts/01_data_prep.R`.

That script should:

- Read the raw file with the correct delimiter.
- Rename raw columns to model columns.
- Convert variable types intentionally.
- Handle missing values intentionally.
- Filter invalid rows.
- Save `data/processed/cleaned_data.rds`.

Examples of type decisions:

- `Y` should usually be numeric.
- `Subject` should usually be an ID or factor.
- `Group` should usually be categorical.
- `Time` should be numeric if modeling a trend.
- `Time` should be a factor if comparing time points as categories.

## 5. Test Locally With a Small Subset

Before using HPC, test the pipeline locally with a small file.

For CSV:

```bash
mkdir -p data/raw
head -100 data/raw/your_file.csv > data/raw/dataset.csv
```

Then run:

```bash
conda activate glmm-hpc
Rscript scripts/01_data_prep.R
Rscript scripts/02_run_glmm.R
Rscript scripts/03_results.R
```

If the model is too complex for a small test file, temporarily simplify the formula in `config/config.yml`:

```yaml
formula: "Y ~ Time * Group + (1 | Subject)"
```

Restore or expand the full model after the small test passes.

## 6. Inspect the Processed Data

After `scripts/01_data_prep.R` runs, inspect the cleaned data before fitting the model.

```r
df <- readRDS("data/processed/cleaned_data.rds")
str(df)
summary(df)
head(df)
table(df$Group)
table(df$Time)
length(unique(df$Subject))
```

Confirm these points:

- `Y` is numeric.
- `Subject` identifies repeated observations or random-effect groups.
- `Group` has the expected categories.
- `Time` has the expected scale and type.
- There are enough observations per subject, group, and time point.
- The cleaning step did not remove too many rows unexpectedly.

## 7. Run the Pipeline Locally

Once the processed data looks correct, run the full pipeline locally:

```bash
make all
```

Expected outputs:

```text
data/processed/cleaned_data.rds
results/model_fitted.rds
results/tables/fixed_effects.csv
results/tables/random_effects.csv
results/figures/residuals_vs_fitted.png
results/logs/session_info.txt
```

## 8. Commit Before HPC

When the local test passes, commit the data-preparation and config changes.

```bash
git add scripts/01_data_prep.R config/config.yml
git commit -m "Adapt data prep for real dataset"
git push origin main
```

Do not commit raw data or generated results unless there is a deliberate reason. The repository is configured to ignore those paths.

## 9. Run on HPC

On HPC:

```bash
cd ~/projects/sagarmatha-microbiome-analysis
git pull --ff-only
conda run -n glmm-hpc Rscript scripts/01_data_prep.R
qsub -v RUN_MODE=dryrun hpc/job_per_mag_lmertest.pbs
qsub hpc/job_per_mag_lmertest.pbs
```

Monitor the job:

```bash
qstat
qstat -u <your_hpc_username>
```

After the job finishes, check:

```bash
ls results/hpc_per_mag_lmertest/tables/
ls results/hpc_per_mag_lmertest/logs/
ls results/hpc_per_mag_lmertest/session/
```

## Practical Rule

If `scripts/01_data_prep.R` does not produce a clean, well-understood `data/processed/cleaned_data.rds` locally, do not run `scripts/02_run_glmm.R` on HPC yet.
