# Diversity and per-MAG manuscript synthesis

## Purpose

This post-processing workflow creates a reader-first publication package from
three completed analyses:

- sample-level alpha diversity and its longitudinal random-intercept LMMs;
- longitudinal MAG beta diversity with restricted-permutation PERMANOVA; and
- the primary feature-wise per-MAG `lmerTest` screen.

It does not modify raw data, rerun the alpha- or beta-diversity pipelines, or
replace the canonical per-MAG result table. The only model fitting performed by
the script is a display-only refit of up to three qualifying per-MAG models per direction so their
documented fixed-effect trends and confidence bands can be drawn consistently.

## Entrypoint

Run from the repository root:

```bash
Rscript scripts/manuscript_diversity_per_mag/01_build_main_figures_tables.R
```

On Nebula, use the tracked OpenPBS wrapper:

```bash
qsub hpc/job_manuscript_diversity_per_mag.pbs
```

For the first run after pulling this workflow, update the existing environment
so the high-resolution `ragg` renderer is installed:

```bash
qsub hpc/update_glmm_env.pbs
```

Wait for that job to finish successfully before submitting the manuscript job.
The update uses `environment.yml`, validates the shared beta-diversity runtime,
and installs `r-ragg` through Conda; no interactive R package installation is
required on a compute node.

The job reads completed artifacts under `results/alpha_diversity*`,
`results/beta_diversity/`, and `results/hpc_per_mag_lmertest/`, plus
`data/processed/cleaned_data.rds` and the exact `data/raw/GDTB.csv` taxonomy
mapping. Required files and their column structures are checked before any output is created.
Sites that keep ignored data outside the checkout can set
`BLACKICE_RAW_DATA_DIR` and `BLACKICE_PROCESSED_DATA_DIR` before submission.

## Outputs

Outputs are written under `results/manuscript_diversity_per_mag/`:

```text
alpha_diversity/
  main_figure_alpha_diversity.{pdf,svg,png,tiff}
  main_table_alpha_diversity.{csv,tsv,md,pdf,svg,png,tiff}
beta_diversity/
  main_figure_beta_diversity.{pdf,svg,png,tiff}
  main_table_beta_diversity.{csv,tsv,md,pdf,svg,png,tiff}
per_mag_lmm/
  main_figure_per_mag_lmm.{pdf,svg,png,tiff}
  main_table_per_mag_lmm.{csv,tsv,md,pdf,svg,png,tiff}
supplementary/
  supplementary_evidence_manifest.{csv,tsv,md}
session_info.txt
```

PDF and SVG are the editable/vector sources. PNG and LZW-compressed TIFF are
exported at 600 dpi for manuscript systems that require raster artwork. The
table CSV/TSV files retain unrounded values; presentation files use
reader-facing formatting.

## Design and scientific boundaries

The color roots are transcribed from the user-specified Ordnance Survey
GeoDataViz Toolkit palette sheet, version 0.9.1. Time and direction also use
ordering, trajectories, point outlines, and labels so meaning does not depend
on color alone. Figures use a 7.2-inch double-column footprint, sans-serif
typography, white backgrounds, quiet axes, vector-first output, and explicit
panel letters.

The main beta-diversity inference is limited to the three pre-specified primary
metrics in `config/beta_diversity.yml`: Bray-Curtis, binary Jaccard, and
Aitchison. The per-MAG table preserves exact MAG identifiers and joins taxonomy
only by `MAG = GDTB.csv Genome`. Nominal per-MAG candidates are always labeled
exploratory unless supported by the configured primary BH-FDR criterion.

## Safe reruns

Existing main figures are protected by default. After reviewing a prior run,
replace only this post-processing output directory with:

```bash
Rscript scripts/manuscript_diversity_per_mag/01_build_main_figures_tables.R --force
```

Do not use broad cleanup commands against `results/`. Review `git status`
before and after every run because the source result directories contain
selected tracked evidence.


## Public replacement and external bundles

The public canonical entrypoint above is assembled from the approved
`01_build_main_figures_tables_public.R` replacement; the private original is
preserved and excluded. The PBS entrypoint and canonical output filenames remain
compatible. No analytical model or selection threshold is changed.

`PUBLIC_DIVERSITY_INPUT_ROOT` may select a separately authorized bundle with the
same `results/`, `config/`, and default `data/` layout. Existing
`BLACKICE_RAW_DATA_DIR` and `BLACKICE_PROCESSED_DATA_DIR` overrides still apply.
`--output-dir` takes precedence over `PUBLIC_DIVERSITY_OUTPUT_DIR`, then the
normal output default. `PUBLIC_DIVERSITY_RASTER_DPI` defaults to 600; synthetic
smoke tests use a smaller value for rendering only, without changing analysis.

The required inputs are the four alpha tables, eight beta statistics/ordination
tables, per-MAG model table, model-ready RDS, exact taxonomy and model-decision
configuration listed by the script. Missing required inputs stop execution.
Captions derive counts, available adjusted-p summaries and displayed selections
from supplied results. An unavailable adjusted value is not interpreted as evidence
of no discovery. The renderer does not establish why omitted metrics were unavailable
or independently verify upstream permutation restrictions or correction families.

Supplementary source directories are optional only for additional evidence
indexing; they cannot substitute for required input files. The original inventory
formats and columns are preserved. `supplementary_inventory_status.csv` records
included, missing or empty source directories, with readable status in the existing
inventory Markdown. Empty selections retain required figures and table schemas,
with an explicit explanatory panel instead of invented data or substitute models.
