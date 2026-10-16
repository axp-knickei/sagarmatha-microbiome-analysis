# Synthetic beta-diversity smoke-test data

Generated with seed 20260711 using only R random draws and fictional metadata.
Contains 16 samples, four fictional subjects, four timepoints and 60 fictional MAGs.
No study input, participant record, real MAG identifier or taxonomy is copied.
The synthetic dates, ages, sex labels and altitudes have no scientific meaning.
No tree is supplied, allowing explicit UniFrac-skip behavior to be tested.

Regenerate from the repository root:

```bash
Rscript analysis/beta_diversity/09_create_smoke_test_data.R --config config/beta_diversity_test.yml
```

Earlier Git history contains real-data-derived fixtures. Replacing current files
does not remove those historical disclosures; review history before publication.
