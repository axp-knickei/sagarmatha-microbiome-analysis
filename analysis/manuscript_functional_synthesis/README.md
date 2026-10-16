# Manuscript functional synthesis

This workflow creates a reader-first, manuscript-facing synthesis from two
completed production workflows:

- `results/mag_functional_annotation/` (Stages 0-10)
- `results/functional_community_dynamics/` (Stages 0-4)

It does not rerun discovery analyses, change source results, or read raw data.
It creates one five-panel main figure, two compact main tables, exact panel
data, publication-sized PDF/PNG/SVG exports, and a completion marker under
`results/manuscript_functional_synthesis/`.

Run locally from the repository root:

```bash
bash scripts/run_manuscript_functional_synthesis.sh \
  --config config/manuscript_functional_synthesis.yml
```

The OpenPBS entry point is `hpc/job_manuscript_functional_synthesis.pbs`.
Machine-specific source-result locations can be supplied through
`MANSYN_MAGFUNC_RESULTS_DIR` and `MANSYN_FUNCDYN_RESULTS_DIR`.

The figure uses exact swatches from GDV colour palettes v0.9.1, including its
colour-blind-safe blue-lilac-pink bivariate family, explicit panel letters,
Arial-compatible typography, vector exports, and a 183-mm page width. Panel c
adapts the outer labels, translucent fills, matching outlines, and bold region
counts demonstrated by the [R Graph Gallery Venn guide](https://r-graph-gallery.com/14-venn-diagramm),
while explicitly treating equal-circle areas as schematic rather than
proportional. Panel b shows mean +/- standard error; panel d is a numerical driver decomposition and
does not test individual genera; panel e uses separately standardized distance
scales and cannot be interpreted as a percentage of biological activity.
