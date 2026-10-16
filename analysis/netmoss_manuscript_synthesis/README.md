# NetMoss manuscript synthesis

This post-processing workflow builds one reader-first, four-panel main figure
and two compact main tables from the preserved historical exploratory NetMoss
run under `results/netmoss/`. It does not recalculate SparCC correlations,
NetMoss scores, p-values, q-values, or scientific thresholds.

The main figure contains:

- participant support before and after within-period aggregation;
- the ranked NetMoss score landscape;
- abundance change versus NetMoss score, emphasizing that the two quantities
  have different meanings;
- paired enhanced association networks using the same visual grammar as
  `scripts/netmoss/create_enhanced_network_figure.R`.

Exact colors are drawn from Ordnance Survey GeoDataViz colour palettes v0.9.1.
Solid network edges retain the original `|r| >= 0.40` rule. Faint dashed edges
between `0.25 <= |r| < 0.40` are visual context only and are explicitly
identified as below threshold.

Run locally from the repository root:

```bash
bash scripts/run_netmoss_manuscript_synthesis.sh --force
```

The OpenPBS entrypoint is `hpc/job_netmoss_manuscript_synthesis.pbs`. Generated
outputs are written under `results/netmoss/manuscript_synthesis/` and remain
ignored under the repository result policy.

This workflow intentionally represents the preserved `T1/T2/Ta` versus
`T8/T9/T10` result. It must not be relabeled as the current Ta/Tb-excluded
configuration. The output remains exploratory and cannot establish microbial
interactions, stable network rewiring, or validated biomarkers.
