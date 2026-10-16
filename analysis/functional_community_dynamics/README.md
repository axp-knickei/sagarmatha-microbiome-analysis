# Functional Community Dynamics Workflow

This workflow asks whether the overall predicted functional profile changes
across T1–T10, whether it moves less than the 441-MAG taxonomic profile, whether
matched taxonomic and functional transitions move together, and whether the
few participants with T1 observations show descriptive return toward their
initial profile at T10.

It is independent of the completed MAG functional candidate workflow. It does
not require GFF or protein FASTA files and does not perform locus analysis.

## Stages

1. `00_validate_inputs.R` validates paths, schemas, checksums, exclusions, and
   the expected sample count.
2. `01_prepare_profiles.R` aligns the taxonomic matrix and four functional
   layers, filters terms, and writes Table 26.
3. `02_run_multivariate_analysis.R` computes Hellinger, Bray–Curtis, and
   Aitchison sensitivities; runs participant-restricted PERMANOVA and
   dispersion tests; and writes Tables 27–28.
4. `03_analyze_longitudinal_stability.R` calculates reference and adjacent
   distances, matched taxonomic–functional models, participant-level
   sensitivity checks, and Tables 29–31.
5. `04_make_figures_and_report.R` creates Figures 16–19 and the dated
   result-dependent report.

`05_create_smoke_test_data.R` is a test fixture generator and is not a
production stage.

## Run

From the repository root:

```bash
bash scripts/run_functional_community_dynamics.sh \
  --config config/functional_community_dynamics.yml
```

On OpenPBS:

```bash
qsub hpc/job_functional_community_dynamics.pbs
```

The production configuration expects the synchronized MAG functional results.
Machine-specific raw-data paths can be passed through `FUNCDYN_ABUNDANCE_FILE`,
`FUNCDYN_METADATA_FILE`, and `FUNCDYN_MAGFUNC_RESULTS_DIR`.

## Interpretation contract

- KO/Hellinger is primary; KEGG-module and EC are sensitivities; COG is broad
  and descriptive.
- Functional values are predicted genomic potential, not expression, enzyme
  activity, pathway completeness, or metabolite production.
- Functional and taxonomic distances are linked because the functional values
  are derived from carrier-MAG abundance.
- Smaller functional movement can support smoothing or redundancy as an
  interpretation, but cannot prove resilience, buffering, or adaptation.
- The participant-level sensitivity treats each participant as one unit and
  checks that transition conclusions are not driven only by repeated
  observations or an unstable random-intercept fit.
- T1-to-T10 recovery is exploratory because paired support may be limited; calculate it from authorized inputs.

Generated outputs remain under the ignored
`results/functional_community_dynamics/` directory.
