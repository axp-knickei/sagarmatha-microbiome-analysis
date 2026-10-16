# MAG Functional Annotation Workflow

This staged R workflow reconstructs functional annotations from a global
MMseqs protein catalogue, maps member proteins back to MAGs, links selected gut
MAGs to taxonomy and abundance, and runs the five analyses prioritized on
2026-08-08.

## Scientific questions

1. Functional convergence and redundancy across sequence clusters.
2. Unknown proteins in small versus large clusters.
3. Taxonomic distribution of functions across the selected gut MAGs.
4. Longitudinal predicted functional potential across `T1`-`T10`.
5. Hysteresis-like equal-recorded-altitude paired contrasts.

The final item is deliberately not labelled proof of a continuous hysteresis
loop. The available metadata support prespecified paired comparisons at 5,200
m, but recorded sampling altitude does not encode the complete intervening
exposure trajectory.

## Stages

| Stage | Script | Role |
| ---: | --- | --- |
| 0 | `00_validate_inputs.R` | Validate paths, preview structures, canonical joins, timepoints, and abundance scale. |
| 1 | `01_prepare_catalogue.R` | Full-file uniqueness, annotation, ID parsing, MAG mapping, and normalized term tables. |
| 2 | `02_catalogue_analyses.R` | Convergence, unknown-cluster, and taxonomic-distribution analyses. |
| 3 | `03_longitudinal_analyses.R` | Unnormalized sample functional potential, mixed models, and paired contrasts. |
| 4 | `04_make_publication_outputs.R` | Numbered Nature-oriented figures and chart map. |
| 5 | `05_build_results_index.R` | Human- and AI-readable index for the original Tables 01-07 and Figures 01-06. |
| 6 | `06_validate_selected_candidates.R` | Same-data validation of `K04028` and EC `1.2.1.87`: carrier attribution, participant trajectories, influence checks, post-hoc contrasts, and captured-fraction sensitivity. |
| 7 | `07_explore_pathway_unknowns.R` | Exploratory functional co-carriage, unknown-singleton audit, novelty/fragmentation readiness, and the follow-up results index. |
| 8 | `08_decompose_shared_unique_carriers.R` | Exact shared-versus-candidate-specific carrier decomposition, component temporal models, and updated follow-up index. |
| 9 | `09_analyze_shared_component_drivers.R` | BH-controlled post-hoc contrasts for the shared component, participant-influence diagnostics, and exact MAG/genus change decomposition. |
| 10 | `10_audit_candidate_annotations.R` | Offline protein-cluster curation of K04028 and EC 1.2.1.87, exact carrier reconciliation, same-representative overlap, and claim-readiness boundaries. |
| test | `06_create_smoke_test_data.R` | Deterministic synthetic fixtures only. |

Shared functions are in `R/common.R`; workflow-local tests are under `tests/`.
The unit-test entrypoint is `run_tests.R`.

## Fixed analysis decisions

- Complete processing and plotting run on Nebula.
- `Ta` and `Tb` are excluded.
- `gut.csv` values are used as provided and are never row-renormalized.
- Full MAG strings are parsed from member protein IDs and matched exactly.
- `GDTB.csv` `Genome` is the taxonomy key; the raw `ordeer` spelling is only
  aliased after reading.
- Protein-cluster member counts are not abundance.
- Longitudinal values are predicted functional potential, not direct gene
  abundance, expression, activity, or metabolite production.
- BH correction is applied within each ontology/model family.

## Run

Production on Nebula:

```bash
qsub hpc/job_mag_functional_annotation.pbs
```

Validate only the full input contract:

```bash
qsub -v MAGFUNC_TO_STAGE=0 hpc/job_mag_functional_annotation.pbs
```

Restart after a completed stage without overwriting it:

```bash
qsub -v MAGFUNC_FROM_STAGE=3 hpc/job_mag_functional_annotation.pbs
```

Run only the focused follow-up after synchronizing completed Stages 0-5:

```bash
qsub -v MAGFUNC_FROM_STAGE=6,MAGFUNC_TO_STAGE=10 hpc/job_mag_functional_annotation.pbs
```

Stage 6 treats `K04028` and EC `1.2.1.87` as candidates selected from the same
dataset. Its leave-one-participant-out and post-hoc p-values are diagnostics,
not independent replication. The unnormalized model remains primary;
row-sum-normalized values are sensitivity-only.

Stage 7 reports MAG-level functional co-carriage but never calls it pathway
completeness. It can describe unknown singletons with current catalogue fields,
but resolving novelty versus fragmentation requires optional protein metrics
(`representative_id`, `protein_length_aa`, `partial_gene`, and
`contig_edge_distance_bp`) and MAG quality metrics (`MAG`, `completeness`, and
`contamination`). Supply those paths through `MAGFUNC_PROTEIN_METRICS_FILE` and
`MAGFUNC_MAG_QUALITY_FILE` when available.

After a completed Stage 7, Stage 8 can be submitted independently with
`MAGFUNC_FROM_STAGE=8,MAGFUNC_TO_STAGE=8`. It partitions the union of the two
candidate carrier sets into shared, K04028-specific, and EC 1.2.1.87-specific
MAGs. Shared membership is defined at MAG level; it does not imply that both
terms occur on the same gene. Shared and candidate-specific components must
reconstruct the Stage-06 candidate potential exactly before models or figures
are produced.

After a completed Stage 8, Stage 9 can be submitted independently with
`MAGFUNC_FROM_STAGE=9,MAGFUNC_TO_STAGE=9`. It tests all 45 T1-T10 contrasts for
the already selected shared component and applies BH correction across
that complete contrast family. It then performs leave-one-participant-out
influence checks and decomposes every paired change exactly into MAG and genus
contributions. These are selection-aware explanations of a same-data result,
not independent confirmation or evidence of functional activity.

After a completed Stage 9, Stage 10 can be submitted independently with
`MAGFUNC_FROM_STAGE=10,MAGFUNC_TO_STAGE=10`. It audits the representative
protein clusters and gut-MAG carriers behind the two selected candidates using
the reviewed metadata stored in `config/mag_functional_annotation.yml`.
Stage 10 requires no network access and performs no new statistical screen. It
reports whether both terms ever occur on the same representative cluster, but
does not infer genomic neighborhood, operon organization, pathway completeness,
expression, enzyme activity, metabolite production, or adaptation.

Use `MAGFUNC_FORCE=1` only when intentionally replacing outputs from stages
being rerun. Private Nebula input paths belong in `MAGFUNC_CLUSTER_FILE` and
`MAGFUNC_ANNOTATION_FILE` environment overrides, not tracked configuration.

Synthetic local validation:

```bash
bash analysis/mag_functional_annotation/tests/run_smoke_test.sh
```

## Output policy

Generated files belong under `results/mag_functional_annotation/`, which is
ignored by Git. `RESULTS_INDEX.md` records the original numbered analysis, and
`FOLLOWUP_RESULTS_INDEX.md` records Tables 08-25 and Figures 07-15 with exact
paths, meanings, and interpretation order. Complete
raw files, RDS intermediates, result matrices, figures, scheduler logs, and
execution records must remain uncommitted.
