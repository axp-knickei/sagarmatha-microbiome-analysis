# NetMoss methods, readiness and companions

Network scores remain exploratory. Internal bootstrap/leave-one-out stability is
not independent validation and does not establish microbial interactions or a
causal altitude effect. Input profiles are generated from separately authorized
abundance, metadata and taxonomy; exact joins and subject aggregation are required.
The main configuration excludes Ta/Tb and compares T1/T2/T3 with T8/T9/T10;
the distinct historical configuration retains its Ta-inclusive design. Never
relabel one scenario's output as another. Pinned dependencies are documented in
[the environment guide](public-environment.md).

## Files and roles

| File | Role |
|---|---|
| `config/netmoss_analysis.yml` | Contrast, aggregation, filtering, network, NetMoss2, and figure settings |
| `scripts/netmoss/install_netmoss2.R` | Install/check the two pinned GitHub packages and rare dependencies |
| `scripts/netmoss/validate_environment.R` | Fail-fast environment validator |
| `scripts/netmoss/run_netmoss_pipeline.R` | Prepare participant-aggregated inputs, build networks, and run NetMoss2 |
| `scripts/netmoss/create_publication_outputs.R` | Produce source data, tables, figures, legends, and Methods text from NetMoss output |
| `hpc/update_netmoss_env.pbs` | Update the existing Conda environment and install pinned packages |
| `hpc/job_netmoss_readiness.pbs` | Run the read-only suitability assessment |
| `hpc/job_netmoss_analysis.pbs` | Run exploratory NetMoss2 and publication stages |
| `hpc/job_netmoss_publication.pbs` | Rebuild publication outputs without rerunning NetMoss2 |
| `hpc/job_netmoss_smoke_test.pbs` | Scheduler-level test using test/synthetic data only |
| `hpc/tests/test_netmoss_pbs.sh` | Static and executable PBS smoke-test harness |

Separate companion files are indexed in the matched-timepoint guide; they do
not alter or replace the NetMoss2 files in this table.

## 1. Update the environment on Nebula

After pulling the repository, submit:

```bash
cd ~/projects/sagarmatha-microbiome-analysis
qsub hpc/update_netmoss_env.pbs
```

The job loads `intelpython/2022.2.0`, updates the existing `glmm-hpc`
environment from `environment.yml`, installs the pinned GitHub revisions, and
verifies every required package. The update stage requires internet access; the
analysis jobs do not download software.

Monitor it with:

```bash
qstat -u "$USER"
tail -f logs/blackice_netmoss_env.pbs.log
```

If the cluster uses different names:

```bash
qsub -v CONDA_MODULE_NAME=intelpython/2022.2.0,CONDA_ENV_NAME=glmm-hpc \
  hpc/update_netmoss_env.pbs
```

## 2. Run the HPC smoke test

Always run this before the real-data analysis:

```bash
qsub hpc/job_netmoss_smoke_test.pbs
```

The smoke job checks PBS syntax and directives, rejects Slurm commands, parses
all R scripts, runs readiness diagnostics on the test microbiome data, prepares
paired test networks, and generates publication outputs from deterministic
synthetic NetMoss-like results. Synthetic outputs are strictly software tests
and must never be interpreted biologically.

Expected final log line:

```text
All NetMoss OpenPBS smoke tests passed. # Expected harness completion message
```

The same harness can be run locally:

```bash
bash hpc/tests/test_netmoss_pbs.sh /tmp/blackice_netmoss_smoke
Rscript scripts/netmoss/tests/testthat.R
```

## 3. Confirm the scientific contrast

Review `config/netmoss_analysis.yml` before submission. The conservative default:

- aggregates MAGs to genus;
- excludes Ta and Tb before analysis;
- compares T1/T2/T3 with T8/T9/T10;
- summarizes repeated samples within participant using the median;
- uses a 50% prevalence screen;
- closes each retained profile to relative abundance;
- constructs separate SparCC networks;
- uses a fixed random seed and pinned NetMoss2 revision.

The contrast is an inferred early-versus-late expedition scenario. It is not a
case–control study, and it must be biologically approved before use. Participant
aggregation prevents pseudoreplication but leaves very few independent profiles.

## 4. Run exploratory NetMoss2

The full job deliberately requires explicit acknowledgement:

```bash
qsub -v NETMOSS_ACKNOWLEDGE_EXPLORATORY=YES hpc/job_netmoss_analysis.pbs
```

Without that variable the job exits before analysis. The job writes:

```text
results/netmoss/prepared/
results/netmoss/networks/
results/netmoss/tables/netmoss_scores.csv
results/netmoss/objects/netmoss_result.rds
results/netmoss/logs/netmoss_run_summary.txt
results/netmoss/reproducibility/netmoss_manifest.yml
results/netmoss/reproducibility/sessionInfo.txt
```

The wrapper constructs its own validated condition networks. It does not call
the NetMoss2 `netBuild(..., method = "pearson")` helper because inspection of
the upstream source found that its single-file control branch calculates the
control correlation from the case object.

## 5. Create publication outputs

The full job runs this automatically. To regenerate formatting only:

```bash
qsub hpc/job_netmoss_publication.pbs
```

Outputs under `results/netmoss/publication/` include:

- a composite NetMoss score, abundance, degree-change, and paired-network figure
  in PDF, SVG, and LZW-compressed 600-dpi TIFF;
- a source-data table for the main figure;
- complete taxon results and condition-network edge tables;
- a figure legend;
- a reproducible Methods draft;
- a reporting checklist, manifest, and `sessionInfo()`.

These files are formatted to support manuscript preparation, but journal
suitability still requires author review, correct biological labels, accessibility
checks, independent network-stability analyses, and compliance with the exact
Nature-family journal instructions current at submission.

## Interpretation safeguards

- Association edges are not direct ecological interactions.
- NetMoss score and FDR thresholds must be declared before inspecting outcomes.
- Network/module stability should be evaluated with participant-level bootstrap
  or leave-one-participant-out sensitivity analysis before biological claims.
- Concordance with LMM/GLMM, differential abundance, beta diversity, and Random
  Forest results should be reported without treating agreement as independent
  replication.
- The 6,500 m sampling point is not a summit sample.
- Generated results and scheduler logs remain ignored by Git; scripts,
  configuration, documentation, and tests are tracked.


## Readiness execution

From the repository root:

```bash
Rscript scripts/check_netmoss_readiness.R
```

With the project Conda environment:

```bash
conda run -n glmm-hpc Rscript scripts/check_netmoss_readiness.R
```

Paths and thresholds are configured in `config/netmoss_readiness.yml`. Command
line overrides are also available:

```bash
Rscript scripts/check_netmoss_readiness.R \
  --abundance data/raw/gut.csv \
  --metadata data/raw/sample_info.csv \
  --taxonomy data/raw/GDTB.csv \
  --output-dir results/netmoss_readiness
```

The script refuses to replace its named outputs unless `--force` is supplied.

For pinned NetMoss2 installation, OpenPBS jobs, smoke tests, exploratory
execution, and manuscript-oriented result rendering, see
[`docs/netmoss-methods.md`](netmoss-methods.md).

All results are written to `results/netmoss_readiness/`. The principal report is
`netmoss_readiness_summary.txt`; machine-readable tables record the dataset,
metadata, sample matching, longitudinal completeness, sparsity, prevalence,
filtering sensitivity, scenarios, and assumption checklist. Every requested
figure is exported as both PDF and SVG, and `sessionInfo.txt` records the R
environment.

## Low/initial versus high-altitude contrast

Configuration `config/netmoss_altitude_analysis.yml` uses T1/T2/T9/T10 versus
T3–T8, excluding Ta/Tb. Collapse repeated observations to a median profile per
participant and group. This controls repeated-sample counting but does not fit a
paired longitudinal network model. The low/initial group mixes initial and
recovery states; expedition time and other exposures remain confounded.

- Configuration: `config/netmoss_altitude_analysis.yml`
- Shared analysis: `scripts/netmoss/run_netmoss_pipeline.R`
- OpenPBS job: `hpc/job_netmoss_altitude_analysis.pbs`
- Generated outputs: `results/netmoss_altitude/`

Run the existing NetMoss smoke test first:

```bash
qsub hpc/job_netmoss_smoke_test.pbs
```

Then submit the third comparison with explicit acknowledgement:

```bash
qsub -v NETMOSS_ALTITUDE_ACKNOWLEDGE_EXPLORATORY=YES \
  hpc/job_netmoss_altitude_analysis.pbs
```

The new configuration and output directory are separate from the preserved
`T1/T2/T3` versus `T8/T9/T10` workflow, so this analysis does not relabel or
overwrite the existing NetMoss results.

## Stability: Predeclared analysis

- Group A: T1, T2, T9, and T10, described as low altitude and initial reach.
- Group B: T3 through T8, described as the high-altitude expedition period.
- Ta and Tb are excluded.
- Repeated observations are collapsed to a median profile per participant and
  group before network construction.
- Genus primary filter: detected in at least 80% of participants in both groups.
- Family primary filter: detected in at least 50% of participants in both groups.
- Matched-cohort prevalence sensitivities: 50%, 70%, and 80% at both ranks.
- The existing NetMoss minimum module size remains 20. Scenarios retaining
  fewer than 20 taxa are recorded as not estimable.

Record the retained-taxon count for every scenario from the supplied inputs.
Scenarios below the existing minimum module size are recorded as not estimable;
do not infer estimability from a scenario label alone.

## Stability: Stability assessment

For every estimable scenario, the workflow constructs a reference SparCC
network pair and NetMoss result, then performs 1,000 deterministic participant
bootstrap replicates and every leave-one-participant-out reconstruction. Common
participants are sampled as paired identities in both groups. In the secondary
cohort, high-only participants are resampled as a separate stratum.

The reference taxon set is fixed before resampling. Failed fits remain in the
failure table and reduce the success fraction; they are never replaced with
zeroes or favorable estimates.

An operationally stable network scenario requires at least 90% successful
bootstrap fits and median edge-vector similarity of at least 0.70 in both
conditions. A reference edge additionally requires the same sign and absolute
correlation of at least 0.40 in at least 80% of successful bootstraps and the configured minimum number of matched-cohort leave-one-out runs. A reference top-10 NetMoss taxon uses the
same resampling-frequency requirements. Module adjusted Rand indices and
degenerate-module warnings are reported separately.

These thresholds are operational diagnostics, not universal guarantees of
network validity.

## Stability: Entrypoints and outputs

- Configuration: `config/netmoss_altitude_stability.yml`
- Shared functions: `scripts/netmoss/altitude_stability_helpers.R`
- Analysis: `scripts/netmoss/run_altitude_stability.R`
- Report: `scripts/netmoss/summarize_altitude_stability.R`
- OpenPBS analysis: `hpc/job_netmoss_altitude_stability.pbs`
- OpenPBS smoke test: `hpc/job_netmoss_altitude_stability_smoke_test.pbs`
- Generated output: `results/netmoss_altitude_stability/`

Principal tables include complete scenario accounting, resampling metrics,
failed fits, edge stability, module stability, NetMoss rank stability,
participant influence, and genus-to-family concordance. The generated report is
`report/netmoss_altitude_stability_report.md`.

Run the smoke test first:

```bash
qsub hpc/job_netmoss_altitude_stability_smoke_test.pbs
```

Then submit the full internal stability analysis:

```bash
qsub -v NETMOSS_STABILITY_ACKNOWLEDGE_EXPLORATORY=YES \
  hpc/job_netmoss_altitude_stability.pbs
```

If either primary scenario fails its operational thresholds, the correct
conclusion is that the corresponding network is unstable and unsuitable for
biological claims. If it passes, it remains an exploratory high-altitude
expedition-period association rather than a causal altitude effect.

## Matched-timepoint abundance companion

The workflow:

1. joins abundance observations to metadata by exact `sampleID`;
2. joins all 441 MAGs to taxonomy by exact `Genome` and sums them to genus;
3. closes each sample to relative abundance and applies a documented
   zero-replacement and CLR transformation;
4. excludes `Ta` and `Tb`, then enumerates all 45 pairs from `T1` through `T10`
   before examining results;
5. retains only participants observed at both timepoints within each pair;
6. records every participant–genus change so individual trajectories remain
   inspectable;
7. uses an exact paired sign-flip test only for pairs with at least six matched
   participants and genera detected in at least three matched participants;
8. applies one Benjamini–Hochberg correction across every testable
   genus-by-timepoint-pair hypothesis;
9. requires a global FDR result, a bootstrap interval excluding zero,
   75% direction concordance, and 80% leave-one-participant-out sign stability
   before assigning the predeclared `candidate` label.

Pairs below the participant threshold remain in the outputs as descriptive
records. No pair is removed because its result looks uninteresting, and a high
effect or low raw p-value alone is not a candidate.

- Configuration: `config/netmoss_matched_timepoint_companion.yml`
- Shared functions: `scripts/netmoss/matched_timepoint_helpers.R`
- Analysis: `scripts/netmoss/run_matched_timepoint_companion.R`
- Figures/report: `scripts/netmoss/summarize_matched_timepoint_companion.R`
- Advanced OpenPBS job: `hpc/job_netmoss_matched_timepoint_companion.pbs`
- OpenPBS smoke test: `hpc/job_netmoss_matched_timepoint_smoke_test.pbs`

Generated outputs are written under
`results/netmoss_matched_timepoint_companion/`:

```text
tables/timepoint_support.csv
tables/timepoint_pair_support.csv
tables/all_genus_pair_results.csv
tables/all_participant_genus_changes.csv
tables/predeclared_candidates.csv
tables/candidate_participant_changes.csv
figures/
report/matched_timepoint_companion_report.md
logs/analysis_summary.txt
reproducibility/analysis_manifest.yml
reproducibility/sessionInfo.txt
```

The complete tables are authoritative. Top-ranked diagnostic figures are not
permission to bypass the candidate rules.

The supplied launcher uses OpenPBS CPU resources and does not request a GPU.
Confirm current site queues and resource availability before submission.

Run the test job first:

```bash
qsub hpc/job_netmoss_matched_timepoint_smoke_test.pbs
```

After reviewing this document, submit the advanced analysis:

```bash
qsub -v NETMOSS_MATCHED_ACKNOWLEDGE_EXPLORATORY=YES \
  hpc/job_netmoss_matched_timepoint_companion.pbs
```

The advanced job uses `workq`, 8 CPU workers, 16 GB RAM, and 10,000
participant-bootstrap replicates per pair. Work is parallelized across the 45
approved timepoint pairs. BLAS/OpenMP threads are held at one per worker to prevent CPU
oversubscription.

Use the launcher's resource requests as documented defaults, subject to site
policy. Assess runtime and memory with the supplied authorized inputs before
increasing parallelism; do not change scientific settings to improve runtime.

This workflow is exploratory and abundance-based. If no comparison passes all
predeclared safeguards, the correct conclusion is that the cohort did not
provide robust all-pairs evidence—not that a more favorable threshold or
timepoint subset should be searched.

For public execution and input requirements, see the
[workflow guide](public-workflow.md).
