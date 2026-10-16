# Configuration Guide

Configuration is split by workflow so scientific settings and output locations
can be reviewed without editing analysis code. YAML uses two-space indentation.

| File | Owner / consumer | Purpose |
|---|---|---|
| `config.yml` | `scripts/01_data_prep.R`, legacy `02_run_glmm.R` | Raw-data paths, processed-data path, and legacy single-model settings |
| `model_decision_rules.yml` | Primary per-MAG workflow | Model status, sparsity, FDR, and interpretation rules |
| `response_scale_decision.yml` | Primary per-MAG workflow | Working primary and sensitivity response-scale decision |
| `beta_diversity.yml` | Beta-diversity pipeline | Production input, distance, permutation, resource, and output settings |
| `beta_diversity_test.yml` | Beta-diversity smoke test | Test-data override of the production beta-diversity configuration |
| `netmoss_readiness.yml` | NetMoss readiness assessment | Input discovery, Ta/Tb exclusion, and readiness thresholds |
| `netmoss_analysis.yml` | NetMoss analysis and publication scripts | Exploratory T1/T2/T3 versus T8/T9/T10 scenario, Ta/Tb exclusion, filtering, pinned package revision, and output settings |
| `netmoss_historical_publication.yml` | Historical NetMoss publication reproduction | Preserved T1/T2/Ta versus T8/T9/T10 settings; use an isolated output override |
| `netmoss_altitude_analysis.yml` | NetMoss analysis and publication scripts | Separate exploratory T1/T2/T9/T10 versus T3–T8 low/initial-to-high-altitude scenario and output path |
| `netmoss_altitude_stability.yml` | NetMoss altitude stability scripts | Matched-cohort genus/family scenarios, prevalence sensitivities, deterministic resampling, operational stability thresholds, and isolated outputs |
| `netmoss_manuscript_synthesis.yml` | NetMoss manuscript post-processing | Preserved historical result inputs, unchanged score/edge thresholds, reader-first figure settings, and generated manuscript output path |
| `netmoss_matched_timepoint_companion.yml` | Matched-timepoint companion scripts | Ta/Tb exclusion, all 45 T1–T10 pairs, CLR transformation, global correction, stability gates, and advanced execution settings; does not configure NetMoss2 |
| `random_forest_config.yml` | `scripts/run_random_forest_model.R` | Exploratory Random Forest data, feature, model, and output settings |
| `random_forest_subject_validation.yml` | Subject-aware validation pipeline | Validation mode, grouped resampling, feature filtering, and output settings |
| `mag_functional_annotation.yml` | MAG functional annotation pipeline | Full-file QC, ontology parsing, T1-T10 decisions, unnormalized functional potential, paired contrasts, selected-candidate follow-up, reviewed offline candidate curation and claim boundaries, shared/unique carrier decomposition, shared-component post-hoc correction family, optional novelty/fragmentation evidence, threading, and output paths |
| `functional_community_dynamics.yml` | Functional community dynamics pipeline | Identical T1–T10 sample contract, KO-primary functional profile analysis, taxonomic comparator, distance sensitivities, restricted permutations, filtering, FDR, and output paths |
| `manuscript_functional_synthesis.yml` | Manuscript functional synthesis | Completed source-result roots, selected candidates, ordered timepoints, figure dimensions, export resolution, and generated synthesis-output paths |

## Authority and duplication

`model_decision_rules.yml` is the detailed operational source for the per-MAG
modeling rules. `response_scale_decision.yml` is the concise response-scale
decision record. Values repeated between them must remain consistent, but
maintenance work must not change either scientific decision without explicit
approval.

`beta_diversity_test.yml` intentionally repeats production keys so it can be
passed directly to the same scripts. When the production schema changes,
review the test configuration in the same change.

The `project.author` placeholder in `config.yml` is retained because its final
value has not been confirmed. It does not authorize guessing contributor
identity.

## Environment files

- `environment.yml` is the primary cross-language Conda environment used by
  local and HPC workflows.
- `renv.lock` records a narrower R package snapshot. It does not currently
  replace the Conda environment or enumerate every Conda dependency.

Do not regenerate either file as cosmetic maintenance: dependency resolution
can change package versions and scientific reproducibility.

## Machine-specific values

Prefer documented environment-variable overrides or ignored `hpc/*.local.*`
notes for machine-specific paths. Do not commit credentials, private HPC paths,
or raw-data locations. The public subject-validation configuration has no private snapshot fallback;
supply authorized processed input at its configured path.
