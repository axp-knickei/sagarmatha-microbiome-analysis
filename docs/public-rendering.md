# Public rendering and validation

These public replacements format separately supplied analytical outputs. The public
archive contains neither study outputs nor the research inputs required to produce
them. Generated figures, tables, captions, panel data and runtime provenance are
research outputs; keep them outside the public source archive.

## Per-MAG figures and tables

Run the `_public.R` entrypoints listed in the
[replacement provenance map](releases/v1.0.0-public-replacements.csv) from the
repository root. Existing data-path defaults remain supported. To use external
inputs and keep output separate, set these environment variables:

```bash
export PUBLIC_LMM_RESULTS=/path/to/authorized/per_mag_lmertest_hpc_results.csv
export PUBLIC_CLEANED_DATA=/path/to/authorized/cleaned_data.rds
export PUBLIC_TAXONOMY=/path/to/authorized/GDTB.csv
export PUBLIC_RENDER_OUTPUT_DIR=/path/to/private/generated-renderings
Rscript scripts/plot_figure_4_hpc_per_mag_lmertest_results_public.R
```

Each entrypoint uses only the inputs required by its original implementation.
Figure 6 needs the taxonomy `label` column as well as the documented taxonomic
ranks; `ordeer` remains the expected order-column spelling. Table 1 calculates
aggregate participant metadata from authorized cleaned data. Figures 5 and 8
refit the existing model for selected MAGs, so they also require cleaned data and
`lme4`. Their nominal candidate-selection rules remain exploratory even when
other MAGs pass FDR correction. Figure 10 retains its existing FDR-first selection
and nominal fallback. No script recomputes or changes the primary FDR family.

Counts, captions and status summaries describe the supplied results. Historical
run timestamps were removed from per-MAG output stems; stable stems are used.
Existing output files at the chosen destination may be replaced, as in the
original renderers. Use a separate destination for each run.

## Functional synthesis

```bash
export MANSYN_MAGFUNC_RESULTS_DIR=/path/to/authorized/mag-functional-results
export MANSYN_FUNCDYN_RESULTS_DIR=/path/to/authorized/functional-dynamics-results
export MANSYN_OUTPUT_DIR=/path/to/private/functional-synthesis
bash scripts/run_manuscript_functional_synthesis.sh --config config/manuscript_functional_synthesis.yml
```

The public renderer and validator read completed stage outputs, including the
Stage 10 Table 23 protein evidence, Table 24 carrier membership and Table 25
curation summary. Carrier counts and protein overlap must agree across those
files. Protein overlap may be zero or nonzero; neither outcome proves complete
pathways. Cohort counts come from supplied candidate trajectories. Transition
panel data must agree with the supplied primary-layer transition table.

The configured candidates and `strongest_shared_contrast` are retained analytical
choices. The renderer does not silently select a different contrast when inputs
change. Genus contributions for that configured contrast are derived from Table
21, using the existing eight-largest-absolute-contributions plus remaining-genera
aggregation. This prevents a pre-rendered Figure 15 for a different contrast from
being mislabeled. The leading genus and before/after labels come from supplied
rows. Reported p/q values and model status replace predetermined interpretation;
inspect upstream convergence and singularity diagnostics before scientific claims.
A reported p-value of zero is printed as zero, without inventing an unreported
numerical upper bound. Return summaries remain descriptive, with no cohort-level
recovery inference.

## Historical exploratory NetMoss synthesis

```bash
export NETMOSS_MANUSCRIPT_INPUT_DIR=/path/to/authorized/netmoss-results
export NETMOSS_MANUSCRIPT_GUT_FILE=/path/to/authorized/gut.csv
export NETMOSS_MANUSCRIPT_TAXONOMY_FILE=/path/to/authorized/GDTB.csv
export NETMOSS_MANUSCRIPT_OUTPUT_DIR=/path/to/private/netmoss-synthesis
bash scripts/run_netmoss_manuscript_synthesis.sh --config config/netmoss_manuscript_synthesis.yml
```

The historical T1/T2/Ta versus T8/T9/T10 grouping, aggregation, thresholds and
exploratory scientific boundary remain unchanged. Counts come from the supplied
manifest and score table. The validator compares scores, taxon identities, ranked
table membership and every displayed network edge against source files, checking
threshold flags rather than expecting particular genera or associations. Palette
checks apply to categories actually displayed. The renderer supports fewer than
five taxonomic phyla without inventing categories.

## Synthetic checks

```bash
Rscript tests/public_release/test-public-renderers.R
Rscript tests/public_release/test-public-synthesis.R
```

Fixtures are created from artificial values in temporary directories. No study
file is read and no repository result is regenerated. The first check renders all
nine per-MAG replacements with and without FDR discoveries. The second renders
both synthesis workflows with different artificial counts, gene overlap and
contrast direction, including zero/nonzero protein overlap and empty networks,
then checks that inconsistent outputs are rejected. Raster
resolution is reduced for smoke checks only. These checks do not certify a full
scientific rerun with authorized research inputs or the future frozen release.

## Provenance

The map records relative original/replacement paths and SHA-256 content hashes.
Hashes identify the reviewed bytes without distributing original content, private
repository addresses or Git history. Original files remain unchanged privately;
only approved public replacements may enter the future clean repository. Refresh
replacement hashes whenever a replacement changes, and verify the final archive
against its approved manifest before release.


## Combined diversity and per-MAG package

The public canonical command is:

```bash
export PUBLIC_DIVERSITY_INPUT_ROOT=/path/to/authorized/analysis-bundle
export PUBLIC_DIVERSITY_OUTPUT_DIR=/path/to/private/diversity-rendering
Rscript scripts/manuscript_diversity_per_mag/01_build_main_figures_tables.R
```

See the [workflow guide](../scripts/manuscript_diversity_per_mag/README.md#public-replacement-and-external-bundles)
for required inputs, path precedence, optional supplementary status and output
formats. The private original remains excluded; the approved public replacement
occupies its canonical destination. The combined workflow is tested within the
existing public-renderer suite with artificial result scenarios and missing-input
checks. Authorized-output numerical comparisons are retained privately. Neither
successful rendering nor a nominal p-value establishes biological causation.


## Random Forest feature stability

The public canonical entrypoint is assembled from
`plot_random_forest_feature_stability_public.R`; its private original is excluded.
Existing CLI arguments and output filenames/formats remain compatible:

```bash
Rscript scripts/plot_random_forest_feature_stability.R \
  --feature-results-dir /path/to/authorized/run/features \
  --taxonomy-file /path/to/authorized/GDTB.csv \
  --output-dir /path/to/private/rf-feature-rendering
```

The five feature CSVs, exact taxonomy and resolved configuration remain required.
The existing renderer contract requires 441 features, 13 LOSO folds, a complete
eligibility grid, disabled `top_n`, nondegenerate rank-IQR tertiles and a fold-variable
feature for the eligibility panel. These checks are preserved; the renderer is
not claimed to support every possible cohort/feature layout. Feature importance,
rankings, quantiles, grouping and numerical agreement tolerance are unchanged.

The pooled overall and per-class performance CSVs are optional reporting context.
Missing files, columns or values are labeled unavailable; no metric is reconstructed
from another file. Known zero-recall counts use only classes with available recall,
with missing classes identified separately. Invalid probabilities, duplicate class
identifiers and ambiguous multiple pooled rows cause explicit errors. No qualitative
performance category is assigned. Predictive importance remains exploratory and
is not an effect size, significance test or causal conclusion.

The candidate-ineligibility statement uses displayed input rows. Validation checks
negative-value preservation against supplied values, including when none is negative.
Default raster quality remains 400 dpi PNG and 600 dpi combined TIFF; optional
`PUBLIC_RF_RASTER_DPI` changes rendering resolution only for smaller smoke checks.
Synthetic fixtures use artificial identifiers, taxonomy, importance values and class
labels under the existing dimension contract, without reading research data.
