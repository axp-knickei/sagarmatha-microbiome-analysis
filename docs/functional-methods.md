# Functional analysis methods

Use [the functional workflow README](../analysis/mag_functional_annotation/README.md)
for inputs, stage order, configuration, execution and outputs. Inputs are supplied
separately with authorization. Read [exact identifier rules](data/input-interpretation.md).
The configured candidate terms are scientific settings, not distributed findings.
The follow-up is selection-aware same-data analysis, not independent replication.

## Fixed primary contract

- `Ta` and `Tb` remain excluded.
- `gut.csv` values remain unnormalized in the primary model.
- Participant is retained as a random intercept.
- A predicted functional potential is the sum of the provided values for MAGs
  carrying a function; it is not expression, activity, flux, or metabolite
  production.
- Candidate contributions count each carrier MAG value once and do not weight
  by protein-copy number.
- BH correction is applied across all 45 T1–T10 pairwise contrasts separately
  within each selected candidate.
- Stage-9 shared-component post-hoc tests form one prespecified family of all
  45 T1–T10 contrasts and receive one BH correction across that family.

## Stage 6 — selected-candidate validation

### Carrier MAG and taxonomic attribution

For each candidate, the workflow reads the corresponding Stage-02 MAG-function
table and lists every exact carrier MAG. It joins `GDTB.csv` using the complete
MAG identifier and calculates each MAG's descriptive contribution to candidate
potential across samples. Contributions are also aggregated by genus and
phylum.

This identifies which taxa numerically construct the candidate signal. It does
not show that those taxa express the function.

### Participant trajectories and leave-one-participant-out refits

The original mixed model is refitted after excluding each participant in turn.
The output records the omnibus raw p-value, singular-fit status, and maximum
change in the estimated timepoint profile.

These raw p-values diagnose influence. They are not reinterpreted as new
confirmatory p-values because the two functions were selected from the
original same-data screen and the remaining thousands of ontology terms are
not refitted in every leave-one-out dataset.

### Controlled post-hoc contrasts

The workflow estimates all 45 fixed-effect T1–T10 differences for each
candidate from the full participant-aware model. Contrasts are expressed as
`after − before` on the log10 scale, with 95% Wald intervals, ratios on the
original scale, and BH adjustment within candidate.

The post-hoc analysis answers which timepoint pairs may drive the already
supported omnibus result. It should not be used to claim a steady altitude
effect, because expedition timepoint also contains other changing exposures.

### Captured-MAG-fraction sensitivity

Three model scenarios are reported side by side:

1. the approved unnormalized primary-model refit;
2. the same response with centered log10 `gut.csv` row sum as a covariate;
3. an explicitly secondary response divided by the row sum of the 441 retained
   MAGs.

The secondary normalized model does not replace the primary model. If the
signal weakens materially, it should be described as sensitive to the captured
MAG fraction.

## Stage 7 — pathway co-carriage and unknown-singleton audit

### Functional co-carriage

For each candidate carrier set, the workflow compares the presence of KO,
KEGG-module, KEGG-pathway, EC, and CAZy terms in carrier versus noncarrier gut
MAGs using Fisher tests and BH correction within candidate and comparison
ontology.

This is an exploratory **MAG-level co-carriage analysis**. It does not establish
pathway completeness, biochemical activity, or causal coordination. Carrier
MAGs may also be taxonomically related, so enrichment can reflect lineage.

### Unknown singleton proteins

The available catalogue supports descriptive unknown/unannotated singleton
fractions by:

- global versus gut-selected catalogue scope;
- construction-source specimen code;
- retained bin label (`strict`, `permissive`, or `orig`);
- GTDB phylum and genus.

These counts are catalogue entries, not participant-level replicates.

The current two TSVs do **not** contain the evidence needed to decide whether
unknown singleton proteins are novel rather than fragmented predictions. The
workflow creates an explicit readiness table and can use optional inputs when
available:

- protein metrics: `representative_id`, `protein_length_aa`, `partial_gene`,
  and `contig_edge_distance_bp`;
- MAG quality: `MAG`, `completeness`, and `contamination`.

The `_length_` component in many `NODE` identifiers is contig length, not
protein length, and must not be substituted for amino-acid length.

## Stage 8 — shared-versus-candidate-specific carrier decomposition

Stage 8
partitions the union into three mutually exclusive classes: MAGs carrying both
selected terms, MAGs carrying only K04028, and MAGs carrying only EC 1.2.1.87.
This is MAG-set overlap; it does not demonstrate that the terms occur on the
same gene or belong to a complete biochemical pathway.

For every candidate-sample, Stage 8 sums the provided `gut.csv` values once
for shared carriers and once for candidate-specific carriers. The two
components must exactly reconstruct the Stage-06 candidate potential. It then
fits participant-aware exploratory timepoint models for the common shared
component, each candidate-specific component, and the within-candidate log10
balance between shared and candidate-specific components. BH correction is
applied within the absolute-component and balance-model families.

The component models are selection-aware explanations of the original
same-data signals. They are not independent replication, gene expression, or
evidence that a carrier MAG uses the encoded function in vivo.

## Stage 9 — shared-component post-hoc and driver analysis

Stage 9 analyzes the component shared by both configured candidates. It refits the participant-aware
shared-component model, calculates every one of the 45 T1–T10 differences, and
applies BH correction across the complete contrast family. A contrast is
treated as supported only when this adjusted value is below the configured
threshold; unadjusted rankings remain exploratory.

The full shared-component model is then refitted after excluding each
participant. These leave-one-participant-out raw p-values and profile changes
are influence diagnostics, not new replication.

For every timepoint pair, the analysis compares only participants observed at
both timepoints. Their paired shared-component change is decomposed into exact
MAG contributions and exact genus contributions. Positive values mean that a
MAG or genus contributes to a higher shared potential at the later timepoint;
negative values mean it contributes to a lower value. The contributions add
back to the observed mean paired change. Large opposing positive and negative
contributions are recorded as cancellation rather than hidden.

This decomposition is numerical attribution. It does not show gene expression,
activity, adaptation, or a causal effect of altitude.
