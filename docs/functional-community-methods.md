# Functional community methods

Use the [workflow README](../analysis/functional_community_dynamics/README.md) for execution and required stage inputs.

## Purpose

This whole-profile workflow consumes generated functional profiles and requires
neither GFF nor protein FASTA inputs.

The workflow asks four questions:

1. Does predicted functional profile shape differ across T1–T10?
2. Is functional movement smaller than movement in the 441-MAG taxonomic
   profile, consistent with redundancy or smoothing?
3. Do matched taxonomic and functional transitions move together within
   participants?
4. Do participants with both T1 and T10 show descriptive return
   toward their T1 profile?

## Inputs and cohort

The taxonomic comparator is rebuilt from `gut.csv` on the exact same included
sample order as every functional layer. Observation metadata comes from the
exact `sampleID` join to `sample_info.csv`. Ta and Tb are excluded. The
production contract expects 86 T1–T10 observations.

KO is primary. KEGG-module and EC-labelled profiles are sensitivities. COG is a
broad descriptive layer. KEGG-module labels do not demonstrate complete
modules.

## Distances and inference

Within each layer, profiles are closed across retained features. Hellinger
distance is primary, Bray–Curtis is a sensitivity, and Aitchison is an
additional sensitivity after validated multiplicative zero replacement.
Filtering removes zero-variance terms and terms present in fewer than 5% of
included samples.

Global timepoint evidence uses PERMANOVA with participant entered before
timepoint and permutations restricted within participant. A dispersion test is
reported beside each global test. Functional global p-values are BH-adjusted
across the functional layer–distance family. Pairwise tests are generated only
for supported functional global results, use subjects observed at both
timepoints, and are BH-adjusted across the 45 T1–T10 contrasts within each
layer–distance family.

Participant-specific distance from earliest observation and adjacent observed
transition distance are calculated for every layer and metric. Matched
taxonomic–functional Hellinger transitions are divided by each layer's median
all-sample pairwise distance before participant-random-intercept models test
coupling and relative movement. A functional/taxonomic ratio below one is
compatible with smoothing; it is not proof of ecological resilience.

A simple participant-level sensitivity check summarizes each participant's
transition coupling with a Spearman correlation and relative movement with a
median log ratio, then tests those participant summaries against zero. This
uses each participant as one unit and provides separate exploratory evidence
about consistency across participants. It is not a direct statistical replacement
for the mixed model: the correlation is unadjusted for transition duration,
whereas the coupling model includes scaled transition duration as a fixed effect.

For the correlation sensitivity, a participant must have at least three matched
transitions and nonzero variation in both standardized distance variables.
Otherwise the correlation is recorded as missing. Only finite participant
correlations enter the two-sided signed-rank test against zero (`exact = FALSE`);
at least three finite participant summaries are required. The existing BH family
includes non-descriptive functional layers within each sensitivity family.
These eligibility and calculation rules remain unchanged.

## Diagnostic reporting boundaries

The current fitting helper records `converged` from the absence of any lme4
convergence message. A boundary/singularity message can therefore produce
`converged = FALSE` even when the optimizer completes with code `0`. Inspect
optimizer completion, singularity and diagnostic messages separately; this
legacy flag alone does not distinguish optimizer failure from a singular fit.
Successful optimizer completion does not certify the model's statistical
reliability. A BH-significance flag likewise does not certify diagnostic quality.

Participant-level sensitivity describes an association with a different estimand
from the duration-adjusted regression slope. Functional profiles are constructed
from carrier-MAG abundances, so their association with taxonomic movement does
not establish independent biological causation. Interpret any participant
omission checks within the examined scope; they do not establish robustness to
every possible subset or influence on the pooled mixed-model coefficient.

## Recovery boundary

T1-to-T10 return is descriptive and exploratory. Determine actual paired support
from authorized inputs, display exact paired counts, and do not relabel later
visits as baselines or infer cohort-wide recovery from limited paired support.

## Outputs and safety

Tables 26–31, Figures 16–19, distance objects, session information, and the
result-dependent interpretation report are generated beneath the ignored
`results/functional_community_dynamics/` directory. Raw data, RDS files,
distance objects, figures, tables, scheduler logs, and other generated outputs
must not be committed.
