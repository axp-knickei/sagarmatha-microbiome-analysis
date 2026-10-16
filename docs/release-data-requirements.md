# Release input contract and data availability

## Intended access policy — confirmed contacts and conditional permissions

The owner confirmed controlled access / available upon reasonable request as the
intended public policy. The public `sagarmatha-microbiome-analysis` repository and
Zenodo code archive exclude underlying study data, metadata and restricted
processed data. The author Alex Prima will receive and coordinate data-access
requests at either owner-supplied public contact address:
[alexprima@student.ub.ac.id](mailto:alexprima@student.ub.ac.id) or
[alex.prima@uq.net.au](mailto:alex.prima@uq.net.au).

### Approved statement for use at public release

> The analytical code and supporting technical documentation are publicly available
> through the `sagarmatha-microbiome-analysis` repository. The underlying study data,
> metadata, and restricted processed data are not included in the public code
> archive. These materials may be made available upon reasonable request to the
> authors, subject to applicable permissions and data-use conditions.

This is the approved intended GitHub/Zenodo data-availability statement. Until
publication is verified, describe the code as being prepared for public release;
the first sentence above is not a claim that the repository is already public.
The study data are not described as openly available or publicly downloadable.
"May be made available" does not promise approval or provision of every requested
component. No public accession or contact address is invented.

### Request contact and permission authority are separate

**Request contact established:** Alex Prima, by explicit owner confirmation, can
receive and coordinate requests. **Public contact details established:** the
owner supplied both addresses above for public data-access requests. CFF's single
author email field records the first listed address; both remain available for
request coordination. No additional affiliation is inferred from an email domain.

**Permission authority remains conditional:** receiving a request does not establish
unrestricted authority to redistribute all data. Access to particular components
may depend on collaborator, institutional, ethical, ownership or data-use
requirements. The responsible approving parties and component-specific conditions
have not been verified in the repository. Confirm applicable requirements before
providing any component; do not imply automatic approval of reasonable requests.

Earlier repository documentation did not establish a request contact; the owner's
confirmation now establishes the coordination role. It does not establish new
redistribution permissions, a data-use agreement, a study accession, or a named
collaborator/institutional approval role. Software authorship, affiliations and
the MIT software license alone are not recorded as study-data permission evidence.
Do not substitute accessions or availability statements from cited research.

### Components require separate confirmation

| Component | Repository evidence and access/reuse status |
|---|---|
| Study abundance inputs and other underlying research data | Excluded; intended controlled access. Applicable sharing authority, permitted recipients and reuse/redistribution conditions are pending. |
| Participant metadata | Excluded; request coordination and public contacts are established, but component-specific sharing conditions remain pending. Do not infer that permission for abundance files also covers metadata. |
| Processed biological data and metadata-derived exports, including the original window CSV | Excluded; processing or aggregation does not establish permission. Any external provision and reuse conditions require separate confirmation. |
| Study protein catalogue and representative/member mappings | External functional-workflow inputs; authority to share and applicable reuse conditions remain unconfirmed. |
| Project taxonomy and functional-annotation tables | External inputs derived using GTDB/eggNOG and associated annotations. Project-data permission and applicable third-party terms need separate verification; exact database versions and terms are not established by this repository. |

Different components may require different conditions, but the repository does not
establish a verified component-by-component set of terms. The existing MIT
license covers repository software and associated documentation; it is not
recorded as a study-data sharing authorization. No data-use agreement is invented
or presumed to exist.

### Code release and later data provision

The confirmed public contacts and conditional request policy are sufficient for
the data-availability documentation of this code-only archive. Underlying study
data are excluded, so publishing the approved code does not require asserting
blanket redistribution permission for those data. Applicable component-specific
permissions and reuse conditions must be verified before any future provision.
Permission to distribute the included code and documentation is separate; retain
the MIT notice and the third-party palette attribution described in the
[README](../README.md#third-party-palette-material). No additional permission
authority, institutional approval or data-use agreement is inferred here.

## Required external input schemas

| Local filename | Format and required content | Consumer |
|---|---|---|
| `data/raw/gut.csv` | CSV; unique `sampleID`, then 441 numeric nonnegative full MAG-name columns; original input has 103 rows | Preprocessing, diversity, RF via cleaned data, networks, functional potential |
| `data/raw/sample_info.csv` | CSV; `sampleID`, `group`, `Timepoint`, `date`, `gender`, `age`, `altitude`; `group` supplies Subject; dates are YYYYMMDD; metadata has 119 rows | Sample-level join and longitudinal covariates |
| `data/raw/GDTB.csv` | CSV; unique `Genome` exactly matching full MAG names; `label`, `domain`, `phylum`, `class`, `ordeer`, `family`, `genus`, `species`, `species_note` used as taxonomy annotations | Taxonomy joins and aggregation; preserve the original `ordeer` spelling on input |
| `data/raw/MAG_gene_original_final_0805_mmseqs_cluster.tsv` | Two-column, headerless representative/member protein mapping; member IDs must resolve to full MAG IDs | Functional annotation |
| `data/raw/MAG_gene_original_final_0805_mmseqs_rep_seq.emapper.annotations.tsv` | eggNOG TSV with detected query header and ontology columns, including seed_ortholog, evalue, score, Description, Preferred_name, KEGG_ko, KEGG_Module, EC, COG_category and PFAMs | Functional annotation and follow-up |

The functional workflow validates the exact annotation contract at stage 0;
see [functional methods](functional-methods.md). Optional protein
length/partial-gene/contig-edge and MAG-quality tables support explicitly labeled
readiness analyses and are not a license to infer unavailable biological claims.
No FASTQ/BAM input is consumed by the documented downstream workflows. Upstream
sequencing, assembly, MAG quality/dereplication, abundance normalization, MMseqs
identity/coverage and annotation database versions still require methods records.

For schemas with fictional values, see [synthetic beta fixtures](../data/test_beta_diversity/README.md).
They are test data, not substitutes for reproduction of study results.

## Joins and cleaning

Join abundance observations to metadata by `sampleID`; join MAGs to taxonomy by
exact `Genome`, retaining complete identifiers. Read [gut.csv provenance](data/input-interpretation.md)
before interpreting MAG names. Embedded specimen/date/altitude fields are
construction provenance, not observation covariates.

Primary preprocessing types Timepoint as T1–T10 and drops missing model-critical
values, excluding Ta/Tb. The subject baseline is calculated before that final
filter; document this existing behavior rather than changing it during packaging.
Output is `data/processed/cleaned_data.rds`, with sampleID, MAG, Y, Subject,
Timepoint, date_converted and days_from_baseline_subject plus joined metadata.
Do not commit this object. Alpha and production beta have separate raw-data paths
and do not inherit the cleaned cohort automatically.

## External locations

### Excluded window metadata and regeneration

The original `results/metadata_expedition_windows/proposed_biological_windows.csv`
is excluded from v1.0.0. It includes aggregate dates, altitude summaries and
participant counts derived from study metadata. The public repository/Zenodo
archive supplies neither this research metadata nor the processed biological
dataset. The six predefined windows and rationale remain in
`scripts/evaluate_expedition_windows.R`.

Follow [window CSV regeneration](public-workflow.md#window-csv-regeneration)
in an isolated authorized checkout. The original table used metadata extracted
from the 86-sample processed dataset, not the unfiltered 119-row raw metadata or
103-sample diversity cohort. Supply a matching authorized processed RDS at
`config/config.yml`'s `data.processed` path, normally
`data/processed/cleaned_data.rds`; the instructions explain how to prevent raw
metadata from taking precedence. This is an external-input requirement, not an
operational data-access mechanism or a claim of verified historical reproduction.

### Workflow-specific input overrides

For preprocessing and alpha preparation, use `BLACKICE_RAW_DIR` when inputs are
outside the checkout. Other workflows have separate configuration paths and
named overrides; this variable is not universal. Functional overrides include
`MAGFUNC_CLUSTER_FILE`, `MAGFUNC_ANNOTATION_FILE`, `MAGFUNC_ABUNDANCE_FILE`,
`MAGFUNC_METADATA_FILE` and `MAGFUNC_TAXONOMY_FILE`. Keep site paths in ignored
local configuration, not tracked YAML or public logs.

## Public/private boundary

The public repository supplies code and synthetic fixtures only. Inputs, metadata,
processed datasets and generated outputs are obtained separately or created
locally with authorization. The software license does not authorize data reuse.
The intended controlled-access model and public request contacts are established.
Component-specific permission/data-use conditions remain pending verification
before any provision. No public study accession is documented.
Preserve private provenance; no private Git history is imported into the public
repository. Review exact public staged files before publication.
