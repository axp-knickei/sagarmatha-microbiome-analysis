# sagarmatha-microbiome-analysis

This journal-associated analytical and reproducibility resource for the Sagarmatha
(Mount Everest) microbiome study provides reviewed analytical code, configuration,
environments, validation materials and supporting technical documentation. The active
primary workflow is a per-MAG LMM; the older single-model GLMM template remains
available for compatibility. Diversity, subject-aware Random Forest, exploratory
NetMoss and functional workflows have separate input and validation contracts.

## Inputs and execution

Study inputs, participant metadata, processed biological datasets and generated
scientific outputs are not distributed. Obtain authorized inputs separately and
follow the [input contract](docs/release-data-requirements.md). The original window
CSV is excluded; its six predefined definitions and rationale remain in the
generator. Full regeneration requires authorized metadata matching the original
processed selection; see [the instructions](docs/public-workflow.md#window-csv-regeneration).

Create the Conda environment from `environment.yml` and use the
[environment guide](docs/public-environment.md). Run scripts from the repository
root. Follow the [workflow guide](docs/public-workflow.md) and
[documentation index](docs/README.md) for execution and validation.

## Repository structure

| Location | Role |
|---|---|
| [scripts](scripts/README.md), [analysis](analysis/README.md) | Analytical code, helpers and rendering code |
| [config](config/README.md) | Required workflow configuration and scientific settings |
| [hpc](hpc/README.md) | Optional site-specific OpenPBS and Slurm launchers |
| `tests/` and workflow-local tests | Validation code; some integration checks require authorized external inputs |
| [Synthetic beta fixtures](data/test_beta_diversity/README.md) | Fictional smoke-test inputs |
| `data/raw/`, `data/processed/`, `results/` | Empty placeholders; authorized local inputs and generated outputs are excluded |
| `docs/` | Public technical methods and release guidance |

## Validation and interpretation

Use the documented synthetic tests and input/environment validators before an
authorized run. Preserve diagnostics, status flags, restricted resampling designs,
multiple-testing families and interpretation limits. A synthetic smoke test or
syntax check does not establish complete study reproduction. The historical
`make report` target is omitted from the public Makefile; it is not a supported
release command. The private original Makefile remains unchanged.

## Release, citation and reuse

The private research/development repository and its history remain private. A
separate public repository will be assembled only from the final approved manifest
and start with fresh history. Its approved name is `sagarmatha-microbiome-analysis`.
The GitHub owner and exact URL will be recorded after repository creation and
verification; no clone URL, release URL or badge is claimed yet. See
[release policy](docs/public-release-policy.md)
and [open gates](docs/public-release-checklist.md).

Software license: [MIT](LICENSE). The approved intended public data-access policy
is controlled access / available upon reasonable request:

> The analytical code and supporting technical documentation are publicly available
> through the `sagarmatha-microbiome-analysis` repository. The underlying study data,
> metadata, and restricted processed data are not included in the public code
> archive. These materials may be made available upon reasonable request to the
> authors, subject to applicable permissions and data-use conditions.

This statement is for use at publication; the code is currently being prepared
for public release. Alex Prima will receive and coordinate requests at
[alexprima@student.ub.ac.id](mailto:alexprima@student.ub.ac.id) or
[alex.prima@uq.net.au](mailto:alex.prima@uq.net.au).
Component-specific permissions and reuse conditions remain pending verification
before any data provision.
Request coordination does not imply unrestricted redistribution authority or
automatic request approval. See the
[data-access statement and conditions](docs/release-data-requirements.md#intended-access-policy--confirmed-contacts-and-conditional-permissions).
No public study accession is established.

[Software citation](CITATION.cff) records the approved name, scope and authorship.
Repository URL, release date and Zenodo DOI remain pending verified identifiers.
No release has been published.

### Third-party palette material

Palette values in the diversity/per-MAG and subject-aware Random Forest manuscript
renderers are credited to the Ordnance Survey
[GeoDataViz Toolkit](https://github.com/OrdnanceSurvey/GeoDataViz-Toolkit), palette
sheet version 0.9.1, and retain the upstream
[Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/)
terms. Contains public sector information licensed under the Open Government
Licence v3.0. MIT applies to original project code and documentation; it does not
replace the terms for this third-party palette material or external software and
data.
