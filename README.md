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

The private research/development repository and its history remain private. The
separate release repository is
[axp-knickei/sagarmatha-microbiome-analysis](https://github.com/axp-knickei/sagarmatha-microbiome-analysis).
The release source uses fresh Git history and contains only reviewed code and
technical materials. Private research history and excluded materials are not
part of this resource. Verify version-specific releases and citation identifiers
at the repository; the repository URL alone does not establish public access or
an archived release. See
[release policy](docs/public-release-policy.md)
and [open gates](docs/public-release-checklist.md).

Software license: [MIT](LICENSE). The approved intended public data-access policy
is controlled access / available upon reasonable request:

> The analytical code and supporting technical documentation are publicly available
> through the `sagarmatha-microbiome-analysis` repository. The underlying study data,
> metadata, and restricted processed data are not included in the public code
> archive. These materials may be made available upon reasonable request to the
> authors, subject to applicable permissions and data-use conditions.

This is the approved statement for public distribution of the code; it does not
itself establish repository visibility or the existence of an archived release.
Alex Prima receives and coordinates requests at
[alexprima@student.ub.ac.id](mailto:alexprima@student.ub.ac.id) or
[alex.prima@uq.net.au](mailto:alex.prima@uq.net.au).
Component-specific permissions and reuse conditions remain pending verification
before any data provision.
Request coordination does not imply unrestricted redistribution authority or
automatic request approval. See the
[data-access statement and conditions](docs/release-data-requirements.md#intended-access-policy--confirmed-contacts-and-conditional-permissions).
No public study accession is supplied in this code archive.

[Software citation](CITATION.cff) records the approved name, scope and authorship.
Version `v1.0.0` was released on 2026-10-20 and is available from the
[GitHub Release](https://github.com/axp-knickei/sagarmatha-microbiome-analysis/releases/tag/v1.0.0)
and [Zenodo record](https://zenodo.org/records/23497840).

Recommended citation:

Prima, A. (2026). *sagarmatha-microbiome-analysis: analytical and reproducibility
code for the Sagarmatha (Mount Everest) microbiome study* (v1.0.0) [Software].
Zenodo. https://doi.org/10.5281/zenodo.23497840

Use the [version DOI](https://doi.org/10.5281/zenodo.23497840) to cite this exact
release. The [concept DOI](https://doi.org/10.5281/zenodo.23497839) represents all
versions.

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
