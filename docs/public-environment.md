# Public software environment

## Supported setup versus historical reproduction

`environment.yml` is the primary cross-language dependency declaration. It pins
R to 4.3 and Python to 3.11, with many unpinned R and Python packages. Recorded
primary per-MAG and historical NetMoss sessions report R 4.3.3 on Linux.
`renv.lock` records R 4.3.3 and 13 packages, chiefly a narrower beta-diversity
snapshot; it does not enumerate all workflows or replace Conda.

Create an environment with `conda env create -f environment.yml`, then activate
`glmm-hpc`. Do not update an environment used for historical reproduction in
place. Capture its actual package/session evidence first. Never generate a
purported historic lockfile from today's package resolution.

## Setup interruption recovery

If an SSH session disconnects during environment creation, reconnect and inspect
the environment before starting an analysis. Use a durable `tmux` or `screen`
session for permitted interactive setup, or the tracked
`hpc/update_glmm_env.pbs` compute-node job when updating an existing environment.
Follow the [Nebula update procedure](nebula-repository-update.md) for cluster
requirements; run-specific historical environments must not be updated in place.

```bash
conda run -n glmm-hpc R --version
bash scripts/validate_beta_diversity_environment.sh
Rscript scripts/netmoss/install_netmoss2.R --check-only
```

Use the R checks in the activated `glmm-hpc` environment and the validators
relevant to the intended workflow. A disconnected setup session does not establish
that the environment is usable; inspect the setup log and missing dependencies
before resubmitting scientific work. This is general recovery guidance, not a
record of a historical failure or evidence of exact-result reproduction.

## Dependency groups

| Workflow | Principal dependencies; full declarations remain authoritative |
|---|---|
| Preparation and shared tables | dplyr, tidyr, readr, lubridate, yaml, here, data.table |
| Mixed models | lme4, lmerTest, Matrix; the retained legacy fit script also invokes lme4. Additional declared packages do not prove that an execution path uses them. |
| Diversity | vegan, ape, permute, GUniFrac/phyloseq for available tree methods; Python pandas, numpy, scikit-bio for alpha comparison |
| Random Forest | randomForest, pROC and shared manipulation/plotting packages |
| Rendering | ggplot2, ggrepel, scales, patchwork, gridExtra, svglite, ragg; system fonts affect rendering |
| Functional workflows | data.table, mixed-model/ecology packages, jsonlite and shared plotting packages |
| NetMoss | NetMoss2, rsparcc, WGCNA, igraph and runtime/build dependencies in the installer and Conda file |
| Validation/provenance | testthat, digest, sessioninfo, yaml; Python standard library for release validation |

R base/recommended packages are supplied with R. `tibble` is used directly and
is normally installed through dplyr; this does not establish a historical version.
Shell runners require Bash; Git and Make support reproduction/maintenance.
OpenPBS or Slurm is needed only for its matching scheduler wrappers; rsync/SSH
are optional transfer utilities. No scheduler is needed for local smoke tests.
MMseqs/eggNOG outputs are upstream inputs, not tools invoked by these downstream
functional workflows; their generating software/database versions remain to verify.

## NetMoss installation

After Conda setup, install and check pinned packages:

```bash
Rscript scripts/netmoss/install_netmoss2.R
Rscript scripts/netmoss/install_netmoss2.R --check-only
```

The installer pins NetMoss2 to `79da905efccfc8bec9c08a4f447bf1ddb12223f6`
and r-sparcc to `4ac4c6c6e83219d0e970b2b8b28a9fef8c0eb30b`.
WeightedROC and plotROC are installed from CRAN if missing; exact historical
versions are not pinned by that script. Installation requires network access;
scientific jobs should use the prepared environment without live installations.


## Reproduction limits

Capture actual package versions, session information, seeds, configuration and
input hashes for each authorized run. The Conda declaration contains unpinned
packages; the narrower renv snapshot does not certify the historical environment.
Alpha fitting uses lmerTest/Satterthwaite when available and an explicitly labeled
normal-approximation lme4 fallback otherwise. Do not treat these inference branches
as interchangeable or alter installed packages simply to force a desired result.
Historical NetMoss settings remain separate in `config/netmoss_historical_publication.yml`;
current NetMoss source explicitly loads its required runtime packages. No research
outputs or historical runtime objects are supplied by this archive.

## Runtime validation and historical evidence

Check actual dependencies per workflow; an installed R executable alone does not
establish a complete runtime. The narrow renv snapshot, the Conda declaration and
a historical session record serve different purposes. Record the selected
executable, package versions and input/configuration hashes for each validation.
A fresh installation supports execution checks; it does not recreate an
unrecorded historical build.

Use separate validation environments when needed. Preserve existing scientific
runtimes, particularly the alpha inference branches described above. Run the
aggregate public test runner only in the clean tree with approved mappings;
the private tree also contains excluded research-dependent originals.
Keep logs and generated synthetic artifacts outside the release-source inventory.

## Verified local runtime scope

Release preparation checked the required local dependency groups for all 19
registered workflow branches against the assembled source. Six workflow
environment validators passed, including the pinned NetMoss revision check.
Package inventories, actual library search paths and session information are
retained privately. This closes the local dependency gate; it does not certify
every declared installation recipe or establish scientific reproducibility by
package availability alone.

| Runtime profile | Verified use and boundary |
|---|---|
| Preserved model environment | R 4.3.3, lme4 1.1-37, lmerTest 3.1-3 and Matrix 1.6-5 remain available for the accepted primary/Class 1 production evidence. An artificial mixed-model API check passed. This environment lacks svglite, so it is not certified for complete inline SVG output. |
| Default system R | R 4.3.3 with ecology, Random Forest, NetMoss, functional-analysis and graphics packages passed the workflow preflights. Its absence of lmerTest preserves the selected alpha normal-approximation branch. |
| Separate software-compatibility library | System R with isolated readr 2.2.0, patchwork 1.3.2 and lmerTest 3.2-1 supplies the combined model/rendering dependency path. All directly referenced R namespaces except conditional GUniFrac load; artificial table I/O, mixed-model fitting and SVG/PNG/PDF exports passed. This newer combination is not claimed to reproduce the preserved model environment or study results. |
| Python virtual environment | Invoke its own `bin/python`, or activate it before using the alpha comparison. numpy 2.5.1, pandas 3.0.3 and scikit-bio 0.7.3 load, including the alpha functions imported by the released script. The checked interpreter is Python 3.12; the declared Conda recipe requests 3.11 and was not freshly solved. Resolving the virtual-environment executable's symlink and invoking the system Python loses the virtual-environment package context. |

The isolated library is validation evidence, not a distributed runtime. Reproduce
its setup in a separate user-controlled library and retain its actual library
search paths. Prepending a validation library must preserve existing library paths;
replacing the user-library setting can hide already installed packages. Do not
enable the lmerTest-containing compatibility profile globally for a selected alpha
fallback reproduction: doing so changes the code's inference branch. Existing
production installations were not updated.

The current no-tree configuration is supported. GUniFrac/tree-enabled execution
remains unverified and requires an authorized matching tree and its dependencies.
Live OpenPBS/Slurm execution was not tested; schedulers are optional local-workflow
alternatives. A fresh complete `glmm-hpc` Conda solve remains unverified. `renv`
and `snakecase` are not loadable in the checked system runtime; neither is invoked
by the retained analytical execution paths checked here, so no unnecessary
installation was performed. The legacy report interface remains omitted.

The diversity/per-MAG and RF reporting issues are resolved in approved replacements
without changing scientific calculations. Their originals remain private. See the
[current release gates](public-release-checklist.md#current-evidence-and-remaining-release-gates).
Dependency verification does not authorize freezing or publication.
