# Public v1.0.0 release policy and decisions

Current scope: 299 approved inclusions, 293 exclusions and 17 replacement pairs.
Downstream workflow evidence is reconciled; approved replacements address the
previously identified renderer content holds.
See the [current release gates](public-release-checklist.md). Dated entries below
preserve earlier decisions and validation boundaries.

Current public request contacts: Alex Prima at
[alexprima@student.ub.ac.id](mailto:alexprima@student.ub.ac.id) or
[alex.prima@uq.net.au](mailto:alex.prima@uq.net.au), supplied by the owner on
2026-10-15. Earlier dated statements that contact details were pending are
historical; component-specific permissions remain conditional before provision.

## 2026-10-01 — separate public repository and code-only scope

The owner approved a separate public code repository with fresh initial history.
The private research repository, full Git history, existing outputs, drafts and
excluded documentation remain preserved. Local staging is authorized and assembled for inspection; no public Git repository
has been initialized.

Include analytical code, required configuration, environments, tests, synthetic
fixtures and technical documentation needed to understand, run and validate the
workflows. Exclude study inputs, metadata, processed biological data, generated
scientific outputs, manuscript drafts, participant-linked reports, scientific
handoffs, planning/status records, prompts and internal research notes.

Mixed originals are excluded unchanged; necessary execution and method contracts
are preserved in the public guides indexed in [docs/README.md](README.md). Real
cohort identity/count assertions remain in a private NetMoss test file; public
NetMoss tests retain its hand-constructed resampling and configuration checks.
Unnecessary hard-coded participant facts were removed from generated-report
wording without changing analysis calculations. The RF private-snapshot fallback
is unset; supply authorized processed data explicitly.

## Window definitions and external inputs

The original window CSV is excluded. Six predefined window definitions and their
rationale remain unchanged in the generator. Regeneration requires separately
obtained authorized metadata matching the original processed selection. The public
archive provides neither study metadata nor processed biological data.

## Provenance and approval

Before assembly, approve every path in the public manifest. Record the private
source revision and a private mapping of local changes, then hash the actual
approved public files. A private commit identifier alone cannot describe a dirty
working tree. Publish only a reviewed provenance summary and included-file hashes;
do not copy private Git objects, branches, tags, logs or excluded evidence.

The approved public repository name is `sagarmatha-microbiome-analysis`.
Citation omits `repository-code` until the GitHub owner and actual URL are verified. The setup helper requires an explicit approved URL when
cloning; it can still operate within an existing checkout. Retained files no
longer present the private repository as the public release destination.

File-scope approval does not certify external data access/reuse, complete
scientific reruns, the frozen snapshot, or publication. See [release gates](public-release-checklist.md).

## 2026-10-02 — public rendering and validation replacements

The owner approved replacing the twelve unresolved result-bearing source files
while preserving private originals. Review found the related NetMoss synthesis
renderer also embeds findings; it receives the same treatment. All thirteen
originals remain unchanged and excluded. Public replacements consume separately
authorized outputs, calculate counts/labels and check source-file agreement without
requiring a particular biological finding. Models, selection rules, thresholds,
window definitions and exploratory boundaries are retained. See the
[rendering guide](public-rendering.md) and
[relative-path/hash correspondence map](releases/v1.0.0-public-replacements.csv).
The nine per-MAG renderers use stable output stems and optional external path
variables. The functional renderer describes reported evidence without inserting
predetermined scientific interpretations; its configured contrast is preserved.

## Historical verification — original replacement batch

At that checkpoint the manifest had 299 approved inclusions, 291 approved exclusions and no
unresolved file-scope paths. Citation identity and software scope are approved;
future URL, DOI and release-date metadata are publication dependencies. Inclusion applies to reviewed current bytes,
not an older Git commit. The new synthetic test drivers exercise every replacement
and reject deliberate cross-file inconsistencies. All thirteen originals retained
their pre-review hashes. No private Git history or research output is included.
Full scientific reruns and final frozen snapshot verification remain open.

Relevant synthetic NetMoss and functional helper tests passed. R files, shell/
scheduler files and documented Bash templates passed parsing, and Make dry runs
and diff whitespace checks passed. No authorized full scientific workflow was run
and no clean public repository was assembled or initialized.

## 2026-10-02 — approved public repository identity

The owner selected `sagarmatha-microbiome-analysis` as the separate journal-associated
analytical and reproducibility resource for the Sagarmatha (Mount Everest)
microbiome study. The software citation, README, release notes, descriptive project
labels, example checkout paths and setup directory default use this identity.
No GitHub owner, repository URL, DOI, release URL or publication record is invented.
CITATION.cff is approved for inclusion with URL-dependent fields omitted until
verified. This closes file-scope decisions; it does not freeze current file bytes
or authorize assembly, Git initialization or publication.

Scientifically meaningful Everest/BlackICE study terminology, input environment
variable names, scheduler identifiers and log filenames remain unchanged. Log
messages used by parsers are execution contracts, not private repository links.
The manifest and replacement map retain relative excluded-original filenames
for provenance; these contain neither private repository addresses nor Git history.

## 2026-10-02 — approved public data-access statement and request contact

The owner approved controlled access / available upon reasonable request as the
intended public policy. [The input contract](release-data-requirements.md) and
README contain the approved statement, explicitly staged for publication because
the code repository is not public yet. Study data, metadata and restricted
processed data remain outside the code archive.

Alex Prima can receive and coordinate requests, by explicit owner confirmation.
This coordination role is separate from authority to release individual data
components. Specific contact details remain pending until explicitly supplied;
no email, accession or data-use agreement is invented. Provision may depend on
applicable collaborator, institutional, ethical, ownership or data-use conditions;
no automatic request approval or unrestricted redistribution authority is claimed.
Component-specific permission and reuse conditions remain pending verification.

File scope remains complete: 299 inclusions, 291 exclusions. The legacy report
interface is resolved by the owner-approved public-only omission; scientific validation, contact details,
permission conditions, exact snapshot verification and publication gates remain
separate. No assembly or publication is authorized by this policy.

## 2026-10-02 — omit obsolete public report command

The owner approved removing `make report` solely from the public-release Makefile.
The private original remains unchanged and excluded; `Makefile.public` is approved
for assembly as `Makefile`. The manifest and replacement map record this destination.
No report.Rmd is created. All retained recipes, scientific scripts, thresholds and
outputs are unchanged. This closes the legacy-interface gate; it does not establish
full workflow validation or authorize assembly, Git initialization or publication.

## 2026-10-03 — synthetic public Random Forest tests and staging refresh

The owner approved replacing the research-input-dependent RF test with a fictional
fixture test while preserving the private original and scientific implementation.
The manifest maps the public replacement to the canonical test path; the provenance
map records original/replacement hashes. A public runner executes synthetic/unit
checks only. Two obsolete private-name backup patterns were removed from .gitignore.
An excluded historical filename remains in manifest provenance; it reveals no
private repository location, URL or history. Local staging assembly/refresh is
now authorized for inspection only; no Git initialization or publication is authorized.
Scope: 299 approved inclusions, 291 exclusions, zero unresolved path decisions.
Production scientific validation and final frozen-snapshot approval remain open.

The refreshed 299-file staging checkpoint passed all nine public synthetic/unit
suites. The private RF test and scientific implementation were hash-verified as
unchanged. Staging inventory hashes remain a local validation checkpoint rather
than an approved frozen v1.0.0 release. No production study-data run was performed.

## 2026-10-04 — accepted Class 1 numerical reproducibility limitation

The owner approved **passed with numerically negligible differences** for the
Class 1 production validation. The failed strict historical comparison is retained
separately from this acceptance; no tolerance, analytical setting, historical
result or input was changed. See the [Class 1 limitation](class1-methods.md#accepted-production-reproducibility-limitation)
for deterministic rerun evidence and the unproven historical cause. Detailed
research-level diagnostics remain private. The existing local staging checkpoint
is unchanged; documentation updates require a later authorized refresh and final
snapshot validation. Other production validations and publication gates remain open.

## 2026-10-15 — accelerated technical completion and evidence retention

The owner authorized documentation/provenance consolidation, targeted authorized-
input validation and reconstruction of inspection-only staging. File scope
remains 299 approved inclusions and 291 approved exclusions. This authorization
does not extend to a frozen candidate, Git initialization, commits, tags, remotes,
push, public visibility or archival publication.

Distinguish directly inspectable artifacts, retained validation reports and
owner-accepted limitations. Earlier temporary staging and some detailed diagnostic
files are unavailable; their successful historical checks do not certify current
staged bytes. Recheck replacements, mappings, disclosure and synthetic tests in
the reconstructed clean tree. Preserve durable detailed evidence privately and
publish only technical methodology and public-safe validation information.
Controlled access, conditional permission authority and pending contact details
remain unchanged. No full upstream-method or exact historical-build claim is made.

The reconstructed inspection checkpoint contains exactly 299 approved files and
passed all nine public synthetic/unit suites. Replacement mappings/hashes, local
links and anchors, syntax, Make dry runs, citation schema, license and disclosure
checks passed. Generated synthetic test artifacts were archived privately and
removed from the source checkpoint. Scientific input files remain excluded.
Runtime availability was recorded by workflow; the broad all-workflow dependency
and scientific validation gates remain open. Final hashes/owner approval and all
GitHub/Zenodo actions are deferred.

## 2026-10-15 — dependency closure and retained-renderer content hold

Required local dependency evidence now covers the registered workflow branches,
with preserved production runtimes and a separate compatibility library; see
[verified runtime scope](public-environment.md#verified-local-runtime-scope).
Conditional tree methods, fresh complete environment installation and live
scheduler execution are not claimed as verified.

The retained diversity/per-MAG synthesis builder was found to contain predetermined
research findings, and its isolated execution failed while indexing missing
optional supplementary evidence. Its original is preserved. An additional
public-safe replacement needs owner approval; neither scientific calculations nor
the approved manifest dispositions were changed. The renderer gate is reopened
and full-workflow closure remains pending. New detailed production comparisons
and rendered research outputs are retained only in the ignored private audit.

The 299-file staging checkpoint and its source-only archive remain unchanged.
They are inspection evidence, not disclosure-cleared or frozen release files.
Approved documentation changes require a later staging refresh and final checks.


## 2026-10-15 — authorized diversity replacement and current gates

The owner approved an additional public-safe diversity/per-MAG renderer. Preserve
its private original and install the replacement at the canonical destination;
record both hashes and the mapping. Counts are 299 approved inclusions and 292
exclusions. Captions now derive reporting from supplied outputs; scientific models,
thresholds, selection rules and numerical calculations remain unchanged. Missing
optional supplementary directories receive explicit status; required input failures
remain errors. Detailed numerical comparisons and research outputs stay private.

The scoped scientific-workflow gate is satisfied by the reconciled retained and
targeted evidence; accepted limitations are preserved. A focused scan identified
unconditional interpretation in an additional RF methods note. Public-renderer
clearance remains open for that specific source pending owner approval. This is
not approval of a frozen candidate, Git operations or publication. Older entries
above remain dated records and do not override these current dispositions.


## 2026-10-15 — approved RF public reporting replacement

The owner approved preserving the RF feature-stability original privately and
using a public replacement at the canonical destination. Both performance judgments
are replaced by neutral supplied metrics and explicit missingness. Scientific
calculations, feature rankings, input contracts and statistical procedures remain
unchanged. The negative-value preservation check no longer requires a particular
observed sign outcome. Synthetic and authorized-output checks passed; a focused
review found no further predetermined renderer statement requiring remediation.
The scope is 299 inclusions, 293 exclusions and 17 provenance pairs. Current
inspection checks and later frozen approval remain distinct from publication.

The refreshed inspection tree passed all nine public suites and the final source,
link, syntax, metadata and disclosure checks. Renderer gate #5 is complete;
scientific-workflow gate #20 remains complete within its accepted scope and
limitations. Owner approval and exact final-snapshot verification remain separate.

## 2026-10-15 — public data-request email addresses supplied

The owner supplied `alexprima@student.ub.ac.id` or `alex.prima@uq.net.au` as public
request contacts for Alex Prima. Both are documented without assigning a preferred
address. CITATION.cff uses the first listed address in its single email field.
Email domains do not establish a new affiliation, data ownership or authority to
redistribute study data. The intended controlled-access policy is unchanged;
applicable component-specific permissions and reuse conditions must be verified
before data provision. No mail delivery/ownership test was performed and no
message was sent. The contact requirement is satisfied, but combined gate #12
remains open for its conditional permission/reuse requirements.

## 2026-10-16 — approved palette attribution and assembly wording

The owner approved a README notice identifying the Ordnance Survey GeoDataViz
Toolkit palette sheet v0.9.1 and its Open Government Licence v3.0 terms. MIT remains
the license for original project code and documentation; third-party palette and
external software/data terms are not replaced. No palette values or analytical
calculations change.

The release procedure now instructs assembly by each included row's explicit
`public_destination_path` and verification of all 17 replacement mappings. Its
examples include the public Makefile, synthetic RF test and canonical destinations
for the diversity/per-MAG and RF feature-stability renderers. No manifest
disposition, destination, source mapping or approved file count changes. This
approval permits documentation changes and inspection checks only; final snapshot
approval, freezing, Git operations and publication remain separate.

## 2026-10-16 — final code-only permission review

Public contacts and the owner-approved conditional reasonable-request statement
are established. Underlying study data remain excluded. Component-specific
sharing authorities and reuse conditions are unverified and must be checked
before any future data provision; they are not a prerequisite for describing
conditional access in this code-only release. No permission is invented or
approved by this review. The code's existing MIT notice and approved GeoDataViz
attribution remain in place. Final owner approval should confirm the right to
publish the included original work and retain applicable third-party terms.

Checklist #12 combines completed public-contact documentation with future
before-provision verification. Its checkbox remains open without treating the
excluded data's pending permissions as a code-publication blocker. #9 and #19
explicitly require a later frozen snapshot; current candidate technical checks
do not establish that state. No publication or snapshot approval is implied.

The final technical review removed study-specific retained-taxon counts, an
observed timing claim and an undated cluster-status assertion from the public
NetMoss methods guide, preserving a byte-identical copy privately. The guide now
asks users to record counts and resource needs from supplied inputs. Prevalence,
module-size, resampling and interpretation rules are unchanged. Launcher help
messages and two source comments now point to retained public instructions;
parsed R expressions and all job commands/settings are preserved. One public
replacement hash is refreshed for its corrected usage comment, without changing
any source/destination mapping or file-scope disposition. These documentation
fixes are part of the authorized final technical review, not new analyses.
