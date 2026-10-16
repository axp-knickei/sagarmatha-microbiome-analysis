# Public v1.0.0 release gates

Status: local inspection-only staging and refresh are authorized. No Git
initialization, commit, tag or publication is authorized. The historical private `RELEASE_CHECKLIST.md` remains
unchanged as a research preparation record and is excluded from the public archive.

Current inspection status (2026-10-15): 299 files are approved for inclusion and
293 for exclusion, with 17 original/replacement mappings. The RF feature-stability
replacement passed synthetic and authorized-output checks and all nine public suites
in refreshed inspection staging. The local dependency and scoped scientific-workflow
gates remain supported with documented
limitations. Renderer gate #5 is complete for the approved public source mappings.
Inspection staging is not a frozen or owner-approved release candidate.

## Scope and disclosure

- [x] Owner selected a separate clean public repository; private history preserved.
- [x] Code-only boundary applied to research-bearing documentation; mixed originals preserved privately and essential instructions retained in public guides.
- [x] Approve `sagarmatha-microbiome-analysis` as public repository identity; finalize URL-independent software citation and setup directory name.
- [ ] After repository creation, verify GitHub owner/URL and add real repository links; add version DOI and release date only after publication.
- [x] Replace result-bearing renderers/validators with input-derived public versions; preserve originals and scientific calculations (see public-rendering.md and replacement map).
- [x] Complete path-level inclusion/exclusion decisions, including CITATION.cff (299 inclusions; 293 exclusions; no unresolved file-scope rows).
- [x] Verify the local staging checkpoint against approved source/destination mappings; remove synthetic runtime artifacts from the staged source inventory.
- [x] Rebuild the inspection-only staging checkpoint from approved mappings and revalidate file scope, hashes, disclosure, dependency availability, links and public tests; remaining all-workflow dependency validation is tracked separately below.
- [ ] Reverify the exact final frozen snapshot and exclusion boundary after all approved changes.
- [x] Owner approved the intended controlled-access / reasonable-request public statement; study data, metadata and restricted processed data remain excluded.
- [x] Establish Alex Prima as request contact for receiving/coordinating requests, separately from permission authority.
- [ ] Record explicitly supplied public contact details; verify applicable component-specific permission authorities and data/annotation reuse conditions before data provision.

The public-contact portion of #12 is complete: Alex Prima receives requests at
[alexprima@student.ub.ac.id](mailto:alexprima@student.ub.ac.id) or
[alex.prima@uq.net.au](mailto:alex.prima@uq.net.au). Component-specific permission
authorities and reuse conditions remain conditional and must be verified before
data provision; the combined checklist item remains unchecked.

For this code-only release, the confirmed contacts and conditional-access statement
satisfy the public data-availability documentation requirement. No study data are
being provided. The remaining component-specific permission checks apply before
future data provision; they are not a blocker to publishing the approved code.
Keeping #12 unchecked does not imply that excluded data must be made distributable
before the code release, or that data-sharing permission has been granted.

## Reproducibility and frozen snapshot

- [x] Owner approved omitting the obsolete `make report` target from the public Makefile; private original preserved, public replacement and destination recorded.
- [x] Local staging installed Makefile.public as Makefile; all 11 retained target dry runs match the private original.
- [x] Prepare synthetic RF test replacement; preserve the research-dependent private test and scientific implementation, and map the public replacement to the canonical test path.
- [x] Remove obsolete private-name backup patterns from public .gitignore; retain only excluded-relative-filename provenance, without private location/history.
- [x] Validate all required execution dependencies in the assembled clean tree.
- [x] Run all nine public synthetic/unit suites in the local staged tree, including synthetic RF; validate links, syntax, Make dry runs and replacement provenance.
- [ ] Repeat required checks on the exact final frozen snapshot after remaining release changes.
- [x] Validate documented full workflows/window regeneration with separately authorized inputs where required; record limitations honestly.
- [x] Owner accepted Class 1 production validation as passed with numerically negligible differences; retain the documented failed strict historical comparison and unchanged tolerance, explicitly noting that detailed temporary comparison artifacts are unavailable, with the [reproducibility limitation](class1-methods.md#accepted-production-reproducibility-limitation) documented. This does not complete the full-workflow gate above.
- [x] Record the local staging checkpoint inventory, approved source/destination mappings and SHA-256 hashes outside the staged source tree.
- [ ] Record final frozen-file hashes and public-safe provenance after release-candidate approval.
- [ ] Approve the complete frozen snapshot and remaining disclosure checks.

## Publication — separate explicit authorization required

- [ ] Create and review release commits in the separate public repository.
- [ ] Authorize public visibility, v1.0.0 tagging, GitHub Release and Zenodo integration/archive.
- [ ] Verify DOI/version/archive contents and update citation metadata.

## Current evidence and remaining release gates

Dependency gate #17 covers the required local paths with the current no-tree
configuration. Checks used staged source, workflow-specific package loading,
environment validators and artificial model/rendering API checks, with logs and
runtime inventories retained privately. See [runtime scope](public-environment.md#verified-local-runtime-scope).
It does not certify a fresh Conda solve, optional tree-enabled execution or live
scheduler jobs, and does not establish historical result identity.

The diversity/per-MAG and RF feature-stability replacements preserve their private
originals and occupy canonical public destinations. The RF replacement reports
supplied metrics neutrally and labels missing values unavailable without inferring
them. Candidate eligibility wording derives from supplied rows. Its negative-value
check verifies preservation rather than requiring a particular observed finding.
Targeted synthetic and authorized-output checks and all nine public suites passed
in the refreshed inspection tree. Required-input failure cases, missing performance
summaries and numeric precision are covered. No additional renderer content hold
was found by the focused review; gate #5 is complete.

Gate #20 is complete for the approved downstream code-release scope. The private
19-workflow register reconciles accepted production reports, new targeted alpha/
functional sensitivity checks, window regeneration, functional-community checks
and successful authorized-output rendering. Optional exploratory branches are not
promoted to fully reproduced study evidence. Upstream annotation provenance,
missing older temporary bundles, Class 1 numerical differences and EC diagnostic
limits remain disclosed; closure does not certify every model's inference or
full upstream reconstruction. Final manuscript approval remains separate.

There are currently 18 checked and 9 unchecked gates. No unresolved file-scope
decisions remain. Public contact details are recorded; component-specific permission
terms before data provision, final snapshot checks, owner approval
and publication are separate requirements.
