# MAG identifiers and input interpretation

## Required keys

Join observations from `gut.csv` to `sample_info.csv` by the complete `sampleID`.
Raw metadata `group` supplies `Subject`; it is not a treatment group. A MAG column's
complete name is its stable feature key. Join it exactly to taxonomy `GDTB.csv`
column `Genome`. Keep the supplied spelling of `GDTB.csv`.

## Construction provenance versus observation metadata

Embedded `B`, `C` and `D` denote the specimen type associated with the retained
representative bin: fecal, saliva and skin, respectively. They do not override
later gut-specific selection represented by membership in `gut.csv`.
Embedded dates and altitude-like components describe construction provenance and
must not supply a sample's date, altitude, participant or specimen type. Use the
row-level metadata join for observation covariates. Do not infer undocumented
meanings for `EVR`, `QML`, `strict`, `permissive` or `orig`.

## Functional joins and scale

For functional annotation, propagate representative-protein annotations to
members, then recover MAG ownership from each member protein identifier. Retain
the full MAG prefix ending in `.strict`, `.permissive` or `.orig` and validate it
against the exact taxonomy and selected abundance columns. A `NODE` contig-length
component is not amino-acid protein length.

Do not infer upstream abundance normalization, MAG assembly/quality thresholds,
MMseqs identity/coverage settings or annotation database versions from identifiers.
These require separately obtained provenance. Values are described as
relative-abundance-like; this archive does not establish their upstream denominator.
Primary functional workflows use supplied values without row renormalization;
explicit normalized sensitivity branches do not replace that primary calculation.

See [input schemas and access requirements](../release-data-requirements.md).
No real participant identifiers, sample examples or biological measurements are
provided in this technical guide.
