# Version and reference identity

The software reports `MycoWES 1.0.0`. Reference resources retain version 0.28.
The engine records its source digest, actual input/reference identities, preset,
quality thresholds, risk/MAPQ policies and execution settings in every run.

Two canonical modes are available. Standard (the default, also accepted as
manuscript) uses 67-target competition/anchor, 100 bp / MAPQ20 / NM0.08 and
gap10/LC180. Fast includes the former relax2 configuration: 18-target
competition, the shared 67-target anchor pool, 60 bp / MAPQ5 / NM0.12,
native MAPQ and legacy risk rules. Relax2 and legacy simple are compatibility
aliases for fast; legacy normal/hard select standard. Both modes retain the
complete competitor set. PRIMARY and STRICT do not select reference multiplicity.

Source scripts under `reproduce/historical/` retain their original scientific
design and commit labels. Root-parameter adaptations are listed in
`reproduce/SOURCE_INDEX.tsv`. Current figure calculations use the separate
source-data archive. Software-version naming does not replace a recorded analysis
identity.
