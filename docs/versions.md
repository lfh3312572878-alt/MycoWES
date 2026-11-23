# Version and reference identity

The software reports `MycoWES 1.0.0`. Reference resources retain version 0.28.
The engine records its source digest, actual input/reference identities, preset,
quality thresholds, risk/MAPQ policies and execution settings in every run.

The manuscript preset uses 67-target competition and 67-target anchor assessment.
The default standard and fast profiles use 18-target competition and the shared
67-target anchor pool. All retain the complete competitor set. PRIMARY and STRICT
do not select reference multiplicity.

Source scripts under `reproduce/historical/` retain their original scientific
design and commit labels. Root-parameter adaptations are listed in
`reproduce/SOURCE_INDEX.tsv`. Current figure calculations use the separate
source-data archive. Software-version naming does not replace a recorded analysis
identity.
