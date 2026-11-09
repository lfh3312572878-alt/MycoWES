# Versions and manuscript correspondence

`VERSION` and the local Git commit identify this revision. Starting files are the 20-file MycoWES 0.3.0-rc1 snapshot at repository commit ae1f1d4b9f0e33baf05944540b8c5b9bfd59268a, confirmed as the authenticated remote main baseline on 2026-11-09. The source Git objects were unavailable at local import: 573ec03 is an explicit file import, not a reconstructed parent/history. The local remediation history and any remote review branch must be identified separately; no synchronization is inferred from their names.

The reviewed corrective engine is fixed at [remote commit 802200ba32e01fbb0a8e86ccb2fd304d97b31ca8](https://github.com/lfh3312572878-alt/MycoWES/commit/802200ba32e01fbb0a8e86ccb2fd304d97b31ca8) on the private `remediation/20261109` branch. The remote main baseline remains ae1f1d4b9f0e33baf05944540b8c5b9bfd59268a; the review branch is not a published release.

Historical source hashes are in `reproduce/SOURCE_INDEX.tsv`. A 77c4e85 label, including full commit strings inside preserved scripts, is an execution-lineage record, not evidence that its Git object has been recovered or is an ancestor here. New synthetic validation and saved-summary code belong to this revision, not the original paper executable.

| Context | Competitive target | Length bp | MAPQ | NM rate |
| --- | --- | ---: | ---: | ---: |
| normal default | one representative per species | 100 | 20 | 0.08 |
| simple matching recorded RELAX2 choice | one representative per species | 60 | 5 | 0.12 |
| hard matching multi-target baseline choice | complete target set | 100 | 20 | 0.08 |

Effective length is min(configured length, nominal PE length). Paper configurations must come from each experiment's commands, configs and outputs, not this table alone. Fig.2B combines target recovery at 60/5/0.12 with non-target routing at 100/20/0.08; it is not a single-configuration ablation.

`paper_mapping.tsv` separates historical evidence, current engine behavior and new validation. Final Fig.2C, Fig.4A–E and some composite plot inputs still require author-provided bindings. No old table has been relabelled as their source.

In that table, historical output directories without an absolute prefix are under `/data/h011/project/analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1`; `E` denotes its `Fig3_Fig4_TXT_COMPLETION_EXEC_20261010_v2` subdirectory. These are provenance locators, not public download URLs or instructions to write into the original project. Delivered code paths are relative to this repository. External material unavailable in the release is explicitly noted.

This revision addresses input validation, effective parameter delivery and run/output integrity. Default references, ordered exact-sequence grouping and scientific grade boundaries are retained. Historical analyses were not rerun or silently updated. Missing historical input-integrity evidence or source/output bindings limit conclusions about their impact.

The repository remains private with no published software release. Authenticated Zenodo inspection on 2026-11-09 confirmed the reference dataset version 0.28 as an unpublished Draft: DOI 10.5281/zenodo.23776490 is reserved and will be registered only on publication. Copyright license, review access for third parties, a published software DOI and final manuscript-to-release acceptance require author confirmation. The reserved reference DOI is not a software DOI.
