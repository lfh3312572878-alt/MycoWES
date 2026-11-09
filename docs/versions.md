# Versions and analysis profiles

`VERSION` identifies the software build; the Git commit identifies the exact repository contents. The release candidate is based on the MycoWES 0.3.0-rc1 source snapshot at ae1f1d4b9f0e33baf05944540b8c5b9bfd59268a. Local commit 573ec03 records the explicit snapshot import. Local and remote commit histories are identified separately.

The runtime engine is fixed at [remote commit 802200ba32e01fbb0a8e86ccb2fd304d97b31ca8](https://github.com/lfh3312572878-alt/MycoWES/commit/802200ba32e01fbb0a8e86ccb2fd304d97b31ca8) on `remediation/20261109`. Subsequent documentation-only commits do not change the runtime engine or default settings.

Historical source hashes are in `reproduce/SOURCE_INDEX.tsv`. Labels such as 77c4e85 inside preserved scripts identify their execution lineage, separately from this repository's Git history. The synthetic example and saved-summary validator validate the current software and supplied saved data respectively.

| Context | Competitive target | Length bp | MAPQ | NM rate |
| --- | --- | ---: | ---: | ---: |
| normal default | one representative per species | 100 | 20 | 0.08 |
| simple matching recorded RELAX2 choice | one representative per species | 60 | 5 | 0.12 |
| hard matching multi-target baseline choice | complete target set | 100 | 20 | 0.08 |

Effective length is min(configured length, nominal PE length). These presets select reference and quality settings; per-experiment commands and input preparation are recorded in the analysis source map.

`paper_mapping.tsv` lists the included historical sources, current engine entries and software validation, with their respective configurations, inputs and output locations. Historical and current validation outputs remain distinct.

In that table, historical output directories without an absolute prefix are under `/data/h011/project/analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1`; `E` denotes its `Fig3_Fig4_TXT_COMPLETION_EXEC_20261010_v2` subdirectory. These are provenance locators, not public download URLs or instructions to write into the original project. Included code paths are relative to this repository.

This revision implements input validation, effective parameter delivery and run/output integrity. Default references, ordered exact-sequence grouping and scientific grade boundaries are retained. Historical analyses and their saved outputs are unchanged.

The release candidate remains on a private branch. The reference dataset v0.28 is a Zenodo draft with reserved DOI 10.5281/zenodo.23776490, not a software DOI. Copyright license, publication and third-party access are managed separately by the authors.
