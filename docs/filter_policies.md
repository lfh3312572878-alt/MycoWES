# MAPQ and risk policies

The default quality gates are 100 aligned bases / MAPQ 20 / NM rate 0.08. Defaults are `--mapq-policy native` and `--risk-policy legacy`. Policy selection is recorded with each run. Adapter/DF2 exclusions, deduplication, anchor retention (0.90) and grading remain separate.

## Multi-assembly MAPQ

Native MAPQ is never overwritten. `assembly-ties-only` can exempt a low-MAPQ mate only when its reported alternatives all belong to one target species, the primary is highest-scoring, at least two assemblies share that score and each assembly has one reported locus.

`assembly-top-ties-gap10` uses the same-species top-score assembly tie, allows lower-scoring same-species loci and requires every reported other-species or competitor hit to have an alignment-score deficit of at least 10. Both primary mates and their unique best-target species must agree. Supplementary hits, missing identity/score and reporting-cap cases (at least 101 hits) cannot qualify. Aligned-length, NM and paired-support gates still apply.

## Risk decisions

| Policy | Treatment of existing annotations |
| --- | --- |
| `legacy` | Exclude any shared-risk overlap. |
| `low-complexity-advisory` | Retain low-complexity-only annotations; exclude other or mixed categories. |
| `low-complexity-unique-pair` | Retain low-complexity-only evidence only when both native end classes are target and both unique best-target species agree. |
| `low-complexity-nonrisk100-pair` | Apply unique-pair and require at least 100 actual non-low-complexity aligned bases on each mate. |
| `low-complexity-nonrisk200-pair` | Apply unique-pair, positive non-low-complexity support on both mates and at least 200 bp combined. |
| `low-complexity-nonrisk180-pair` | As above, with at least 180 bp combined. |
| `paired-high-cde-only` | Exclude an original-core-passing target pair only if each mate has at least 100 actual aligned bases and every such base is inside the high-CDE union. Other annotations are advisory. |

Support counts reference-aligned CIGAR M/=/X bases outside the frozen low-complexity union. Insertions, clipping, D/N and the unsequenced mate gap do not contribute. Low-complexity exceptions leave other/mixed exclusions unchanged. STRICT additionally applies DF2 capture intervals.

The paired high-CDE policy requires the v0.28 shared BED with SHA256 `aec94df33f252c76777300f45c0e5c3bbf121aed3f9192009ba95592cde8ff8e`. Original annotations and policy decisions are separate in `06_filter/risk_policy_audit.tsv`. See `resources/risk_dictionary.tsv` for resource categories.

## Selection and interpretation

The gap-10 MAPQ exception and 180-bp support rule were selected after inspecting development data. They change scientific filtering and are not independent validation or timing optimizations. Unique concordant alignments can still be wrong. Examples with no wrong candidates before risk filtering cannot measure its error-prevention benefit. Compare correct recovery and wrong assignments on complete fixed inputs; keep unavailable and anchor-failed results distinct from zero.

Use a new output directory when changing policies. Mocked wiring tests and the native synthetic example cover different scopes; wiring tests alone do not validate native mapping or biological performance.
