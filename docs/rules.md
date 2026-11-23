# Executed evidence rules

## Competition and quality

The parser compares each mate's best target and competitor scores. Target support requires a target alignment and either no competitor alignment or a strictly larger best target score; a tie is ambiguous. A one-point advantage can pass the ordinary comparison: no five-point ΔAS gate is applied. Initial candidates may be single-ended, but quality filtering requires `paired_target_consistent`: both mates support the target compartment and the parser finds a shared species/genus/family/order lowest common ancestor. This is not a global same-unique-species guarantee or a SAM proper-pair requirement.

The risk annotator consumes support class, both mates' aligned length/MAPQ and the maximum normalized edit distance NM/L_align (query-aligned CIGAR M/I/=/X bases). Existing NM-rate CLI/column names are retained. In standard mode, each mate must pass native MAPQ>=20 or its eligible assembly-top-ties-gap10 exception, and both mates must pass. Only the exception requires matching unique primary/best-target species across mates plus the recorded top-assembly-tie and >=10-AS alternative-hit gap checks. This is not a universal ΔAS10 gate. Proper-pair and soft-clipping values are recorded but are not independent exclusion gates. PRIMARY applies adapter/shared-risk exclusions with the selected low-complexity policy; the standard low-complexity-only exception additionally requires matching unique best-target species and paired non-low-complexity support >=180 bp. STRICT adds DF2 capture intervals. Product-specific capture candidates are separate sensitivity-analysis resources.

## Deduplication and anchor assessment

Sequence groups use the ordered uppercase (R1,R2) key. Mate swapping and reverse-complement normalization are not performed; later coordinate grouping does not imply sequence canonicalization.

After quality/risk filtering and exact-sequence deduplication, species assignment ranks original target alignments by number of mapped mates, then summed best per-mate AS; a tie in both values is unassigned. There is no additional unconditional same-unique-species gate. Assigned fragments are remapped to candidate target assemblies; retained pairs have both mates on a common anchor contig. Anchor retention is retained exact-deduplicated fragments divided by that species' input exact-deduplicated fragments, not the risk-mask survival fraction. Anchors rank by retained paired fragments, total alignment score, fewer contigs, larger assembly and accession. Even single-representative competition uses the full target FASTA here.

Retention below 0.90 gives ANCHOR_FAIL for that candidate-species detection; its breadth-dependent grades are unavailable, not those of the whole sample or other species. NA means unavailable/not-called, not zero. None means evaluated absence of retained evidence, not a failed run. A positive final count alone is not Supportive/Strong.

## Spatial grades

Coordinate families share selected assembly, contig, fragment start/end and orientation. Bins use coordinate midpoints; default width is 250 bp. C_ind (code: `c_ind`) counts independent coordinate families; c_max counts those in the most occupied bin; Fmax=c_max/C_ind. These are not raw-read counts/fractions. Changing width changes bin counts and maximum-bin fraction. Historical descriptive rebinning did not re-estimate the formal 250-bp grade.

| Grade | Criteria after successful anchor QC |
| --- | --- |
| Strong | ≥5 coordinate families; ≥3 bins; physical fragment-span union ≥300 bp; Fmax≤0.50 |
| Supportive | ≥3 coordinate families; ≥2 bins; Fmax<0.80; Strong not satisfied |
| Weak | retained evidence not satisfying higher grades |
| None | evaluated zero evidence |

Fmax=0.50 can be Strong if its other criteria pass. Fmax=2/3 can be Supportive but not Strong. There is no universal Fmax≥0.50 Weak override.

Legacy covered_bp is the union of paired-fragment outer spans on anchor contigs. It includes the insert between mates and is not CIGAR aligned-base coverage. Its name is retained for output compatibility; scientific interpretation must use physical fragment-span coverage.
