# Executed evidence rules

## Competition and quality

The parser compares each mate's best target and competitor scores. Target support requires a strictly larger target score; a tie is ambiguous. A one-point advantage can pass: no five-point ΔAS gate is applied. Paired support requires the appropriate unambiguous target-species agreement, not a SAM proper-pair flag.

The risk annotator consumes support class, both mates' aligned length/MAPQ and the maximum normalized NM rate. Proper-pair and soft-clipping values are recorded but are not independent exclusion gates. PRIMARY applies adapter/shared-risk exclusions; STRICT adds DF2 capture intervals. Product-specific capture candidates are separate sensitivity-analysis resources.

## Deduplication and anchor assessment

Sequence groups use the ordered uppercase (R1,R2) key. Mate swapping and reverse-complement normalization are not performed; later coordinate grouping does not imply sequence canonicalization.

After quality/risk filtering and exact-sequence deduplication, fragments assigned to each species are remapped to candidate target assemblies. Anchor retention is retained exact-deduplicated fragments divided by that species' input exact-deduplicated fragments, not the risk-mask survival fraction. Anchors rank by retained paired fragments, total alignment score, fewer contigs, larger assembly and accession. Even single-representative competition uses the full target FASTA here.

Retention below 0.90 gives ANCHOR_FAIL; breadth-dependent grades are unavailable. NA means unavailable/not-called, not zero. None means evaluated absence of retained evidence, not a failed run. A positive final count alone is not Supportive/Strong.

## Spatial grades

Coordinate families share selected assembly, contig, fragment start/end and orientation. Bins use coordinate midpoints; default width is 250 bp. Changing width changes bin counts and maximum-bin fraction. Historical descriptive rebinning did not re-estimate the formal 250-bp grade.

| Grade | Criteria after successful anchor QC |
| --- | --- |
| Strong | ≥5 coordinate families; ≥3 bins; physical fragment-span union ≥300 bp; Fmax≤0.50 |
| Supportive | ≥3 coordinate families; ≥2 bins; Fmax<0.80; Strong not satisfied |
| Weak | retained evidence not satisfying higher grades |
| None | evaluated zero evidence |

Fmax=0.50 can be Strong if its other criteria pass. Fmax=2/3 can be Supportive but not Strong. There is no universal Fmax≥0.50 Weak override.

Legacy covered_bp is the union of paired-fragment outer spans on anchor contigs. It includes the insert between mates and is not CIGAR aligned-base coverage. Its name is retained for output compatibility; scientific interpretation must use physical fragment-span coverage.
