# Methods and parameters

Two analysis modes are available. `standard` (the default) selects 67 target
assemblies, the shared competitor library, 67 candidate anchor assemblies,
100 bp / MAPQ20 / NM0.08, assembly-top-ties-gap10 and
low-complexity-nonrisk180-pair. `fast` selects 18 target representatives for
competition, the same 67-assembly anchor pool and risk resources,
60 bp / MAPQ5 / NM0.12, native MAPQ and legacy risk rules.
`manuscript` aliases standard; `relax2` and legacy `--difficulty simple` alias
fast. Legacy normal/hard alias standard. Canonical standard/fast names are
recorded per run. Explicit CLI quality thresholds override the selected mode.
Fast and its aliases require nominal paired-end read length at least 60 bp.

## Competition and quality

Target and competitor sequences share one joint index and one competitive
alignment. For each mate, target support requires a target alignment and either
no competitor alignment or a strictly greater best target AS than competitor AS;
a target/competitor tie is ambiguous. Initial candidates can include one target
mate with an uninformative mate. Quality filtering subsequently requires
`paired_target_consistent`: both mates support the target compartment and their
best target hits have a shared lowest common ancestor at species, genus, family
or order rank. This support class is not a global guarantee that both mates
uniquely identify the same species.

Aligned length applies to both primary mates. The maximum of the two primary
normalized edit distances, NM/L_align, is at most 0.08 in standard or 0.12 in fast,
unless explicitly overridden. L_align counts query bases
in CIGAR M, I, = and X operations. CLI/output names containing `nm-rate` or
`NM_rate` are retained for compatibility. Effective minimum aligned length is the
smaller of the configured threshold and nominal read length. Proper-pair and
soft-clipping are recorded, not additional gates.

The standard `assembly-top-ties-gap10` policy also admits low native MAPQ
when top-scoring hits belong to one target species, tie across at least two
assemblies with one top locus per assembly, and all reported competing-species
or competitor hits are at least 10 AS below the best hit. Supplementary hits,
missing score/identity, a non-best primary hit, or a reporting cap of 101 hits
prevent this exception. Each mate must satisfy native MAPQ >=20 or its eligible
exception; both mates must pass. The exception additionally requires both primary
alignments and both unique best-target labels to identify the same target species.
It is not an initial-candidate rule or an unconditional species-uniqueness gate on
the native-MAPQ path. There is no universal 10-AS gap requirement for that path.

## Risk and deduplication

The standard `low-complexity-nonrisk180-pair` rule retains an otherwise
quality-passing pair with low-complexity-only risk when both mates' unique best
target species agree, both actually aligned mates have positive
non-low-complexity support and their sum is at least
180 bp. Other or mixed risk classes remain excluded. Risk overlap uses aligned
blocks, not the unsequenced inter-mate gap. PRIMARY uses the shared mask;
STRICT additionally applies the shared design-family-support >=2 capture mask.
Adapter rules are unchanged. Exact deduplication uses ordered R1/R2 sequences;
reverse-complement pairs are not normalized together.

## Anchor selection and evidence tiers

After exact deduplication, species assignment uses the original joint alignments:
species rank first by the number of mates with target alignments, then by summed
best per-mate AS. A tie in both values is unassigned. This step does not add an
unconditional requirement that both mates uniquely identify the same species.
Anchor remapping retains pairs with both mates mapping to a common anchor contig.

Anchor assemblies are ordered by retained exact-deduplicated fragments, summed
pair AS, fewer contigs, greater assembly length, then accession. Retention is
anchor-retained / species-assigned exact-deduplicated fragments and must be at
least 0.90. A failed anchor yields ANCHOR_FAIL with unavailable coordinate metrics
and no breadth-dependent grade for that candidate-species detection, not for the
entire sample. Other species are assessed independently.

Coordinate groups use contig, paired outer start/end and orientation. Midpoints
define occupied bins (250 bp by default). C_ind (code: `c_ind`) is the independent
coordinate count, c_max is the independent-coordinate count in the most occupied
bin, and F_max=c_max/C_ind is the largest-bin fraction, not a raw-read fraction.
Physical fragment-span coverage is the
union of paired outer intervals and includes the unsequenced insert gap; it is
not aligned-base coverage. Strong requires >=5 independent coordinates, >=3
occupied bins, physical-span union >=300 bp and largest-bin fraction <=0.50.
Supportive requires >=3 coordinates, >=2 bins and largest-bin fraction <0.80.
Other retained evidence is Weak; absence is None. PRIMARY and STRICT are branches,
not quality presets or independent observations.

## Exact spatial null model

The selected anchor assembly length L determines M=max(2,ceil(L/bin_bp)). Expected
occupied bins under uniform placement are M[1-((M-1)/M)^c_ind]. Under a focal
artifact null, each coordinate escapes the largest-occupancy bin with probability
epsilon=0.01. With k_out=c_ind-c_max, the reported probability is the exact upper
binomial tail sum from k_out to c_ind. The three added columns are descriptive
statistics and do not change tier assignments. They are NA for absent evidence
or failed anchors. The implementation is `exact_spatial_dispersion_stats`.
The calculation uses floating-point arithmetic; a tail reported as 0.0 can be
below the floating-point range rather than mathematically zero.

An optional `MYCOWES_ANCHOR_CACHE` directory can supply MANIFEST.json and
byte-identical read-only per-assembly FASTA files. Without a cache the same
anchor implementation writes per-run FASTA files. Cached file identities are
checked before use.
