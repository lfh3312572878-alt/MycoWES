# Human reference identity

The public engine uses a user-configured GRCh38 minimap2 index and retains both-unmapped pairs. It does not infer a human reference patch from a filename, download a human assembly or add CHM13 automatically.

An available historical source identity, `MycoWES_GRCh38_v1`, reports:

| Field | Recorded value |
| --- | --- |
| FASTA sequence count | 455 |
| Total bases | 3,209,286,105 |
| FASTA bytes | 3,273,481,150 |
| FASTA file SHA256 | 5be01555d98347fdb3714dc84c6f77c9d8bc774adcf32c6f7a8fa06f5baf5e51 |
| Sequence identity SHA256 | 283229c8ad1bb1844f9419a9f5667827002405a8d557bf04d3c9615d4c08f66c |
| Source index SHA256 | 701f1a449a88fe4bf4b6c265bb663abfd6f0330761ab16d22b8aa53eb19b3587 |

Identifiers include UCSC-style primary chromosomes and alternate/unplaced records (for example chr1, chr11_KI270721v1_random and chr14_GL000009v2_random); these examples are not an exhaustive contig list. The identity is recorded in `reference_metadata/host_ref/identity.tsv` in the preserved source export. Bind a run to its recorded FASTA and local index identity rather than inferring a patch level from the filename or substituting a current default reference.

When building a local index, retain the FASTA acquisition/accession and `.fai`, minimap2 version/command, and the actual resulting index identity. A historical source-index hash is not proof that a locally rebuilt `.mmi` is byte-identical.
