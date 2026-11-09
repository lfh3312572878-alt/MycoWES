# Reference v0.28 resource description

Software: https://github.com/lfh3312572878-alt/MycoWES

Reference DOI: https://doi.org/10.5281/zenodo.23776490

This separate sidecar supersedes the download/software instructions and placeholders in the reference package's old README; it does not edit or repackage the original archive. Authenticated Zenodo inspection on 2026-11-09 confirmed the existing reference dataset version 0.28 as an unpublished Draft. DOI 10.5281/zenodo.23776490 is reserved and will be registered only on publication. It does not yet identify a published reference archive and is not a software version DOI. Third-party review access must be arranged by the authors. Repository inspection confirmed Private visibility and unchanged main baseline ae1f1d4b9f0e33baf05944540b8c5b9bfd59268a. The reviewed corrective engine is fixed at [commit 802200ba32e01fbb0a8e86ccb2fd304d97b31ca8](https://github.com/lfh3312572878-alt/MycoWES/commit/802200ba32e01fbb0a8e86ccb2fd304d97b31ca8) on the private `remediation/20261109` branch, not a published software release.

Original archive: `MycoWES_Reference_v0.28.tar.gz`, 648,693,043 bytes; recorded SHA256 `6f4749795c9455156a6cbcd33aeff41f4c9bdf498df971c131e142c36a9b2936`. That archive was not modified or newly rehashed here. Its `SOURCE_RESOURCES.tsv` and `FILES.tsv` identify supplied FASTA/taxonomy/BED assets; `.mmi` indexes are generated locally and are not distributed. Source-index provenance hashes must not be presented as hashes of a user's rebuilt local index.

The normal/simple competitive reference has 18 target representatives (18 species, 118 contigs), plus 560 competitor assemblies (553 taxa, 7823 contigs). The hard reference has 67 target assemblies (18 species, 783 contigs), with the same competitors. Both use the full target FASTA for anchor assessment. Internal paths and accessions remain unchanged; descriptions use stable public terms rather than changing resource keys.

`risk_dictionary.tsv` distinguishes shared risk annotations, STRICT-only DF2 capture intervals, product-specific capture candidates and historical human-screen windows. Records may overlap; record count and interval-union length are not interchangeable. No BED was reduced to obtain a preferred coverage percentage. Runtime resources remain in the separate original reference package rather than being duplicated into the code package.

`selection.md`, `representatives.tsv` and the preserved selection manifest document actual membership and the available ranking/tie-break evidence. A source script's prespecified/validated labels are not independent proof of timing or ablation validity. Missing ranking-generation and exclusion-validation records are explicitly identified.

## Sources and redistribution

Microbial assembly and contig accessions in the taxonomy/FASTA identify NCBI genome resources. Source accessions can be resolved through `https://www.ncbi.nlm.nih.gov/datasets/genome/<assembly_accession>/`. Missing original download dates or URLs have not been inferred. Human-reference metadata are described separately in `host_reference.md`.

Commercial capture BED provenance includes Agilent S07604514 Regions (V6), Agilent S33266436 Regions (V8), IDT xGen Exome v2 target intervals and Twist Exome 2.0 v2.0.2 covered-target intervals. Target-region and covered-target strata are not interchangeable. The analysis measures sequence similarity to BED-region sequences, not bait hybridization or melting temperature. These source BEDs are not redistributed in this code revision; authors must confirm source-specific access and redistribution permissions before a public reference release.

The reference record's CC BY metadata does not override third-party sequence/design terms or license the software. Code licensing remains a copyright-holder decision. No private preview credential, patient data, raw sequencing input or large intermediate index is included in this code package.

## Available construction sources

Selected original builders are preserved under `reproduce/historical/reference/` and `reproduce/historical/masks/`, with hashes in `reproduce/SOURCE_INDEX.tsv`. The one-representative builder is parameterized and uses the supplied manifest. The v0.28 builder filters an earlier reference and reuses an existing ablation joint FASTA/index; it is not a full de novo database acquisition pipeline. The capture builder requires the separately obtained human FASTA/BED designs and bedtools/minimap2. Complete shared-mask category construction and a chain from the old 31 human windows through a 50-bp buffer to current runtime intervals are unavailable; these are missing historical materials, not new methods reconstructed in this revision.
