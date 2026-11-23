# Reference v0.28 resource description

Software: MycoWES v1.0.0 (MIT License)

Reference DOI: https://doi.org/10.5281/zenodo.23776490

Reference resource record: 23776490; version 0.28. The resource version and
software version are tracked separately.

Original archive: `MycoWES_Reference_v0.28.tar.gz`, 648,693,043 bytes; recorded SHA256 `6f4749795c9455156a6cbcd33aeff41f4c9bdf498df971c131e142c36a9b2936`. Its `SOURCE_RESOURCES.tsv` and `FILES.tsv` identify supplied FASTA/taxonomy/BED assets; `.mmi` indexes are generated locally and are not distributed. Source-index provenance hashes must not be presented as hashes of a user's rebuilt local index.

The two analysis modes use these competitive references: standard/default has 67 target assemblies (18 species, 783 contigs), plus 560 competitor assemblies (553 taxa, 7823 contigs), with 100 bp / MAPQ20 / NM0.08 and gap10/LC180. Fast has 18 target representatives (18 species, 118 contigs), with the same competitors, 60 bp / MAPQ5 / NM0.12 and native MAPQ/legacy risk rules. Manuscript aliases standard; relax2 and legacy simple alias fast; legacy normal/hard alias standard. In either mode, selected targets and the complete competitor set form one joint FASTA/index for a single competitive alignment, not two independent alignments. Both use the full 67-target FASTA for anchor assessment. Internal paths and accessions remain unchanged; descriptions use stable public terms rather than changing resource keys.

`risk_dictionary.tsv` distinguishes shared risk annotations, STRICT-only DF2 capture intervals, product-specific capture candidates and historical human-screen windows. Records may overlap; record count and interval-union length are not interchangeable. Runtime resources remain in the separate original reference package rather than being duplicated into the code package.

`selection.md`, `representatives.tsv` and the preserved selection manifest document the supplied membership and the ranking/tie-break keys used by the representative builder.

## Sources and redistribution

Microbial assembly and contig accessions in the taxonomy/FASTA identify NCBI genome resources. Source accessions can be resolved through `https://www.ncbi.nlm.nih.gov/datasets/genome/<assembly_accession>/`. Human-reference metadata are described separately in `host_reference.md`.

Commercial capture BED provenance includes Agilent S07604514 Regions (V6), Agilent S33266436 Regions (V8), IDT xGen Exome v2 target intervals and Twist Exome 2.0 v2.0.2 covered-target intervals. Target-region and covered-target strata are not interchangeable. The analysis measures sequence similarity to BED-region sequences, not bait hybridization or melting temperature. These source BEDs are obtained separately under their source-specific access and redistribution terms and are not redistributed in this code revision.

The reference record's CC BY metadata does not override third-party sequence/design terms or license the software. The software is distributed under the MIT License. No private preview credential, patient data, raw sequencing input or large intermediate index is included in this code package.

## Available construction sources

Selected original builders are preserved under `reproduce/historical/reference/` and `reproduce/historical/masks/`, with hashes in `reproduce/SOURCE_INDEX.tsv`. The one-representative builder is parameterized and uses the supplied manifest. The v0.28 builder takes an earlier reference as input, filters its members and reuses the recorded joint FASTA/index. The capture builder requires separately obtained human FASTA/BED designs and bedtools/minimap2. Runtime risk filtering uses the supplied BED boundaries without additional padding. Historical human-screen windows are described separately in `risk_dictionary.tsv` and are not substituted for the runtime BED.
