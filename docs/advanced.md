# Presets and input preparation

## Competitive and anchor references

The default standard mode uses 67-target competition/anchor with gap10/LC180.
Fast explicitly selects 18-target competition with the shared 67-target anchor
pool, native MAPQ and legacy risk rules. Mapper/trimming threads can be set
with `--threads`; additional samtools threads with `--samtools-threads`.
The defaults remain 4 and 1.

| Mode | Competitive target | Competitors | Anchor target | Quality gates |
| --- | --- | --- | --- | --- |
| standard (default) | 67 assemblies / 18 species / 783 contigs | 560 assemblies / 553 taxa / 7823 contigs | complete 67-assembly target FASTA | 100 bp / MAPQ20 / NM0.08; gap10/LC180 |
| fast | 18 representatives / 18 species / 118 contigs | same complete set | same 67-assembly target FASTA | 60 bp / MAPQ5 / NM0.12; native/legacy |

`manuscript` is a compatibility alias for standard; `relax2` and legacy
`--difficulty simple` select fast. Legacy normal/hard select standard. These
aliases are accepted explicitly and recorded as canonical standard/fast modes.
PRIMARY and STRICT are produced for either mode; STRICT is not a difficulty
level. Actual assemblies per species follow metadata.

Use `--preset fast` for the 18-representative, relaxed quality configuration,
using a new output directory. Fast and its aliases require nominal PE length
at least 60. Threshold priority is explicit CLI value, then selected mode;
inherited `CORE_*` variables are ignored by the public CLI. Explicit
`--min-aligned-bp 100 --min-mapq 20 --max-nm-rate 0.08` overrides quality gates
but does not change the selected reference, MAPQ or risk policy. Length must
be positive, MAPQ 0..255, and normalized edit distance NM/L_align finite and
nonnegative; L_align counts query-aligned M/I/=/X CIGAR bases. Existing NM-rate
CLI names are unchanged. Length is still capped at nominal PE length.
Compatible legacy aliases may be supplied together with a preset; conflicting
canonical modes are rejected. Internal profile keys `ONE_REP_TARGET_DEVELOPMENT`
and `MULTI_TARGET_V028` retain file compatibility; their paths are resolved in
`workflow/reference_profiles.py`.

Download and extract reference v0.28 as described in the main README. Set `REF`
to the root containing both `baseline_v028/` and `ONE_REP_TARGET_DEVELOPMENT/`.
The package does not include `.mmi` files. Build the joint competitive index
once for the selected mode; users do not need both indexes.

For the default `standard` mode:

```bash
REF="/path/to/reference"
minimap2 -x sr -d \
  "$REF/baseline_v028/01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.sr.mmi" \
  "$REF/baseline_v028/01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.fa"
python3 mycowes.py reference --reference-root "$REF" --reference-mode MULTI_TARGET_V028 --validate
export MYCOWES_HOST_MMI="/path/to/hg38.sr.mmi"
python3 mycowes.py run --sample sample_multi --r1 R1.fastq.gz --r2 R2.fastq.gz \
  --out output_multi --reference-root "$REF" --nominal-pe-length 150 --preset standard
```

For `fast`, use this index instead:

```bash
REF="/path/to/reference"
minimap2 -x sr -d \
  "$REF/ONE_REP_TARGET_DEVELOPMENT/ONE_REP_TARGET.AB_CDE_JOINT.sr.mmi" \
  "$REF/ONE_REP_TARGET_DEVELOPMENT/ONE_REP_TARGET.AB_CDE_JOINT.fa"
python3 mycowes.py reference --reference-root "$REF" --reference-mode ONE_REP_TARGET_DEVELOPMENT --validate
export MYCOWES_HOST_MMI="/path/to/hg38.sr.mmi"
python3 mycowes.py run --sample sample_fast --r1 R1.fastq.gz --r2 R2.fastq.gz \
  --out output_fast --reference-root "$REF" --nominal-pe-length 150 --preset fast
```

`MYCOWES_HOST_MMI` must point to a separately prepared human index for the full
FASTQ workflow. It is independent of either microbial joint index.

An explicit reference mode must agree with the preset; missing resources do not trigger fallback. These choices reproduce resource/quality selections, not every historical input preparation or saved-output rescore.

## Preprocessing and host subtraction

The full entry uses Cutadapt for paired adapter trimming and a 50-base minimum read length; either short mate rejects the pair. It does not implement fastp, poly-G removal or Q20 trimming. Host subtraction uses minimap2 `-ax sr` against the configured human index, then samtools `-f 12 -F 2304` to retain both-unmapped pairs and exclude secondary/supplementary records. It does not silently add CHM13 or use BWA/85%-query-coverage rules.

A historical host identity records 455 UCSC-style contigs and 3,209,286,105 bases, including primary chromosomes and alternate/unplaced records. See `resources/host_reference.md`; each new run records its actual local index identity.

## BAM conversion

The analysis CLI accepts paired FASTQ. Recover paired reads externally; retain singleton/unpaired outputs separately:

```bash
samtools collate -u -O input.bam | \
 samtools fastq -F 0x900 -n -1 recovered_R1.fastq -2 recovered_R2.fastq \
 -0 unpaired.fastq -s singletons.fastq -
```

Record the source BAM identity and samtools version. Samtools restores sequencer orientation; do not reverse-complement again. Do not treat singleton/unpaired reads as paired inputs or negative results. Use the full pipeline unless the historical experiment explicitly started with pre-depleted/post-trimming inputs. This command does not recover sequence absent from the source BAM.

## Historical resource construction

`reproduce/historical/reference/build_one_rep_reference.py` accepts selection manifest, taxonomy, FASTA and output directory; it refuses existing output files. The supplied selection manifest's SHA256 matches the recorded profile hash and selects the same 18 accessions as its distributed taxonomy.

Historical drivers accept --project-root, --ref-root and --out-dir; see
reproduce/README.md for their input layout and fixed historical designs.
