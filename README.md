# MycoWES v1.0.0

MycoWES evaluates paired-end Mollicutes evidence in human exome libraries through
trimming, host subtraction, competitive alignment, fragment quality/risk filtering,
deduplication, anchor assessment and spatial grading.

## Install and run

Linux or WSL requires Python 3, Bash, Cutadapt, minimap2, samtools, GNU coreutils,
awk and util-linux. Runtime Python uses only the standard library.

```bash
conda env create -f environment.yml
conda activate mycowes
```

Download and extract `MycoWES_Reference_v0.28.tar.gz` from the
[v0.28 reference record](https://doi.org/10.5281/zenodo.23776490).
Set `REF` to the extracted reference root containing both `baseline_v028/` and
`ONE_REP_TARGET_DEVELOPMENT/`. The reference package supplies FASTA and resource
files; `.mmi` indexes are built locally, once per selected competitive profile.
For the default `standard` mode, build the joint index and validate resources
before running:

```bash
REF="/path/to/reference"
minimap2 -x sr -d \
  "$REF/baseline_v028/01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.sr.mmi" \
  "$REF/baseline_v028/01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.fa"
python3 mycowes.py reference --reference-root "$REF" --reference-mode MULTI_TARGET_V028 --validate
export MYCOWES_HOST_MMI="/path/to/hg38.sr.mmi"
python3 mycowes.py run --preset standard --reference-root "$REF" \
  --sample sample01 --r1 R1.fastq.gz --r2 R2.fastq.gz --out output/sample01 \
  --nominal-pe-length 150 --threads 4 --samtools-threads 1
```

Build only the competitive index for the mode you use. For `fast`, follow its
indexing/validation commands in [docs/advanced.md](docs/advanced.md) and run with
`--preset fast`. The full FASTQ workflow also requires a separately prepared
human index at `MYCOWES_HOST_MMI`; the microbial joint index does not replace it.

There are two analysis modes. `standard` (the default) uses 67 target assemblies
for competition and anchor assessment, 100 bp / MAPQ20 / NM0.08,
assembly-top-ties-gap10 and low-complexity-nonrisk180-pair. `fast` uses 18 target
representatives for competition, the shared 67-assembly anchor pool,
60 bp / MAPQ5 / NM0.12, native MAPQ and legacy risk rules.
`manuscript` is a compatibility alias of standard; `relax2` and the legacy
`--difficulty simple` select fast. Legacy normal/hard select standard.
Run records store the canonical mode standard or fast. Both modes use the same
complete competitor library, jointly indexed with the selected targets for a
single competitive alignment. Explicit settings are saved per run.

Initial candidates and quality-filtered evidence are distinct stages. Paired
target support uses a shared taxonomic lowest common ancestor, not a global
guarantee of unique species agreement between mates; the low-MAPQ and
low-complexity exceptions have additional species-agreement requirements. See
`docs/methods_and_parameters.md` for the exact executed rules. NM0.08 denotes a
maximum per-mate normalized edit distance NM/L_align of 0.08.

Use `--host-depleted` only for paired reads already trimmed and retained as both
human-unmapped. PRIMARY and STRICT are separate branches. The final species
summary reports independent coordinates, anchor retention, tiers and spatial-null
statistics. See `docs/methods_and_parameters.md` for definitions and BAM conversion
in `docs/advanced.md`.

## Verify and reproduce

```bash
python3 mycowes.py --version
python3 -B -m unittest discover -s tests -v
python3 examples/native_smoke.py --out /tmp/mycowes-native-smoke
python3 reproduce/summarize_saved.py --self-test
```

Reference builders and historical analysis sources are described in
`reproduce/README.md`. Figure workbooks and calculation scripts are distributed
in the separate source-data archive. Code is MIT licensed; reference and
third-party resources retain their own terms.
