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
export MYCOWES_HOST_MMI=/path/to/hg38.sr.mmi
python3 mycowes.py run --preset standard --reference-root /path/to/reference \
  --sample sample01 --r1 R1.fastq.gz --r2 R2.fastq.gz --out output/sample01 \
  --nominal-pe-length 150 --threads 4 --samtools-threads 1
```

Reference resources: v0.28, https://doi.org/10.5281/zenodo.23776490.
Retain the `baseline_v028/` and `ONE_REP_TARGET_DEVELOPMENT/` directory layout.

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
