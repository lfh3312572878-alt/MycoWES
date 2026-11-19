# MycoWES v1.0.0-beta.1

MycoWES assesses paired-end microbial sequence evidence in human whole-exome sequencing data through trimming, human-read subtraction, competitive alignment, fragment filtering, exact-sequence deduplication, anchor assessment and spatial grading. PRIMARY and STRICT are separate output branches.

## Installation

Use Linux or Ubuntu/WSL with Python, Bash, Cutadapt, minimap2, samtools and standard Unix utilities.

```bash
conda env create -f environment.yml
conda activate mycowes
python3 mycowes.py --version
```

Obtain the [v0.28 reference resources](https://doi.org/10.5281/zenodo.23776490) separately. Retain internal paths, set `REF` to the extracted `reference` directory and `MYCOWES_HOST_MMI` to your GRCh38 minimap2 index. See [reference preparation and BAM conversion](docs/advanced.md).

## Run

The default uses 18 target representatives for competition, the complete competitor set and all 67 target assemblies for anchors. Quality gates are 100 aligned bases, MAPQ 20 and maximum NM rate 0.08.

```bash
python3 mycowes.py run --sample sample01 --r1 R1.fastq.gz --r2 R2.fastq.gz \
  --out output/sample01 --reference-root "$REF" --nominal-pe-length 150
```

For the complete-target analysis used in the current study:

```bash
python3 mycowes.py run --sample sample01 --r1 R1.fastq.gz --r2 R2.fastq.gz \
  --out output/sample01_complete --reference-root "$REF" --nominal-pe-length 150 \
  --difficulty hard --min-aligned-bp 100 --min-mapq 20 --max-nm-rate 0.08 \
  --mapq-policy assembly-top-ties-gap10 \
  --risk-policy low-complexity-nonrisk180-pair
```

`--host-depleted` accepts already trimmed, paired host-unmapped input. Use a new output directory when inputs, references or settings change. Parameter scans retain their per-point gates. See [profiles](docs/versions.md), [evidence rules](docs/rules.md) and [filter policies](docs/filter_policies.md).

## CPU profiles

The default uses 8 / 32 / 8 / 4 threads for Cutadapt / host minimap2 / microbial minimap2 / samtools additional threads. Add `--cpu-profile light` for 2 / 2 / 2 / 1. Anchor minimap2 uses 2 threads. These are per-program settings; overlapping processes and I/O add threads. Light mode reduces concurrency, not index memory requirements. CPU availability produces a warning, without changing scientific settings or automatically retrying.

## Results and verification

Initial assignment, quality-pass fragments, filtered fragments and final independent coordinates are different measurements. Retention is not itself a Supportive/Strong call; unavailable or anchor-failed evidence is not a negative result.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -B -m unittest discover -s tests -v
python3 examples/native_smoke.py --out /tmp/mycowes-native-example
python3 reproduce/summarize_saved.py
```

The [synthetic example](examples/README.md) exercises the native pipeline. [Analysis sources](reproduce/README.md), [figure bindings](docs/current_analysis_profiles.tsv) and [source identities](reproduce/SOURCE_INDEX.tsv) are supplied separately from patient data, reads, BAMs and large indexes.

## Resources and licensing

- Software: [lfh3312572878-alt/MycoWES](https://github.com/lfh3312572878-alt/MycoWES)
- Reference v0.28: [10.5281/zenodo.23776490](https://doi.org/10.5281/zenodo.23776490), reserved for the reference draft, not a software DOI.

Code licensing remains to be specified by the copyright holder. Reference CC BY metadata does not license the software or override third-party resource terms.
