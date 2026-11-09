# MycoWES

MycoWES evaluates paired-end microbial sequence evidence in human exome libraries through adapter trimming, host subtraction, competitive alignment, fragment quality/risk filtering, deduplication, anchor assessment and spatial grading. PRIMARY and STRICT results remain separate. Evidence grades are not clinical diagnoses.

The default remains one representative per target species for competition, all competitor genomes, and the complete target assembly set for anchor assessment. Software version: `python3 mycowes.py --version`. Reference release: v0.28.

## Access and licensing

Repository: https://github.com/lfh3312572878-alt/MycoWES

Reference identifier: https://doi.org/10.5281/zenodo.23776490

Authenticated repository inspection on 2026-11-09 confirmed Private visibility, main commit ae1f1d4b9f0e33baf05944540b8c5b9bfd59268a and no published release. The reference DOI was reserved for a draft in the last verified reference record; its current publication state has not been authenticated. It is not a published software archive DOI. Review access must be arranged by the authors. This revision does not publish resources or change visibility.

No code license has been selected. The copyright holder must supply the applicable license before public distribution; the reference record's CC BY metadata does not license this code. See [resource provenance and redistribution](resources/README_v028.md).

## Installation

Use Linux or Ubuntu/WSL with Python, Bash, Cutadapt, minimap2, samtools, GNU coreutils, awk and util-linux. Runtime Python uses the standard library. The dependency environment is in `environment.yml`; historical plotting scripts require their own stated libraries.

```bash
conda env create -f environment.yml
conda activate mycowes
python3 mycowes.py --version
```

## Reference preparation

Retain the archive's internal paths. The original reference archive is separate and is not modified by this revision.

```bash
mkdir -p resources
tar -xzf MycoWES_Reference_v0.28.tar.gz -C resources
export REF="$PWD/resources/reference"
python3 - "$REF" <<'PY'
import subprocess, sys
from workflow.reference_profiles import bindings
a = bindings(sys.argv[1])
subprocess.run(['minimap2', '-x', 'sr', '-d', a['MYCOWES_JOINT_MMI'],
                a['MYCOWES_JOINT_FASTA']], check=True)
PY
```

The human reference is not supplied. Build from the exact selected GRCh38 FASTA, preserve its source/accession and contig inventory, and record its local index identity. An unverified index must not be labelled GRCh38.p14 or CHM13.

```bash
minimap2 -x sr -d /path/to/hg38.mmi /path/to/hg38.fa
export MYCOWES_HOST_MMI=/path/to/hg38.mmi
python3 mycowes.py reference --reference-root "$REF" --validate
python3 mycowes.py run --sample sample01 \
  --r1 /path/to/R1.fastq.gz --r2 /path/to/R2.fastq.gz \
  --out /path/to/output/sample01 --reference-root "$REF" --nominal-pe-length 150
```

Use a new output directory when changing inputs, references or scientific configuration. Default gates are 100 aligned bases per mate, MAPQ20 and maximum NM rate 0.08; the effective length gate is the smaller of configured length and nominal PE length. See [rules](docs/rules.md), [presets and BAM conversion](docs/advanced.md), and [versions](docs/versions.md).

Historical compartment/profile abbreviations are defined in [terms](docs/terms.md); their resource paths are kept stable.

## Validation and paper sources

Run the small [native synthetic example](examples/README.md). It validates software, not manuscript performance. Unit tests alone do not establish scientific reproducibility.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -B -m unittest discover -s tests -v
for f in workflow/*.sh config/default.env; do bash -n "$f"; done
python3 examples/native_smoke.py --out /tmp/mycowes-native-example
python3 reproduce/summarize_saved.py
```

The shared saved-output entry validates and summarizes the supplied historical spike-in summary without running a mapper. [Paper mapping](docs/paper_mapping.tsv) separates confirmed sources from missing final-figure bindings. Necessary historical sources are in `reproduce/historical/`, with byte-level provenance in `reproduce/SOURCE_INDEX.tsv`. Several retain their original project paths: they are source records, not turnkey launchers. Do not run them against the original project tree.

The code package contains the engine, environment, synthetic example, selected historical scripts, small non-patient summaries and resource/version dictionaries. It does not contain patient-level tables, clinical correspondence, original reads, BAMs, private access URLs, large indexes or a complete database-construction history. Missing mask-construction and final-figure inputs are identified rather than replaced with current defaults.
