# MycoWES

MycoWES evaluates paired-end microbial sequence evidence in human whole-exome sequencing data. It performs adapter trimming, human-host subtraction, competitive microbial alignment, fragment quality/risk filtering, deduplication, anchor checks and spatial evidence grading. PRIMARY and STRICT outputs are retained separately.

The default uses one representative genome per target species and retains the full competitor compartment. Version: **0.3.0-rc1**. Compatible reference release: **v0.28**.

## Installation

Use Linux or an existing Ubuntu/WSL environment. Required software: Python >=3.9, Bash, Cutadapt, minimap2, samtools, gzip, GNU coreutils, awk and util-linux (`column`). Python modules use only the standard library.

```bash
conda env create -f environment.yml
conda activate mycowes
python3 mycowes.py --version
```

Code acquisition location: [GitHub repository](https://github.com/lfh3312572878-alt/MycoWES).

## Reference download and preparation

Reference database (v0.28): [Zenodo DOI: 10.5281/zenodo.23776490](https://doi.org/10.5281/zenodo.23776490) (reserved DOI; private preview available for peer review).

Verify the reference archive against its supplied SHA256, then extract it into `resources/`:

```bash
mkdir -p resources
tar -xzf MycoWES_Reference_v0.28.tar.gz -C resources
export REF="$PWD/resources/reference"
```

The internal reference layout is preserved. Do not rename its directories or accessions. Reference sequences, annotations and indices are not included in the software ZIP. The separate reference package supplies the sequences and annotations, but not `.mmi` files.

Build the default microbial index once, using the same paths resolved by the program:

```bash
python3 - "$REF" <<'PY'
import subprocess, sys
from workflow.reference_profiles import bindings
a = bindings(sys.argv[1])
subprocess.run(['minimap2', '-x', 'sr', '-d', a['MYCOWES_JOINT_MMI'],
                a['MYCOWES_JOINT_FASTA']], check=True)
PY
```

Use an existing GRCh38 host index, or prepare one separately with the appropriate GRCh38 FASTA:

```bash
minimap2 -x sr -d /path/to/hg38.mmi /path/to/hg38.fa
export MYCOWES_HOST_MMI=/path/to/hg38.mmi
python3 mycowes.py reference --reference-root "$REF" --validate
```

Do not use a different human assembly while labelling it GRCh38. The host-subtraction strategy is unchanged: retain pairs where both mates are unmapped, using `-f 12 -F 2304` before name sorting. No additional host reference is silently added.

## Default run

```bash
python3 mycowes.py run --sample sample01 \
  --r1 /path/to/sample01_R1.fastq.gz --r2 /path/to/sample01_R2.fastq.gz \
  --out /path/to/output/sample01 --reference-root "$REF" --nominal-pe-length 150
```

All example paths are placeholders. The default quality gates are 100 aligned bases, MAPQ 20 and maximum NM rate 0.08. As in the existing implementation, the effective aligned-length gate is the smaller of the configured gate and nominal PE length. Anchor retention is 0.90; deduplication and spatial grading remain unchanged. Reference choice alone does not change these gates or PRIMARY/STRICT logic.

## Outputs

Each `--out` is one sample's independent directory. Numbered stages contain trimmed pairs, host-depleted pairs, joint alignments, fragment evidence, risk annotations, PRIMARY/STRICT selections, exact-sequence deduplication, independent coordinates and evidence grades. The pipeline can generate BAM/FASTQ/log files while running; none are shipped in this release.

`run_manifest.tsv` records native completion and parameters. `reference_selection.json` records software version, reference identity, selected settings and input paths. A selection record alone is not proof of completion. Output reuse across different inputs, settings or reference profiles is rejected. Use a new output directory when changing them.

None/Weak/Supportive/Strong describe evidence grades, not clinical diagnosis. A positive final count is not automatically Supportive or Strong. Preserve NA and failure states rather than treating them as zero.

## Tests

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -B -m unittest discover -s tests -v
for f in workflow/*.sh config/default.env; do bash -n "$f"; done
```

Tests generate synthetic temporary fixtures; no clinical data or reference database is included. Release verification results and temporary native-test outputs are distributed separately from the public software/reference archives.

Additional reference configurations are available for advanced use. [Advanced use](docs/advanced.md)
