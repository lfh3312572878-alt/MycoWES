# Native validation

The deterministic synthetic example runs real minimap2, samtools and Cutadapt
on Linux. It covers normal evidence, zero candidates, anchor failure and spatial
grade boundaries. Its synthetic BED requires an explicit legacy risk override;
it therefore checks default/standard/manuscript reference, MAPQ and quality
wiring, but does not validate the full LC180 default. A separate fast run checks
the single-representative/native binding and 60 bp / MAPQ5 / NM0.12 gates. CLI wiring unit tests verify that
implicit standard, explicit standard and manuscript select gap10/LC180 without
an override. LC180 interval loading retains the unchanged v0.28 BED identity gate.

Run:
```bash
python3 -B -m unittest discover -s tests -v
python3 examples/native_smoke.py --out /tmp/mycowes-native-smoke
```

The output contains commands, dependency versions, generated input digests,
expected/observed counts, quality settings and branch summaries. It is separate
from the source package.
