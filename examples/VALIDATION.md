# Native validation

MycoWES 1.0.0 was tested on Linux using real minimap2, samtools and Cutadapt.
The deterministic synthetic example covers normal evidence, zero candidates,
anchor failure, and spatial-grade boundaries. Default and explicit standard
settings produce identical results. A multi-profile example tests manuscript
reference/MAPQ wiring with an explicit legacy risk rule for its synthetic BED;
the manuscript LC180 policy remains bound to the v0.28 BED.

Run:
```bash
python3 -B -m unittest discover -s tests -v
python3 examples/native_smoke.py --out /tmp/mycowes-native-smoke
```

The output contains commands, dependency versions, generated input digests,
expected/observed counts, quality settings and branch summaries. It is separate
from the source package.
