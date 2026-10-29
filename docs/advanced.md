# Advanced use

The default is `--difficulty normal`: one representative per target species, quality gates 100 bp / MAPQ20 / NM rate 0.08. Complete competitor references remain in the joint reference.

Optional existing presets:

| Preset | Target reference | Quality gates |
| --- | --- | --- |
| `simple` | one representative per species | frozen RELAX2 60 / 5 / 0.12 |
| `normal` (default) | one representative per species | 100 / 20 / 0.08 |
| `hard` | original multi-genome reference | 100 / 20 / 0.08 |

Here `hard` means a multi-genome reference, not newly tightened quality thresholds. Actual genomes per species follow the supplied reference metadata, not an assumption that every species has more than one genome. Both PRIMARY and STRICT are produced in every preset; STRICT is not a difficulty level.

Use the default run command and add `--difficulty simple` or `--difficulty hard`, using a new `--out`. Simple requires nominal PE length >=60 to preserve its frozen 60 bp gate. Host strategy, adapter handling, deduplication, anchor retention and spatial grades are not changed by the preset.

Internal bindings retained for compatibility:

- `ONE_REP_TARGET_DEVELOPMENT`: default/single-representative profile; identity `MycoWES_AB_CDE_v0.28_ONE_REP_TARGET_DEV_20261015`.
- `MULTI_TARGET_V028`: multi-genome profile; identity `MycoWES_AB_CDE_v0.28_reference_qc_beta`.

Inspect or validate the multi-genome resources:

```bash
python3 mycowes.py reference --reference-root "$REF" --reference-mode MULTI_TARGET_V028 --validate
```

Before first multi-genome use, build its index (not distributed):

```bash
minimap2 -x sr -d "$REF/baseline_v028/01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.sr.mmi" \
  "$REF/baseline_v028/01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.fa"
```

`run --reference-mode` may also be supplied, but must agree with the selected preset; conflicts fail instead of overriding it. Missing reference resources do not trigger fallback to a different profile. Program help retains these options and each run records the resolved profile and parameters.

The historical single-representative configuration uses the shared baseline target FASTA for anchor/assembly steps, while joint alignment uses the one-representative target compartment plus all competitors. This resource binding has been preserved rather than silently pruning its target FASTA.
