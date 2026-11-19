# Explicit-use modes

These modes are not listed in the normal command help and are never enabled automatically. They use the common MycoWES v1.0.0-beta.1 engine and record the selected mode in both run identities and native manifests. They are configurations, not separate software versions.

## Validation62 Strict

```bash
python3 mycowes.py run --special-mode strict62 --sample sample01 \
  --r1 R1.fastq.gz --r2 R2.fastq.gz --out output/sample01_strict62 \
  --reference-root "$REF" --nominal-pe-length 150
```

Fixed binding: v0.28 complete67-target competition plus the frozen complete competitor set, complete67-target anchor reference, aligned length100, native per-mate MAPQ20 and maximum NM/aligned-length0.08. PRIMARY excludes any shared-risk overlap and adapter-core evidence; STRICT_DF2 additionally excludes DF2 capture intervals. Anchor retention is0.90 and dispersion bins are250bp. Use the STRICT_DF2 branch for the designated Strict analysis; PRIMARY is also retained. This is not the gap10/LC180 configuration, and BroadExclusion is not enabled.

The binding was checked against the retained62-library execution configuration and native PRIMARY/STRICT rules. This package update does not rerun the62 libraries or assert bytewise equivalence of outputs from different engine revisions. Retained measurements keep their actual source identities.

## Simple exploratory analysis

```bash
python3 mycowes.py run --special-mode explore --sample sample01 \
  --r1 R1.fastq.gz --r2 R2.fastq.gz --out output/sample01_explore \
  --reference-root "$REF" --nominal-pe-length 150
```

Fixed binding:18 single target representatives plus the complete frozen competitor set, complete67-target anchors, RELAX2 gates60bp / native MAPQ5 / NM0.12, any-shared-risk exclusion, anchor retention0.90 and250bp bins. PRIMARY and STRICT_DF2 remain separate outputs. This is an exploratory configuration, not an independent validation result.

Both modes default to the benchmark CPU profile8/32/8/4, with anchor2. Add `--cpu-profile light` for2/2/2/1 without changing either mode's scientific settings. Add `--host-depleted` only for already trimmed, paired host-unmapped input. Conflicting scientific overrides are rejected. Nominal PE length must be at least100 for strict62 or60 for explore. Always select a new output directory when changing modes.
