# Changelog

## v1.0.0-beta.1

- Standard 100 aligned bases / MAPQ 20 / NM rate 0.08 gates, explicit overrides and an optional 60 / 5 / 0.12 preset.
- Single-representative competition by default, optional complete-target competition and complete target anchors.
- Opt-in multi-assembly MAPQ and low-complexity risk policies, retaining native scores and original annotations.
- Benchmark (8 / 32 / 8 / 4) and light (2 / 2 / 2 / 1) CPU profiles.
- Paired FASTQ validation, input/reference identities and protection against incomplete or incompatible output reuse.
- Native synthetic example, boundary tests, resource dictionaries and analysis mappings.

The complete-target MAPQ and risk options were selected during iterative development. Their rules and selection context are in `docs/filter_policies.md`; they are not the default or independent validation results.
