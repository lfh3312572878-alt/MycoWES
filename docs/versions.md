# Versions and analysis profiles

Software: **MycoWES v1.0.0-beta.1**. Reference resources: **v0.28**.

| Profile | Competitive targets | Competitors | Anchor targets | Length / MAPQ / NM |
| --- | --- | --- | --- | --- |
| Default standard | 18 representatives | Complete set | 67 assemblies | 100 / 20 / 0.08 |
| Explicit relax2 | 18 representatives | Complete set | 67 assemblies | 60 / 5 / 0.12 |
| Current study | 67 assemblies | Complete set | 67 assemblies | 100 / 20 / 0.08 |
| Parameter scans | Recorded binding | Recorded set | Recorded target set | Per-point settings |

The current study selects `--difficulty hard --mapq-policy assembly-top-ties-gap10 --risk-policy low-complexity-nonrisk180-pair`. PRIMARY and STRICT are output branches, not quality presets. See [filter policies](filter_policies.md) for rules and development context.

The default CPU profile is 8 / 32 / 8 / 4; `--cpu-profile light` selects 2 / 2 / 2 / 1. Anchor minimap2 uses 2 threads. CPU profiles do not change scientific settings.

`current_analysis_profiles.tsv` maps analyses to configurations and frozen source tables. `paper_mapping.tsv` and `reproduce/SOURCE_INDEX.tsv` identify earlier sources separately. Saved results retain their actual execution identities; a version-label or packaging change is not a new analysis or timing measurement.

`PACKAGE_IDENTITY.json` and `SHA256SUMS` identify this distribution. The reference DOI and licensing are separate from the software version.
