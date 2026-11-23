# Changelog

## v1.0.0

- Documentation clarification, 2026-11-23: distinguish initial candidates, paired
  taxonomic-LCA support, per-mate quality/rescue and later species assignment;
  define the joint index, species-local anchor failure, c_max and NM/L_align.
  Executable code, parameters, resource data and version are unchanged.
- Add manuscript and fast presets with explicit reference and quality/policy binding.
- Preserve 100 bp / MAPQ20 / NM0.08 manuscript gates, assembly-top-ties-gap10,
  paired low-complexity non-risk support >=180 bp, and all evidence tiers.
- Append exact focal-artifact-null statistics using selected anchor assembly length.
- Expose mapper/trimming and additional samtools thread controls.
- Include reference builders, parameterized analysis sources and native examples.
- Release the code under the MIT License.
