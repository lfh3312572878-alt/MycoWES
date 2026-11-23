# Changelog

## v1.0.0

- Default configuration correction, 2026-11-24: align implicit/default standard,
  explicit standard and manuscript with 67-target competition, 67-target anchor
  assessment and gap10/LC180. Keep fast and relax2 as explicit 18-target options.
  Merge relax2 into fast: two canonical modes, standard and fast (60/5/0.12).
  Keep manuscript/relax2 as documented aliases for standard/fast; legacy
  normal/hard and simple resolve to standard and fast, respectively.
  Workflow implementations, reference data and evidence tiers are unchanged.

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
