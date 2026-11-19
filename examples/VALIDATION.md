# Native validation

## 19 November 2026: configuration checks only

The CPU-profile and explicit-use mode revision passed58 unit tests, including argument forwarding, fixed scientific bindings, conflict rejection, CPU-profile run identity and omission of the special modes from regular help. Workflow/config shell syntax was checked in Ubuntu/WSL. The benchmark profile is8/32/8/4, light is2/2/2/1, and anchor stays2. No native pipeline analysis or performance measurement was run for this revision, as requested. The native evidence below belongs to its dated implementations and is not validation of the changed package.

## 16 November 2026: explicit standard/relax2 package

Status: PASS on actual Linux dependencies for `0.3.0-rc2+standard.20261116`, implementation commit `bf273d56a481ccd31f3313873c27098b1373ecbb`. All 36 unit tests passed in both Windows and Linux; all workflow/config shell syntax checks passed. The six existing synthetic cases below were run without mocking any pipeline subprocess, plus one additional run of the same normal input with explicit `--preset standard --min-aligned-bp 100 --min-mapq 20 --max-nm-rate 0.08`.

Implicit standard and explicit 100/20/0.08 produced exactly equal filter counts and PRIMARY/STRICT scientific records. Native manifests confirmed effective length 100, MAPQ 20 and NM rate 0.08; runtime settings remained THREADS=4 and SAMTOOLS_THREADS=1. The normal case yielded PRIMARY 8 / STRICT 7 retained pairs, coordinate counts 7 / 6, both Strong. Fmax 0.50 was Strong, Fmax 2/3 was Supportive, anchor failure remained ANCHOR_FAIL with coordinate metrics NA, and the zero-candidate case remained None.

Test environment: Python 3.13.15, minimap2 2.30-r1287, samtools 1.21, Cutadapt 5.2. Synthetic random-DNA seed: 20261109. The actual current engine/config/VERSION hashes were recorded and checked unchanged during the run. Native validation summary SHA256: `5345aeaedcb70628fe553b60ef9d291b831d862c5182c4f3e90da129f6b73b46`.

This validates only the current small synthetic example and default/explicit parameter equivalence. It is not a rerun of historical manuscript analyses, a comparison against relax2, or validation of the complete reference database. Generated inputs, indices, BAM and logs remain outside the release tree. No software or database was publicly released.

## Historical validation on 9 November 2026

Historical evidence for the November 9 revision only. The values below are preserved as that earlier run, not re-labelled as the November 16 validation above.

Status: PASS for all six cases using the public CLI and actual dependencies in Ubuntu/WSL. No pipeline subprocess was mocked. This validates the remediation code, not the manuscript's historical performance or the complete reference database.

Tested software: `MycoWES 0.3.0-rc1+remediation.20261109`; Python 3.13.15, minimap2 2.30-r1287, samtools 1.21, Cutadapt 5.2. Random DNA seed: `20261109`. Input generation, commands and assertions are in `native_smoke.py` and `expected_summary.json`.

| Case | Observed native result |
| --- | --- |
| Normal | 9 target candidates; 1 shared-risk exclusion; PRIMARY 8, STRICT 7; exact/coordinate dedup 7 and 6; both Strong; anchor retention 1.0. One DF2 fragment is removed only in STRICT. |
| Fmax 0.50 | Six independent coordinates in four 250-bp bins, physical span union 1220 bp; Strong in both branches. |
| Fmax 2/3 | Three coordinates in two 250-bp bins; Supportive in both branches. |
| Bin size 100 | Same coordinate positions as the preceding case yield three occupied bins and Fmax 1/3; Supportive, with actual `bin_size_bp=100`. |
| Anchor failure | Three filtered exact fragments, zero retained by the deliberately unrelated anchor; `ANCHOR_FAIL`, with coordinate/coverage metrics NA. |
| Zero candidate | Competitor-only input yields PRIMARY 0, STRICT 0 and None, not ANCHOR_FAIL. |

Every case reached `run_state.json: COMPLETE`. The normal case also contains a synthetic host pair removed by host depletion and a competitor pair excluded from target evidence. Generated inputs, indices, BAM and native logs were kept outside the source/release tree.

The run recorded all engine/config/VERSION file hashes and asserted they did not change during testing. SHA256 of the sorted compact JSON engine-hash manifest: `20615a8a858c52c2fa3a09f958c9bbfc5dd9297f3fa6fc7a32509ae17653a324`. Full local validation summary SHA256: `c9872eec50c66586a95a3742ef499c188b854ee58bf4e222078bf270500e1c63`.

To repeat, activate the environment and run the command in `README.md` with a new output directory. The resulting `native_validation_summary.json` includes actual commands, tool paths, input hashes, engine hashes and output tables. Any subsequent engine change requires another run before claiming that changed revision passed this validation.
