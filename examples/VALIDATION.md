# Native validation on 9 November 2026

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
