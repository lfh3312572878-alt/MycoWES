# Native synthetic example

Run in Linux or WSL after activating the environment in `environment.yml`:

```sh
python3 examples/native_smoke.py --out /tmp/mycowes-native-example
```

The output directory must be new and outside the source tree. The script generates random DNA with seed `20261109`, paired 150-bp reads, a synthetic host, one target and one competitor. It builds local minimap2 indices and invokes the public `mycowes.py run` entry point using actual Cutadapt, minimap2 and samtools. No mapper or pipeline subprocess is mocked. It does not download data.

This tests the current implementation, not historical manuscript performance. The toy host is not GRCh38; all reference sequences and risk intervals are synthetic. Profile IDs only satisfy the resolver's layout contract and do not certify v0.28 database content. The anchor-failure fixture intentionally supplies an unrelated anchor sequence with the same species label to exercise the failure state; it is not a recommended reference setup.

`expected_summary.json` defines the independently specified assertions: the normal case loses one shared-risk fragment in both branches, one DF2 fragment in STRICT, and an exact duplicate during deduplication; a competitor and host fragment do not become target calls. Two other native cases exercise `Fmax=0.50` (Strong) and `Fmax=2/3` (Supportive). An anchor-retention-zero case yields `ANCHOR_FAIL` with NA coordinate metrics. A competitor-only input yields no target evidence, not anchor failure.

An additional `--bin-size 100` run uses the same three coordinate positions as the 250-bp `Fmax=2/3` case: occupied bins change from two to three and the largest-bin fraction changes from 2/3 to 1/3. This verifies that the parameter reaches the native analysis, not merely its metadata.

`native_validation_summary.json` records dependency versions, full commands, input hashes and observed outputs. FASTQ, FASTA, indices, BAM and logs are generated only under the requested output directory and are not part of the source/release package. Reusing or overwriting an existing example directory is refused.

The six original cases now use the true implicit standard default. The normal-case input is additionally run once in a separate directory with explicit `--preset standard --min-aligned-bp 100 --min-mapq 20 --max-nm-rate 0.08`. Filter counts and PRIMARY/STRICT branch records must match exactly, excluding absolute output paths. The summary records both results and asserts selected and native effective quality values; both configuration files are included in the engine identity. This does not require relax2 to give the same result.
