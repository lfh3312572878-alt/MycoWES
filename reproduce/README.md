# Analysis source records and saved-output entry

Run from the code directory:

```bash
python3 reproduce/summarize_saved.py
python3 reproduce/summarize_saved.py --self-test
python3 reproduce/summarize_saved.py --risk-bed /path/to/shared.bed --capture-bed /path/to/DF2.bed
```

This new entry reads only the supplied non-patient historical summaries and selection manifest. It validates denominators and returns a JSON summary to stdout; it does not map reads, select outcomes or claim to reproduce the current final manuscript plots. Output may be redirected to a new file for review. No original output is overwritten.

Optional existing BED inputs are checked as zero-based half-open intervals and summarized once by category (the label before the first pipe) and across all categories. Overlaps are merged per contig; category unions must not be added to obtain the combined mask union. This reports existing boundaries without rebuilding or buffering a mask. The documented denominator is the target FASTA's 52,101,540 bases, not the human capture-design size.

The six small `historical/saved_summary/Fig2*_source.tsv` files are outputs of the historical spike-in figure revision. The figure number is historical: these are not the current Fig.2 mechanism/scan inputs. Four microbial taxa were injected into three pre-depleted WES backgrounds with fixed seeds 20260930–20260934 and doses 1,2,3,5,10,20,50. The design has 420 nonzero technical conditions and three independent blanks; 60 repeated zero labels in the old 480-label design do not create extra inputs. An any-non-None detection includes Weak and is not equivalent to Supportive/Strong.

## Preserved source scripts

`SOURCE_INDEX.tsv` records each supplied source's original SHA256 and provided-file SHA256. Every source transfer in this revision is byte-identical. `source_archive_20261029/` and `local_execution_snapshot/` are source-artifact identifiers, not download URLs. Remote source paths identify the historical project files; the saved-summary entry operates only on local inputs.

| Subdirectory | Purpose and required external input | Expected original output |
| --- | --- | --- |
| reference | Filter explicit selection manifest/taxonomy/FASTA; preserve the v0.28 filter builder | one-representative FASTA/taxonomy/exclusions and QC; locally built indexes |
| masks | Extract human capture-BED regions, map with minimap2sr secondaryyes N100, retain ≥40 aligned bp and identity≥0.90; bedtools, human FASTA and separately obtained capture BEDs required | product candidate BED and mapping QC; saved-read association/mask summaries |
| experiments/run_spike_grid_pre77c.py | Four public microbial WGS runs ERR14129335, SRR1655810, SRR15402804, SRR35637520; three selected pre-depleted public WES backgrounds and frozen sample selection; no download performed | historical mixed FASTQ/run/status records; main420conditions plus shared blanks |
| experiments/confirmed_branches.py | Saved historical runtimes/alignments, challenge assemblies/map and frozen branch manifests; extra taxa and SF4 rescores | historical truth-challenge JSON and parameter/branch tables |
| experiments/sf3_dose_extension_runner.py | Same frozen branch module, seven-dose design and existing dose50 results; adds four non-target taxa's 1–20-pair conditions | 360 additional technical conditions plus60 existing dose50 conditions, combined with the three extra target taxa separately |
| experiments/f4a_98cell_runner_20261014.zero_safe.py | Frozen98-cell ablation plan, savedalignments/outputs and historical engine | branch-specific deduplication/dispersion results; not a universal ablation launcher |
| plot/plot_Fig2_PRE_BROAD_v3.R | Original Figure2_PRE_BROAD condition/stage/tier tables; R ggplot2, patchwork, scales | historical spike-in figure and small summaries, not current final panels |
| plot/build_fig3_continuous.py | Frozen comparator species scores/status and original423condition design; NumPy/Matplotlib | continuous PR/AP and LOBO tool tables and historical figure |

Historical drivers retain original absolute project paths and may write, lock or launch tools when explicitly executed. They are preserved for source review, not invoked by the common entry. Do not point them at the original project or start an old queue. Their headers may declare a scientific commit; those strings are provenance records, not a recovered Git history.

See `docs/paper_mapping.tsv` for the analysis source map and `resources/` for the runtime resource dictionary and representative list. The package does not contain private clinical records, FASTQ/BAM, large databases, commercial source BEDs or credentials.
