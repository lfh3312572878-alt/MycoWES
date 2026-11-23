# Analysis and reference sources

The saved-summary entry reads the supplied historical tables:
```bash
python3 reproduce/summarize_saved.py --self-test
python3 reproduce/summarize_saved.py --risk-bed shared.bed --capture-bed DF2.bed
```

Historical drivers use explicit roots:
```bash
python3 reproduce/historical/experiments/confirmed_branches.py \
  --project-root /path/to/workspace --ref-root /path/to/reference-workspace \
  --out-dir /path/to/new-results
```

The project workspace contains the original analysis/, software/ and envs/
subdirectories. The reference workspace contains human, capture and microbial
resource subdirectories. Original relative dataset names and fixed seeds are
preserved. Provide those input tables/reads before executing a historical driver;
`--help` lists the common root arguments. Output-producing roots use --out-dir.

The original builders and their adapted copies are indexed by SHA256 in
SOURCE_INDEX.tsv. Reference builders require the selection manifest, source FASTA,
taxonomy and resource metadata. Capture builders additionally require separately
obtained platform BEDs and the human reference. No clinical reads, BAMs, large
indexes or commercial BEDs are included.

The saved spike-in tables describe 420 nonzero technical conditions and three
backgrounds with fixed seeds 20260930–20260934 and doses 1,2,3,5,10,20,50.
They retain their historical design. Current manuscript workbooks and analysis
entries are provided in the separate source-data archive.
