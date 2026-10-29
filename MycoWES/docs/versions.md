# Versions and parameters

Software: MycoWES 0.3.0-rc1. Compatible database: v0.28, with profile identities recorded in configuration and each run. Complete Git commit is unavailable for the preserved runtime snapshot; file SHA256 values identify the release. This packaging version is not claimed to be an unchanged historical manuscript executable.

| Context | Target reference | Minimum aligned bp | Minimum MAPQ | Maximum NM rate |
| --- | --- | ---: | ---: | ---: |
| Release default | one representative per target species | 100 | 20 | 0.08 |
| Existing optional simple preset | one representative per target species | 60 | 5 | 0.12 |
| Existing optional multi-genome preset | original multi-genome baseline | 100 | 20 | 0.08 |
| Frozen manuscript RELAX2 experiments | one-representative development reference | 60 | 5 | 0.12 |
| Historical pre-Broad baseline records (77c4e85 lineage) | baseline reference defined by that run | 100 | 20 | 0.08 |

The RELAX2 values are supported by the preserved exploration configuration and code-change audit; the default release values come from the existing baseline configuration. They must not be merged into one manuscript parameter claim. No manuscript text, completed analyses or scientific thresholds were changed for this packaging task. Nominal-PE length handling remains as documented in the default README.

Shared settings: anchor retention 0.90; 250 bp dispersion bins; Supportive >=3 coordinates in >=2 bins and maximum-bin fraction <0.80; Strong >=5 coordinates in >=3 bins, covered bases >=300 and maximum-bin fraction <=0.50. Exact implementation remains authoritative for other edge states.

Default execution resource settings are 4 mapper/trimmer threads and 1 samtools additional thread. They are independent of reference selection. These are execution settings, not biological evidence thresholds. External scheduling must account for all piped processes and other tasks.

License text must be supplied or confirmed by the copyright holder before public release. No new license has been assigned by this packaging work.
