# Research checklist

## Completed

- [x] Use `VC_SQRT` rather than raw `NONE` Hi-C contacts.
- [x] Record ENCODE accessions, file types, genome build, and MD5 values.
- [x] Split train, validation, and test data by chromosome.
- [x] Compare against a direct genomic-distance baseline.
- [x] Include naive, distance-only, ATAC-only, and combined ablations.
- [x] Report MSE, MAE, R², Pearson, and Spearman.
- [x] Add calibration and distance-stratified residual analysis.
- [x] Save per-pair held-out predictions.
- [x] Reproduce the primary result with modular code.
- [x] Add unit tests, dependency pins, and CI configuration.
- [x] Distinguish executed evidence from planned expansion.

## Next execution stage

- [ ] Extract GM12878 autosomes chr1–22.
- [ ] Run the locked-model confirmatory analysis on chr1–15 and chr22.
- [ ] Inspect chromosome-specific failures and residual plots.
- [ ] Report the chromosome-bootstrap interval.
- [ ] Extract matched K562 autosomes.
- [ ] Run the predefined K562 train/test split.
- [ ] Run and label the zero-shot transfer analysis as exploratory.

## Before public release

- [ ] Run the full notebook from a fresh runtime.
- [ ] Preserve the resulting environment/version log.
- [ ] Review generated figures at full resolution.
- [ ] Add a repository URL and DOI to `CITATION.cff` when available.
- [ ] State all unsuccessful and null analyses.
- [ ] Avoid causal or whole-genome reconstruction claims.
