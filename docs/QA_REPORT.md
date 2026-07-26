# Quality-assurance report

Date: 2026-07-26

## Checks performed

- All Python modules and scripts compiled successfully.
- Eleven automated tests passed, including:
  - processed-data shape, finiteness, non-negativity, and symmetry checks;
  - feature ordering and endpoint-symmetry checks;
  - chromosome-region boundary checks;
  - distance-baseline correctness;
  - calibration edge cases;
  - chromosome-level bootstrap behaviour; and
  - an integration regression test for the validated chr21 result.
- The expansion notebook contains 16 valid notebook cells and every code cell
  passes Python syntax compilation.
- The included processed dataset contains the expected ATAC and Hi-C arrays
  for chr16–chr21 and matches MD5
  `6a1507b409b4f6c825f3ba1f52b12033`.
- The modular pipeline reproduced the primary result:
  - model MSE: 0.5354537584;
  - distance-baseline MSE: 0.7602440858;
  - MSE improvement: 29.568178%;
  - Pearson: 0.8207269; and
  - Spearman: 0.8583626.
- The generated four-panel diagnostic figure was visually inspected.

## Environment limitation

The full-autosome extraction and K562 analysis were not run in this workspace
because their optional extraction dependencies and large ENCODE transfers were
not available here. Those stages are guarded by explicit notebook switches and
must not be represented as completed.
