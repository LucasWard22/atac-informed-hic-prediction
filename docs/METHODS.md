# Methods

## Study objective

The primary estimand is the reduction in held-out mean squared error obtained
by adding local ATAC-seq features to a genomic-distance model for short-range
cis Hi-C contacts.

## Data

The analysis uses GRCh38 ENCODE files listed in `configs/datasets.json`.
Hi-C contacts are read at 100 kb resolution with `VC_SQRT` normalisation.
ATAC-seq signal p-value bigWigs are averaged into the same bins. Missing or
non-finite signal values are set to zero. Source checksums are verified for
downloaded ATAC files and recorded with the expected Hi-C checksum in the
processed-data provenance file.

The processed dataset validator requires:

- a one-dimensional ATAC vector and square Hi-C matrix per chromosome;
- equal ATAC and Hi-C bin counts;
- finite, non-negative values; and
- a symmetric Hi-C matrix within a numerical tolerance of 1e-5.

## Outcome and features

For every within-chromosome pair separated by 1–15 bins (100 kb–1.5 Mb), the
outcome is:

```text
log1p(VC_SQRT-normalised Hi-C contact)
```

The six predictors are:

1. genomic separation in 100 kb bins;
2. lower endpoint `log1p(ATAC)` signal;
3. higher endpoint `log1p(ATAC)` signal;
4. mean endpoint `log1p(ATAC)` signal;
5. absolute endpoint ATAC difference; and
6. endpoint ATAC product.

Pairs whose endpoint bins both have zero total Hi-C coverage are excluded.
Pairs never cross chromosomes.

## Model

The primary model is `HistGradientBoostingRegressor`. The candidate grid varies
maximum leaf nodes and L2 regularisation while fixing learning rate, iterations,
and random seed. Predictions are clipped at zero because the transformed
contact outcome is non-negative.

Baselines:

- **Naive mean:** mean training outcome for every pair.
- **Distance only:** mean training outcome for each discrete separation.
- **ATAC only:** the same model without the distance feature.
- **Distance + ATAC:** the primary model.

The direct distance baseline is the main comparator because Hi-C contact
frequency is strongly distance dependent.

## Validated split

- Training: chr16, chr17, chr18, chr19
- Model selection: chr20
- Primary held-out test: chr21 bins 130 onward (long arm)

The chromosome 21 restriction preserves the original pre-specified test
region. The test set contains 4,935 pairs.

## Confirmatory full-autosome design

Model settings are selected once on the pilot split and then locked. The model
is refitted on chr16–20 and evaluated on chr1–15 and chr22, which were not used
in the pilot's fitting, model selection, or primary test. Each chromosome is
reported separately. The primary confirmatory summary is the median chromosome
MSE improvement over the distance-only baseline with a non-parametric bootstrap
interval resampling chromosomes.

## K562 validation

The K562 within-line analysis reuses the locked GM12878 hyperparameters,
trains on chr1–18, and evaluates chr19–22. This tests whether the modelling
strategy generalises without reusing held-out chromosomes for tuning.

The GM12878 model is also applied to K562 without refitting. This zero-shot
analysis is explicitly exploratory because assay scale, batch, and K562
karyotype can differ. Failure would not invalidate the within-GM12878 result.

## Metrics and diagnostics

- MSE and MAE on the log1p contact scale
- R²
- Pearson and Spearman correlations
- percentage MSE and MAE improvement over distance only
- calibration slope and intercept
- mean and standard deviation of residuals
- residuals by genomic separation and observed-contact decile
- chromosome-level bootstrap interval for confirmatory summaries

No pair-level p-value is used: nearby genomic pairs are correlated and would
produce anti-conservative uncertainty estimates.

## Reproducibility

Dependencies are pinned in `pyproject.toml` and `requirements.txt`. All
randomised procedures use seed 42. Processed data include a provenance JSON
sidecar. The validated test predictions are saved row by row with chromosome,
bin coordinates, distance, observation, baseline, model prediction, and
residual.
