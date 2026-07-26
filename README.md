# ATAC-informed prediction of short-range Hi-C contacts

[![Tests](https://github.com/LucasWard22/atac-informed-hic-prediction/actions/workflows/tests.yml/badge.svg)](https://github.com/LucasWard22/atac-informed-hic-prediction/actions/workflows/tests.yml)

## Abstract

Chromatin accessibility and three-dimensional genome organisation are closely
related, but genomic distance remains the strongest predictor of contact
frequency in Hi-C data. This project therefore asks a focused question:

> Does local ATAC-seq signal improve prediction of short-range cis Hi-C contact
> strength beyond genomic distance alone?

Matched GM12878 ATAC-seq and Hi-C data were obtained from ENCODE and summarised
in 100 kb genomic bins. Hi-C contacts were balanced using `VC_SQRT`
normalisation, and only intrachromosomal pairs separated by 100 kbâ€“1.5 Mb were
included. A distance-dependent baseline was compared with a gradient-boosted
model containing genomic distance and ATAC-derived features. Chromosomes were
kept separate during training, model selection and testing to reduce genomic
data leakage.

On the independent chromosome 21 test region, adding ATAC features reduced mean
squared error by **29.57%** and mean absolute error by **29.73%** relative to
the distance-only baseline. The combined model also improved on every
chromosome in a nested chromosome-level validation across chromosomes 16â€“21.
These results support a limited conclusion: local chromatin accessibility
contains useful predictive information about short-range GM12878 contact
strength beyond genomic distance alone.

## Study design

The main analysis used a strict chromosome-level split:

| Dataset partition | Chromosomes | Purpose |
|---|---|---|
| Training | chr16â€“chr19 | Fit candidate models |
| Validation | chr20 | Select model complexity and regularisation |
| Independent test | chr21 long arm | Final performance assessment |

The held-out chromosome 21 region was not used during training or model
selection. The primary test contained **4,935 genomic pairs**.

For each pair, the model used:

- genomic separation in 100 kb bins;
- the lower and higher ATAC signal at the two loci;
- mean ATAC signal;
- absolute difference in ATAC signal; and
- the product of the two ATAC signals.

ATAC values and Hi-C contact strengths were transformed using `log1p`. The
final model was a `HistGradientBoostingRegressor` with 7 maximum leaf nodes,
L2 regularisation of 0.1, a learning rate of 0.05 and 250 boosting iterations.
These settings were selected on chromosome 20 before testing on chromosome 21.

## Results

### ATAC improves prediction on the independent chromosome 21 test

The distance-only baseline captured much of the expected distance-decay
relationship. However, adding ATAC-derived features improved every reported
test metric.

| Held-out chr21 metric | Distance only | Distance + ATAC |
|---|---:|---:|
| MSE | 0.7602 | **0.5355** |
| MAE | 0.6351 | **0.4463** |
| RÂ² | 0.5346 | **0.6722** |
| Pearson correlation | 0.7663 | **0.8207** |
| Spearman correlation | 0.7878 | **0.8584** |

This corresponds to a **29.57% reduction in MSE** and a **29.73% reduction in
MAE**. The calibration slope was 0.964 and the mean residual was âˆ’0.027
`log1p` contact units, indicating little overall bias on the held-out region.

![GM12878 held-out chromosome 21 validation](results/modular_reproduction/chr21_validation_diagnostics.png)

**Figure 1 | Validation of ATAC-informed contact prediction on held-out
chromosome 21.** Top left, observed and predicted `log1p` contact strengths for
the distance-plus-ATAC model; colour intensity represents the density of
genomic pairs and the dashed line shows perfect agreement. Top right, mean
residual and interquartile range at each genomic separation. The mean residual
remained close to zero across most distances, although prediction was slightly
high at the shortest separations. Bottom left, model ablation showing mean
squared error for the naive mean, distance-only, ATAC-only and combined models.
Distance was the dominant individual predictor, but the lowest error was
obtained when distance and ATAC were combined. Bottom right, decile calibration
for the distance-only and combined models. The distance-plus-ATAC predictions
followed the ideal calibration line more closely across most of the observed
contact range.

### The model recovers broad contact structure within the held-out region

To examine the predictions in their genomic context, observed and predicted
contact maps were compared across a 1.6 Mb window of the held-out chromosome 21
region. The predicted map reproduced the main distance-dependent pattern:
contacts were strongest close to the diagonal and generally weakened as
genomic separation increased. It also retained broad differences in predicted
contact strength across the window.

The prediction was smoother than the observed Hi-C map and did not reproduce
every local fluctuation. This is expected from a model based only on genomic
distance and two-locus ATAC signal. The contact map therefore provides a
qualitative view of model behaviour rather than evidence of exact recovery of
fine-scale chromatin structure. Quantitative performance was assessed across
the complete held-out test region using the metrics reported above.

![Observed and predicted Hi-C contacts across the held-out chromosome 21 region](results/modular_reproduction/chr21_observed_vs_predicted_VC_SQRT.png)

**Figure 2 | Observed and predicted contact structure within a held-out
chromosome 21 window.** Left, observed `VC_SQRT`-balanced Hi-C contact
strengths across chr21:34.5â€“36.1 Mb. Right, contact strengths predicted by the
fixed distance-plus-ATAC model for the same genomic pairs. Values are shown as
`log1p` contact strength in 100 kb bins using a shared colour scale, allowing
direct visual comparison between panels. Grey diagonal cells represent
self-interactions, which were excluded from model fitting and evaluation. The
model recovered the broad distance-dependent contact pattern but produced a
smoother map than the observed data.

### Distance and ATAC provide complementary information

The ablation analysis showed that ATAC signal alone did not explain contact
strength as effectively as genomic distance. The combined model nevertheless
performed substantially better than either information source alone.

| Model | MSE | MAE | RÂ² | Pearson |
|---|---:|---:|---:|---:|
| Naive training mean | 1.7145 | 1.0082 | âˆ’0.0496 | â€” |
| Distance only | 0.7602 | 0.6351 | 0.5346 | 0.7663 |
| ATAC only | 1.3828 | 0.8893 | 0.1535 | 0.3933 |
| Distance + ATAC | **0.5355** | **0.4463** | **0.6722** | **0.8207** |

This is biologically plausible: genomic distance describes the strong average
decay in contact frequency, whereas local accessibility adds information about
differences between pairs at similar separations.

### The result is consistent across chromosomes 16â€“21

To test whether the result depended on a favourable chromosome 21 split, a
nested chromosome-level analysis was performed. Each chromosome was held out
in turn, while model settings were selected using a separate validation
chromosome. The combined model improved on the distance baseline for all six
test chromosomes.

- Chromosomes improved: **6 of 6**
- Median MSE reduction: **45.80%**
- Median MAE reduction: **29.98%**
- Lowest MSE reduction: **9.02%** on chromosome 19

The magnitude of improvement differed between chromosomes, so the median
should not be interpreted as a genome-wide estimate. The analysis does,
however, show that the main finding is not restricted to chromosome 21.

## Interpretation

The results indicate that chromatin accessibility contributes information that
is not fully captured by genomic distance. Distance remained the strongest
single predictor, while ATAC alone had limited predictive value. The improved
performance of the combined model therefore reflects complementary rather than
substitutive information.

This project does not show that ATAC-seq can replace Hi-C, nor does it establish
that accessibility directly causes stronger chromatin contacts. The model was
trained and tested within one cell line, at one resolution and over a defined
short-range interval. The conclusion should therefore remain specific to
prediction of short-range GM12878 cis contacts under the conditions tested.

## Methods summary

### Data processing

Matched GRCh38 GM12878 datasets were downloaded from ENCODE. Hi-C matrices were
extracted at 100 kb resolution using `VC_SQRT` balancing. ATAC-seq signal was
averaged into the corresponding genomic bins. Chromosomes with missing or
non-finite values were rejected during validation, and the processed Hi-C
matrices were checked for compatible dimensions, non-negative values and
symmetry.

### Pairwise feature construction

Pairs one to fifteen bins apart were generated separately for each chromosome.
Pairs involving bins without Hi-C coverage were excluded. Feature construction
was symmetric with respect to the two loci, preventing the arbitrary order of
the pair from changing its representation.

### Baseline and model selection

The distance baseline predicted the mean training-set contact strength at each
discrete genomic separation. Candidate gradient-boosted models were trained on
chromosomes 16â€“19 and ranked by MSE on chromosome 20. The selected model was
then evaluated once on the chromosome 21 long arm. Predictions were clipped at
zero because the transformed contact target was non-negative.

### Evaluation

Performance was assessed using MSE, MAE, RÂ², Pearson correlation and Spearman
correlation. Calibration, residuals by genomic separation and model ablation
were examined in addition to headline error metrics. Nested validation kept
the outer test chromosome separate from the chromosome used for model
selection.

Full implementation details are provided in
[`docs/METHODS.md`](docs/METHODS.md), with interpretation and quality checks in
[`docs/VALIDATED_REPORT.md`](docs/VALIDATED_REPORT.md) and
[`docs/QA_REPORT.md`](docs/QA_REPORT.md).

## Source datasets

All source files are GRCh38 ENCODE releases.

| Cell line | Assay | ENCODE file | Output |
|---|---|---|---|
| GM12878 | Hi-C | [ENCFF916GWS](https://www.encodeproject.org/files/ENCFF916GWS/) | Mapping-quality-thresholded contact matrix |
| GM12878 | ATAC-seq | [ENCFF667MDI](https://www.encodeproject.org/files/ENCFF667MDI/) | Signal p-value bigWig |
| K562 | Hi-C | [ENCFF080DPJ](https://www.encodeproject.org/files/ENCFF080DPJ/) | Mapping-quality-thresholded contact matrix |
| K562 | ATAC-seq | [ENCFF600FDO](https://www.encodeproject.org/files/ENCFF600FDO/) | Signal p-value bigWig |

The K562 files are included in the data manifest for the planned external
validation; they were not used to produce the validated results reported
above. File accessions and MD5 checksums are recorded in
[`configs/datasets.json`](configs/datasets.json). Source data remain subject to
ENCODE's data-use policies.

## Repository structure

```text
.
â”œâ”€â”€ configs/                  # Dataset accessions and checksums
â”œâ”€â”€ data/
â”‚   â”œâ”€â”€ processed/            # Validated chr16â€“21 processed dataset
â”‚   â””â”€â”€ README.md             # Data provenance and licensing notes
â”œâ”€â”€ docs/                     # Methods, validation and quality reports
â”œâ”€â”€ notebooks/
â”‚   â”œâ”€â”€ 01_GM12878_validated_chr16_21.ipynb
â”‚   â””â”€â”€ 02_genomewide_external_validation.ipynb
â”œâ”€â”€ results/
â”‚   â”œâ”€â”€ validated/            # Original validated result tables and model
â”‚   â””â”€â”€ modular_reproduction/ # Reproduced results, predictions and figures
â”œâ”€â”€ scripts/                  # Extraction and analysis entry points
â”œâ”€â”€ src/atac_hic/             # Reusable analysis package
â””â”€â”€ tests/                    # Automated data, feature and model tests
```

## Reproducing the validated analysis

Create an isolated environment and install the project:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Reproduce the held-out chromosome 21 analysis:

```bash
python scripts/run_validated.py \
  --data data/processed/GM12878_VC_SQRT_chr16_21.npz \
  --output-directory results/my_reproduction
```

Run the automated tests:

```bash
pytest
```

The continuous-integration workflow also runs the test suite following each
push or pull request.

## Planned validation

The repository includes a separate workflow for two extensions:

1. evaluation on the remaining GM12878 autosomes using model settings fixed
   before those chromosomes are examined; and
2. matched K562 evaluation, followed by an exploratory GM12878-to-K562 transfer
   test.

These analyses are implemented in
[`notebooks/02_genomewide_external_validation.ipynb`](notebooks/02_genomewide_external_validation.ipynb)
and [`scripts/run_expansion.py`](scripts/run_expansion.py), but they have not yet
been executed. They should not be described as successful validation until the
required outputs have been generated and reviewed.

## Limitations

- The validated dataset contains only chromosomes 16â€“21 from GM12878.
- The analysis is restricted to 100 kb bins and separations up to 1.5 Mb.
- Nearby genomic pairs are correlated, so pair-level observations are not
  independent biological replicates.
- ATAC-derived features may capture technical or genomic covariates as well as
  biological accessibility.
- Predictive improvement does not demonstrate a causal relationship between
  accessibility and chromatin contact formation.
- Full-autosome and cross-cell-line confirmation remain outstanding.

## Citation and licence

If you use this repository, please cite the software and data sources described
in [`CITATION.cff`](CITATION.cff). The project code is released under the
[MIT Licence](LICENSE).
