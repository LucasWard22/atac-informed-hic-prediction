# Validated result report

## Question

Does endpoint ATAC-seq accessibility improve prediction of short-range
GM12878 cis Hi-C contact strength beyond genomic separation?

## Primary held-out result

The model was trained on chromosomes 16–19, selected on chromosome 20, and
tested on 4,935 pairs from the long arm of chromosome 21.

| Metric | Distance only | Distance + ATAC | Change |
|---|---:|---:|---:|
| MSE | 0.7602 | 0.5355 | 29.57% lower |
| MAE | 0.6351 | 0.4463 | 29.73% lower |
| R² | 0.5346 | 0.6722 | +0.1376 |
| Pearson | 0.7663 | 0.8207 | +0.0545 |
| Spearman | 0.7878 | 0.8584 | +0.0706 |

The calibration slope was 0.964 with an intercept of 0.177. The mean residual
(observed − predicted) was −0.027, indicating little average bias on the
transformed scale.

## Ablation

| Model | MSE | MAE | R² |
|---|---:|---:|---:|
| Naive mean | 1.7145 | 1.0082 | −0.0496 |
| Distance only | 0.7602 | 0.6351 | 0.5346 |
| ATAC only | 1.3828 | 0.8893 | 0.1535 |
| Distance + ATAC | 0.5355 | 0.4463 | 0.6722 |

Distance explains most of the predictable structure, as expected. ATAC alone
is substantially weaker than distance, but adding ATAC to distance improves
the held-out result. This is the scientifically useful finding.

## Chromosome-level robustness

The nested analysis improved on all six tested chromosomes:

| Chromosome | MSE improvement |
|---|---:|
| chr16 | 30.40% |
| chr17 | 61.21% |
| chr18 | 83.02% |
| chr19 | 9.02% |
| chr20 | 67.73% |
| chr21 | 28.54% |

Median improvement was 45.80%; median MAE improvement was 29.98%. The range,
especially the 9.02% result on chr19, shows why individual chromosomes and
uncertainty should be reported rather than only a pooled pair-level score.

## Interpretation

The result supports a predictive association between local accessibility and
short-range chromatin contacts in GM12878 at 100 kb resolution. It does not
show causality, does not reconstruct a full Hi-C map from ATAC alone, and does
not yet demonstrate external cell-line generalisation.

## Most important next evidence

The highest-value next test is the locked-model confirmation on previously
unused GM12878 chromosomes. The second is within-K562 held-out chromosome
validation. Zero-shot GM12878-to-K562 transfer should be interpreted only as
an exploratory stress test.
