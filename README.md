# ATAC-informed prediction of short-range Hi-C contacts

This repository tests a focused biological question:

> Does local chromatin accessibility add predictive information about
> short-range cis chromatin contact strength beyond genomic distance alone?

The analysis uses matched ENCODE ATAC-seq and Hi-C data, chromosome-level data
splits, a direct distance-decay baseline, model ablations, calibration and
residual diagnostics, and a pre-specified expansion to previously unused
chromosomes and a second cell line.

## What has been validated

The completed GM12878 pilot uses 100 kb bins and cis pairs separated by
100 kb–1.5 Mb. Chromosomes 16–19 are used for training, chromosome 20 for model
selection, and the long arm of chromosome 21 for the primary held-out test.
No pair from the test chromosome is used during training or tuning.

| Held-out chr21 metric | Distance only | Distance + ATAC |
|---|---:|---:|
| MSE | 0.7602 | 0.5355 |
| MAE | 0.6351 | 0.4463 |
| R² | 0.5346 | 0.6722 |
| Pearson | 0.7663 | 0.8207 |
| Spearman | 0.7878 | 0.8584 |

This is a **29.57% reduction in MSE** and **29.73% reduction in MAE** relative
to the distance-only baseline. The held-out calibration slope is 0.964 and
the mean residual is −0.027 log1p contact units.

### Held-out chromosome 21 diagnostics

![Held-out chromosome 21 validation diagnostics](results/modular_reproduction/chr21_validation_diagnostics.png)

The combined distance-and-ATAC model outperformed the distance-only baseline on the independent chromosome 21 test region. The diagnostic panels show observed-versus-predicted contacts, residual behaviour and performance across genomic distances.

Nested chromosome-level validation across chromosomes 16–21 improved over
the distance baseline on all 6 chromosomes, with a median MSE improvement of
45.80%. These values are stored in `results/validated/`.

The conclusion supported by these results is deliberately narrow:
endpoint ATAC accessibility contains useful information for predicting
short-range GM12878 Hi-C contact strength beyond genomic distance. The current
evidence does not establish a universal relationship across all chromosomes,
cell types, resolutions, or long-range contacts.

## Research-grade additions

- Exact ENCODE accessions, output types, genome build, and source checksums.
- `VC_SQRT` Hi-C normalisation rather than raw `NONE` contacts.
- Chromosome-level train/validation/test separation.
- Direct distance-decay, naive-mean, and ATAC-only baselines.
- Locked hyperparameters before confirmatory chromosome evaluation.
- Per-pair held-out predictions for audit and reanalysis.
- Calibration, residual-by-distance, contact-quantile, and ablation outputs.
- Chromosome-level bootstrap intervals rather than treating correlated pairs
  as independent samples.
- A full-autosome GM12878 confirmation workflow.
- A matched K562 within-line evaluation and clearly labelled exploratory
  GM12878-to-K562 transfer test.
- Modular source code, unit tests, pinned dependencies, and continuous
  integration.

## Repository layout

```text
.
├── configs/                 # Dataset manifest
├── data/
│   ├── processed/           # Small validated six-chromosome dataset
│   └── README.md            # Data provenance and licensing notes
├── docs/
│   ├── METHODS.md
│   ├── QA_REPORT.md
│   ├── RESEARCH_CHECKLIST.md
│   └── VALIDATED_REPORT.md
├── notebooks/
│   ├── 01_GM12878_validated_chr16_21.ipynb
│   └── 02_genomewide_external_validation.ipynb
├── results/
│   ├── validated/           # Original validated outputs
│   └── modular_reproduction/# Exact modular reproduction + diagnostics
├── scripts/
│   ├── extract_encode.py
│   ├── run_validated.py
│   └── run_expansion.py
├── src/atac_hic/            # Reusable analysis package
└── tests/                   # Automated correctness checks
```

## Quick start

Create an isolated environment and install the project:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Reproduce the completed held-out chromosome 21 analysis:

```bash
python scripts/run_validated.py \
  --data data/processed/GM12878_VC_SQRT_chr16_21.npz \
  --output-directory results/my_reproduction
```

Run the tests:

```bash
pytest
```

For Google Colab, open
`notebooks/02_genomewide_external_validation.ipynb`. Its first section
reproduces the validated result immediately; the later switches control the
large ENCODE downloads and expansion experiments.

## Confirmatory expansion

The expansion is intentionally separated from the validated evidence:

1. Reproduce the chr16–21 pilot and lock the selected model settings.
2. Refit on chr16–20.
3. Evaluate once on previously unused GM12878 chromosomes
   chr1–15 and chr22.
4. Summarise performance by chromosome and compute a chromosome bootstrap
   interval for median MSE improvement.
5. With K562 data, train a K562 model using the locked settings and evaluate
   held-out chr19–22.
6. Apply the GM12878 model directly to K562 as an **exploratory** transfer
   analysis. K562's aneuploid cancer karyotype makes this a stress test, not a
   clean biological replication.

After extracting full-autosome datasets, run:

```bash
python scripts/run_expansion.py \
  --gm12878-data data/processed/GM12878_VC_SQRT_autosomes.npz \
  --k562-data data/processed/K562_VC_SQRT_autosomes.npz \
  --output-directory results/expansion
```

Do not describe the expansion as successful until those output files have
actually been generated and reviewed.

## Source datasets

All files are GRCh38 ENCODE releases.

| Cell line | Assay | ENCODE file | Output |
|---|---|---|---|
| GM12878 | Hi-C | [ENCFF916GWS](https://www.encodeproject.org/files/ENCFF916GWS/) | Mapping-quality-thresholded contact matrix |
| GM12878 | ATAC-seq | [ENCFF667MDI](https://www.encodeproject.org/files/ENCFF667MDI/) | Signal p-value bigWig |
| K562 | Hi-C | [ENCFF080DPJ](https://www.encodeproject.org/files/ENCFF080DPJ/) | Mapping-quality-thresholded contact matrix |
| K562 | ATAC-seq | [ENCFF600FDO](https://www.encodeproject.org/files/ENCFF600FDO/) | Signal p-value bigWig |

The file IDs and MD5 values are also machine-readable in
`configs/datasets.json`. Source data remain subject to ENCODE's data-use
policies.

## Reproducibility status

- **Executed and verified here:** package import, data validation, feature
  construction, source compilation, and exact reproduction of the primary
  chr21 metrics.
- **Included but not executed here:** full-autosome ENCODE extraction,
  confirmatory GM12878 analysis, and K562 analysis. Those stages require large
  downloads and remote `.hic` access and are designed to run in Colab or a
  suitable workstation.

See `docs/VALIDATED_REPORT.md` for interpretation and `docs/METHODS.md` for
the full protocol.
