# ATAC-informed prediction of Hi-C contacts

An exploratory computational-biology project testing whether chromatin accessibility can improve prediction of short-range chromatin contacts in the GM12878 cell line. A gradient-boosted model combines ATAC-seq features with genomic distance and is evaluated on chromosomes not used for training.

![Project summary](figures/ATAC_HiC_final_summary_figure.png)

## Main result

On the held-out long arm of chromosome 21, adding ATAC-derived features improved prediction over a distance-only baseline:

| Metric | Result |
|---|---:|
| MSE improvement | 23.59% |
| MAE improvement | 27.92% |
| Pearson correlation | 0.739 |
| Spearman correlation | 0.816 |
| ATAC permutation test | p = 0.0099 |

The improvement was strongest at longer separations within the modelled 100 kb–1.5 Mb range. In an exploratory leave-one-chromosome-out analysis, the ATAC model improved MSE on all six analysed chromosomes (mean improvement 46.06%).

## Approach

- **Cell line:** GM12878
- **Inputs:** ENCODE ATAC-seq signal and Hi-C contacts
- **Resolution:** 100 kb
- **Modelled pairs:** 100 kb to 1.5 Mb apart
- **Features:** genomic distance, accessibility at both loci, mean accessibility, absolute difference, and product
- **Model:** `HistGradientBoostingRegressor`
- **Primary split:** chromosomes 16–19 for training, chromosome 20 for model selection, and the chromosome 21 long arm for the final test
- **Baseline:** mean contact intensity at each genomic distance, learned from training chromosomes only

![Observed and predicted chromosome 21 contact window](figures/chr21_observed_vs_predicted_1p6Mb.png)

## Repository contents

- `notebooks/ATAC_to_HiC_GM12878.ipynb` — cleaned Colab workflow
- `figures/` — summary and observed-versus-predicted figures
- `results/` — evaluation tables, permutation values, predictions, and JSON summary
- `models/` — fitted pairwise gradient-boosting model

The saved model was created with scikit-learn 1.6.1; that version is pinned in `requirements.txt` for compatible loading.

## Reproduce the analysis

Open the notebook in Google Colab and run the cells from top to bottom. The notebook installs its Python dependencies and downloads the required ENCODE files. The ATAC-seq bigWig is approximately 1.3 GB, so allow time and sufficient Colab disk space.

Data accessions:

- ATAC-seq: [ENCFF667MDI](https://www.encodeproject.org/files/ENCFF667MDI/) from [ENCSR637XSC](https://www.encodeproject.org/experiments/ENCSR637XSC/)
- Hi-C: [ENCFF916GWS](https://www.encodeproject.org/files/ENCFF916GWS/) from [ENCSR830NVY](https://www.encodeproject.org/experiments/ENCSR830NVY/)

## Interpretation and limitations

This is a proof-of-concept association study, not a method for replacing Hi-C experiments or establishing causality. It uses one cell line, six chromosomes, 100 kb bins, short-range contacts, and unnormalised Hi-C counts where KR normalisation was unavailable. The leave-one-chromosome-out analysis is exploratory because chromosome 20 was used during model selection. Broader validation should use additional chromosomes, biological replicates, cell types, normalisation strategies, and external datasets.

## Licence

Code is released under the MIT License. ENCODE data remain subject to their original data-use terms.
