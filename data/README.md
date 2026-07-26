# Data

The compact processed dataset for the validated six-chromosome experiment is
included under `processed/` so the primary result can be reproduced without
downloading the 1.43 GB ATAC bigWig again. Full-autosome datasets and original
source files are intentionally excluded.

## Validated GM12878 inputs

- Hi-C: ENCODE `ENCFF916GWS`
- ATAC-seq: ENCODE `ENCFF667MDI`
- Assembly: GRCh38
- Resolution: 100 kb
- Hi-C processing: observed, VC_SQRT balanced
- ATAC processing: mean signal p-value in matching 100-kb bins

## External K562 inputs

- Hi-C: ENCODE `ENCFF080DPJ`
- ATAC-seq: ENCODE `ENCFF600FDO`
- Assembly: GRCh38
- Resolution: 100 kb
- Hi-C processing: observed, VC_SQRT balanced
- ATAC processing: mean signal p-value in matching 100-kb bins

Both cell lines therefore use GRCh38 files with matching ENCODE output types.
K562 is a biologically stringent external target rather than a replicate:
it is an aneuploid chronic-myeloid-leukaemia cell line.

## Creating processed datasets

From the repository root:

```bash
python scripts/extract_encode.py \
  --dataset gm12878 \
  --output data/processed/GM12878_autosomes_100kb_VC_SQRT.npz

python scripts/extract_encode.py \
  --dataset k562 \
  --output data/processed/K562_autosomes_100kb_VC_SQRT.npz
```

The ATAC bigWigs are approximately 1.4 GB for GM12878 and 0.86 GB for K562.
The `.hic` files are accessed remotely by `hic-straw`; they are not downloaded
in full.

Each extraction also writes a `.provenance.json` file recording accessions,
processing settings, chromosomes and the local ATAC-source checksum.

## Repository policy

Do not commit `.bigWig`, `.hic`, full-autosome `.npz`, or other large generated
files. The validated six-chromosome `.npz` and small fitted model are deliberate
reproducibility exceptions. Source data remain subject to ENCODE's data-use
policy.
