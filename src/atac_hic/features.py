from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

FEATURE_NAMES = (
    "distance_bins",
    "lower_ATAC",
    "higher_ATAC",
    "mean_ATAC",
    "ATAC_difference",
    "ATAC_product",
)


@dataclass(frozen=True)
class PairwiseDataset:
    """Pairwise feature matrix, distance column, target and genomic metadata."""

    features: np.ndarray
    distances: np.ndarray
    target: np.ndarray
    metadata: pd.DataFrame

    def __post_init__(self):
        n_rows = len(self.target)
        if self.features.shape != (n_rows, len(FEATURE_NAMES)):
            raise ValueError("Unexpected feature-matrix dimensions")
        if self.distances.shape != (n_rows, 1):
            raise ValueError("Unexpected distance-feature dimensions")
        if len(self.metadata) != n_rows:
            raise ValueError("Metadata and target lengths differ")


def _chromosome_pairs(
    chromosome: str,
    atac: np.ndarray,
    hic: np.ndarray,
    maximum_distance_bins: int,
    minimum_bin: int,
) -> PairwiseDataset:
    atac = np.asarray(atac, dtype=float)
    hic = np.asarray(hic, dtype=float)
    if hic.shape != (len(atac), len(atac)):
        raise ValueError(f"{chromosome}: incompatible ATAC and Hi-C dimensions")

    atac = np.log1p(np.clip(atac, 0, None))
    coverage = np.sum(hic, axis=1) > 0
    feature_blocks = []
    distance_blocks = []
    target_blocks = []
    metadata_blocks = []

    for distance in range(1, maximum_distance_bins + 1):
        bin_a = np.arange(minimum_bin, len(atac) - distance)
        bin_b = bin_a + distance
        valid = coverage[bin_a] & coverage[bin_b]
        bin_a = bin_a[valid]
        bin_b = bin_b[valid]

        accessibility_a = atac[bin_a]
        accessibility_b = atac[bin_b]
        features = np.column_stack([
            np.full(len(bin_a), distance),
            np.minimum(accessibility_a, accessibility_b),
            np.maximum(accessibility_a, accessibility_b),
            (accessibility_a + accessibility_b) / 2,
            np.abs(accessibility_a - accessibility_b),
            accessibility_a * accessibility_b,
        ])

        feature_blocks.append(features)
        distance_blocks.append(np.full((len(bin_a), 1), distance))
        target_blocks.append(np.log1p(hic[bin_a, bin_b]))
        metadata_blocks.append(pd.DataFrame({
            "chromosome": chromosome,
            "bin_a": bin_a.astype(int),
            "bin_b": bin_b.astype(int),
            "distance_bins": distance,
        }))

    return PairwiseDataset(
        features=np.vstack(feature_blocks).astype(np.float32),
        distances=np.vstack(distance_blocks).astype(np.float32),
        target=np.concatenate(target_blocks).astype(np.float32),
        metadata=pd.concat(metadata_blocks, ignore_index=True),
    )


def concatenate_pairwise(datasets: list[PairwiseDataset]) -> PairwiseDataset:
    if not datasets:
        raise ValueError("At least one dataset is required")
    return PairwiseDataset(
        features=np.vstack([dataset.features for dataset in datasets]),
        distances=np.vstack([dataset.distances for dataset in datasets]),
        target=np.concatenate([dataset.target for dataset in datasets]),
        metadata=pd.concat(
            [dataset.metadata for dataset in datasets],
            ignore_index=True,
        ),
    )


def build_pairwise_dataset(
    data: np.lib.npyio.NpzFile | dict[str, np.ndarray],
    chromosomes: list[str],
    maximum_distance_bins: int = 15,
    chromosome_minimum_bins: dict[str, int] | None = None,
) -> PairwiseDataset:
    """Build leakage-safe, short-range cis pair features by chromosome."""

    chromosome_minimum_bins = chromosome_minimum_bins or {}
    datasets = []
    for chromosome in chromosomes:
        datasets.append(_chromosome_pairs(
            chromosome=chromosome,
            atac=data[f"{chromosome}_ATAC"],
            hic=data[f"{chromosome}_HiC"],
            maximum_distance_bins=maximum_distance_bins,
            minimum_bin=chromosome_minimum_bins.get(chromosome, 0),
        ))
    return concatenate_pairwise(datasets)
