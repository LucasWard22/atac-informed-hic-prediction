"""ATAC-informed prediction of short-range cis Hi-C contacts."""

from .config import DatasetConfig, DATASETS
from .features import FEATURE_NAMES, PairwiseDataset, build_pairwise_dataset
from .models import DistanceBaseline, ModelSettings, make_model

__all__ = [
    "DATASETS",
    "DatasetConfig",
    "DistanceBaseline",
    "FEATURE_NAMES",
    "ModelSettings",
    "PairwiseDataset",
    "build_pairwise_dataset",
    "make_model",
]

__version__ = "1.0.0"
