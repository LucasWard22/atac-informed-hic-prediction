from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_squared_error


@dataclass(frozen=True)
class ModelSettings:
    max_leaf_nodes: int = 7
    l2_regularisation: float = 0.1
    learning_rate: float = 0.05
    max_iter: int = 250
    random_state: int = 42


class DistanceBaseline:
    """Mean training-set contact at each discrete genomic separation."""

    def __init__(self):
        self.means_: dict[int, float] | None = None

    def fit(self, distances: np.ndarray, target: np.ndarray):
        distance_values = np.asarray(distances).reshape(-1)
        self.means_ = {
            int(distance): float(target[distance_values == distance].mean())
            for distance in np.unique(distance_values)
        }
        return self

    def predict(self, distances: np.ndarray) -> np.ndarray:
        if self.means_ is None:
            raise RuntimeError("DistanceBaseline must be fitted before prediction")
        values = np.asarray(distances).reshape(-1)
        missing = sorted(set(map(int, np.unique(values))) - set(self.means_))
        if missing:
            raise ValueError(f"Unseen genomic distances: {missing}")
        return np.asarray([self.means_[int(value)] for value in values])


def make_model(settings: ModelSettings) -> HistGradientBoostingRegressor:
    return HistGradientBoostingRegressor(
        learning_rate=settings.learning_rate,
        max_iter=settings.max_iter,
        max_leaf_nodes=settings.max_leaf_nodes,
        l2_regularization=settings.l2_regularisation,
        random_state=settings.random_state,
    )


def select_settings(
    training_features: np.ndarray,
    training_target: np.ndarray,
    validation_features: np.ndarray,
    validation_target: np.ndarray,
    candidates: list[ModelSettings],
) -> tuple[ModelSettings, list[dict]]:
    """Select settings on a chromosome-level validation partition."""

    best_settings = None
    best_mse = np.inf
    rows = []

    for settings in candidates:
        model = make_model(settings)
        model.fit(training_features, training_target)
        prediction = np.clip(model.predict(validation_features), 0, None)
        mse = mean_squared_error(validation_target, prediction)
        rows.append({
            "max_leaf_nodes": settings.max_leaf_nodes,
            "l2_regularisation": settings.l2_regularisation,
            "validation_MSE": mse,
        })
        if mse < best_mse:
            best_mse = mse
            best_settings = settings

    if best_settings is None:
        raise RuntimeError("No candidate model was fitted")
    return best_settings, rows


DEFAULT_CANDIDATES = [
    ModelSettings(max_leaf_nodes=7, l2_regularisation=0.1),
    ModelSettings(max_leaf_nodes=15, l2_regularisation=0.1),
    ModelSettings(max_leaf_nodes=15, l2_regularisation=1.0),
    ModelSettings(max_leaf_nodes=31, l2_regularisation=1.0),
    ModelSettings(max_leaf_nodes=31, l2_regularisation=10.0),
]
