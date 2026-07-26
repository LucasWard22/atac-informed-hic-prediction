from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


def regression_metrics(
    observed: np.ndarray,
    predicted: np.ndarray,
) -> dict[str, float]:
    observed = np.asarray(observed)
    predicted = np.asarray(predicted)
    is_constant = np.ptp(predicted) < 1e-12
    return {
        "MSE": float(mean_squared_error(observed, predicted)),
        "MAE": float(mean_absolute_error(observed, predicted)),
        "R2": float(r2_score(observed, predicted)),
        "Pearson": (
            np.nan
            if is_constant
            else float(pearsonr(observed, predicted).statistic)
        ),
        "Spearman": (
            np.nan
            if is_constant
            else float(spearmanr(observed, predicted).statistic)
        ),
    }


def compare_with_baseline(
    observed: np.ndarray,
    model_prediction: np.ndarray,
    baseline_prediction: np.ndarray,
) -> dict[str, float]:
    model_metrics = regression_metrics(observed, model_prediction)
    baseline_metrics = regression_metrics(observed, baseline_prediction)
    mse_improvement = (
        np.nan
        if baseline_metrics["MSE"] == 0
        else (
            (baseline_metrics["MSE"] - model_metrics["MSE"])
            / baseline_metrics["MSE"]
            * 100
        )
    )
    mae_improvement = (
        np.nan
        if baseline_metrics["MAE"] == 0
        else (
            (baseline_metrics["MAE"] - model_metrics["MAE"])
            / baseline_metrics["MAE"]
            * 100
        )
    )
    return {
        **{f"model_{key}": value for key, value in model_metrics.items()},
        **{f"baseline_{key}": value for key, value in baseline_metrics.items()},
        "MSE_improvement_percent": mse_improvement,
        "MAE_improvement_percent": mae_improvement,
    }


def calibration_metrics(
    observed: np.ndarray,
    predicted: np.ndarray,
) -> dict[str, float]:
    """Fit observed = intercept + slope × predicted."""

    if np.ptp(predicted) < 1e-12:
        slope = np.nan
        intercept = float(np.mean(observed))
    else:
        slope, intercept = np.polyfit(predicted, observed, 1)
    residual = observed - predicted
    return {
        "calibration_slope": float(slope),
        "calibration_intercept": float(intercept),
        "mean_residual": float(np.mean(residual)),
        "residual_sd": float(np.std(residual, ddof=1)),
    }


def stratified_residuals(
    observed: np.ndarray,
    predicted: np.ndarray,
    distances: np.ndarray,
    quantiles: int = 10,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    residual = observed - predicted
    distance_values = np.asarray(distances).reshape(-1)

    by_distance = (
        pd.DataFrame({
            "distance_bins": distance_values.astype(int),
            "residual": residual,
            "absolute_residual": np.abs(residual),
        })
        .groupby("distance_bins", as_index=False)
        .agg(
            pairs=("residual", "size"),
            mean_residual=("residual", "mean"),
            MAE=("absolute_residual", "mean"),
        )
    )

    quantile_bin = pd.qcut(
        observed,
        q=quantiles,
        labels=False,
        duplicates="drop",
    )
    by_contact_quantile = (
        pd.DataFrame({
            "contact_quantile": quantile_bin,
            "observed": observed,
            "predicted": predicted,
            "residual": residual,
            "absolute_residual": np.abs(residual),
        })
        .groupby("contact_quantile", as_index=False)
        .agg(
            pairs=("residual", "size"),
            mean_observed=("observed", "mean"),
            mean_predicted=("predicted", "mean"),
            mean_residual=("residual", "mean"),
            MAE=("absolute_residual", "mean"),
        )
    )
    return by_distance, by_contact_quantile


def chromosome_bootstrap_interval(
    chromosome_results: pd.DataFrame,
    metric: str,
    iterations: int = 10_000,
    confidence: float = 0.95,
    random_state: int = 42,
) -> dict[str, float]:
    """Bootstrap the median across chromosomes, not correlated pairs."""

    values = chromosome_results[metric].dropna().to_numpy(dtype=float)
    if len(values) < 2:
        raise ValueError("At least two chromosomes are required")
    rng = np.random.default_rng(random_state)
    samples = rng.choice(values, size=(iterations, len(values)), replace=True)
    estimates = np.median(samples, axis=1)
    alpha = 1 - confidence
    return {
        "estimate": float(np.median(values)),
        "lower": float(np.quantile(estimates, alpha / 2)),
        "upper": float(np.quantile(estimates, 1 - alpha / 2)),
        "confidence": confidence,
        "chromosomes": len(values),
        "iterations": iterations,
    }
