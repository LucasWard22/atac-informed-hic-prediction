from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = {
    "model": "#176B87",
    "baseline": "#A7C7E7",
    "accent": "#E07A5F",
    "grid": "#D8DEE9",
}


def plot_validation_summary(
    observed: np.ndarray,
    predicted: np.ndarray,
    baseline_predicted: np.ndarray,
    distances: np.ndarray,
    ablation: pd.DataFrame,
    output_path: str | Path,
) -> Path:
    """Create a publication-ready four-panel held-out-test diagnostic."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    observed = np.asarray(observed)
    predicted = np.asarray(predicted)
    baseline_predicted = np.asarray(baseline_predicted)
    distances = np.asarray(distances).reshape(-1)
    residual = observed - predicted

    figure, axes = plt.subplots(2, 2, figsize=(12, 9))
    figure.suptitle(
        "GM12878 held-out chromosome 21 validation",
        fontsize=16,
        fontweight="bold",
    )

    maximum = float(max(observed.max(), predicted.max()))
    axes[0, 0].hexbin(
        observed,
        predicted,
        gridsize=35,
        mincnt=1,
        cmap="Blues",
    )
    axes[0, 0].plot([0, maximum], [0, maximum], "--", color=COLORS["accent"])
    axes[0, 0].set(
        xlabel="Observed log1p contact",
        ylabel="Predicted log1p contact",
        title="Observed versus predicted",
    )

    residual_frame = pd.DataFrame({
        "distance": distances,
        "residual": residual,
    })
    residual_summary = (
        residual_frame.groupby("distance", as_index=False)
        .agg(
            mean=("residual", "mean"),
            lower=("residual", lambda x: np.quantile(x, 0.25)),
            upper=("residual", lambda x: np.quantile(x, 0.75)),
        )
    )
    axes[0, 1].fill_between(
        residual_summary["distance"],
        residual_summary["lower"],
        residual_summary["upper"],
        color=COLORS["baseline"],
        alpha=0.45,
        label="interquartile range",
    )
    axes[0, 1].plot(
        residual_summary["distance"],
        residual_summary["mean"],
        marker="o",
        color=COLORS["model"],
        label="mean residual",
    )
    axes[0, 1].axhline(0, linestyle="--", color=COLORS["accent"])
    axes[0, 1].set(
        xlabel="Genomic separation (100 kb bins)",
        ylabel="Observed − predicted",
        title="Residuals by genomic separation",
    )
    axes[0, 1].legend(frameon=False)

    metric_order = ["Naive mean", "Distance only", "ATAC only", "Distance + ATAC"]
    ablation_plot = ablation.set_index("model").reindex(metric_order)
    axes[1, 0].bar(
        np.arange(len(ablation_plot)),
        ablation_plot["MSE"],
        color=[
            "#CBD5E1",
            COLORS["baseline"],
            "#81B29A",
            COLORS["model"],
        ],
    )
    axes[1, 0].set_xticks(
        np.arange(len(ablation_plot)),
        ["Naive", "Distance", "ATAC", "Distance\n+ ATAC"],
    )
    axes[1, 0].set(
        ylabel="Mean squared error",
        title="Model ablation on held-out chr21",
    )

    calibration = (
        pd.DataFrame({
            "observed": observed,
            "model": predicted,
            "baseline": baseline_predicted,
        })
        .assign(bin=lambda frame: pd.qcut(
            frame["model"],
            q=10,
            labels=False,
            duplicates="drop",
        ))
        .groupby("bin", as_index=False)
        .agg(
            observed=("observed", "mean"),
            model=("model", "mean"),
            baseline=("baseline", "mean"),
        )
    )
    axes[1, 1].plot(
        calibration["observed"],
        calibration["model"],
        "o-",
        color=COLORS["model"],
        label="Distance + ATAC",
    )
    axes[1, 1].plot(
        calibration["observed"],
        calibration["baseline"],
        "o-",
        color=COLORS["baseline"],
        label="Distance only",
    )
    calibration_max = float(calibration["observed"].max())
    axes[1, 1].plot(
        [0, calibration_max],
        [0, calibration_max],
        "--",
        color=COLORS["accent"],
        label="ideal",
    )
    axes[1, 1].set(
        xlabel="Mean observed contact",
        ylabel="Mean predicted contact",
        title="Decile calibration",
    )
    axes[1, 1].legend(frameon=False)

    for axis in axes.flat:
        axis.grid(axis="y", color=COLORS["grid"], linewidth=0.7, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)

    figure.tight_layout()
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(figure)
    return output_path


def plot_chromosome_improvements(
    results: pd.DataFrame,
    output_path: str | Path,
    title: str,
) -> Path:
    """Plot MSE improvement over the distance-only baseline by chromosome."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ordered = results.copy()
    ordered["chromosome_number"] = (
        ordered["chromosome"].str.replace("chr", "", regex=False).astype(int)
    )
    ordered = ordered.sort_values("chromosome_number")
    values = ordered["MSE_improvement_percent"].to_numpy()
    colors = np.where(values >= 0, COLORS["model"], COLORS["accent"])

    figure, axis = plt.subplots(figsize=(11, 5))
    axis.bar(ordered["chromosome"], values, color=colors)
    axis.axhline(0, color="#334155", linewidth=1)
    axis.set(
        xlabel="Held-out chromosome",
        ylabel="MSE improvement over distance-only baseline (%)",
        title=title,
    )
    axis.grid(axis="y", color=COLORS["grid"], linewidth=0.7, alpha=0.7)
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(figure)
    return output_path
