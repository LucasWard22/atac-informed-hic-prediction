#!/usr/bin/env python
import json
from argparse import ArgumentParser
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from atac_hic.evaluation import (
    calibration_metrics,
    compare_with_baseline,
    regression_metrics,
    stratified_residuals,
)
from atac_hic.features import build_pairwise_dataset
from atac_hic.models import (
    DEFAULT_CANDIDATES,
    DistanceBaseline,
    make_model,
    select_settings,
)
from atac_hic.plotting import plot_validation_summary
from atac_hic.provenance import write_run_metadata


def parse_args():
    parser = ArgumentParser(
        description="Reproduce the validated GM12878 chr16–chr21 analysis."
    )
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    data = np.load(args.data, allow_pickle=False)

    training = build_pairwise_dataset(
        data,
        ["chr16", "chr17", "chr18", "chr19"],
    )
    validation = build_pairwise_dataset(data, ["chr20"])
    test = build_pairwise_dataset(
        data,
        ["chr21"],
        chromosome_minimum_bins={"chr21": 130},
    )

    settings, selection_rows = select_settings(
        training.features,
        training.target,
        validation.features,
        validation.target,
        DEFAULT_CANDIDATES,
    )
    model = make_model(settings).fit(training.features, training.target)
    prediction = np.clip(model.predict(test.features), 0, None)

    baseline = DistanceBaseline().fit(
        training.distances,
        training.target,
    )
    baseline_prediction = baseline.predict(test.distances)

    results = compare_with_baseline(
        test.target,
        prediction,
        baseline_prediction,
    )
    results.update(calibration_metrics(test.target, prediction))
    results.update({
        "selected_max_leaf_nodes": settings.max_leaf_nodes,
        "selected_l2_regularisation": settings.l2_regularisation,
        "test_pairs": len(test.target),
    })

    atac_only_settings, _ = select_settings(
        training.features[:, 1:],
        training.target,
        validation.features[:, 1:],
        validation.target,
        DEFAULT_CANDIDATES,
    )
    atac_only_model = make_model(atac_only_settings).fit(
        training.features[:, 1:],
        training.target,
    )
    atac_only_prediction = np.clip(
        atac_only_model.predict(test.features[:, 1:]),
        0,
        None,
    )

    ablation_rows = []
    prediction_sets = {
        "Naive mean": np.full(len(test.target), training.target.mean()),
        "Distance only": baseline_prediction,
        "ATAC only": atac_only_prediction,
        "Distance + ATAC": prediction,
    }
    for name, values in prediction_sets.items():
        ablation_rows.append({
            "model": name,
            **regression_metrics(test.target, values),
        })

    by_distance, by_quantile = stratified_residuals(
        test.target,
        prediction,
        test.distances,
    )

    pd.DataFrame(selection_rows).to_csv(
        args.output_directory / "model_selection.csv",
        index=False,
    )
    pd.DataFrame(ablation_rows).to_csv(
        args.output_directory / "model_ablation.csv",
        index=False,
    )
    by_distance.to_csv(
        args.output_directory / "residuals_by_distance.csv",
        index=False,
    )
    by_quantile.to_csv(
        args.output_directory / "residuals_by_contact_quantile.csv",
        index=False,
    )
    prediction_table = test.metadata.copy()
    prediction_table["observed_log1p_contact"] = test.target
    prediction_table["model_prediction"] = prediction
    prediction_table["distance_baseline_prediction"] = baseline_prediction
    prediction_table["residual"] = test.target - prediction
    prediction_table.to_csv(
        args.output_directory / "chr21_test_predictions.csv.gz",
        index=False,
        compression="gzip",
    )
    plot_validation_summary(
        test.target,
        prediction,
        baseline_prediction,
        test.distances,
        pd.DataFrame(ablation_rows),
        args.output_directory / "chr21_validation_diagnostics.png",
    )
    joblib.dump(model, args.output_directory / "primary_model.joblib")
    with open(args.output_directory / "primary_results.json", "w") as handle:
        json.dump(results, handle, indent=2)
    write_run_metadata(
        args.output_directory / "run_metadata.json",
        analysis="validated GM12878 chr16-chr21 reproduction",
        input_paths={"processed_dataset": args.data},
        parameters={
            "training_chromosomes": ["chr16", "chr17", "chr18", "chr19"],
            "validation_chromosome": "chr20",
            "test_chromosome": "chr21",
            "test_minimum_bin": 130,
            "resolution_bp": 100_000,
            "maximum_distance_bins": 15,
            "hic_normalisation": "VC_SQRT",
            "random_state": 42,
        },
    )

    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
