#!/usr/bin/env python
import json
from argparse import ArgumentParser
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from atac_hic.config import AUTOSOMES
from atac_hic.evaluation import (
    calibration_metrics,
    chromosome_bootstrap_interval,
    compare_with_baseline,
    stratified_residuals,
)
from atac_hic.features import (
    build_pairwise_dataset,
    concatenate_pairwise,
)
from atac_hic.models import (
    DEFAULT_CANDIDATES,
    DistanceBaseline,
    make_model,
    select_settings,
)
from atac_hic.plotting import plot_chromosome_improvements
from atac_hic.provenance import write_run_metadata

PILOT_TRAIN = ["chr16", "chr17", "chr18", "chr19"]
PILOT_VALIDATION = ["chr20"]
PRIMARY_TEST = ["chr21"]
CONFIRMATORY_CHROMOSOMES = [
    chromosome
    for chromosome in AUTOSOMES
    if chromosome not in PILOT_TRAIN + PILOT_VALIDATION + PRIMARY_TEST
]


def parse_args():
    parser = ArgumentParser(
        description=(
            "Run confirmatory whole-autosome GM12878 evaluation and optional "
            "K562 within- and cross-cell-line validation."
        )
    )
    parser.add_argument("--gm12878-data", type=Path, required=True)
    parser.add_argument("--k562-data", type=Path)
    parser.add_argument("--output-directory", type=Path, required=True)
    return parser.parse_args()


def chromosome_cache(data, chromosomes):
    cache = {}
    for chromosome in chromosomes:
        minimum = {"chr21": 130} if chromosome == "chr21" else None
        cache[chromosome] = build_pairwise_dataset(
            data,
            [chromosome],
            chromosome_minimum_bins=minimum,
        )
    return cache


def evaluate_by_chromosome(
    model,
    baseline,
    cache,
    chromosomes,
    evaluation_name,
):
    rows = []
    predictions = {}
    for chromosome in chromosomes:
        dataset = cache[chromosome]
        prediction = np.clip(model.predict(dataset.features), 0, None)
        baseline_prediction = baseline.predict(dataset.distances)
        rows.append({
            "evaluation": evaluation_name,
            "chromosome": chromosome,
            "pairs": len(dataset.target),
            **compare_with_baseline(
                dataset.target,
                prediction,
                baseline_prediction,
            ),
            **calibration_metrics(dataset.target, prediction),
        })
        predictions[chromosome] = prediction
    return pd.DataFrame(rows), predictions


def save_diagnostics(
    output_directory,
    prefix,
    cache,
    predictions,
):
    distance_tables = []
    quantile_tables = []
    for chromosome, prediction in predictions.items():
        by_distance, by_quantile = stratified_residuals(
            cache[chromosome].target,
            prediction,
            cache[chromosome].distances,
        )
        by_distance.insert(0, "chromosome", chromosome)
        by_quantile.insert(0, "chromosome", chromosome)
        distance_tables.append(by_distance)
        quantile_tables.append(by_quantile)

    pd.concat(distance_tables, ignore_index=True).to_csv(
        output_directory / f"{prefix}_residuals_by_distance.csv",
        index=False,
    )
    pd.concat(quantile_tables, ignore_index=True).to_csv(
        output_directory / f"{prefix}_residuals_by_contact_quantile.csv",
        index=False,
    )


def main():
    args = parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    gm_data = np.load(args.gm12878_data, allow_pickle=False)
    gm_cache = chromosome_cache(gm_data, AUTOSOMES)

    pilot_training = concatenate_pairwise([
        gm_cache[chromosome] for chromosome in PILOT_TRAIN
    ])
    pilot_validation = gm_cache["chr20"]
    settings, selection_rows = select_settings(
        pilot_training.features,
        pilot_training.target,
        pilot_validation.features,
        pilot_validation.target,
        DEFAULT_CANDIDATES,
    )

    # Lock settings after pilot selection, then refit on chr16–chr20.
    pilot_refit = concatenate_pairwise([
        gm_cache[chromosome]
        for chromosome in PILOT_TRAIN + PILOT_VALIDATION
    ])
    gm_model = make_model(settings).fit(
        pilot_refit.features,
        pilot_refit.target,
    )
    gm_baseline = DistanceBaseline().fit(
        pilot_refit.distances,
        pilot_refit.target,
    )

    gm_results, gm_predictions = evaluate_by_chromosome(
        gm_model,
        gm_baseline,
        gm_cache,
        CONFIRMATORY_CHROMOSOMES,
        "GM12878 locked-model confirmatory",
    )
    gm_results.to_csv(
        args.output_directory / "GM12878_confirmatory_chromosomes.csv",
        index=False,
    )
    plot_chromosome_improvements(
        gm_results,
        args.output_directory / "GM12878_confirmatory_improvements.png",
        "Locked GM12878 model on previously unused chromosomes",
    )
    pd.DataFrame(selection_rows).to_csv(
        args.output_directory / "pilot_model_selection.csv",
        index=False,
    )
    save_diagnostics(
        args.output_directory,
        "GM12878_confirmatory",
        gm_cache,
        gm_predictions,
    )
    gm_interval = chromosome_bootstrap_interval(
        gm_results,
        "MSE_improvement_percent",
    )

    summary = {
        "locked_settings": {
            "max_leaf_nodes": settings.max_leaf_nodes,
            "l2_regularisation": settings.l2_regularisation,
        },
        "GM12878_confirmatory_chromosomes": CONFIRMATORY_CHROMOSOMES,
        "GM12878_chromosomes_improved": int(
            (gm_results["MSE_improvement_percent"] > 0).sum()
        ),
        "GM12878_median_MSE_improvement_percent": float(
            gm_results["MSE_improvement_percent"].median()
        ),
        "GM12878_chromosome_bootstrap": gm_interval,
    }

    joblib.dump(
        gm_model,
        args.output_directory / "GM12878_locked_confirmatory_model.joblib",
    )

    if args.k562_data:
        k562_data = np.load(args.k562_data, allow_pickle=False)
        k562_cache = chromosome_cache(k562_data, AUTOSOMES)

        # Within-K562 evaluation uses locked GM12878 settings and a predefined
        # chromosome split. No K562 test chromosome is used for tuning.
        k562_train_chromosomes = [f"chr{x}" for x in range(1, 19)]
        k562_test_chromosomes = ["chr19", "chr20", "chr21", "chr22"]
        k562_training = concatenate_pairwise([
            k562_cache[chromosome]
            for chromosome in k562_train_chromosomes
        ])
        k562_model = make_model(settings).fit(
            k562_training.features,
            k562_training.target,
        )
        k562_baseline = DistanceBaseline().fit(
            k562_training.distances,
            k562_training.target,
        )
        k562_results, k562_predictions = evaluate_by_chromosome(
            k562_model,
            k562_baseline,
            k562_cache,
            k562_test_chromosomes,
            "K562 within-cell-line",
        )
        k562_results.to_csv(
            args.output_directory / "K562_within_cell_validation.csv",
            index=False,
        )
        plot_chromosome_improvements(
            k562_results,
            args.output_directory / "K562_within_cell_improvements.png",
            "K562 within-cell-line validation",
        )
        save_diagnostics(
            args.output_directory,
            "K562_within_cell",
            k562_cache,
            k562_predictions,
        )
        joblib.dump(
            k562_model,
            args.output_directory / "K562_locked_settings_model.joblib",
        )

        # Exploratory zero-shot transfer: the GM12878 model and GM12878
        # distance baseline are applied directly to K562.
        transfer_results, transfer_predictions = evaluate_by_chromosome(
            gm_model,
            gm_baseline,
            k562_cache,
            k562_test_chromosomes,
            "GM12878-to-K562 zero-shot transfer",
        )
        transfer_results.to_csv(
            args.output_directory / "GM12878_to_K562_transfer.csv",
            index=False,
        )
        plot_chromosome_improvements(
            transfer_results,
            args.output_directory / "GM12878_to_K562_improvements.png",
            "Exploratory GM12878-to-K562 transfer",
        )
        save_diagnostics(
            args.output_directory,
            "GM12878_to_K562",
            k562_cache,
            transfer_predictions,
        )

        summary.update({
            "K562_test_chromosomes": k562_test_chromosomes,
            "K562_within_chromosomes_improved": int(
                (k562_results["MSE_improvement_percent"] > 0).sum()
            ),
            "K562_within_median_MSE_improvement_percent": float(
                k562_results["MSE_improvement_percent"].median()
            ),
            "K562_transfer_chromosomes_improved": int(
                (transfer_results["MSE_improvement_percent"] > 0).sum()
            ),
            "K562_transfer_median_MSE_improvement_percent": float(
                transfer_results["MSE_improvement_percent"].median()
            ),
            "K562_transfer_status": "exploratory",
        })

    with open(
        args.output_directory / "expansion_summary.json",
        "w",
    ) as handle:
        json.dump(summary, handle, indent=2)
    input_paths = {"GM12878_processed_dataset": args.gm12878_data}
    if args.k562_data:
        input_paths["K562_processed_dataset"] = args.k562_data
    write_run_metadata(
        args.output_directory / "run_metadata.json",
        analysis="locked-model GM12878 confirmation and K562 validation",
        input_paths=input_paths,
        parameters={
            "pilot_training_chromosomes": PILOT_TRAIN,
            "pilot_validation_chromosomes": PILOT_VALIDATION,
            "original_primary_test_chromosomes": PRIMARY_TEST,
            "GM12878_confirmatory_chromosomes": CONFIRMATORY_CHROMOSOMES,
            "resolution_bp": 100_000,
            "maximum_distance_bins": 15,
            "hic_normalisation": "VC_SQRT",
            "random_state": 42,
        },
    )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
