import unittest
from pathlib import Path

import numpy as np

from atac_hic.evaluation import compare_with_baseline
from atac_hic.features import build_pairwise_dataset
from atac_hic.models import (
    DEFAULT_CANDIDATES,
    DistanceBaseline,
    make_model,
    select_settings,
)


class ValidatedResultRegressionTest(unittest.TestCase):
    def test_primary_chr21_result_is_reproduced(self):
        root = Path(__file__).resolve().parents[1]
        data_path = (
            root
            / "data"
            / "processed"
            / "GM12878_VC_SQRT_chr16_21.npz"
        )
        with np.load(data_path, allow_pickle=False) as data:
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

        settings, _ = select_settings(
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
        metrics = compare_with_baseline(
            test.target,
            prediction,
            baseline.predict(test.distances),
        )

        self.assertEqual(settings.max_leaf_nodes, 7)
        self.assertAlmostEqual(settings.l2_regularisation, 0.1)
        self.assertAlmostEqual(metrics["model_MSE"], 0.5354537584, places=5)
        self.assertAlmostEqual(
            metrics["MSE_improvement_percent"],
            29.5681784,
            places=4,
        )
        self.assertAlmostEqual(metrics["model_Pearson"], 0.8207269, places=5)
        self.assertAlmostEqual(metrics["model_Spearman"], 0.8583626, places=5)


if __name__ == "__main__":
    unittest.main()
