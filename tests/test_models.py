import numpy as np
import unittest

from atac_hic.evaluation import (
    calibration_metrics,
    chromosome_bootstrap_interval,
    compare_with_baseline,
)
from atac_hic.models import DistanceBaseline


class ModelAndEvaluationTests(unittest.TestCase):
    def test_distance_baseline_uses_training_means(self):
        distance = np.array([[1], [1], [2], [2]])
        target = np.array([2.0, 4.0, 10.0, 14.0])
        baseline = DistanceBaseline().fit(distance, target)
        prediction = baseline.predict(np.array([[1], [2]]))
        self.assertTrue(np.allclose(prediction, [3.0, 12.0]))

    def test_comparison_reports_positive_improvement(self):
        observed = np.array([0.0, 1.0, 2.0, 3.0])
        baseline = np.array([1.5, 1.5, 1.5, 1.5])
        model = np.array([0.1, 0.9, 2.1, 2.9])
        results = compare_with_baseline(observed, model, baseline)
        self.assertGreater(results["MSE_improvement_percent"], 0)
        self.assertGreater(results["MAE_improvement_percent"], 0)

    def test_perfect_calibration(self):
        observed = np.arange(1, 6, dtype=float)
        metrics = calibration_metrics(observed, observed)
        self.assertTrue(np.isclose(metrics["calibration_slope"], 1))
        self.assertTrue(np.isclose(metrics["calibration_intercept"], 0))

    def test_constant_calibration_is_handled(self):
        observed = np.arange(1, 6, dtype=float)
        predicted = np.ones(5)
        metrics = calibration_metrics(observed, predicted)
        self.assertTrue(np.isnan(metrics["calibration_slope"]))
        self.assertEqual(metrics["calibration_intercept"], observed.mean())

    def test_chromosome_bootstrap_uses_chromosome_rows(self):
        import pandas as pd

        table = pd.DataFrame({"improvement": [10, 20, 30, 40]})
        interval = chromosome_bootstrap_interval(
            table,
            "improvement",
            iterations=500,
        )
        self.assertEqual(interval["chromosomes"], 4)
        self.assertLessEqual(interval["lower"], interval["estimate"])
        self.assertLessEqual(interval["estimate"], interval["upper"])


if __name__ == "__main__":
    unittest.main()
