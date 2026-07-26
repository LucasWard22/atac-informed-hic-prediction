import unittest

import numpy as np

from atac_hic.features import FEATURE_NAMES, build_pairwise_dataset


def synthetic_data():
    atac = np.array([0, 1, 2, 3, 4, 5], dtype=np.float32)
    hic = np.ones((6, 6), dtype=np.float32)
    hic = (hic + hic.T) / 2
    return {
        "chr1_ATAC": atac,
        "chr1_HiC": hic,
    }


class FeatureTests(unittest.TestCase):
    def test_pairwise_dimensions_and_order(self):
        dataset = build_pairwise_dataset(
            synthetic_data(),
            ["chr1"],
            maximum_distance_bins=2,
        )
        self.assertEqual(dataset.features.shape[1], len(FEATURE_NAMES))
        self.assertEqual(
            dataset.distances.shape,
            (len(dataset.target), 1),
        )
        self.assertEqual(len(dataset.metadata), len(dataset.target))
        self.assertEqual(set(dataset.distances[:, 0]), {1, 2})

    def test_feature_symmetry_under_anchor_swap(self):
        dataset = build_pairwise_dataset(
            synthetic_data(),
            ["chr1"],
            maximum_distance_bins=1,
        )
        lower = dataset.features[:, 1]
        higher = dataset.features[:, 2]
        self.assertTrue(np.all(lower <= higher))
        self.assertTrue(
            np.allclose(dataset.features[:, 3], (lower + higher) / 2)
        )
        self.assertTrue(np.allclose(dataset.features[:, 4], higher - lower))
        self.assertTrue(np.allclose(dataset.features[:, 5], lower * higher))

    def test_minimum_bin_is_respected(self):
        dataset = build_pairwise_dataset(
            synthetic_data(),
            ["chr1"],
            maximum_distance_bins=1,
            chromosome_minimum_bins={"chr1": 3},
        )
        self.assertGreaterEqual(dataset.metadata["bin_a"].min(), 3)


if __name__ == "__main__":
    unittest.main()
