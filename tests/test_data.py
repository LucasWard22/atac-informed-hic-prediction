import unittest

import numpy as np

from atac_hic.data import validate_processed_dataset


class ProcessedDataTests(unittest.TestCase):
    def test_processed_dataset_validation(self):
        data = {
            "chr1_ATAC": np.ones(4),
            "chr1_HiC": np.eye(4),
        }
        rows = validate_processed_dataset(data, ["chr1"])
        self.assertEqual(rows[0]["bins"], 4)
        self.assertEqual(rows[0]["symmetry_error"], 0)

    def test_asymmetric_hic_is_rejected(self):
        hic = np.eye(4)
        hic[0, 1] = 2
        data = {
            "chr1_ATAC": np.ones(4),
            "chr1_HiC": hic,
        }
        with self.assertRaisesRegex(ValueError, "symmetry"):
            validate_processed_dataset(data, ["chr1"])


if __name__ == "__main__":
    unittest.main()
