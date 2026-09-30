"""Two-objective feedback and second-round integration checks."""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import qLogNEHVI as bo
import simu_experiment
from bo_utils import load_experiment_data
from test_support import small_pool


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        return reader.fieldnames, list(reader)


class MultiRoundtripTests(unittest.TestCase):
    def test_example_input_and_second_round(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            experiment, pool, _, pool_rows = small_pool(root)
            data = load_experiment_data(experiment, pool, label="multi")
            self.assertEqual(data.train_X.shape, (8, 7))
            self.assertEqual(data.train_Y.shape, (8, 2))
            self.assertEqual(len(data.candidate_indices), 5)
            output = root / "run"
            argv = ["qLogNEHVI.py", "--experiment", str(experiment), "--pool", str(pool),
                    "--output", str(output), "--batch-size", "2", "--mc-samples", "16",
                    "--fit-maxiter", "30", "--directions", "max", "min"]
            with patch.object(sys, "argv", argv):
                args = bo.parse_args()
            bo.run(args)
            self.assertTrue((output / "candidate_predictions.csv").exists())
            fields, first = read_rows(output / "recommendation_1.csv")
            self.assertEqual(fields, list(pool_rows[0]) + ["logCE", "Conductivity"])
            self.assertEqual(len(first), 2)
            self.assertTrue(all(row["logCE"] == row["Conductivity"] == ""
                                for row in first))
            first_ref = json.loads((output / "training_summary.json").read_text())["ref_point"]
            _, observations = read_rows(output / "observation.csv")
            self.assertEqual(len(observations), 9)
            self.assertEqual(sum(bool(row["logCE"] and row["Conductivity"])
                                 for row in observations), 8)
            self.assertIsNone(bo.run(args))
            self.assertFalse((output / "recommendation_2.csv").exists())

            simu_experiment.run(SimpleNamespace(output=output, recommendations=None,
                                                seed=27))
            bo.run(args)
            _, second = read_rows(output / "recommendation_2.csv")
            self.assertEqual(len(second), 2)
            self.assertTrue((output / "candidate_predictions_2.csv").exists())
            _, observations = read_rows(output / "observation.csv")
            self.assertEqual(sum(bool(row["logCE"] and row["Conductivity"])
                                 for row in observations), 10)
            self.assertEqual(len({row["sample_id"] for row in observations}), len(observations))
            self.assertTrue(set(row["sample_id"] for row in first).isdisjoint(
                row["sample_id"] for row in second))
            summary = json.loads((output / "training_summary.json").read_text())
            self.assertEqual(summary["round"], 2)
            self.assertEqual(summary["n_observed"], 10)
            self.assertEqual(summary["ref_point"], first_ref)


if __name__ == "__main__":
    unittest.main()
