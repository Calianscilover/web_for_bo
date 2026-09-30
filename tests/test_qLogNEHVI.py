"""Multi-objective feedback and second-round integration checks."""
import csv
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from test_support import small_pool
import qLogNEHVI as bo
import simu_experiment
from bo_utils import load_experiment_data, save_csv
from visualize_bo import build_visualization


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

    def test_four_objectives_roundtrip(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, pool, experiments, pool_rows = small_pool(root)
            fields = list(experiments[0]) + ["Viscosity", "Cost"]
            rows = []
            for i, row in enumerate(experiments):
                complete = row["logCE"].strip() and row["Conductivity"].strip()
                rows.append(dict(row, Viscosity=str(5 + (i * 7) % 5) if complete else "",
                                 Cost=str(10 - i) if complete else ""))
            experiment = root / "experiment.csv"
            save_csv(experiment, fields, rows)
            targets = ["logCE", "Conductivity", "Viscosity", "Cost"]
            data = load_experiment_data(experiment, pool, label="multi", target_names=targets)
            self.assertEqual(data.train_Y.shape, (8, 4))

            output = root / "run"
            argv = ["qLogNEHVI.py", "--experiment", str(experiment), "--pool", str(pool),
                    "--output", str(output), "--batch-size", "2", "--mc-samples", "16",
                    "--fit-maxiter", "30", "--targets", *targets,
                    "--directions", "max", "max", "min", "min"]
            with patch.object(sys, "argv", argv):
                args = bo.parse_args()
            bo.run(args)
            fields, first = read_rows(output / "recommendation_1.csv")
            self.assertEqual(fields, list(pool_rows[0]) + targets)
            self.assertEqual(len(first), 2)
            prediction_fields, _ = read_rows(output / "candidate_predictions.csv")
            self.assertEqual(len(prediction_fields), 1 + 2 * len(targets))
            summary = json.loads((output / "training_summary.json").read_text())
            self.assertEqual(summary["target_names"], targets)
            self.assertEqual(len(summary["ref_point"]), 4)

            simu_experiment.run(SimpleNamespace(output=output, recommendations=None, seed=5))
            bo.run(args)
            _, second = read_rows(output / "recommendation_2.csv")
            self.assertEqual(len(second), 2)
            self.assertTrue({row["sample_id"] for row in first}.isdisjoint(
                row["sample_id"] for row in second))
            visual = json.loads(build_visualization(output, experiment).read_text())
            self.assertEqual(visual["mode"], "multi")
            self.assertEqual(visual["targets"], targets)
            self.assertGreaterEqual(len(visual["pareto"]), 1)
            self.assertTrue(all(point["value"] >= 0 for point in visual["progress"]))

    def test_rejects_invalid_objective_counts(self):
        base = ["qLogNEHVI.py", "--targets"]
        for extra in (["A"], ["A", "B", "C", "D", "E"], ["A", "A"],
                      ["A", "B", "C", "--directions", "max", "min"],
                      ["A", "B", "C", "--ref-point", "1", "2"],
                      ["A", "B", "--noise-std", "0.1", "0.2", "0.3"]):
            with (patch.object(sys, "argv", base + extra), patch("sys.stderr", io.StringIO()),
                  self.assertRaises(SystemExit)):
                bo.parse_args()
        with patch.object(sys, "argv", base + ["A", "B", "C"]):
            self.assertEqual(bo.parse_args().directions, ["max", "max", "max"])
        with tempfile.TemporaryDirectory() as directory:
            experiment, pool, _, _ = small_pool(Path(directory))
            with self.assertRaises(ValueError):
                load_experiment_data(experiment, pool, label="multi",
                                     target_names=["logCE", "Conductivity", "a", "b", "c"])


if __name__ == "__main__":
    unittest.main()
