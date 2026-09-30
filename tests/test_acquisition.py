"""Batch acquisition checks on a small pool with unmeasured candidates."""
import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import test_support  # noqa: F401  (puts algorithms/ on sys.path)
from bo_utils import (fit_model_data, fit_multi_model_data, initial_reference,
                      load_experiment_data, multi_recommendation, predict_candidates,
                      recommendation, save_csv, save_recommendations)


class AcquisitionTests(unittest.TestCase):
    def test_single_and_multi_batches_are_distinct_and_export_recipes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pool = root / "pool.csv"
            experiment = root / "experiment.csv"
            fields = ["sample_id", "compound_0", "mass_ratio_0", "mass_ratio_1",
                      "logCE", "Conductivity"]
            rows = [dict(sample_id=f"recipe_{i}", compound_0=f"salt_{i}",
                         mass_ratio_0=i / 6, mass_ratio_1=1 - i / 6,
                         logCE=1 + i / 5, Conductivity=2 + i / 4)
                    for i in range(7)]
            save_csv(pool, fields, rows)
            measured = [dict(row) for row in rows]
            for row in measured[3:]:
                row["logCE"] = ""
                row["Conductivity"] = ""
            save_csv(experiment, fields, measured)
            common = dict(kernel="default", ard=True, matern_nu=2.5,
                          lengthscale_init=0.5, fit_maxiter=30,
                          batch_size=3, mc_samples=16, pool_batch_size=4,
                          seed=2026)

            single = load_experiment_data(experiment, pool, label="single")
            single_args = SimpleNamespace(**common, minimize=False, noise_std=0.1)
            single_model = fit_model_data(single.train_X, single.train_Y, single_args)
            single_mean, single_std = predict_candidates(
                single, single_model, directions=["max"], batch_size=2)
            self.assertEqual(single_mean.shape, (4, 1))
            self.assertTrue(np.isfinite(single_std).all())
            single_ids = recommendation(single, single_model, single_args)
            self.assertEqual(len(single_ids), 3)
            self.assertEqual(len(set(single_ids)), 3)
            self.assertTrue(set(single_ids).issubset(set(single.candidate_indices)))

            multi = load_experiment_data(experiment, pool, label="multi")
            multi_args = SimpleNamespace(**common, directions=["max", "min"],
                                         noise_std=[0.1, 0.1])
            multi_model = fit_multi_model_data(multi.train_X, multi.train_Y, multi_args)
            multi_mean, multi_std = predict_candidates(
                multi, multi_model, directions=multi_args.directions, batch_size=2)
            self.assertEqual(multi_mean.shape, (4, 2))
            self.assertTrue(np.isfinite(multi_std).all())
            ref_point = initial_reference(multi.train_Y.numpy(), multi_args.directions)
            multi_ids = multi_recommendation(multi, multi_model, multi_args, ref_point)
            self.assertEqual(len(multi_ids), 3)
            self.assertEqual(len(set(multi_ids)), 3)
            self.assertTrue(set(multi_ids).issubset(set(multi.candidate_indices)))

            output = root / "recommendations.csv"
            save_recommendations(output, pool, multi, multi_ids)
            with output.open(encoding="utf-8-sig", newline="") as f:
                reader = csv.DictReader(f)
                self.assertEqual(reader.fieldnames, fields)
                recommendations = list(reader)
            for record in recommendations:
                i = int(record["sample_id"].split("_")[1])
                self.assertEqual(record["compound_0"], rows[i]["compound_0"])
                self.assertEqual(record["logCE"], "")
                self.assertEqual(record["Conductivity"], "")


if __name__ == "__main__":
    unittest.main()
