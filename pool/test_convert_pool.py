"""Tests for deterministic pool conversion and dual ratio outputs."""
from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

import convert_pool as app


class PoolConversionTests(unittest.TestCase):
    def read(self, path):
        with path.open(encoding="utf-8-sig", newline="") as stream:
            return list(csv.DictReader(stream))

    def test_full_pool_has_dual_ratios_and_mole_features(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            args = app.parse_args(["--output-dir", str(output)])
            features_path, catalog_path, manifest_path = app.convert(args)
            features = self.read(features_path)
            catalog = self.read(catalog_path)
            self.assertEqual(len(features), 9656)
            self.assertEqual(len(catalog), 9656)
            self.assertEqual(list(features[0]), [f"mole_ratio_{i}" for i in range(7)])
            for index in range(7):
                self.assertIn(f"compound_{index}", catalog[0])
                self.assertIn(f"smiles_{index}", catalog[0])
                self.assertIn(f"mass_ratio_{index}", catalog[0])
                self.assertIn(f"mole_ratio_{index}", catalog[0])
            self.assertEqual(catalog[0]["compound_0"], "LiDFOB")
            self.assertEqual(catalog[0]["compound_1"], "NDFA")
            self.assertEqual(len({row["sample_id"] for row in catalog}), len(catalog))
            self.assertTrue(all(abs(sum(map(float, row.values())) - 1) < 1e-8
                                for row in features))
            manifest = json.loads(manifest_path.read_text())
            self.assertEqual(manifest["settings"]["feature_basis"], "mole")
            self.assertEqual(manifest["counts"]["samples"], 9656)

    def test_ids_are_stable_and_mass_features_are_supported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first_args = app.parse_args(["--output-dir", str(root / "first")])
            second_args = app.parse_args([
                "--output-dir", str(root / "second"), "--feature-basis", "mass"
            ])
            _, first_catalog, _ = app.convert(first_args)
            second_features, second_catalog, _ = app.convert(second_args)
            first_ids = [row["sample_id"] for row in self.read(first_catalog)]
            second_ids = [row["sample_id"] for row in self.read(second_catalog)]
            self.assertEqual(first_ids, second_ids)
            self.assertTrue(all(name.startswith("mass_ratio_")
                                for name in self.read(second_features)[0]))


if __name__ == "__main__":
    unittest.main()
