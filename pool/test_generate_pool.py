"""Tests for presence-aware candidate sampling in generate_pool.py."""
from __future__ import annotations

import copy
import csv
import tempfile
import unittest
from pathlib import Path

import generate_pool as app

CONFIG = {
    "name": "presence test",
    "batch_mass_g": 5,
    "mass_step_g": 0.001,
    "salt": {"name": "LiPF6", "smiles": "[Li+].F[P-](F)(F)(F)(F)F",
             "final_mass_fraction_bounds": [0.08, 0.16]},
    "solvents": [
        {"name": "DME", "smiles": "COCCOC", "balance": True,
         "pool_mass_fraction_bounds": [0.3, 1.0]},
        {"name": "TTE", "smiles": "FC(F)C(F)(F)COC(F)(F)C(F)F",
         "pool_mass_fraction_bounds": [0.0, 0.6]},
    ],
    "additives": [
        {"name": "FEC", "smiles": "O=C1OCC(F)O1", "final_mass_fraction_bounds": [0.0, 0.05]},
        {"name": "VC", "smiles": "O=C1OC=CO1", "final_mass_fraction_bounds": [0.01, 0.05]},
    ],
    "sampling": {"candidate_power": 10, "seed": 7},
}


def run(config):
    with tempfile.TemporaryDirectory() as directory:
        stats = app.generate(config, directory)
        with (Path(directory) / "pool.csv").open(encoding="utf-8-sig", newline="") as stream:
            return stats, list(csv.DictReader(stream))


def grams(row, name):
    return float(row[f"{name}_mass_g"])


class PresenceSamplingTests(unittest.TestCase):
    def test_optional_components_are_sometimes_absent(self):
        stats, rows = run(CONFIG)
        self.assertEqual(stats["presence_patterns"], 4)
        for name in ("TTE", "FEC"):
            absent = sum(grams(row, name) == 0 for row in rows) / len(rows)
            self.assertGreater(absent, 0.15, name)
            self.assertLess(absent, 0.35, name)
            self.assertTrue(all(grams(row, name) >= 0.001 - 1e-12
                                for row in rows if grams(row, name) > 0))
        for name in ("LiPF6", "DME", "VC"):
            self.assertTrue(all(grams(row, name) > 0 for row in rows), name)
        patterns = {row["presence_pattern"] for row in rows}
        self.assertIn("LiPF6;DME;VC", patterns)
        for row in rows:
            self.assertAlmostEqual(sum(grams(row, n) for n in ("LiPF6", "DME", "TTE", "FEC", "VC")), 5)
            self.assertTrue(0.08 - 1e-12 <= float(row["LiPF6_mass_fraction"]) <= 0.16 + 1e-12)
            self.assertLessEqual(float(row["FEC_mass_fraction"]), 0.05 + 1e-12)
            solvent = grams(row, "DME") + grams(row, "TTE")
            self.assertGreaterEqual(grams(row, "DME") / solvent, 0.3 - 1e-9)
            self.assertLessEqual(grams(row, "TTE") / solvent, 0.6 + 1e-9)

    def test_zero_absence_keeps_nearly_every_component(self):
        config = copy.deepcopy(CONFIG)
        config["sampling"]["absence_probability"] = 0
        _, rows = run(config)
        full = sum(row["presence_pattern"] == "LiPF6;DME;TTE;FEC;VC" for row in rows)
        self.assertGreater(full / len(rows), 0.98)

    def test_cosolvent_is_restored_when_balance_cannot_fill_pool(self):
        config = copy.deepcopy(CONFIG)
        config["solvents"][0]["pool_mass_fraction_bounds"] = [0.3, 0.7]
        config["sampling"]["absence_probability"] = 0.5
        _, rows = run(config)
        self.assertTrue(all(grams(row, "TTE") > 0 for row in rows))
        self.assertTrue(any(grams(row, "FEC") == 0 for row in rows))

    def test_minimum_nonzero_mass_is_respected(self):
        config = copy.deepcopy(CONFIG)
        config["additives"][0]["minimum_nonzero_mass_g"] = 0.05
        config["solvents"][1]["minimum_nonzero_mass_g"] = 0.2
        _, rows = run(config)
        fec = [grams(row, "FEC") for row in rows if grams(row, "FEC") > 0]
        tte = [grams(row, "TTE") for row in rows if grams(row, "TTE") > 0]
        self.assertTrue(fec and min(fec) >= 0.05 - 1e-12)
        self.assertTrue(tte and min(tte) >= 0.2 - 1e-12)

    def test_rejects_invalid_absence_probability(self):
        for value in (-0.1, 1, float("nan")):
            config = copy.deepcopy(CONFIG)
            config["sampling"]["absence_probability"] = value
            with self.assertRaises(ValueError):
                app.validate(config)


if __name__ == "__main__":
    unittest.main()
