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

    def test_salt_molality_counts_only_solvents(self):
        _, rows = run(CONFIG)
        span = app.salt_molality_range(CONFIG)
        weight, = span["molecular_weights"]
        self.assertAlmostEqual(weight, 151.905, places=2)
        for row in rows:
            expected = 1000 * grams(row, "LiPF6") / weight / (grams(row, "DME") + grams(row, "TTE"))
            value = float(row["salt_molality_mol_kg"])
            self.assertAlmostEqual(value, expected, places=8)
            self.assertTrue(span["low"] - 1e-9 <= value <= span["high"] + 1e-9)

    def test_salt_list_with_one_salt_matches_salt_object(self):
        config = copy.deepcopy(CONFIG)
        config["salts"] = [config.pop("salt")]
        self.assertEqual(run(config), run(CONFIG))

    def test_rejects_invalid_absence_probability(self):
        for value in (-0.1, 1, float("nan")):
            config = copy.deepcopy(CONFIG)
            config["sampling"]["absence_probability"] = value
            with self.assertRaises(ValueError):
                app.validate(config)


def two_salt_config(lipf6=(0.0, 0.12), lidfob=(0.0, 0.06)):
    config = copy.deepcopy(CONFIG)
    config["salts"] = [
        dict(config.pop("salt"), final_mass_fraction_bounds=list(lipf6)),
        {"name": "LiDFOB", "smiles": "F[B-]1(OC(C(O1)=O)=O)F.[Li+]",
         "final_mass_fraction_bounds": list(lidfob), "minimum_nonzero_mass_g": 0.01},
    ]
    return config


class TwoSaltSamplingTests(unittest.TestCase):
    def test_salts_are_optional_alike_and_never_both_absent(self):
        _, rows = run(two_salt_config())
        present = [(grams(row, "LiPF6") > 0, grams(row, "LiDFOB") > 0) for row in rows]
        self.assertNotIn((False, False), present)
        for j in (0, 1):
            absent = sum(not pair[j] for pair in present) / len(rows)
            self.assertGreater(absent, 0.15)
            self.assertLess(absent, 0.30)
        for row in rows:
            self.assertAlmostEqual(sum(grams(row, n) for n in
                                       ("LiPF6", "LiDFOB", "DME", "TTE", "FEC", "VC")), 5)
            self.assertLessEqual(float(row["LiPF6_mass_fraction"]), 0.12 + 1e-12)
            self.assertLessEqual(float(row["LiDFOB_mass_fraction"]), 0.06 + 1e-12)
            for name, minimum in (("LiPF6", 0.001), ("LiDFOB", 0.01)):
                self.assertTrue(grams(row, name) == 0 or grams(row, name) >= minimum - 1e-12)

    def test_salt_with_positive_lower_bound_is_always_present(self):
        _, rows = run(two_salt_config(lidfob=(0.01, 0.06)))
        self.assertTrue(all(grams(row, "LiDFOB") > 0 for row in rows))
        self.assertTrue(any(grams(row, "LiPF6") == 0 for row in rows))

    def test_zero_absence_keeps_both_salts(self):
        config = two_salt_config()
        config["sampling"]["absence_probability"] = 0
        _, rows = run(config)
        both = sum(grams(row, "LiPF6") > 0 and grams(row, "LiDFOB") > 0 for row in rows)
        self.assertGreater(both / len(rows), 0.98)

    def test_molality_sums_both_salts(self):
        config = two_salt_config()
        _, rows = run(config)
        span = app.salt_molality_range(config)
        lipf6, lidfob = span["molecular_weights"]
        for row in rows:
            expected = 1000 * (grams(row, "LiPF6") / lipf6 + grams(row, "LiDFOB") / lidfob) / (
                grams(row, "DME") + grams(row, "TTE"))
            value = float(row["salt_molality_mol_kg"])
            self.assertAlmostEqual(value, expected, places=8)
            self.assertTrue(span["low"] - 1e-9 <= value <= span["high"] + 1e-9)

    def test_rejects_invalid_salt_lists(self):
        three = two_salt_config()
        three["salts"].append({"name": "LiFSI", "smiles": "[Li+].O=S(=O)(F)[N-]S(=O)(=O)F",
                               "final_mass_fraction_bounds": [0, 0.05]})
        both = two_salt_config()
        both["salt"] = copy.deepcopy(CONFIG["salt"])
        crowded = two_salt_config(lipf6=(0.0, 0.6), lidfob=(0.0, 0.4))
        duplicate = two_salt_config()
        duplicate["salts"][1]["name"] = "LiPF6"
        tiny = two_salt_config(lidfob=(0.0, 0.001))
        for config in (three, both, crowded, duplicate, tiny):
            with self.assertRaises(ValueError):
                app.validate(config)


if __name__ == "__main__":
    unittest.main()
