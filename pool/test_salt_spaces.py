"""Fixed salt-space sampling and compatibility with the existing pool pipeline."""
import csv
import io
import tempfile
import unittest
from collections import Counter
from pathlib import Path

import generate_pool as app
from convert_pool import convert, parse_args
from import_experiments import RecipeError, import_recipes, template


def five_salt_config():
    salts = [
        ("LiPF6", "[Li+].F[P-](F)(F)(F)(F)F"),
        ("LiFSI", "[Li+].O=S(=O)(F)[N-]S(=O)(=O)F"),
        ("LiTFSI", "[Li+].O=S(=O)(C(F)(F)F)[N-]S(=O)(=O)C(F)(F)F"),
        ("LiBF4", "[Li+].F[B-](F)(F)F"),
        ("LiDFOB", "F[B-]1(OC(C(O1)=O)=O)F.[Li+]"),
    ]
    return {
        "name": "five salts, three solvents", "batch_mass_g": 5, "mass_step_g": .001,
        "salts": [{"name": name, "smiles": smiles, "final_mass_fraction_bounds": [.02, .25]}
                  for name, smiles in salts],
        "salt_selection": {"sizes": [1, 2]},
        "solvents": [
            {"name": "DME", "smiles": "COCCOC", "balance": True,
             "pool_mass_fraction_bounds": [.1, .8]},
            {"name": "TTE", "smiles": "FC(F)(OCC(F)(F)C(F)F)C(F)F",
             "pool_mass_fraction_bounds": [.1, .8]},
            {"name": "EC", "smiles": "O=C1OCCO1", "pool_mass_fraction_bounds": [.1, .8]},
        ],
        "sampling": {"candidate_power": 10, "seed": 2026, "require_all_solvents": True},
    }


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        return reader.fieldnames, list(reader)


def generate(config):
    with tempfile.TemporaryDirectory() as directory:
        stats = app.generate(config, directory)
        return stats, read_rows(Path(directory) / "pool.csv")[1]


class SaltSpaceTests(unittest.TestCase):
    def test_five_salts_cover_exactly_fifteen_spaces_with_total_budget(self):
        config = five_salt_config()
        stats, rows = generate(config)
        self.assertEqual((stats["generated"], stats["components"]), (1024, 8))
        self.assertEqual(len(stats["salt_spaces"]), 15)
        self.assertEqual(Counter(s["salt_count"] for s in stats["salt_spaces"]), {1: 5, 2: 10})
        self.assertEqual(sum(s["generated"] for s in stats["salt_spaces"]), 1024)
        quotas = [s["generated"] for s in stats["salt_spaces"]]
        self.assertLessEqual(max(quotas) - min(quotas), 1)
        observed = Counter()
        names = [item["name"] for item in config["salts"] + config["solvents"]]
        for row in rows:
            active = tuple(name for name in names[:5] if float(row[f"{name}_mass_g"]) > 0)
            self.assertIn(len(active), (1, 2))
            observed[active] += 1
            self.assertAlmostEqual(sum(float(row[f"{name}_mass_g"]) for name in names), 5)
            for name in names[:5]:
                fraction = float(row[f"{name}_mass_fraction"])
                self.assertTrue(fraction == 0 or .02 - 1e-12 <= fraction <= .25 + 1e-12)
            self.assertTrue(all(float(row[f"{name}_mass_g"]) > 0 for name in names[5:]))
        self.assertEqual(len(observed), 15)
        for space in stats["salt_spaces"]:
            self.assertEqual(space["feasible"], observed[tuple(space["salts"])])
        self.assertEqual(stats["feasible"], sum(s["feasible"] for s in stats["salt_spaces"]))
        self.assertEqual(stats["feasible"], stats["single_salt_candidates"] + stats["dual_salt_candidates"])

    def test_single_and_dual_only_modes_and_zero_absence(self):
        for size in (1, 2):
            config = five_salt_config()
            config["salt_selection"]["sizes"] = [size]
            config["sampling"]["absence_probability"] = 0
            stats, rows = generate(config)
            self.assertEqual(len(stats["salt_spaces"]), 5 if size == 1 else 10)
            self.assertTrue(all(sum(float(row[f"{s['name']}_mass_g"]) > 0
                                    for s in config["salts"]) == size for row in rows))

    def test_required_zero_lower_bound_solvents_and_selected_salt_minimum(self):
        config = five_salt_config()
        for solvent in config["solvents"]:
            solvent["pool_mass_fraction_bounds"] = [0, 1]
            solvent["minimum_nonzero_mass_g"] = .05
        for salt in config["salts"]:
            salt["final_mass_fraction_bounds"] = [0, .25]
            salt["minimum_nonzero_mass_g"] = .02
        _, rows = generate(config)
        for row in rows:
            for solvent in config["solvents"]:
                self.assertGreaterEqual(float(row[f"{solvent['name']}_mass_g"]), .05 - 1e-12)
            for salt in config["salts"]:
                mass = float(row[f"{salt['name']}_mass_g"])
                self.assertTrue(mass == 0 or mass >= .02 - 1e-12)

    def test_seed_reproducibility(self):
        config = five_salt_config()
        first = generate(config)
        self.assertEqual(first, generate(config))
        config["sampling"]["seed"] += 1
        self.assertNotEqual(first[1], generate(config)[1])

    def test_optional_solvents_and_additives_still_work_in_fixed_spaces(self):
        config = five_salt_config()
        config["sampling"]["require_all_solvents"] = False
        config["solvents"][0]["pool_mass_fraction_bounds"] = [.1, 1]
        for solvent in config["solvents"][1:]:
            solvent["pool_mass_fraction_bounds"] = [0, .8]
        config["additives"] = [{"name": "VC", "smiles": "O=C1OC=CO1",
                                "final_mass_fraction_bounds": [0, .05]}]
        _, rows = generate(config)
        for name in ("TTE", "EC", "VC"):
            self.assertTrue(any(float(row[f"{name}_mass_g"]) == 0 for row in rows))
            self.assertTrue(any(float(row[f"{name}_mass_g"]) > 0 for row in rows))
        self.assertTrue(all(sum(float(row[f"{salt['name']}_mass_g"]) > 0
                                for salt in config["salts"]) in (1, 2) for row in rows))

    def test_molality_range_uses_legal_combinations(self):
        config = five_salt_config()
        span = app.salt_molality_range(config)
        self.assertGreater(span["low"], 0)
        _, rows = generate(config)
        for row in rows:
            self.assertTrue(span["low"] - 1e-9 <= float(row["salt_molality_mol_kg"])
                            <= span["high"] + 1e-9)

    def test_rejects_invalid_modes_budget_and_unweighable_selected_salt(self):
        for selection in ({"sizes": []}, {"sizes": [3]}, {"sizes": [True]},
                          {"sizes": [1, 1]}, "invalid"):
            config = five_salt_config()
            config["salt_selection"] = selection
            with self.assertRaises(ValueError):
                app.validate(config)
        config = five_salt_config()
        config["sampling"]["candidate_power"] = 3
        with self.assertRaisesRegex(ValueError, "预算不足"):
            app.validate(config)
        config = five_salt_config()
        config["salts"][0]["final_mass_fraction_bounds"] = [.00021, .00039]
        with self.assertRaises(ValueError):
            app.validate(config)
        config = five_salt_config()
        config["salts"] = config["salts"][:1]
        config["salt_selection"]["sizes"] = [2]
        with self.assertRaises(ValueError):
            app.validate(config)

    def test_conversion_templates_and_import_preserve_eight_features(self):
        config = five_salt_config()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            app.generate(config, root)
            convert(parse_args(["--pool", str(root / "pool.csv"), "--config",
                                str(root / "config_snapshot.json"), "--output-dir",
                                str(root / "converted"), "--feature-basis", "mass"]))
            feature_fields, features = read_rows(root / "converted" / "pool_features.csv")
            self.assertEqual(feature_fields, [f"mass_ratio_{i}" for i in range(8)])
            self.assertTrue(any(float(row["mass_ratio_0"]) == 0 for row in features))
            fields, rows = template(root, ["Conductivity"], source="pool")
            stream = io.StringIO()
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(dict(row, Conductivity=1+i) for i, row in enumerate(rows[:3]))
            result = import_recipes(root, stream.getvalue(), ["Conductivity"])
            self.assertEqual((result["samples"], result["added"], result["out_of_bounds"]), (3, 0, 0))
            self.assertEqual(result["complete"], 3)
            # Disallowed three-salt recipes must not enter the candidate pool on upload.
            bad = {field: "" for field in fields}
            bad.update(experiment_id="invalid", Conductivity=1)
            for salt in config["salts"][:3]:
                bad[f"{salt['name']}_mass_g"] = .3
            for solvent in config["solvents"]:
                bad[f"{solvent['name']}_mass_g"] = 1
            before = (root / "pool.csv").read_bytes()
            stream = io.StringIO()
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerow(bad)
            with self.assertRaises(RecipeError):
                import_recipes(root, stream.getvalue(), ["Conductivity"])
            self.assertEqual(before, (root / "pool.csv").read_bytes())


if __name__ == "__main__":
    unittest.main()
