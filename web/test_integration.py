"""Small local end-to-end checks for the UI/API and BO scripts."""
import csv
import io
import tempfile
import time
import unittest
import sys
from pathlib import Path

import server
sys.path.insert(0, str(server.ROOT / "algorithms"))
from visualize_bo import hypervolume


CONFIG = {
    "name": "integration-test", "temperature_C": 25,
    "batch_mass_g": 5, "mass_step_g": 0.0001,
    "salt": {"name": "LiDFOB", "smiles": "F[B-]1(OC(C(O1)=O)=O)F.[Li+]",
             "final_mass_fraction_bounds": [0.1, 0.3]},
    "solvents": [
        {"name": "NDFA", "smiles": "CN(C)C(=O)C(F)(F)F", "balance": True,
         "pool_mass_fraction_bounds": [0.5, 1]},
        {"name": "TTE", "smiles": "FC(F)(OCC(F)(F)C(F)F)C(F)F",
         "pool_mass_fraction_bounds": [0, 0.5]},
    ], "sampling": {"candidate_power": 5, "seed": 42},
}


class LocalFlow(unittest.TestCase):
    def test_exact_two_objective_hypervolume(self):
        records = [{"values": [2, 1]}, {"values": [1, 2]}, {"values": [0.5, 0.5]}]
        self.assertAlmostEqual(hypervolume(records, ["max", "max"], [0, 0]), 3)
        self.assertAlmostEqual(hypervolume(records[:2], ["min", "min"], [3, 3]), 3)

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        server.OUTPUT = Path(self.temp.name)
        (server.OUTPUT / "designs").mkdir()
        (server.OUTPUT / "runs").mkdir()
        self.client = server.app.test_client()

    def tearDown(self):
        self.temp.cleanup()

    def wait(self, url):
        for _ in range(180):
            result = self.client.get(url).get_json()
            if result["status"] in {"succeeded", "failed"}:
                self.assertEqual(result["status"], "succeeded", result.get("error"))
                return result
            time.sleep(0.5)
        self.fail(f"Timed out waiting for {url}")

    def test_generation_upload_single_and_second_round(self):
        result = self.client.post("/api/v1/designs", json=CONFIG)
        self.assertEqual(result.status_code, 202, result.get_json())
        design_id = result.get_json()["design_id"]
        self.wait(f"/api/v1/designs/{design_id}")
        catalog = server.OUTPUT / "designs" / design_id / "converted" / "pool_catalog.csv"
        fields, rows = server.table(catalog)
        self.assertGreaterEqual(len(rows), 8)
        for row in rows:
            salt = float(row["mass_ratio_0"])
            solvent_total = sum(float(row[f"mass_ratio_{i}"]) for i in (1, 2))
            self.assertAlmostEqual(solvent_total, 1 - salt)
            self.assertAlmostEqual(float(row["mass_ratio_1"]) / solvent_total
                                   + float(row["mass_ratio_2"]) / solvent_total, 1)
        stream = io.StringIO()
        writer = csv.DictWriter(stream, fieldnames=fields + ["Conductivity"])
        writer.writeheader()
        for i, row in enumerate(rows[:4]):
            writer.writerow(dict(row, Conductivity=str(3 + i / 3)))
        upload = self.client.post(f"/api/v1/designs/{design_id}/observations",
                                  json={"targets": ["Conductivity"], "csv": stream.getvalue()})
        self.assertEqual(upload.status_code, 200, upload.get_json())
        self.assertEqual(upload.get_json()["complete"], 4)
        run = self.client.post(f"/api/v1/designs/{design_id}/optimization-runs",
                               json={"mode": "single", "targets": ["Conductivity"],
                                     "directions": ["max"], "batch_size": 2,
                                     "mc_samples": 16, "fit_maxiter": 10})
        self.assertEqual(run.status_code, 202, run.get_json())
        run_id = run.get_json()["run_id"]
        info = self.wait(f"/api/v1/optimization-runs/{run_id}")
        self.assertEqual(info["round"], 1)
        first_visual = self.client.get(f"/api/v1/optimization-runs/{run_id}/visualization").get_json()
        self.assertEqual(first_visual["mode"], "single")
        self.assertEqual(len(first_visual["fit"]["Conductivity"]["points"]), 4)
        self.assertEqual(len(first_visual["progress"]), 1)
        rec = self.client.get(f"/api/v1/optimization-runs/{run_id}/recommendations/1").get_json()
        self.assertEqual(len(rec["rows"]), 2)
        self.assertTrue(all(row["Conductivity"] == "" for row in rec["rows"]))
        simulation = self.client.post(f"/api/v1/optimization-runs/{run_id}/simulate", json={})
        self.assertEqual(simulation.status_code, 200, simulation.get_json())
        after_feedback = self.client.get(f"/api/v1/optimization-runs/{run_id}/visualization").get_json()
        self.assertEqual(len(after_feedback["progress"]), 2)
        self.assertEqual(after_feedback["counts"]["simulated_feedback"], 2)
        self.assertEqual(len(server.table(server.OUTPUT / "runs" / run_id / "observation.csv")[1]), 6)
        second = self.client.post(f"/api/v1/optimization-runs/{run_id}/next-round", json={})
        self.assertEqual(second.status_code, 202, second.get_json())
        info = self.wait(f"/api/v1/optimization-runs/{run_id}")
        self.assertEqual(info["round"], 2)
        second_visual = self.client.get(f"/api/v1/optimization-runs/{run_id}/visualization").get_json()
        self.assertEqual(len(second_visual["fit"]["Conductivity"]["points"]), 6)
        self.assertTrue((server.OUTPUT / "runs" / run_id / "observation.csv").exists())
        self.assertTrue((server.OUTPUT / "runs" / run_id / "recommendation_2.csv").exists())
        second_rec = self.client.get(f"/api/v1/optimization-runs/{run_id}/recommendations/2").get_json()
        entered = [{"sample_id": row["sample_id"], "values": {"Conductivity": 5 + i}}
                   for i, row in enumerate(second_rec["rows"])]
        uploaded = self.client.post(f"/api/v1/optimization-runs/{run_id}/feedback",
                                    json={"measurements": entered})
        self.assertEqual(uploaded.status_code, 200, uploaded.get_json())
        self.assertEqual(len(server.table(server.OUTPUT / "runs" / run_id / "observation.csv")[1]), 8)
        third = self.client.post(f"/api/v1/optimization-runs/{run_id}/next-round", json={})
        self.assertEqual(third.status_code, 202, third.get_json())
        info = self.wait(f"/api/v1/optimization-runs/{run_id}")
        self.assertEqual(info["round"], 3)
        third_visual = self.client.get(f"/api/v1/optimization-runs/{run_id}/visualization").get_json()
        self.assertEqual(third_visual["counts"]["simulated_feedback"], 2)
        self.assertEqual(third_visual["counts"]["real_feedback"], 2)
        self.assertEqual(len(third_visual["progress"]), 3)

    def test_multi_objective_recommendation(self):
        response = self.client.post("/api/v1/designs", json=CONFIG)
        design_id = response.get_json()["design_id"]
        self.wait(f"/api/v1/designs/{design_id}")
        fields, rows = server.table(server.OUTPUT / "designs" / design_id /
                                    "converted" / "pool_catalog.csv")
        stream = io.StringIO()
        writer = csv.DictWriter(stream, fieldnames=fields + ["Conductivity", "LCE"])
        writer.writeheader()
        for i, row in enumerate(rows[:5]):
            writer.writerow(dict(row, Conductivity=str(3 + i / 3), LCE=str(60 - i * 2)))
        upload = self.client.post(f"/api/v1/designs/{design_id}/observations",
                                  json={"targets": ["Conductivity", "LCE"],
                                        "csv": stream.getvalue()})
        self.assertEqual(upload.status_code, 200, upload.get_json())
        run = self.client.post(f"/api/v1/designs/{design_id}/optimization-runs",
                               json={"mode": "multi", "targets": ["Conductivity", "LCE"],
                                     "directions": ["max", "max"], "batch_size": 2,
                                     "mc_samples": 16, "fit_maxiter": 10})
        self.assertEqual(run.status_code, 202, run.get_json())
        run_id = run.get_json()["run_id"]
        info = self.wait(f"/api/v1/optimization-runs/{run_id}")
        self.assertEqual(info["round"], 1)
        first_visual = self.client.get(f"/api/v1/optimization-runs/{run_id}/visualization").get_json()
        self.assertEqual(first_visual["mode"], "multi")
        self.assertGreaterEqual(first_visual["progress"][0]["value"], 0)
        self.assertEqual(first_visual["reference_point"], info["summary"]["ref_point"])
        rec = self.client.get(f"/api/v1/optimization-runs/{run_id}/recommendations/1").get_json()
        self.assertEqual(len(rec["rows"]), 2)
        self.assertTrue(all(row["Conductivity"] == row["LCE"] == "" for row in rec["rows"]))
        blocked = self.client.post(f"/api/v1/optimization-runs/{run_id}/next-round", json={})
        self.assertEqual(blocked.status_code, 400)
        feedback = io.StringIO()
        writer = csv.DictWriter(feedback, fieldnames=rec["fields"])
        writer.writeheader()
        for i, row in enumerate(rec["rows"]):
            writer.writerow(dict(row, Conductivity=str(4 + i), LCE=str(55 + i)))
        uploaded = self.client.post(f"/api/v1/optimization-runs/{run_id}/feedback",
                                    json={"csv": feedback.getvalue()})
        self.assertEqual(uploaded.status_code, 200, uploaded.get_json())
        after_feedback = self.client.get(f"/api/v1/optimization-runs/{run_id}/visualization").get_json()
        self.assertEqual(after_feedback["counts"]["real_feedback"], 2)
        self.assertEqual(len(server.table(server.OUTPUT / "runs" / run_id / "observation.csv")[1]), 7)
        self.assertGreaterEqual(after_feedback["progress"][1]["value"],
                                after_feedback["progress"][0]["value"])
        second = self.client.post(f"/api/v1/optimization-runs/{run_id}/next-round", json={})
        self.assertEqual(second.status_code, 202, second.get_json())
        info = self.wait(f"/api/v1/optimization-runs/{run_id}")
        self.assertEqual(info["round"], 2)
        self.assertGreaterEqual(len(self.client.get(f"/api/v1/optimization-runs/{run_id}/visualization").get_json()["pareto"]), 1)


if __name__ == "__main__":
    unittest.main()
