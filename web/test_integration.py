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
                                     "mc_samples": 16, "fit_maxiter": 10,
                                     "kernel": "matern", "matern_nu": 1.5, "ard": False,
                                     "lengthscale_init": "0.3", "noise_std": ["0.05"],
                                     "feature_basis": "mole", "pool_batch_size": 64})
        self.assertEqual(run.status_code, 202, run.get_json())
        run_id = run.get_json()["run_id"]
        info = self.wait(f"/api/v1/optimization-runs/{run_id}")
        self.assertEqual(info["round"], 1)
        summary = info["summary"]
        self.assertEqual((summary["kernel"], summary["matern_nu"], summary["ard"],
                          summary["lengthscale_init"], summary["noise_std"],
                          summary["pool_batch_size"]), ("matern", 1.5, False, 0.3, 0.05, 64))
        self.assertTrue(all(name.startswith("mole_ratio_") for name in summary["feature_names"]))
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
        self.assertEqual(self.client.post(f"/api/v1/optimization-runs/{run_id}/simulate",
                                          json={}).status_code, 200)
        simulated = [record["values"] for record in self.client.get(
            f"/api/v1/optimization-runs/{run_id}/visualization").get_json()["measurements"]
            if record["source"] == "simulated"]
        self.assertEqual(len(simulated), 4)
        self.assertNotEqual(sorted(simulated[:2]), sorted(simulated[2:]))

    def test_additive_pool_and_recipe_import(self):
        config = dict(CONFIG, additives=[{"name": "VC", "smiles": "O=C1OC=CO1",
                                          "final_mass_fraction_bounds": [0.01, 0.05]}])
        design_id = self.client.post("/api/v1/designs", json=config).get_json()["design_id"]
        self.wait(f"/api/v1/designs/{design_id}")
        directory = server.OUTPUT / "designs" / design_id
        fields, rows = server.table(directory / "converted" / "pool_catalog.csv")
        self.assertEqual(rows[0]["compound_3"], "VC")
        for row in rows:
            self.assertEqual(row["functional_additives"], "VC")
            self.assertTrue(0.01 - 1e-12 <= float(row["mass_ratio_3"]) <= 0.05 + 1e-12)
            solvent_total = float(row["mass_ratio_1"]) + float(row["mass_ratio_2"])
            self.assertTrue(0.5 - 1e-12 <= float(row["mass_ratio_1"]) / solvent_total <= 1)

        template = self.client.get(f"/api/v1/designs/{design_id}/recipe-template?target=Conductivity")
        template_text = template.get_data(as_text=True).lstrip("\ufeff")
        header, *examples = template_text.strip().splitlines()
        self.assertEqual(header, "experiment_id,LiDFOB_mass_g,NDFA_mass_g,TTE_mass_g,VC_mass_g,Conductivity")
        self.assertEqual([line.split(",")[0] for line in examples], ["EXAMPLE-01", "EXAMPLE-02"])
        untouched = self.client.post(f"/api/v1/designs/{design_id}/experiment-recipes",
                                     json={"targets": ["Conductivity"], "csv": template_text})
        self.assertEqual(untouched.status_code, 400)
        self.assertEqual(untouched.get_json()["message"], "没有可导入的实验数据")
        self.assertIn("2 行 EXAMPLE", untouched.get_json()["details"][0])
        _, pool = server.table(directory / "pool.csv")
        existing = [pool[0][f"{name}_mass_g"] for name in ("LiDFOB", "NDFA", "TTE", "VC")]
        recipes = "\n".join([
            template_text.strip(),
            "E1," + ",".join(existing) + ",3.5",
            "E2,1.0,2.8,1.0,0.2,4",
            "E2b,1.0,2.8,1.0,0.2,6",
            "E3,2.0,3.0,,,2.5",
            ",,,,,",
            "",
        ])
        unknown = self.client.post(f"/api/v1/designs/{design_id}/experiment-recipes",
                                   json={"targets": ["Conductivity"],
                                         "csv": "experiment_id,LiDFOB_mass_g,FEC_mass_g,Conductivity\nX,1,1,1\n"})
        self.assertEqual(unknown.status_code, 400)
        self.assertIn("FEC_mass_g", unknown.get_json()["details"][0])
        self.assertIn("LiDFOB, NDFA, TTE, VC", unknown.get_json()["details"][1])
        invalid = self.client.post(f"/api/v1/designs/{design_id}/experiment-recipes", json={
            "targets": ["Conductivity"],
            "csv": "experiment_id,LiDFOB_mass_g,NDFA_mass_g,TTE_mass_g,VC_mass_g,Conductivity\n"
                   "A,1,abc,1,0,3\nB,1,-2,1,0,3\nC,,,,,3\nD,1,1,1,0,3,9\nE,1,1,1,0,3\n"}).get_json()
        self.assertEqual(invalid["message"], "有 4 处数据不符合规范，本次未导入任何数据")
        self.assertEqual([detail.split("：")[0] for detail in invalid["details"][:4]],
                         ["第 2 行 NDFA_mass_g", "第 3 行 NDFA_mass_g", "第 4 行", "第 5 行"])
        self.assertFalse((directory / "experiment.csv").exists())
        no_target = self.client.post(f"/api/v1/designs/{design_id}/experiment-recipes", json={
            "targets": ["logCE"], "csv": recipes}).get_json()
        self.assertEqual(no_target["message"], "缺少目标列：logCE")
        self.assertIn("Conductivity", no_target["details"][1])
        imported = self.client.post(f"/api/v1/designs/{design_id}/experiment-recipes",
                                    json={"targets": ["Conductivity"], "csv": recipes})
        self.assertEqual(imported.status_code, 200, imported.get_json())
        summary = imported.get_json()
        self.assertEqual((summary["recipes"], summary["samples"], summary["matched_existing"],
                          summary["added"], summary["out_of_bounds"], summary["replicates_merged"]),
                         (4, 3, 1, 2, 1, 1))
        self.assertEqual(summary["skipped_examples"], 2)
        self.assertEqual(summary["pool_size"], len(rows) + 2)
        _, mapping = server.table(directory / "experiment_mapping.csv")
        by_id = {row["experiment_id"]: row for row in mapping}
        self.assertEqual(by_id["E1"]["sample_id"], rows[0]["sample_id"])
        self.assertEqual(by_id["E1"]["pool_status"], "matched_existing")
        self.assertEqual(by_id["E2"]["sample_id"], by_id["E2b"]["sample_id"])
        self.assertEqual(by_id["E2"]["within_design_bounds"], "true")
        self.assertEqual(by_id["E3"]["within_design_bounds"], "false")
        self.assertNotEqual(by_id["E3"]["chem_group_id"], by_id["E2"]["chem_group_id"])
        _, experiment = server.table(directory / "experiment.csv")
        values = {row["sample_id"]: float(row["Conductivity"]) for row in experiment}
        self.assertEqual(values[by_id["E2"]["sample_id"]], 5)
        again = self.client.post(f"/api/v1/designs/{design_id}/experiment-recipes",
                                 json={"targets": ["Conductivity"], "csv": recipes}).get_json()
        self.assertEqual((again["added"], again["pool_size"]), (0, summary["pool_size"]))
        _, catalog = server.table(directory / "converted" / "pool_catalog.csv")
        self.assertEqual([row["sample_id"] for row in catalog[:len(rows)]],
                         [row["sample_id"] for row in rows])
        run = self.client.post(f"/api/v1/designs/{design_id}/optimization-runs",
                               json={"mode": "single", "targets": ["Conductivity"],
                                     "directions": ["max"], "batch_size": 2,
                                     "mc_samples": 16, "fit_maxiter": 10})
        self.assertEqual(run.status_code, 202, run.get_json())
        info = self.wait(f"/api/v1/optimization-runs/{run.get_json()['run_id']}")
        self.assertEqual(info["round"], 1)
        blocked = self.client.post(f"/api/v1/designs/{design_id}/experiment-recipes",
                                   json={"targets": ["Conductivity"], "csv": recipes})
        self.assertEqual(blocked.status_code, 400)

    def test_pool_template_uses_recipe_layout(self):
        design_id = self.client.post("/api/v1/designs", json=CONFIG).get_json()["design_id"]
        self.wait(f"/api/v1/designs/{design_id}")
        directory = server.OUTPUT / "designs" / design_id
        text = self.client.get(f"/api/v1/designs/{design_id}/experiment-template?target=Conductivity"
                               ).get_data(as_text=True).lstrip("\ufeff")
        rows = list(csv.DictReader(io.StringIO(text)))
        self.assertEqual(list(rows[0]), ["experiment_id", "LiDFOB_mass_g", "NDFA_mass_g",
                                         "TTE_mass_g", "Conductivity"])
        _, catalog = server.table(directory / "converted" / "pool_catalog.csv")
        self.assertEqual([row["experiment_id"] for row in rows], [row["sample_id"] for row in catalog])
        for i, row in enumerate(rows[:3]):
            row["Conductivity"] = str(3 + i)
        stream = io.StringIO()
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
        imported = self.client.post(f"/api/v1/designs/{design_id}/experiment-recipes",
                                    json={"targets": ["Conductivity"], "csv": stream.getvalue()})
        self.assertEqual(imported.status_code, 200, imported.get_json())
        summary = imported.get_json()
        self.assertEqual((summary["matched_existing"], summary["added"], summary["complete"],
                          summary["skipped_unmeasured"], summary["pool_size"]),
                         (3, 0, 3, len(rows) - 3, len(rows)))
        _, experiment = server.table(directory / "experiment.csv")
        self.assertEqual([row["sample_id"] for row in experiment],
                         [row["sample_id"] for row in catalog[:3]])
        seven = (server.ROOT / "pool" / "converted" / "experiment.csv").read_text(encoding="utf-8-sig")
        foreign = self.client.post(f"/api/v1/designs/{design_id}/experiment-recipes",
                                   json={"targets": ["Conductivity"], "csv": seven})
        self.assertEqual(foreign.status_code, 400)
        self.assertIn("FEC", foreign.get_json()["details"][0])

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
        base = {"mode": "multi", "targets": ["Conductivity", "LCE"],
                "directions": ["max", "max"], "batch_size": 2, "mc_samples": 16, "fit_maxiter": 10}
        for invalid, message in (({"ard": False}, "关闭 ARD"),
                                 ({"noise_std": ["0.1", ""]}, "观测噪声标准差需为每个目标"),
                                 ({"noise_std": ["0.1", "-1"]}, "观测噪声标准差必须是正数"),
                                 ({"lengthscale_init": 0}, "初始长度尺度"),
                                 ({"kernel": "linear"}, "核函数"),
                                 ({"ref_point": ["abc", "1"]}, "超体积参考点必须是数字"),
                                 ({"pool_batch_size": 8}, "候选评分批大小")):
            rejected = self.client.post(f"/api/v1/designs/{design_id}/optimization-runs",
                                        json={**base, **invalid})
            self.assertEqual(rejected.status_code, 400)
            self.assertIn(message, rejected.get_json()["message"])
        run = self.client.post(f"/api/v1/designs/{design_id}/optimization-runs",
                               json={**base, "kernel": "rbf", "noise_std": ["0.1", "0.5"],
                                     "ref_point": ["2", "-1.5"]})
        self.assertEqual(run.status_code, 202, run.get_json())
        run_id = run.get_json()["run_id"]
        info = self.wait(f"/api/v1/optimization-runs/{run_id}")
        self.assertEqual(info["round"], 1)
        first_visual = self.client.get(f"/api/v1/optimization-runs/{run_id}/visualization").get_json()
        self.assertEqual(first_visual["mode"], "multi")
        self.assertGreaterEqual(first_visual["progress"][0]["value"], 0)
        self.assertEqual(first_visual["reference_point"], info["summary"]["ref_point"])
        self.assertEqual(info["summary"]["ref_point"], [2.0, -1.5])
        self.assertEqual((info["summary"]["kernel"], info["summary"]["noise_std"]), ("rbf", [0.1, 0.5]))
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
        self.client.post(f"/api/v1/optimization-runs/{run_id}/simulate", json={})
        visual = self.client.get(f"/api/v1/optimization-runs/{run_id}/visualization").get_json()
        points = {record["sample_id"]: record["values"] for record in visual["measurements"]}
        self.assertEqual(len(points), 9)
        front = {sample_id for sample_id, (x, y) in points.items()
                 if not any(a >= x and b >= y and (a, b) != (x, y) for a, b in points.values())}
        self.assertEqual({record["sample_id"] for record in visual["pareto"]}, front)


if __name__ == "__main__":
    unittest.main()
