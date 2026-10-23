"""Serial benchmark worker; only revealed observations are passed to existing BO code."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import shutil
import sys
import time
import traceback
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from offline_benchmark.dataset import read_json, save_json, table, write_csv, verify_frozen


def seed_for(root, repeat, purpose, round_id=0):
    return int.from_bytes(hashlib.sha256(f"{root}:{repeat}:{purpose}:{round_id}".encode()).digest()[:4], "big")


def code_fingerprint():
    files = [*sorted((ROOT / "offline_benchmark").glob("*.py")),
             ROOT / "algorithms" / "bo_utils.py", ROOT / "algorithms" / "acquisitions.py",
             ROOT / "pool" / "convert_pool.py", ROOT / "pool" / "import_experiments.py",
             ROOT / "bayesian_optimization" / "step_03_experiment" / "actions.py"]
    return hashlib.sha256(b"".join(path.read_bytes() for path in files)).hexdigest()


def aggregate(config, trials):
    import numpy as np
    summary, curves = [], []
    planned = config["repeats"]
    for method in config["methods"]:
        rows = [x for x in trials if x["method"] == method]
        complete = [x for x in rows if x["status"] == "completed"]
        hits = [x for x in rows if x.get("first_hit_round") is not None]
        hit_budgets = sorted(x["first_hit_evaluations"] for x in hits)
        median_hit = hit_budgets[math.ceil(planned / 2)-1] if len(hits) >= math.ceil(planned / 2) else None
        item = {"method": method, "planned": planned, "completed": len(complete),
                "failed": sum(x["status"] == "failed" for x in rows), "hits": len(hits),
                "hit_rate": len(hits)/planned, "median_hit_evaluations": median_hit}
        for key in ("final_regret", "auc", "seconds"):
            vals = [x[key] for x in complete]
            item[key] = float(np.median(vals)) if vals else None
        summary.append(item)
        for budget in range(config["n_init"], config["total_budget"]+1):
            vals, best = [], []
            for trial in complete:
                point = [x for x in trial["curve"] if x["evaluations"] <= budget][-1]
                vals.append(point["regret"])
                best.append(point["best_value"])
            curves.append({"method": method, "evaluations": budget, "valid_repeats": len(vals),
                           "hit_rate": sum(x["first_hit_evaluations"] <= budget for x in hits)/planned,
                           "regret": float(np.median(vals)) if vals else None,
                           "q25": float(np.quantile(vals,.25)) if vals else None,
                           "q75": float(np.quantile(vals,.75)) if vals else None,
                           "best_value": float(np.median(best)) if best else None})
    return {"summary": summary, "curves": curves, "trials": trials,
            "completed": sum(x["status"] == "completed" for x in trials),
            "failed": sum(x["status"] == "failed" for x in trials),
            "planned": len(config["methods"])*planned}


def run_trial(benchmark, dataset, config, method, repeat, initial, progress):
    import numpy as np
    import torch
    sys.path.insert(0, str(ROOT / "algorithms"))
    from bo_utils import (load_experiment_data, fit_model_data, recommendation,
                          save_recommendations, sync_observations)
    from bayesian_optimization.step_03_experiment.actions import apply_feedback, validate_upload

    started = time.monotonic()
    child = benchmark / "methods" / method / f"repeat_{repeat:03d}"
    child.mkdir(parents=True, exist_ok=True)
    # A physical copy of X, with no target, is the sole candidate input to BO.
    catalog_path = child / "pool_catalog.csv"
    if not catalog_path.exists():
        shutil.copy2(dataset / "pool_catalog.csv", catalog_path)
    elif catalog_path.read_bytes() != (dataset / "pool_catalog.csv").read_bytes():
        raise ValueError("运行用候选池已改变，不能继续回放")
    fields, catalog = table(catalog_path)
    by_id = {row["sample_id"]: row for row in catalog}
    _, targets = table(dataset / "truth.csv")
    truth = {row["sample_id"]: float(row["target_value"]) for row in targets}
    alias = "offline_target"
    sign = 1 if config["target"]["direction"] == "max" else -1
    best_known = max(sign*x for x in truth.values())
    span = max(truth.values())-min(truth.values())
    initial_path = child / "experiment.csv"
    initial_rows = [{**by_id[sid], alias: truth[sid]} for sid in initial]
    if not initial_path.exists():
        write_csv(initial_path, fields+[alias], initial_rows)
    observed = list(initial)
    trace = [{"round": 0, "sample_id": sid, "value": truth[sid], "evaluations": len(initial),
              "source": "initial", "batch_id": "initial"} for sid in initial]
    curve = []
    previous = read_json(child / "trial.json") if (child / "trial.json").exists() else {}
    trial = {"method": method, "repeat": repeat, "status": "running", "curve": curve,
             "initial_ids": initial, "first_hit_round": None, "first_hit_evaluations": None,
             "trace": trace, "final_regret": None, "auc": None, "seconds": 0,
             "previous_attempts": previous.get("previous_attempts", []) +
                 ([{"status": previous["status"], "seconds": previous.get("seconds",0), "error":previous.get("error")}] if previous else [])}

    def checkpoint(round_id):
        best = max(sign*truth[sid] for sid in observed)
        gap = max(0.0, best_known-best)
        if gap <= 1e-12 and trial["first_hit_round"] is None:
            trial.update(first_hit_round=round_id, first_hit_evaluations=len(observed))
        curve.append({"round": round_id, "evaluations": len(observed), "best_value": best*sign,
                      "regret": gap, "normalized_regret": gap/span})
        write_csv(child / "trace.csv", ["round", "sample_id", "value", "evaluations", "source", "batch_id"], trace)
        save_json(child / "trial.json", trial)
        progress({"method": method, "repeat": repeat, "round": round_id, "evaluations": len(observed)})

    checkpoint(0)
    try:
        round_id = 1
        while len(observed) < config["total_budget"]:
            if (benchmark / "cancel.request").exists():
                trial["status"] = "cancelled"
                break
            q = min(config["batch_size"], config["total_budget"]-len(observed), len(catalog)-len(observed))
            path = child / f"recommendation_{round_id}.csv"
            if not path.exists():
                if method == "random":
                    remaining = sorted(set(by_id)-set(observed))
                    rng = random.Random(seed_for(config["seed"], repeat, "random", round_id))
                    selected = rng.sample(remaining, q)
                    write_csv(path, fields+[alias], [{**by_id[sid], alias: ""} for sid in selected])
                else:
                    # Reconstruct exactly the revealed prefix before fitting (including after resume).
                    write_csv(child / "observation.csv", fields+[alias],
                              [{**by_id[sid], alias: truth[sid]} for sid in observed])
                    data = load_experiment_data(child / "observation.csv", catalog_path,
                                                label="single", target_names=[alias], feature_basis=config["feature_basis"])
                    torch.manual_seed(seed_for(config["seed"], repeat, "fit", round_id))
                    args = SimpleNamespace(**{**config["model"], "batch_size": q, "minimize": sign == -1,
                                              "acquisition": method, "seed": seed_for(config["seed"], repeat, "acquisition", round_id)})
                    model = fit_model_data(data.train_X, data.train_Y, args)
                    indices = recommendation(data, model, args)
                    save_recommendations(path, catalog_path, data, indices)
                    save_json(child / f"training_audit_{round_id}.json",
                              {"observed_ids": observed, "n_observed": len(observed), "input_dim": data.input_dim,
                               "feature_names": data.feature_names, "seed": args.seed})
            _, batch = table(path)
            ids = [row["sample_id"] for row in batch]
            if len(ids) != q or len(set(ids)) != q or set(ids) & set(observed) or not set(ids) <= set(truth):
                raise ValueError("回放推荐存在重复、未知 ID 或超出批次预算")
            filled = [bool(row[alias].strip()) for row in batch]
            if any(filled) and not all(filled):
                raise ValueError("恢复时遇到部分提交的回传批次")
            if all(filled):
                validate_upload(path.read_text(encoding="utf-8-sig"), catalog_path, [alias], set(ids))
                if any(not math.isclose(float(row[alias]), truth[row["sample_id"]], rel_tol=0, abs_tol=1e-12) for row in batch):
                    raise ValueError("已提交回传与冻结目标表不一致")
            else:
                apply_feedback(child, catalog_path, [alias], round_id,
                               measurements=[{"sample_id": sid, "values": {alias: truth[sid]}} for sid in ids])
            sync_observations(initial_path, catalog_path, child, [alias])
            observed.extend(ids)
            trace.extend({"round": round_id, "sample_id": sid, "value": truth[sid],
                          "evaluations": len(observed), "source": "offline_lookup", "batch_id": f"round_{round_id}"} for sid in ids)
            checkpoint(round_id)
            round_id += 1
        else:
            trial["status"] = "completed"
        if trial["status"] == "completed":
            trial["final_regret"] = curve[-1]["regret"]
            trial["auc"] = sum(a["normalized_regret"]*(b["evaluations"]-a["evaluations"])
                               for a,b in zip(curve,curve[1:]))/(config["total_budget"]-config["n_init"])
    except Exception as exc:
        trial.update(status="failed", error=str(exc))
        (child / "error.log").write_text(traceback.format_exc())
    trial["seconds"] = previous.get("seconds", 0) + time.monotonic()-started
    save_json(child / "trial.json", trial)
    return trial


def run(benchmark):
    cache = benchmark / "cache"
    cache.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import torch
    torch.set_num_threads(4)
    config = read_json(benchmark / "config.json")
    dataset = benchmark.parents[1] / "datasets" / config["dataset_id"]
    state = read_json(benchmark / "status.json")
    trials = []
    try:
        if config.get("code_fingerprint", code_fingerprint()) != code_fingerprint():
            raise ValueError("回放代码已变化，请创建新比较任务以保持可重现性")
        frozen = verify_frozen(dataset)
        if frozen["fingerprint"] != config["dataset_fingerprint"]:
            raise ValueError("数据集指纹改变")
        _, catalog = table(dataset / "pool_catalog.csv")
        ids = sorted(row["sample_id"] for row in catalog)
        state.update(status="running", pid=os.getpid(), error=None)
        save_json(benchmark / "status.json", state)

        def progress(current):
            state["current"] = current
            save_json(benchmark / "status.json", state)

        for repeat in range(config["repeats"]):
            rng = random.Random(seed_for(config["seed"], repeat, "initialization"))
            initial = sorted(rng.sample(ids, config["n_init"]))
            init_dir = benchmark / "initializations"
            init_dir.mkdir(exist_ok=True)
            save_json(init_dir / f"repeat_{repeat:03d}.json", {"sample_ids": initial})
            for method in config["methods"]:
                previous = benchmark / "methods" / method / f"repeat_{repeat:03d}" / "trial.json"
                old = read_json(previous) if previous.exists() else None
                trial = old if old and old["status"] == "completed" else run_trial(benchmark, dataset, config, method, repeat, initial, progress)
                trials.append(trial)
                result = aggregate(config, trials)
                save_json(benchmark / "results.json", result)
                state.update(completed=result["completed"], failed=result["failed"])
                save_json(benchmark / "status.json", state)
                if (benchmark / "cancel.request").exists() and len(trials) < result["planned"]:
                    state["status"] = "cancelled"
                    save_json(benchmark / "status.json", state)
                    return
        result = aggregate(config, trials)
        write_csv(benchmark / "summary.csv", list(result["summary"][0]), result["summary"])
        compact = [{k:v for k,v in t.items() if k not in ("curve", "trace", "initial_ids", "previous_attempts")} for t in trials]
        keys = list(dict.fromkeys(k for row in compact for k in row))
        write_csv(benchmark / "trials.csv", keys, compact)
        state.update(status="completed" if not result["failed"] else "completed_with_failures", current=None)
    except Exception as exc:
        state.update(status="failed", error=str(exc))
        (benchmark / "error.log").write_text(traceback.format_exc())
    save_json(benchmark / "status.json", state)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", type=Path, required=True)
    run(parser.parse_args().benchmark.resolve())
