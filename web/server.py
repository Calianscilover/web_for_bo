"""Local UI/API for electrolyte pool generation and Bayesian optimization."""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Lock

from flask import Flask, Response, jsonify, request, send_file
from werkzeug.exceptions import HTTPException


ROOT = Path(__file__).resolve().parents[1]
WEB = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from bayesian_optimization.step_01_formulation.actions import generate_formulations  # noqa: E402
from bayesian_optimization.step_02_pool.actions import convert_formulations, import_demo  # noqa: E402
from bayesian_optimization.step_03_experiment.actions import (  # noqa: E402
    import_recipes, recipe_template, simulate_feedback as simulate_experiment,
    synchronize_feedback, validate_upload, write_table)
from bayesian_optimization.step_04_optimization.actions import (  # noqa: E402
    execute_round, refresh_visualization)

app = Flask(__name__, static_folder=str(WEB / "static"), static_url_path="/static")
app.json.ensure_ascii = False
executor = ThreadPoolExecutor(max_workers=2)
OUTPUT = ROOT / "bo_test"
ID = re.compile(r"^[a-f0-9]{12}$")
TARGET = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,63}$")
CATALOG = ROOT / "pool" / "converted" / "pool_catalog.csv"
EXPERIMENT = ROOT / "pool" / "converted" / "experiment.csv"
DESIGN_FILES = {"config_snapshot.json", "components.csv", "feasible_candidates.csv",
                "pool.csv", "pool_catalog.csv", "pool_features.csv", "pool_manifest.json",
                "experiment.csv", "experiment_mapping.csv"}
RUN_FILE = re.compile(r"^(recommendation_[1-9]\d*|candidate_predictions(?:_[1-9]\d*)?|observation|training_summary|visualization_[1-9]\d*)\.(csv|json)$")
run_locks = {}
locks_lock = Lock()


def lock_for(run_id):
    with locks_lock:
        return run_locks.setdefault(run_id, Lock())


def uid():
    return uuid.uuid4().hex[:12]


def path_for(kind, identifier):
    if not ID.fullmatch(identifier):
        raise ValueError("Invalid ID")
    path = OUTPUT / kind / identifier
    if not path.is_dir():
        raise FileNotFoundError("Task not found")
    return path


def save_json(path, value):
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temp.replace(path)


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def table(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("CSV has missing or duplicate headers")
        return list(reader.fieldnames), list(reader)


def csv_response(fields, rows, name):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    return Response("\ufeff" + stream.getvalue(), mimetype="text/csv",
                    headers={"Content-Disposition": f'attachment; filename="{name}"'})


def send_csv_file(path):
    return send_file(path, as_attachment=True, download_name=path.name)


def target_names(raw, mode=None):
    names = raw if isinstance(raw, list) else []
    if len(names) != len(set(names)) or any(not isinstance(x, str) or not TARGET.fullmatch(x) for x in names):
        raise ValueError("Target columns must be distinct simple column names")
    if mode and len(names) != (1 if mode == "single" else 2):
        raise ValueError("Select one target for single mode or two for multi mode")
    if not 1 <= len(names) <= 2:
        raise ValueError("Select one or two targets")
    return names


def job(directory, action):
    state_path = directory / "status.json"
    state = read_json(state_path)
    state.update(status="running", error=None)
    save_json(state_path, state)
    try:
        updates = action() or {}
        state.update(updates)
        state["status"] = "succeeded"
    except Exception as exc:
        state["status"] = "failed"
        state["error"] = str(exc)
        (directory / "error.log").write_text(traceback.format_exc(), encoding="utf-8")
    save_json(state_path, state)


def design_info(directory):
    state = read_json(directory / "status.json")
    state["design_id"] = directory.name
    state["files"] = sorted(name for name in DESIGN_FILES if (directory / name).exists()
                            or (directory / "converted" / name).exists())
    return state


def run_info(directory):
    state = read_json(directory / "status.json")
    state["run_id"] = directory.name
    state["files"] = sorted(p.name for p in directory.iterdir() if RUN_FILE.fullmatch(p.name))
    summary = directory / "training_summary.json"
    if summary.exists():
        state["summary"] = read_json(summary)
    round_id = state.get("round")
    if round_id:
        rec_path = directory / f"recommendation_{round_id}.csv"
        if rec_path.exists():
            _, rows = table(rec_path)
            missing = sum(not row[name].strip() for row in rows for name in state["targets"])
            state["feedback"] = {"status": "complete" if missing == 0 else "pending",
                                 "missing_values": missing, "rows": len(rows),
                                 "source": state.get("feedback_sources", {}).get(str(round_id))}
    return state


@app.errorhandler(Exception)
def error_response(exc):
    if isinstance(exc, HTTPException):
        status = exc.code
    elif isinstance(exc, FileNotFoundError):
        status = 404
    elif isinstance(exc, (ValueError, KeyError)):
        status = 400
    else:
        status = 500
    return jsonify({"code": type(exc).__name__, "message": str(exc)}), status


@app.get("/")
def home():
    return app.send_static_file("index.html")


@app.get("/api/v1/designs")
def list_designs():
    dirs = sorted((OUTPUT / "designs").glob("*/status.json"), reverse=True)
    return jsonify([design_info(x.parent) for x in dirs])


@app.post("/api/v1/designs")
def create_design():
    config = request.get_json(force=True)
    if not isinstance(config, dict):
        raise ValueError("Design must be a JSON object")
    directory = OUTPUT / "designs" / uid()
    directory.mkdir(parents=True)
    save_json(directory / "status.json", {"name": str(config.get("name", "Untitled design")),
                                           "kind": "generated", "status": "queued"})

    def action():
        stats = generate_formulations(config, directory)
        convert_formulations(directory)
        return {"counts": stats, "temperature_C": config.get("temperature_C")}

    executor.submit(job, directory, action)
    return jsonify(design_info(directory)), 202


@app.post("/api/v1/designs/demo")
def create_demo():
    directory = OUTPUT / "designs" / uid()
    directory.mkdir(parents=True)
    count = import_demo(directory, CATALOG, EXPERIMENT)
    save_json(directory / "status.json", {"name": "七元电解液示例", "kind": "demo",
                                           "status": "succeeded", "temperature_C": None,
                                           "counts": {"feasible": count}})
    return jsonify(design_info(directory)), 201


@app.get("/api/v1/designs/<design_id>")
def get_design(design_id):
    return jsonify(design_info(path_for("designs", design_id)))


@app.get("/api/v1/designs/<design_id>/optimization-runs")
def list_runs(design_id):
    path_for("designs", design_id)
    dirs = sorted((OUTPUT / "runs").glob(f"*/design_{design_id}"),
                  key=lambda p: p.stat().st_mtime, reverse=True)
    return jsonify([run_info(x.parent) for x in dirs])


@app.get("/api/v1/designs/<design_id>/candidates")
def candidates(design_id):
    directory = path_for("designs", design_id)
    fields, rows = table(directory / "converted" / "pool_catalog.csv")
    page = max(1, int(request.args.get("page", 1)))
    size = min(50, max(1, int(request.args.get("size", 10))))
    start = (page - 1) * size
    return jsonify({"fields": fields, "rows": rows[start:start + size],
                    "total": len(rows), "page": page, "size": size})


@app.get("/api/v1/designs/<design_id>/files/<name>")
def design_file(design_id, name):
    if name not in DESIGN_FILES:
        raise ValueError("Unknown file")
    directory = path_for("designs", design_id)
    path = directory / name
    if not path.exists():
        path = directory / "converted" / name
    if not path.is_file():
        raise FileNotFoundError("File not found")
    return send_csv_file(path)


@app.get("/api/v1/designs/<design_id>/experiment-template")
def experiment_template(design_id):
    names = target_names(request.args.getlist("target"))
    fields, rows = recipe_template(path_for("designs", design_id), names, source="pool")
    return csv_response(fields, rows, "experiment_template.csv")


@app.post("/api/v1/designs/<design_id>/observations")
def upload_observations(design_id):
    directory = path_for("designs", design_id)
    payload = request.get_json(force=True)
    names = target_names(payload.get("targets"))
    fields, rows, complete = validate_upload(payload.get("csv"),
                                             directory / "converted" / "pool_catalog.csv", names)
    if any((OUTPUT / "runs").glob(f"*/design_{design_id}")):
        raise ValueError("A run already uses this experiment; create a new design to replace it")
    write_table(directory / "experiment.csv", fields, rows)
    return jsonify({"rows": len(rows), "complete": complete, "targets": names})


@app.get("/api/v1/designs/<design_id>/recipe-template")
def recipe_template_csv(design_id):
    names = target_names(request.args.getlist("target"))
    fields, rows = recipe_template(path_for("designs", design_id), names)
    return csv_response(fields, rows, "recipe_template.csv")


@app.post("/api/v1/designs/<design_id>/experiment-recipes")
def upload_recipes(design_id):
    directory = path_for("designs", design_id)
    payload = request.get_json(force=True)
    names = target_names(payload.get("targets"))
    state = read_json(directory / "status.json")
    if state.get("status") != "succeeded":
        raise ValueError("The candidate pool is not ready")
    if any((OUTPUT / "runs").glob(f"*/design_{design_id}")):
        raise ValueError("A run already uses this experiment; create a new design to replace it")
    summary = import_recipes(directory, payload.get("csv"), names)
    state.setdefault("counts", {})["feasible"] = summary["pool_size"]
    save_json(directory / "status.json", state)
    return jsonify(dict(summary, targets=names))


def execute_bo(run_dir):
    state = read_json(run_dir / "status.json")
    design_dir = path_for("designs", state["design_id"])
    return execute_round(run_dir, design_dir, state)


def sync_feedback(run_dir, state):
    design_dir = path_for("designs", state["design_id"])
    synchronize_feedback(run_dir, design_dir, state["targets"], state["round"])


@app.post("/api/v1/designs/<design_id>/optimization-runs")
def create_run(design_id):
    directory = path_for("designs", design_id)
    if not (directory / "experiment.csv").exists():
        raise ValueError("Upload experiment CSV first")
    payload = request.get_json(force=True)
    mode = payload.get("mode")
    if mode not in ("single", "multi"):
        raise ValueError("mode must be single or multi")
    names = target_names(payload.get("targets"), mode)
    experiment_fields, experiments = table(directory / "experiment.csv")
    if any(x not in experiment_fields for x in names):
        raise ValueError("Selected targets are not present in uploaded experiment CSV")
    complete = sum(all(row[name].strip() for name in names) for row in experiments)
    if complete < 2:
        raise ValueError("At least two complete experimental rows are needed")
    directions = payload.get("directions", ["max"] * len(names))
    if len(directions) != len(names) or any(x not in ("max", "min") for x in directions):
        raise ValueError("Directions must be max or min for each target")
    settings = {}
    for key, lower, upper, default in (("batch_size", 1, 20, 3), ("mc_samples", 16, 2048, 128),
                                        ("fit_maxiter", 10, 1000, 100), ("seed", 0, 2**31 - 1, 2026)):
        value = int(payload.get(key, default))
        if not lower <= value <= upper:
            raise ValueError(f"{key} must be between {lower} and {upper}")
        settings[key] = value
    run_dir = OUTPUT / "runs" / uid()
    run_dir.mkdir(parents=True)
    state = {"status": "queued", "design_id": design_id, "mode": mode,
             "targets": names, "directions": directions, **settings}
    save_json(run_dir / "status.json", state)
    (run_dir / f"design_{design_id}").touch()
    executor.submit(job, run_dir, lambda: execute_bo(run_dir))
    return jsonify(run_info(run_dir)), 202


@app.get("/api/v1/optimization-runs/<run_id>")
def get_run(run_id):
    return jsonify(run_info(path_for("runs", run_id)))


@app.get("/api/v1/optimization-runs/<run_id>/recommendations/<int:round_id>")
def recommendation_preview(run_id, round_id):
    path = path_for("runs", run_id) / f"recommendation_{round_id}.csv"
    fields, rows = table(path)
    return jsonify({"fields": fields, "rows": rows, "round": round_id})


@app.get("/api/v1/optimization-runs/<run_id>/visualization")
def visualization(run_id):
    run_dir = path_for("runs", run_id)
    state = read_json(run_dir / "status.json")
    round_id = int(request.args.get("round", state.get("round", 0)))
    if round_id < 1:
        raise ValueError("No completed model round")
    path = run_dir / f"visualization_{round_id}.json"
    if not path.exists():
        if round_id == state.get("round") and (run_dir / "model.pt").exists():
            refresh_visualization(run_dir)
        else:
            raise FileNotFoundError("Visualization has not been generated")
    return jsonify(read_json(path))


@app.get("/api/v1/optimization-runs/<run_id>/files/<name>")
def run_file(run_id, name):
    if not RUN_FILE.fullmatch(name):
        raise ValueError("Unknown file")
    path = path_for("runs", run_id) / name
    if not path.is_file():
        raise FileNotFoundError("File not found")
    return send_csv_file(path)


@app.post("/api/v1/optimization-runs/<run_id>/feedback")
def upload_feedback(run_id):
    run_dir = path_for("runs", run_id)
    with lock_for(run_id):
        state = read_json(run_dir / "status.json")
        if state["status"] != "succeeded":
            raise ValueError("Wait for the current run to finish")
        round_id = state["round"]
        rec_path = run_dir / f"recommendation_{round_id}.csv"
        if not rec_path.exists():
            raise ValueError("No recommendation available")
        _, recommendations = table(rec_path)
        if not recommendations:
            raise ValueError("Recommendation batch is empty")
        if any(any(row[name].strip() for name in state["targets"]) for row in recommendations):
            raise ValueError("This recommendation already contains feedback")
        expected = {row["sample_id"] for row in recommendations}
        payload = request.get_json(force=True)
        if not isinstance(payload, dict):
            raise ValueError("Feedback must be a JSON object")
        csv_text = payload.get("csv")
        if csv_text is None:
            measurements = payload.get("measurements")
            if not isinstance(measurements, list) or len(measurements) != len(recommendations):
                raise ValueError("Provide one measured row for every recommendation")
            entered = {}
            for item in measurements:
                if not isinstance(item, dict) or item.get("sample_id") not in expected:
                    raise ValueError("Unknown recommended sample_id in feedback")
                if item["sample_id"] in entered:
                    raise ValueError("Duplicate sample_id in feedback")
                values = item.get("values")
                if not isinstance(values, dict) or set(values) != set(state["targets"]):
                    raise ValueError("Each recommended row needs all selected target values")
                entered[item["sample_id"]] = values
            stream = io.StringIO()
            rec_fields = list(recommendations[0])
            writer = csv.DictWriter(stream, fieldnames=rec_fields)
            writer.writeheader()
            for row in recommendations:
                writer.writerow({**row, **entered[row["sample_id"]]})
            csv_text = stream.getvalue()
        design_dir = path_for("designs", state["design_id"])
        fields, rows, complete = validate_upload(csv_text,
            design_dir / "converted" / "pool_catalog.csv", state["targets"], expected)
        if complete != len(rows):
            raise ValueError("Every recommended target must be filled")
        if len(rows) != len(recommendations):
            raise ValueError("Feedback row count does not match recommendation")
        audit = run_dir / "feedback_uploads"
        audit.mkdir(exist_ok=True)
        (audit / f"round_{round_id}_uploaded.csv").write_text(csv_text, encoding="utf-8")
        write_table(rec_path, fields, rows)
        state.setdefault("feedback_sources", {})[str(round_id)] = "real"
        save_json(run_dir / "status.json", state)
        sync_feedback(run_dir, state)
        refresh_visualization(run_dir)
        return jsonify({"rows": len(rows), "round": round_id,
                        "feedback": "complete", "next_round_ready": True})


@app.post("/api/v1/optimization-runs/<run_id>/simulate")
def simulate_feedback(run_id):
    run_dir = path_for("runs", run_id)
    with lock_for(run_id):
        state = read_json(run_dir / "status.json")
        if state["status"] != "succeeded":
            raise ValueError("Wait for the current run to finish")
        _, recommendations = table(run_dir / f"recommendation_{state['round']}.csv")
        if any(any(row[name].strip() for name in state["targets"]) for row in recommendations):
            raise ValueError("This recommendation already contains feedback")
        payload = request.get_json(silent=True) or {}
        message = simulate_experiment(run_dir, int(payload.get("seed", 2026)))
        state.setdefault("feedback_sources", {})[str(state["round"])] = "simulated"
        save_json(run_dir / "status.json", state)
        sync_feedback(run_dir, state)
        refresh_visualization(run_dir)
        return jsonify({"message": message, "round": state["round"],
                        "feedback": "complete", "next_round_ready": True})


@app.post("/api/v1/optimization-runs/<run_id>/next-round")
def next_round(run_id):
    run_dir = path_for("runs", run_id)
    with lock_for(run_id):
        state = read_json(run_dir / "status.json")
        if state["status"] != "succeeded":
            raise ValueError("Wait for the current run to finish")
        _, rows = table(run_dir / f"recommendation_{state['round']}.csv")
        if any(not all(row[name].strip() for name in state["targets"]) for row in rows):
            raise ValueError("Fill all recommended target values before the next round")
        state["status"] = "queued"
        save_json(run_dir / "status.json", state)
        executor.submit(job, run_dir, lambda: execute_bo(run_dir))
        return jsonify(run_info(run_dir)), 202


def main():
    global OUTPUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    OUTPUT = args.output.resolve()
    (OUTPUT / "designs").mkdir(parents=True, exist_ok=True)
    (OUTPUT / "runs").mkdir(parents=True, exist_ok=True)
    app.run(host=args.host, port=args.port, threaded=True, use_reloader=False)


if __name__ == "__main__":
    main()
