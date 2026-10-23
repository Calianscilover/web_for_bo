"""Flask endpoints for the independent offline replay workspace."""
import os
from pathlib import Path
import re
import subprocess
import sys
import uuid

from flask import Blueprint, jsonify, request, send_file

from offline_benchmark import dataset as ds
from offline_benchmark.runner import aggregate, code_fingerprint


def register(app, output_dir, model_settings, lock_for):
    bp = Blueprint("offline", __name__)
    root = Path(__file__).resolve().parents[1]
    processes = {}

    def body():
        payload = request.get_json(force=True)
        if not isinstance(payload, dict):
            raise ds.DataError("请求内容必须为 JSON 对象")
        return payload

    def path(kind, identifier):
        if not re.fullmatch(r"[a-f0-9]{12}", identifier):
            raise ds.DataError("ID 无效")
        directory = output_dir() / "offline" / kind / identifier
        if not directory.is_dir():
            raise FileNotFoundError("离线任务或数据集不存在")
        return directory

    def create(text, name):
        directory = output_dir() / "offline" / "datasets" / uuid.uuid4().hex[:12]
        return jsonify(ds.upload(directory, text, name)), 201

    def dispatch(directory):
        with lock_for("offline_worker_launch"):
            _dispatch(directory)

    def _dispatch(directory):
        # One replay worker at a time; online requests use their original executor.
        for prior in (output_dir() / "offline" / "benchmarks").glob("*/status.json"):
            if prior.parent == directory:
                continue
            existing = ds.read_json(prior)
            if existing.get("status") in ("queued", "running") and existing.get("pid"):
                try:
                    os.kill(existing["pid"], 0)
                except ProcessLookupError:
                    existing.update(status="failed", error="之前的 worker 已停止，可恢复")
                    ds.save_json(prior, existing)
                else:
                    state = ds.read_json(directory / "status.json")
                    state.update(status="failed", error="已有离线任务运行；结束后可恢复此任务")
                    ds.save_json(directory / "status.json", state)
                    raise ds.DataError("已有离线比较运行，请先等待结束或取消")
        log = (directory / "worker.log").open("a")
        try:
            process = subprocess.Popen([sys.executable, "-B", str(root / "offline_benchmark" / "runner.py"),
                                        "--benchmark", str(directory)], cwd=root, stdout=log, stderr=log,
                                       start_new_session=True)
        finally:
            log.close()
        processes[directory.name] = process
        state = ds.read_json(directory / "status.json")
        state["pid"] = process.pid
        ds.save_json(directory / "status.json", state)

    @bp.get("/offline")
    def page():
        return app.send_static_file("offline.html")

    @bp.post("/api/v1/offline-datasets")
    def upload():
        payload = body()
        if payload.get("design_id"):
            identifier = payload["design_id"]
            if not re.fullmatch(r"[a-f0-9]{12}", identifier):
                raise ds.DataError("设计 ID 无效")
            source = output_dir() / "designs" / identifier / "experiment.csv"
            return create(source.read_text(encoding="utf-8-sig"), "existing_experiments.csv")
        return create(payload.get("csv"), payload.get("name", "uploaded.csv"))

    @bp.post("/api/v1/offline-datasets/example")
    def example():
        return create((root / "web" / "examples" / "crr_example.csv").read_text(encoding="utf-8-sig"), "CRR_41_formulations.csv")

    @bp.get("/api/v1/offline-datasets/<identifier>")
    def dataset_info(identifier):
        return jsonify(ds.read_json(path("datasets", identifier) / "status.json"))

    @bp.get("/api/v1/offline-datasets/<identifier>/files/<name>")
    def dataset_file(identifier, name):
        if name not in ("formulation_catalog.csv", "formulation_mapping.csv", "formulation_report.json"):
            raise ds.DataError("不支持下载该数据文件")
        directory = path("datasets", identifier)
        state = ds.read_json(directory / "status.json")
        if name != "formulation_report.json" and state["status"] == "uploaded":
            raise ds.DataError("请先通过配方校验")
        return send_file(directory / name, as_attachment=True, download_name=name)

    @bp.post("/api/v1/offline-datasets/<identifier>/validate-formulations")
    def formulations(identifier):
        payload = body()
        if payload.get("design_id"):
            design_id = payload["design_id"]
            if not re.fullmatch(r"[a-f0-9]{12}", design_id):
                raise ds.DataError("设计 ID 无效")
            payload["design_config"] = ds.read_json(output_dir() / "designs" / design_id / "config_snapshot.json")
        with lock_for("offline_dataset_"+identifier):
            return jsonify(ds.validate_formulations(path("datasets", identifier), payload))

    @bp.post("/api/v1/offline-datasets/<identifier>/target")
    def target(identifier):
        with lock_for("offline_dataset_"+identifier):
            return jsonify(ds.validate_target(path("datasets", identifier), body()))

    @bp.post("/api/v1/offline-datasets/<identifier>/freeze")
    def freeze(identifier):
        with lock_for("offline_dataset_"+identifier):
            return jsonify(ds.freeze(path("datasets", identifier)))

    @bp.post("/api/v1/offline-benchmarks")
    def benchmark():
        from bayesian_optimization.step_04_optimization.actions import ACQUISITIONS
        payload = body()
        directory = path("datasets", payload.get("dataset_id", ""))
        frozen = ds.verify_frozen(directory)
        count = frozen["N"]
        n_init = ds.integer(payload.get("n_init", min(8,count-2)), "初始配方数", 2, count-1)
        budget = ds.integer(payload.get("total_budget", min(24,count)), "总预算", n_init+1, count)
        repeats = ds.integer(payload.get("repeats", 20), "重复次数", 1, 100)
        methods = payload.get("methods", ["qlognei", "qlogei", "qucb", "qpi", "qsr", "random"])
        if (not isinstance(methods, list) or not methods or any(not isinstance(x, str) for x in methods) or len(methods) != len(set(methods))
                or any(x not in (*ACQUISITIONS["single"], "random") for x in methods)):
            raise ds.DataError("请选择不重复的受支持单目标策略")
        settings = model_settings({**payload, "acquisition": "qucb"}, "single", ["offline_target"])
        if settings["feature_basis"] not in frozen.get("available_feature_bases", ["mass", "mole"]):
            raise ds.DataError("此数据集只提供质量特征；请选择质量分数，或补充各组分 SMILES 后重新校验")
        settings["noise_std"] = settings["noise_std"][0] if settings["noise_std"] else None
        config = {"dataset_id": directory.name, "dataset_fingerprint": frozen["fingerprint"],
                  "code_fingerprint": code_fingerprint(),
                  "target": frozen["target"], "n_init": n_init, "total_budget": budget,
                  "batch_size": settings["batch_size"], "repeats": repeats, "methods": methods,
                  "seed": settings["seed"], "feature_basis": settings["feature_basis"], "model": settings}
        destination = output_dir() / "offline" / "benchmarks" / uuid.uuid4().hex[:12]
        destination.mkdir(parents=True)
        ds.save_json(destination / "config.json", config)
        state = {"benchmark_id": destination.name, "status": "queued", "completed": 0, "failed": 0,
                 "planned": repeats*len(methods), "config": config, "current": None}
        ds.save_json(destination / "status.json", state)
        dispatch(destination)
        return jsonify(ds.read_json(destination / "status.json")), 202

    @bp.get("/api/v1/offline-benchmarks")
    def list_benchmarks():
        paths = sorted((output_dir() / "offline" / "benchmarks").glob("*/status.json"), key=lambda x:x.stat().st_mtime, reverse=True)
        return jsonify([ds.read_json(p) for p in paths[:30]])

    @bp.get("/api/v1/offline-benchmarks/<identifier>")
    def info(identifier):
        directory = path("benchmarks", identifier)
        state = ds.read_json(directory / "status.json")
        process = processes.get(identifier)
        if state["status"] in ("queued", "running") and process and process.poll() is not None:
            state.update(status="failed", error="worker 异常退出，请查看 worker.log；可恢复")
            ds.save_json(directory / "status.json", state)
        return jsonify(state)

    @bp.get("/api/v1/offline-benchmarks/<identifier>/results")
    def results(identifier):
        directory = path("benchmarks", identifier)
        filename = directory / "results.json"
        return jsonify(ds.read_json(filename) if filename.exists() else aggregate(ds.read_json(directory / "config.json"), []))

    @bp.post("/api/v1/offline-benchmarks/<identifier>/cancel")
    def cancel(identifier):
        directory = path("benchmarks", identifier)
        if ds.read_json(directory / "status.json")["status"] not in ("queued", "running"):
            raise ds.DataError("该任务当前未运行")
        (directory / "cancel.request").touch()
        return jsonify({"status": "cancel_requested", "message": "当前轮提交后停止"})

    @bp.post("/api/v1/offline-benchmarks/<identifier>/resume")
    def resume(identifier):
        directory = path("benchmarks", identifier)
        with lock_for("offline_benchmark_"+identifier):
            state = ds.read_json(directory / "status.json")
            if state["status"] not in ("cancelled", "failed", "completed_with_failures"):
                raise ds.DataError("只允许恢复已停止或失败的任务")
            config = ds.read_json(directory / "config.json")
            if config.get("code_fingerprint", code_fingerprint()) != code_fingerprint():
                raise ds.DataError("回放代码已经变化，请创建新的比较任务")
            ds.verify_frozen(path("datasets", config["dataset_id"]))
            (directory / "cancel.request").unlink(missing_ok=True)
            state["status"] = "queued"
            ds.save_json(directory / "status.json", state)
            dispatch(directory)
            return jsonify(state), 202

    @bp.get("/api/v1/offline-benchmarks/<identifier>/files/<name>")
    def download(identifier, name):
        if name not in {"config.json", "results.json", "summary.csv", "trials.csv", "worker.log", "error.log"}:
            raise FileNotFoundError("文件不在下载列表")
        return send_file(path("benchmarks", identifier) / name, as_attachment=True, download_name=name)

    app.register_blueprint(bp)
