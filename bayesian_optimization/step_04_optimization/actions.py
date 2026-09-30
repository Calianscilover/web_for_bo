"""Subprocess boundary around the existing single/multi-objective BO scripts."""
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ALGORITHMS = ROOT / "algorithms"


def command_for(state: dict, run_dir: Path, design_dir: Path) -> list[str]:
    """Build an argument vector; no shell interpolation is used."""
    mode = state["mode"]
    command = [sys.executable, str(ALGORITHMS / ("qLogNEI.py" if mode == "single" else "qLogNEHVI.py")),
               "--experiment", str(design_dir / "experiment.csv"), "--pool",
               str(design_dir / "converted" / "pool_catalog.csv"), "--output", str(run_dir),
               "--batch-size", str(state["batch_size"]), "--mc-samples", str(state["mc_samples"]),
               "--fit-maxiter", str(state["fit_maxiter"]), "--pool-batch-size", "128",
               "--seed", str(state["seed"])]
    if mode == "single":
        command += ["--target", state["targets"][0]]
        if state["directions"][0] == "min":
            command.append("--minimize")
    else:
        command += ["--targets", *state["targets"], "--directions", *state["directions"]]
    return command


def refresh_visualization(run_dir: Path):
    """Rebuild visualization_<round>.json from the run's model and observations."""
    result = subprocess.run([sys.executable, str(ALGORITHMS / "visualize_bo.py"),
                             "--run-dir", str(run_dir)], cwd=ROOT,
                            capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError("Visualization failed: " + result.stderr[-2500:])


def execute_round(run_dir: Path, design_dir: Path, state: dict) -> dict:
    """Fit, score, recommend, then build the visualization for one BO round."""
    import json
    result = subprocess.run(command_for(state, run_dir, design_dir), cwd=ROOT,
                            capture_output=True, text=True)
    (run_dir / "algorithm.log").write_text(result.stdout + "\n" + result.stderr,
                                            encoding="utf-8")
    if result.returncode:
        raise RuntimeError(result.stderr[-3000:] or result.stdout[-3000:] or "Algorithm failed")
    summary_path = run_dir / "training_summary.json"
    if not summary_path.exists():
        raise RuntimeError("Algorithm completed without training summary")
    refresh_visualization(run_dir)
    new_round = json.loads(summary_path.read_text(encoding="utf-8")).get("round", 1)
    exhausted = not (run_dir / f"recommendation_{new_round}.csv").exists()
    return {"round": new_round, "message": result.stdout.strip(), "exhausted": exhausted}
