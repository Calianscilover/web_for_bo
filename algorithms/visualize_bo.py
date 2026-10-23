"""Build model-fit and measured-optimization history data for one BO run.

The JSON is intentionally chart-library independent so the same values can be
used by the local UI or exported to another plotting program.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np
import torch
from botorch.models import ModelListGP, SingleTaskGP
from botorch.models.transforms.outcome import Standardize
from botorch.utils.multi_objective.hypervolume import Hypervolume
from botorch.utils.multi_objective.pareto import is_non_dominated
from gpytorch.kernels import MaternKernel, RBFKernel, ScaleKernel


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames:
            raise ValueError(f"Empty CSV: {path}")
        return list(reader.fieldnames), list(reader)


def model_from_checkpoint(run_dir, summary):
    checkpoint = torch.load(run_dir / "model.pt", map_location="cpu", weights_only=True)
    train_x = checkpoint["train_X"]
    train_y = checkpoint["train_Y"]
    models = []
    for j, direction in enumerate(summary["directions"]):
        covariance = None
        if summary["kernel"] != "default":
            kwargs = {"ard_num_dims": train_x.shape[-1] if summary["ard"] else None}
            base = (RBFKernel(**kwargs) if summary["kernel"] == "rbf" else
                    MaternKernel(nu=summary["matern_nu"], **kwargs))
            covariance = ScaleKernel(base)
        signed_y = train_y[:, j:j + 1] * (-1 if direction == "min" else 1)
        noise = summary.get("noise_std")
        noise_j = None if noise is None else (noise[j] if isinstance(noise, list) else noise)
        variance = None if noise_j is None else torch.full_like(signed_y, noise_j**2)
        models.append(SingleTaskGP(train_x, signed_y, train_Yvar=variance,
                                   covar_module=covariance, outcome_transform=Standardize(m=1)))
    model = models[0] if len(models) == 1 else ModelListGP(*models)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    with torch.no_grad():
        fitted = model.posterior(train_x).mean.detach().cpu().numpy()
    signs = np.array([-1 if direction == "min" else 1 for direction in summary["directions"]])
    return train_y.detach().cpu().numpy(), fitted * signs


def collected_measurements(run_dir, experiment_path, targets, feedback_sources):
    def composition(row):
        result = []
        slot = 0
        while f"compound_{slot}" in row:
            fraction = float(row[f"mass_ratio_{slot}"])
            if fraction > 0:
                result.append({"name": row[f"compound_{slot}"], "mass_fraction": fraction})
            slot += 1
        return result

    _, experiments = read_csv(experiment_path)
    records = {}
    for row in experiments:
        if all(row.get(name, "").strip() for name in targets):
            records[row["sample_id"]] = {
                "sample_id": row["sample_id"], "values": [float(row[name]) for name in targets],
                "round": 0, "source": "initial", "composition": composition(row),
            }
    rounds = sorted((int(match.group(1)), path) for path in run_dir.glob("recommendation_*.csv")
                    if (match := re.fullmatch(r"recommendation_(\d+)\.csv", path.name)))
    completed = []
    for number, path in rounds:
        _, batch = read_csv(path)
        ready = [row for row in batch if all(row.get(name, "").strip() for name in targets)]
        if len(ready) != len(batch):
            continue
        completed.append(number)
        for row in ready:
            sample_id = row["sample_id"]
            values = [float(row[name]) for name in targets]
            if sample_id in records:
                if records[sample_id]["values"] != values:
                    raise ValueError(f"Conflicting measurements for {sample_id}")
                continue
            records[sample_id] = {"sample_id": sample_id, "values": values,
                                  "round": number,
                                  "source": feedback_sources.get(str(number), "real"),
                                  "composition": composition(row)}
    if not records:
        raise ValueError("No complete measurements are available for visualization")
    return list(records.values()), completed


def hypervolume(records, directions, reference):
    signs = np.array([-1 if direction == "min" else 1 for direction in directions])
    values = torch.as_tensor([record["values"] for record in records], dtype=torch.double)
    transformed = values * torch.as_tensor(signs, dtype=torch.double)
    reference_signed = torch.as_tensor(reference, dtype=torch.double) * torch.as_tensor(signs)
    valid = (transformed > reference_signed).all(dim=1)
    if not valid.any():
        return 0.0
    pareto = transformed[valid]
    pareto = pareto[is_non_dominated(pareto)]
    return float(Hypervolume(ref_point=reference_signed).compute(pareto))


def build_visualization(run_dir, experiment_path=None):
    run_dir = Path(run_dir)
    summary = json.loads((run_dir / "training_summary.json").read_text(encoding="utf-8"))
    experiment_path = Path(experiment_path or summary["experiment"])
    targets = summary["target_names"]
    directions = summary["directions"]
    round_model = summary["round"]
    state_path = run_dir / "status.json"
    state = json.loads(state_path.read_text(encoding="utf-8")) if state_path.exists() else {}
    records, completed = collected_measurements(
        run_dir, experiment_path, targets, state.get("feedback_sources", {}))
    training_y, fitted_y = model_from_checkpoint(run_dir, summary)
    ids = summary["observed_sample_ids"]
    record_by_id = {record["sample_id"]: record for record in records}
    fit = {}
    for j, target in enumerate(targets):
        actual = training_y[:, j]
        predicted = fitted_y[:, j]
        error = predicted - actual
        variance = float(np.sum((actual - actual.mean()) ** 2))
        fit[target] = {
            "points": [{"sample_id": sample_id, "measured": float(actual[i]),
                        "predicted": float(predicted[i]),
                        "source": record_by_id.get(sample_id, {}).get("source", "initial"),
                        "round": record_by_id.get(sample_id, {}).get("round", 0)}
                       for i, sample_id in enumerate(ids)],
            "mae": float(np.abs(error).mean()),
            "r2": None if variance <= 0 else float(1 - np.sum(error**2) / variance),
        }
    progress = []
    for round_id in [0, *completed]:
        observed = [record for record in records if record["round"] <= round_id]
        if not observed:
            continue
        if len(targets) == 1:
            values = [record["values"][0] for record in observed]
            metric = min(values) if directions[0] == "min" else max(values)
        else:
            metric = hypervolume(observed, directions, summary["ref_point"])
        progress.append({"round": round_id, "value": float(metric),
                         "observed_count": len(observed),
                         "new_count": sum(record["round"] == round_id for record in observed),
                         "simulated_count": sum(record["source"] == "simulated"
                                                for record in observed)})
    pareto = []
    if len(targets) >= 2:
        signs = torch.tensor([-1 if d == "min" else 1 for d in directions])
        values = torch.tensor([record["values"] for record in records], dtype=torch.double)
        mask = is_non_dominated(values * signs)
        pareto = [record for record, included in zip(records, mask.tolist()) if included]
    rec_path = run_dir / f"recommendation_{round_model}.csv"
    recommendation_rows = read_csv(rec_path)[1] if rec_path.exists() else []
    prediction_path = run_dir / ("candidate_predictions.csv" if round_model == 1 else
                                 f"candidate_predictions_{round_model}.csv")
    predictions = read_csv(prediction_path)[1] if prediction_path.exists() else []
    predictions_by_id = {row["sample_id"]: row for row in predictions}
    recommendations = []
    for row in recommendation_rows:
        sample_id = row["sample_id"]
        predicted = predictions_by_id.get(sample_id, {})
        recommendations.append({
            "sample_id": sample_id,
            "predicted": {name: {"mean": float(predicted[f"predicted_mean_{name}"]),
                                  "std": float(predicted[f"predicted_std_{name}"])}
                          for name in targets},
            "measured": {name: float(row[name]) if row[name].strip() else None for name in targets},
            "returned": all(row[name].strip() for name in targets),
        })
    payload = {
        "schema_version": 1, "model_round": round_model,
        "mode": "single" if len(targets) == 1 else "multi",
        "targets": targets, "directions": directions,
        "reference_point": summary.get("ref_point"),
        "fit_kind": "in_sample_posterior_mean",
        "fit": fit, "progress": progress, "pareto": pareto,
        "measurements": records,
        "recommendations": recommendations,
        "counts": {"initial": sum(record["source"] == "initial" for record in records),
                   "real_feedback": sum(record["source"] == "real" for record in records),
                   "simulated_feedback": sum(record["source"] == "simulated" for record in records),
                   "pending": sum(not rec["returned"] for rec in recommendations)},
    }
    path = run_dir / f"visualization_{round_model}.json"
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                         encoding="utf-8")
    temporary.replace(path)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--experiment", type=Path, default=None)
    args = parser.parse_args()
    print(build_visualization(args.run_dir, args.experiment))


if __name__ == "__main__":
    main()
