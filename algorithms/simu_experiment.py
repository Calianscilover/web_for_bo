"""Simulate feedback by filling empty target cells in an issued recommendation CSV."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from bo_utils import existing_recommendation_path, save_csv


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        return reader.fieldnames, list(reader)


def run(args):
    if args.recommendations is not None:
        path = args.recommendations
        output = path.parent
    else:
        output = args.output
        round_id = 1
        if not existing_recommendation_path(output, round_id).exists():
            raise ValueError(f"No recommendation file in {output}")
        while existing_recommendation_path(output, round_id + 1).exists():
            round_id += 1
        path = existing_recommendation_path(output, round_id)
    summary = json.loads((output / "training_summary.json").read_text(encoding="utf-8"))
    targets = summary["target_names"]
    fields, recommendations = read_csv(path)
    observation_fields, observations = read_csv(output / "observation.csv")
    if not recommendations or any(name not in fields or name not in observation_fields
                                  for name in targets):
        raise ValueError("Recommendation or observation targets are missing")
    observed = {row["sample_id"]: row for row in observations}
    ranges = {}
    for name in targets:
        values = [float(row[name]) for row in observations if row[name].strip()]
        if not values:
            raise ValueError(f"No measured values for {name}")
        low, high = min(values), max(values)
        if low == high:
            width = max(abs(low) * 0.1, 1.0)
            low, high = low - width, high + width
        ranges[name] = low, high
    rng = np.random.default_rng(args.seed)
    generated = 0
    for row in recommendations:
        prior = observed.get(row["sample_id"], {})
        for name in targets:
            if row[name].strip():
                continue
            if prior.get(name, "").strip():
                row[name] = prior[name]
            else:
                low, high = ranges[name]
                row[name] = f"{rng.uniform(low, high):.10g}"
                generated += 1
    save_csv(path, fields, recommendations)
    print(f"Filled {generated} simulated measurements in {path}")
    return path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--output", type=Path,
                        help="Run directory; simulate its latest recommendation batch")
    source.add_argument("--recommendations", type=Path,
                        help="Specific recommendation CSV to fill")
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
