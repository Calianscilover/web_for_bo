"""Build a small pool from the current experimental inputs for integration tests."""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "algorithms"))
from bo_utils import save_csv  # noqa: E402


SOURCE = ROOT / "pool/converted"


def small_pool(root, extra=4):
    experiment = SOURCE / "experiment.csv"
    with experiment.open(encoding="utf-8-sig", newline="") as f:
        experiments = list(csv.DictReader(f))
    with (SOURCE / "pool_catalog.csv").open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        fields, pool = reader.fieldnames, list(reader)
    by_id = {row["sample_id"]: row for row in pool}
    used = {row["sample_id"] for row in experiments}
    rows = [by_id[row["sample_id"]] for row in experiments]
    rows.extend([row for row in pool if row["sample_id"] not in used][:extra])
    pool_path = root / "pool.csv"
    save_csv(pool_path, fields, rows)
    return experiment, pool_path, experiments, rows
