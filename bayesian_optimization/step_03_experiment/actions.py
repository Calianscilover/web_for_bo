"""CSV contracts and feedback actions used by the web API."""
import csv
import io
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def read_table(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("CSV has missing or duplicate headers")
        return list(reader.fieldnames), list(reader)


def validate_upload(csv_text, catalog_path, names, expected_ids=None):
    """Require exact catalog recipe columns plus selected targets, with finite values."""
    if not isinstance(csv_text, str) or len(csv_text.encode("utf-8")) > 20_000_000:
        raise ValueError("CSV must be text under 20 MB")
    pool_fields, pool = read_table(catalog_path)
    if set(names) & set(pool_fields):
        raise ValueError("Target columns must not collide with candidate pool columns")
    reader = csv.DictReader(io.StringIO(csv_text.lstrip("\ufeff")))
    if reader.fieldnames != pool_fields + names:
        raise ValueError("CSV headers must exactly match the downloaded template")
    by_id = {row["sample_id"]: row for row in pool}
    rows = list(reader)
    if not rows:
        raise ValueError("CSV has no experiment rows")
    if len(rows) > len(pool):
        raise ValueError("CSV has more rows than the candidate pool")
    seen = set()
    complete = 0
    for i, row in enumerate(rows, 2):
        sample_id = row.get("sample_id", "")
        if None in row or None in row.values():
            raise ValueError(f"Row {i}: malformed CSV row")
        if sample_id in seen:
            raise ValueError(f"Row {i}: duplicate sample_id {sample_id}")
        if sample_id not in by_id:
            raise ValueError(f"Row {i}: sample_id absent from candidate pool: {sample_id}")
        if expected_ids is not None and sample_id not in expected_ids:
            raise ValueError(f"Row {i}: sample_id was not in the issued recommendation")
        seen.add(sample_id)
        reference = by_id[sample_id]
        for field in pool_fields:
            if row[field] != reference[field]:
                try:
                    same_numeric = math.isclose(float(row[field]), float(reference[field]),
                                                rel_tol=0, abs_tol=1e-12)
                except (ValueError, TypeError):
                    same_numeric = False
                if not same_numeric:
                    raise ValueError(f"Row {i}: recipe field {field} differs from candidate pool")
            row[field] = reference[field]
        for name in names:
            value = row[name].strip()
            if value and not math.isfinite(float(value)):
                raise ValueError(f"Row {i}: {name} must be finite")
            row[name] = value
        complete += all(row[name] for name in names)
    if expected_ids is not None and seen != expected_ids:
        raise ValueError("Feedback must contain every issued recommendation row")
    return pool_fields + names, rows, complete


def write_table(path, fields, rows):
    """Atomically replace a CSV file after the complete content is written."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def synchronize_feedback(run_dir: Path, design_dir: Path, targets: list[str], current_round: int):
    """Merge completed recommendation rows into observation.csv."""
    cache = Path(tempfile.gettempdir()) / "electrolyte_web_cache"
    (cache / "matplotlib").mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    root = Path(__file__).resolve().parents[2]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    from bo_utils import sync_observations
    _, next_round, pending = sync_observations(
        design_dir / "experiment.csv", design_dir / "converted" / "pool_catalog.csv",
        run_dir, targets)
    if pending or next_round != current_round + 1:
        raise RuntimeError("Returned observations did not synchronize correctly")


def simulate_feedback(run_dir: Path, seed: int = 2026) -> str:
    """Fill the current recommendation using the existing demo simulator."""
    root = Path(__file__).resolve().parents[2]
    result = subprocess.run([sys.executable, str(root / "simu_experiment.py"),
                             "--output", str(run_dir), "--seed", str(seed)],
                            cwd=root, capture_output=True, text=True)
    if result.returncode:
        raise ValueError(result.stderr[-2000:] or result.stdout[-2000:])
    return result.stdout.strip()
