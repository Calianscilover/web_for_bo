"""Convert a sampled pool or import the existing seven-component example."""
import shutil
from pathlib import Path

from pool.convert_pool import convert, parse_args


def convert_formulations(design_dir: Path, feature_basis: str = "mass") -> Path:
    """Write converted/pool_catalog.csv, pool_features.csv and pool_manifest.json."""
    output_dir = design_dir / "converted"
    args = parse_args([
        "--pool", str(design_dir / "pool.csv"),
        "--config", str(design_dir / "config_snapshot.json"),
        "--output-dir", str(output_dir),
        "--feature-basis", feature_basis,
    ])
    convert(args)
    return output_dir / "pool_catalog.csv"


def import_demo(design_dir: Path, catalog: Path, experiment: Path) -> int:
    """Copy the seven-component reference files into an isolated design directory.

    pool.csv and config_snapshot.json are copied too so recipes outside the pool can be
    imported and reconverted with the same component order.
    """
    converted = design_dir / "converted"
    converted.mkdir(parents=True, exist_ok=True)
    shutil.copy2(catalog, converted / "pool_catalog.csv")
    shutil.copy2(experiment, design_dir / "experiment.csv")
    for name in ("pool_features.csv", "pool_manifest.json"):
        source = catalog.parent / name
        if source.exists():
            shutil.copy2(source, converted / name)
    for name in ("pool.csv", "config_snapshot.json"):
        source = catalog.parent.parent / name
        if source.exists():
            shutil.copy2(source, design_dir / name)
    with catalog.open(encoding="utf-8-sig") as stream:
        return sum(1 for _ in stream) - 1
