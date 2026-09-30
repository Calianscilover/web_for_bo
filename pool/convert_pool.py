"""Convert a formulation pool into BO features and traceable sample metadata."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

import rdkit
from rdkit import Chem
from rdkit.Chem import inchi


ROOT = Path(__file__).resolve().parent
SCRIPT_VERSION = "1.1.0"
ROLE_CATEGORIES = {
    "balance_solvent": "solvents",
    "cosolvent": "solvents",
    "primary_salt": "lithium_salts",
    "salt_additive": "lithium_salts",
    "functional_additive": "functional_additives",
}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pool", type=Path, default=ROOT / "pool.csv")
    parser.add_argument("--config", type=Path, default=ROOT / "config_snapshot.json")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "converted")
    parser.add_argument("--feature-basis", choices=["mole", "mass"], default="mole",
                        help="Ratio basis written to pool_features.csv")
    parser.add_argument("--presence-tolerance", type=float, default=1e-12)
    parser.add_argument("--sum-tolerance", type=float, default=1e-8)
    parser.add_argument("--id-decimals", type=int, default=10)
    args = parser.parse_args(argv)
    if args.presence_tolerance < 0 or not math.isfinite(args.presence_tolerance):
        parser.error("presence-tolerance must be finite and nonnegative")
    if args.sum_tolerance <= 0 or not math.isfinite(args.sum_tolerance):
        parser.error("sum-tolerance must be finite and positive")
    if args.id_decimals < 1:
        parser.error("id-decimals must be positive")
    return args


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_id(prefix, signature, length):
    return f"{prefix}_{hashlib.sha256(signature.encode('utf-8')).hexdigest()[:length]}"


def load_components(config_path):
    config = json.loads(config_path.read_text(encoding="utf-8"))
    raw_components = config.get("components")
    if not isinstance(raw_components, list) or not raw_components:
        raise ValueError("config must contain a nonempty components list")
    components = []
    seen_names = set()
    for index, item in enumerate(raw_components, start=1):
        if not isinstance(item, dict):
            raise ValueError(f"component {index} must be an object")
        missing = [key for key in ("name", "role", "smiles") if not item.get(key)]
        if missing:
            raise ValueError(f"component {index} is missing: {', '.join(missing)}")
        name = str(item["name"]).strip()
        role = str(item["role"]).strip()
        smiles = str(item["smiles"]).strip()
        if name in seen_names:
            raise ValueError(f"duplicate component name: {name}")
        if role not in ROLE_CATEGORIES:
            raise ValueError(f"unsupported role for {name}: {role}")
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"RDKit cannot parse SMILES for {name}: {smiles}")
        canonical_smiles = Chem.MolToSmiles(
            molecule, canonical=True, isomericSmiles=True
        )
        try:
            inchi_key = inchi.MolToInchiKey(molecule)
        except Exception:
            inchi_key = ""
        components.append({
            "slot": index - 1,
            "name": name,
            "role": role,
            "category": ROLE_CATEGORIES[role],
            "smiles": smiles,
            "canonical_smiles": canonical_smiles,
            "inchi_key": inchi_key,
            "chemical_key": inchi_key or canonical_smiles,
        })
        seen_names.add(name)
    return config, components


def read_pool(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("pool CSV must have nonempty unique column names")
        rows = list(reader)
    if not rows:
        raise ValueError("pool CSV is empty")
    return reader.fieldnames, rows


def required_columns(components):
    columns = {"presence_pattern"}
    for component in components:
        name = component["name"]
        columns.update({f"{name}_mass_fraction", f"{name}_mole_fraction", f"{name}_mass_g"})
    return columns


def finite_nonnegative(row, column, source_row):
    try:
        value = float(row[column])
    except (TypeError, ValueError) as error:
        raise ValueError(f"row {source_row}: {column} must be numeric") from error
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"row {source_row}: {column} must be finite and nonnegative")
    return value


def transform_rows(rows, components, feature_basis, presence_tolerance,
                   sum_tolerance, id_decimals):
    names = [component["name"] for component in components]
    feature_fields = [f"{feature_basis}_ratio_{index}" for index in range(len(components))]
    catalog_rows = []
    feature_rows = []
    sample_sources = {}
    feature_sources = {}
    chem_groups = set()

    for source_row, row in enumerate(rows, start=1):
        mass_ratios = {
            name: finite_nonnegative(row, f"{name}_mass_fraction", source_row)
            for name in names
        }
        mole_ratios = {
            name: finite_nonnegative(row, f"{name}_mole_fraction", source_row)
            for name in names
        }
        masses = {
            name: finite_nonnegative(row, f"{name}_mass_g", source_row)
            for name in names
        }
        if abs(sum(mass_ratios.values()) - 1.0) > sum_tolerance:
            raise ValueError(f"row {source_row}: mass ratios do not sum to 1")
        if abs(sum(mole_ratios.values()) - 1.0) > sum_tolerance:
            raise ValueError(f"row {source_row}: mole ratios do not sum to 1")

        active_components = []
        categories = {category: [] for category in set(ROLE_CATEGORIES.values())}
        for component in components:
            name = component["name"]
            mass_active = mass_ratios[name] > presence_tolerance
            mole_active = mole_ratios[name] > presence_tolerance
            if mass_active != mole_active:
                raise ValueError(
                    f"row {source_row}: mass/mole presence disagrees for {name}"
                )
            if mass_active and masses[name] <= 0:
                raise ValueError(f"row {source_row}: positive {name} ratio has zero mass")
            if not mass_active and masses[name] > sum_tolerance:
                raise ValueError(f"row {source_row}: zero {name} ratio has positive mass")
            if mass_active:
                active_components.append(component)
                categories[component["category"]].append(name)

        group_signature = "|".join(
            f"{component['role']}:{component['chemical_key']}"
            for component in active_components
        )
        chem_group_id = stable_id("CG", group_signature, 12)
        chem_groups.add(chem_group_id)
        ratio_signature = "|".join(
            f"{component['chemical_key']}:mass={mass_ratios[component['name']]:.{id_decimals}f}:"
            f"mole={mole_ratios[component['name']]:.{id_decimals}f}"
            for component in components
        )
        sample_id = stable_id("SMP", f"ratio_schema=v1|{ratio_signature}", 16)
        if sample_id in sample_sources:
            raise ValueError(
                f"sample ID collision or duplicate formulation at rows "
                f"{sample_sources[sample_id]} and {source_row}"
            )
        sample_sources[sample_id] = source_row

        selected_ratios = mole_ratios if feature_basis == "mole" else mass_ratios
        feature_key = tuple(selected_ratios[name] for name in names)
        if feature_key in feature_sources:
            raise ValueError(
                f"duplicate {feature_basis}-ratio feature vectors at rows "
                f"{feature_sources[feature_key]} and {source_row}"
            )
        feature_sources[feature_key] = source_row
        feature_rows.append({
            f"{feature_basis}_ratio_{index}": selected_ratios[component["name"]]
            for index, component in enumerate(components)
        })

        catalog = {
            "sample_id": sample_id,
            "chem_group_id": chem_group_id,
            "solvents": ";".join(categories["solvents"]),
            "lithium_salts": ";".join(categories["lithium_salts"]),
            "functional_additives": ";".join(categories["functional_additives"]),
            "source_row": source_row,
            "presence_pattern": row["presence_pattern"].strip(),
        }
        for index, component in enumerate(components):
            name = component["name"]
            catalog.update({
                f"compound_{index}": name,
                f"smiles_{index}": component["canonical_smiles"],
                f"mass_ratio_{index}": mass_ratios[name],
                f"mole_ratio_{index}": mole_ratios[name],
            })
        catalog_rows.append(catalog)

    return feature_fields, feature_rows, list(catalog_rows[0]), catalog_rows, chem_groups


def write_csv(path, fields, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def verify_csv(path, fields, expected_rows):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != fields:
            raise ValueError(f"written columns do not match for {path}")
        row_count = sum(1 for _ in reader)
    if row_count != expected_rows:
        raise ValueError(f"written row count does not match for {path}")


def convert(args):
    config, components = load_components(args.config)
    source_fields, rows = read_pool(args.pool)
    missing = sorted(required_columns(components) - set(source_fields))
    if missing:
        raise ValueError(f"pool CSV is missing columns: {', '.join(missing)}")
    feature_fields, feature_rows, catalog_fields, catalog_rows, chem_groups = transform_rows(
        rows=rows,
        components=components,
        feature_basis=args.feature_basis,
        presence_tolerance=args.presence_tolerance,
        sum_tolerance=args.sum_tolerance,
        id_decimals=args.id_decimals,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    features_path = args.output_dir / "pool_features.csv"
    catalog_path = args.output_dir / "pool_catalog.csv"
    manifest_path = args.output_dir / "pool_manifest.json"
    features_temp = features_path.with_suffix(".csv.tmp")
    catalog_temp = catalog_path.with_suffix(".csv.tmp")
    manifest_temp = manifest_path.with_suffix(".json.tmp")

    write_csv(features_temp, feature_fields, feature_rows)
    write_csv(catalog_temp, catalog_fields, catalog_rows)
    verify_csv(features_temp, feature_fields, len(rows))
    verify_csv(catalog_temp, catalog_fields, len(rows))
    manifest = {
        "schema_version": 2,
        "script_version": SCRIPT_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "rdkit_version": rdkit.__version__,
        "inputs": {
            "pool": {"path": str(args.pool.resolve()), "sha256": sha256_file(args.pool)},
            "config": {"path": str(args.config.resolve()), "sha256": sha256_file(args.config)},
        },
        "settings": {
            "feature_basis": args.feature_basis,
            "presence_tolerance": args.presence_tolerance,
            "sum_tolerance": args.sum_tolerance,
            "id_decimals": args.id_decimals,
        },
        "components": components,
        "counts": {
            "samples": len(rows),
            "unique_samples": len({row["sample_id"] for row in catalog_rows}),
            "chemical_groups": len(chem_groups),
        },
        "outputs": {
            "pool_features.csv": {
                "columns": feature_fields,
                "sha256": sha256_file(features_temp),
            },
            "pool_catalog.csv": {
                "columns": catalog_fields,
                "sha256": sha256_file(catalog_temp),
            },
        },
        "design_name": config.get("design_name", ""),
    }
    manifest_temp.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    os.replace(features_temp, features_path)
    os.replace(catalog_temp, catalog_path)
    os.replace(manifest_temp, manifest_path)
    return features_path, catalog_path, manifest_path


def main(argv=None):
    args = parse_args(argv)
    outputs = convert(args)
    print(f"Converted pool to {args.output_dir}: {', '.join(path.name for path in outputs)}")


if __name__ == "__main__":
    main()
