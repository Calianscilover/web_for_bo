"""Merge measured recipes, including ones outside the sampled pool, into a design.

Input is a recipe CSV: optional ``experiment_id``, one ``<component>_mass_g`` (or
``<component>_mass_fraction``) column per design component, then the target columns.
Catalog-layout files (compound_i + mass_ratio_i, e.g. an older experiment.csv) are also
accepted when their components belong to the design.
Each recipe is matched to pool.csv by its mass fractions (10 decimals, the precision of
sample_id). Unmatched recipes are appended to pool.csv and the pool is reconverted, so
sample_id and chem_group_id come from the same hashes as every sampled candidate.
Replicates of one recipe are averaged per target. Rows whose experiment_id starts with
EXAMPLE (the template's worked examples), rows without any target value and fully blank
rows are skipped.

Outputs in the design directory: pool.csv (appended), converted/, experiment.csv
(catalog columns + targets) and experiment_mapping.csv (experiment_id -> sample_id).
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors

try:
    from pool.convert_pool import convert, parse_args
except ImportError:
    from convert_pool import convert, parse_args

MAPPING_FIELDS = ["experiment_id", "recipe_row", "sample_id", "chem_group_id", "pool_status",
                  "within_design_bounds", "replicates"]
KEY_DECIMALS = 10
BOUND_TOLERANCE = 1e-9
EXAMPLE_PREFIX = "EXAMPLE"


def read_csv_text(text, strict=True):
    """Parse CSV text; strict=False keeps ragged rows (as spreadsheets export) for parse_recipes."""
    reader = csv.DictReader(io.StringIO(text.lstrip("\ufeff")))
    fields = [field.strip() for field in reader.fieldnames or []]
    if not fields or len(fields) != len(set(fields)):
        raise ValueError("Recipe CSV has missing or duplicate headers")
    reader.fieldnames = fields
    rows = list(reader)
    if strict and any(None in row or None in row.values() for row in rows):
        raise ValueError("Recipe CSV has malformed rows")
    return fields, rows


def read_csv_file(path):
    return read_csv_text(Path(path).read_text(encoding="utf-8-sig"))


def write_csv(path, fields, rows):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def load_design(design_dir):
    config_path = design_dir / "config_snapshot.json"
    pool_path = design_dir / "pool.csv"
    if not config_path.exists() or not pool_path.exists():
        raise ValueError("Design has no pool.csv/config_snapshot.json; regenerate it or reload the demo")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    components = []
    for item in config.get("components", []):
        molecule = Chem.MolFromSmiles(str(item.get("smiles", "")))
        if molecule is None:
            raise ValueError(f"Invalid SMILES for {item.get('name')}")
        components.append({"name": str(item["name"]).strip(), "role": item["role"],
                           "molecular_weight": Descriptors.MolWt(molecule)})
    if not components:
        raise ValueError("config_snapshot.json lists no components")
    return config, components


def template_fields(design_dir, targets):
    _, components = load_design(Path(design_dir))
    return ["experiment_id", *(f"{item['name']}_mass_g" for item in components), *targets]


def template(design_dir, targets, source="examples"):
    """One recipe CSV layout for every upload; only the prefilled rows differ.

    source="pool" lists every candidate (experiment_id = sample_id, targets blank);
    source="examples" gives EXAMPLE rows: one recipe copied from the pool, one weighed
    freely. Unused components are blank; EXAMPLE and unmeasured rows are skipped on import.
    """
    design_dir = Path(design_dir)
    fields = template_fields(design_dir, targets)
    _, pool = read_csv_file(design_dir / "pool.csv")
    mass_fields = fields[1:len(fields) - len(targets)]
    if source == "pool":
        _, catalog = read_csv_file(design_dir / "converted" / "pool_catalog.csv")
        if len(catalog) != len(pool):
            raise ValueError("Candidate catalog does not match pool.csv")
        return fields, [{
            "experiment_id": entry["sample_id"],
            **{field: row[field] if float(row[field] or 0) > 0 else "" for field in mass_fields},
            **{name: "" for name in targets},
        } for row, entry in zip(pool, catalog)]
    examples = []
    for number, (row, digits) in enumerate([(pool[0], None), (pool[len(pool) // 2], 2)], 1):
        example = {"experiment_id": f"{EXAMPLE_PREFIX}-{number:02d}"}
        for field in mass_fields:
            mass = float(row[field] or 0)
            example[field] = "" if mass <= 0 else (row[field] if digits is None
                                                   else f"{mass:.{digits}f}")
        example.update({name: f"{1.23 * (i + 1) + 0.1 * number:.2f}"
                        for i, name in enumerate(targets)})
        examples.append(example)
    return fields, examples


def parse_number(raw, label):
    try:
        value = float(raw)
    except ValueError as error:
        raise ValueError(f"{label} must be numeric") from error
    if not math.isfinite(value):
        raise ValueError(f"{label} must be finite")
    return value


def catalog_to_fractions(fields, rows, names):
    """Rewrite catalog-layout rows (compound_i + mass_ratio_i) as <name>_mass_fraction rows."""
    slots = sorted(int(field[9:]) for field in fields
                   if field.startswith("compound_") and field[9:].isdigit())
    unknown, converted = set(), []
    for row in rows:
        values = {f"{name}_mass_fraction": "" for name in names}
        for slot in slots:
            name = (row.get(f"compound_{slot}") or "").strip()
            ratio = (row.get(f"mass_ratio_{slot}") or "").strip()
            if not name or not ratio or float(ratio) == 0:
                continue
            if name not in names:
                unknown.add(name)
            values[f"{name}_mass_fraction"] = ratio
        converted.append({"experiment_id": row.get("sample_id") or row.get("experiment_id") or "",
                          **values, **{k: v for k, v in row.items() if k not in values}})
    if unknown:
        raise ValueError(f"CSV contains components outside this design: {', '.join(sorted(unknown))}; "
                         "download the template of the current design or regenerate the pool")
    extra = [field for field in fields if field not in {"experiment_id"}]
    return ["experiment_id", *(f"{name}_mass_fraction" for name in names), *extra], converted


def parse_recipes(fields, rows, names, targets, batch_mass):
    if "compound_0" in fields and "mass_ratio_0" in fields:
        fields, rows = catalog_to_fractions(fields, rows, names)
    by_suffix = {suffix: [f"{name}{suffix}" for name in names]
                 for suffix in ("_mass_g", "_mass_fraction")}
    unknown = [field for field in fields
               if field.endswith(("_mass_g", "_mass_fraction"))
               and field not in by_suffix["_mass_g"] + by_suffix["_mass_fraction"]]
    if unknown:
        raise ValueError("Recipe columns name components outside this design: "
                         f"{', '.join(unknown)}; add them to the formulation space and regenerate")
    uses = [suffix for suffix, columns in by_suffix.items() if set(columns) & set(fields)]
    if len(uses) != 1:
        raise ValueError("Use either <component>_mass_g or <component>_mass_fraction columns")
    suffix = uses[0]
    missing = [name for name in targets if name not in fields]
    if missing:
        raise ValueError(f"Recipe CSV is missing target columns: {', '.join(missing)}")
    recipes = []
    skipped = {"examples": 0, "unmeasured": 0}
    for line, row in enumerate(rows, 2):
        if None in row and any((value or "").strip() for value in row[None]):
            raise ValueError(f"Row {line}: more values than header columns")
        experiment_id = (row.get("experiment_id") or "").strip()
        if experiment_id.upper().startswith(EXAMPLE_PREFIX):
            skipped["examples"] += 1
            continue
        if not any((row.get(field) or "").strip() for field in fields):
            continue
        measured = {}
        for name in targets:
            raw = (row.get(name) or "").strip()
            measured[name] = parse_number(raw, f"Row {line}: {name}") if raw else None
        if all(value is None for value in measured.values()):
            skipped["unmeasured"] += 1
            continue
        values = []
        for column in by_suffix[suffix]:
            raw = (row.get(column) or "").strip()
            value = parse_number(raw, f"Row {line}: {column}") if raw else 0.0
            if value < 0:
                raise ValueError(f"Row {line}: {column} must be nonnegative")
            values.append(value)
        total = sum(values)
        if total <= 0:
            raise ValueError(f"Row {line}: recipe has no components")
        if suffix == "_mass_fraction" and abs(total - 1) > 1e-6:
            raise ValueError(f"Row {line}: mass fractions sum to {total:.6g}, expected 1")
        fractions = [value / total for value in values]
        recipes.append({
            "row": line,
            "experiment_id": experiment_id or f"row_{line}",
            "fractions": fractions,
            "masses": values if suffix == "_mass_g" else [x * batch_mass for x in fractions],
            "targets": measured,
        })
    if not recipes:
        raise ValueError(f"CSV has no measured rows ({skipped['examples']} {EXAMPLE_PREFIX} rows and "
                         f"{skipped['unmeasured']} rows without target values were ignored); "
                         "fill the target columns of the experiments you ran")
    return recipes, skipped


def bounds_checker(config, names):
    """Return a design-bounds test for web-generated designs, or None when unknown."""
    salt, solvents = config.get("salt"), config.get("solvents")
    if not isinstance(salt, dict) or not isinstance(solvents, list):
        return None
    limits = {str(salt["name"]).strip(): salt["final_mass_fraction_bounds"]}
    limits.update({str(item["name"]).strip(): item["final_mass_fraction_bounds"]
                   for item in config.get("additives") or []})
    pool = {str(item["name"]).strip(): item["pool_mass_fraction_bounds"] for item in solvents}

    def inside(value, bounds):
        return bounds[0] - BOUND_TOLERANCE <= value <= bounds[1] + BOUND_TOLERANCE

    def check(fractions):
        value = dict(zip(names, fractions))
        solvent_total = sum(value[name] for name in pool)
        return (solvent_total > 0
                and all(inside(value[name], bounds) for name, bounds in limits.items())
                and all(inside(value[name] / solvent_total, bounds) for name, bounds in pool.items()))

    return check


def pool_key(row, names):
    return tuple(round(float(row[f"{name}_mass_fraction"]), KEY_DECIMALS) for name in names)


def pool_row(fields, components, recipe):
    moles = [x / item["molecular_weight"] for x, item in zip(recipe["fractions"], components)]
    total_moles = sum(moles)
    row = {field: "" for field in fields}
    row["presence_pattern"] = ";".join(
        item["name"] for item, x in zip(components, recipe["fractions"]) if x > 0)
    for item, mass, fraction, mole in zip(components, recipe["masses"], recipe["fractions"], moles):
        name = item["name"]
        row[f"{name}_mass_g"] = f"{mass:.10g}"
        row[f"{name}_mass_fraction"] = f"{fraction:.12g}"
        row[f"{name}_mole_fraction"] = f"{mole / total_moles:.12g}"
    return row


def reconvert(design_dir):
    converted = design_dir / "converted"
    manifest = converted / "pool_manifest.json"
    settings = json.loads(manifest.read_text(encoding="utf-8"))["settings"] if manifest.exists() else {}
    convert(parse_args([
        "--pool", str(design_dir / "pool.csv"),
        "--config", str(design_dir / "config_snapshot.json"),
        "--output-dir", str(converted),
        "--feature-basis", settings.get("feature_basis", "mass"),
        "--id-decimals", str(settings.get("id_decimals", KEY_DECIMALS)),
    ]))
    return converted / "pool_catalog.csv"


def import_recipes(design_dir, csv_text, targets):
    design_dir = Path(design_dir)
    config, components = load_design(design_dir)
    names = [item["name"] for item in components]
    reserved = {"sample_id", "chem_group_id", "solvents", "lithium_salts", "functional_additives",
                "source_row", "presence_pattern"}
    if any(name in reserved or name.startswith(("compound_", "smiles_", "mass_ratio_", "mole_ratio_"))
           for name in targets):
        raise ValueError("Target columns must not collide with candidate pool columns")
    fields, rows = read_csv_text(csv_text, strict=False)
    recipes, skipped = parse_recipes(fields, rows, names, targets,
                                     float(config.get("batch_mass_g", 5)))
    pool_path = design_dir / "pool.csv"
    original_text = pool_path.read_text(encoding="utf-8-sig")
    pool_fields, pool_rows = read_csv_text(original_text)
    original_count = len(pool_rows)
    index = {pool_key(row, names): i for i, row in enumerate(pool_rows)}
    check = bounds_checker(config, names)
    for recipe in recipes:
        row = pool_row(pool_fields, components, recipe)
        key = pool_key(row, names)
        if key not in index:
            index[key] = len(pool_rows)
            pool_rows.append(row)
        recipe["pool_index"] = index[key]
        recipe["within"] = check(recipe["fractions"]) if check else None

    if len(pool_rows) > original_count:
        write_csv(pool_path, pool_fields, pool_rows)
    try:
        catalog_path = reconvert(design_dir)
    except Exception:
        pool_path.write_text(original_text, encoding="utf-8-sig")
        raise
    catalog_fields, catalog = read_csv_file(catalog_path)
    if len(catalog) != len(pool_rows):
        raise RuntimeError("Converted catalog does not match pool.csv")

    groups = {}
    for recipe in recipes:
        groups.setdefault(recipe["pool_index"], []).append(recipe)
    experiment_rows = []
    for pool_index, group in groups.items():
        row = dict(catalog[pool_index])
        for name in targets:
            values = [item["targets"][name] for item in group if item["targets"][name] is not None]
            row[name] = f"{sum(values) / len(values):.12g}" if values else ""
        experiment_rows.append(row)
    mapping_rows = [{
        "experiment_id": recipe["experiment_id"],
        "recipe_row": recipe["row"],
        "sample_id": catalog[recipe["pool_index"]]["sample_id"],
        "chem_group_id": catalog[recipe["pool_index"]]["chem_group_id"],
        "pool_status": "matched_existing" if recipe["pool_index"] < original_count else "added",
        "within_design_bounds": "" if recipe["within"] is None else str(recipe["within"]).lower(),
        "replicates": len(groups[recipe["pool_index"]]),
    } for recipe in recipes]
    write_csv(design_dir / "experiment.csv", catalog_fields + list(targets), experiment_rows)
    write_csv(design_dir / "experiment_mapping.csv", MAPPING_FIELDS, mapping_rows)
    return {
        "recipes": len(recipes),
        "samples": len(groups),
        "matched_existing": sum(index < original_count for index in groups),
        "added": sum(index >= original_count for index in groups),
        "out_of_bounds": sum(recipe["within"] is False for recipe in recipes),
        "replicates_merged": len(recipes) - len(groups),
        "complete": sum(all(row[name] for name in targets) for row in experiment_rows),
        "pool_size": len(catalog),
        "skipped_examples": skipped["examples"],
        "skipped_unmeasured": skipped["unmeasured"],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--design-dir", type=Path, required=True)
    parser.add_argument("--recipes", type=Path, help="recipe CSV to import")
    parser.add_argument("--targets", nargs="+", required=True)
    parser.add_argument("--template", type=Path, help="write a recipe template and exit")
    parser.add_argument("--template-source", choices=["examples", "pool"], default="examples",
                        help="prefill the template with EXAMPLE rows or with every candidate")
    args = parser.parse_args(argv)
    if args.template:
        write_csv(args.template, *template(args.design_dir, args.targets, args.template_source))
        print(args.template)
        return
    if not args.recipes:
        parser.error("--recipes is required unless --template is given")
    summary = import_recipes(args.design_dir, args.recipes.read_text(encoding="utf-8-sig"),
                             args.targets)
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
