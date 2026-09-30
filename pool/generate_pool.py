"""Generate a measurable salt/solvent candidate pool for convert_pool.py."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from scipy.stats import qmc


def _bounds(value, label):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{label} must have two bounds")
    low, high = map(float, value)
    if not (math.isfinite(low) and math.isfinite(high) and 0 <= low <= high <= 1):
        raise ValueError(f"{label} must satisfy 0 ≤ lower ≤ upper ≤ 1")
    return low, high


def validate(config):
    salt = config.get("salt", {})
    solvents = config.get("solvents", [])
    if not isinstance(solvents, list) or len(solvents) < 2:
        raise ValueError("At least two solvents are required")
    if sum(bool(item.get("balance")) for item in solvents) != 1:
        raise ValueError("Mark exactly one balance solvent")
    components = [salt, *solvents]
    names = [str(item.get("name", "")).strip() for item in components]
    if any(not name or not name.replace("_", "").isalnum() for name in names):
        raise ValueError("Component names must contain only letters, digits or underscores")
    if len(set(names)) != len(names):
        raise ValueError("Component names must be unique")
    mols = []
    for name, item in zip(names, components):
        mol = Chem.MolFromSmiles(str(item.get("smiles", "")))
        if mol is None:
            raise ValueError(f"Invalid SMILES for {name}")
        mols.append(mol)
    salt_bounds = _bounds(salt.get("final_mass_fraction_bounds"), "Salt fraction")
    bounds = [_bounds(item.get("pool_mass_fraction_bounds"), f"{name} fraction")
              for name, item in zip(names[1:], solvents)]
    if sum(x[0] for x in bounds) > 1 + 1e-12 or sum(x[1] for x in bounds) < 1 - 1e-12:
        raise ValueError("Solvent bounds cannot sum to one")
    batch = float(config.get("batch_mass_g", 5))
    step = float(config.get("mass_step_g", 0.0001))
    if not (math.isfinite(batch) and math.isfinite(step) and batch > 0 and step > 0):
        raise ValueError("Batch mass and weighing step must be positive")
    units = round(batch / step)
    if units < 100 or abs(units * step - batch) > 1e-8:
        raise ValueError("Batch mass must be an exact multiple of the weighing step (at least 100 steps)")
    sampling = config.get("sampling", {})
    power = int(sampling.get("candidate_power", 10))
    seed = int(sampling.get("seed", 2026))
    if not 3 <= power <= 14:
        raise ValueError("candidate_power must be between 3 and 14")
    temperature = float(config.get("temperature_C", 25))
    if not math.isfinite(temperature):
        raise ValueError("Temperature must be finite")
    ordered = [{"name": names[0], "role": "primary_salt",
                "smiles": Chem.MolToSmiles(mols[0]), "molecular_weight": Descriptors.MolWt(mols[0])}]
    for item, name, mol in zip(solvents, names[1:], mols[1:]):
        ordered.append({"name": name, "role": "balance_solvent" if item.get("balance") else "cosolvent",
                        "smiles": Chem.MolToSmiles(mol), "molecular_weight": Descriptors.MolWt(mol)})
    return ordered, salt_bounds, bounds, batch, step, units, power, seed


def generate(config, output_dir):
    components, salt_bounds, bounds, batch, step, units, power, seed = validate(config)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = qmc.Sobol(d=len(bounds), scramble=True, seed=seed).random_base2(power)
    salt_ratio = salt_bounds[0] + samples[:, 0] * (salt_bounds[1] - salt_bounds[0])
    solvent_ratio = np.zeros((len(samples), len(bounds)))
    remaining = np.ones(len(samples))
    for i in range(len(bounds) - 1):
        low = max(0.0, bounds[i][0])
        lower = np.maximum(low, remaining - sum(x[1] for x in bounds[i + 1:]))
        upper = np.minimum(bounds[i][1], remaining - sum(x[0] for x in bounds[i + 1:]))
        solvent_ratio[:, i] = lower + samples[:, i + 1] * np.maximum(0, upper - lower)
        remaining -= solvent_ratio[:, i]
    solvent_ratio[:, -1] = remaining
    fractions = np.column_stack([salt_ratio, (1 - salt_ratio)[:, None] * solvent_ratio])
    raw_units = fractions * units
    integral = np.floor(raw_units).astype(int)
    remainder = units - integral.sum(axis=1)
    order = np.argsort(-(raw_units - integral), axis=1, kind="stable")
    for row_index, count in enumerate(remainder):
        integral[row_index, order[row_index, :count]] += 1
    names = [item["name"] for item in components]
    weights = np.array([item["molecular_weight"] for item in components])
    fields = ["presence_pattern"]
    for name in names:
        fields.extend([f"{name}_mass_g", f"{name}_mass_fraction", f"{name}_mole_fraction"])
    rows = []
    seen = set()
    for amount in integral:
        signature = tuple(int(x) for x in amount)
        if signature in seen:
            continue
        mass = amount / units
        if not (salt_bounds[0] - 1e-12 <= mass[0] <= salt_bounds[1] + 1e-12):
            continue
        solvent_total = mass[1:].sum()
        if solvent_total <= 0:
            continue
        actual_solvents = mass[1:] / solvent_total
        if any(not (low - 1e-12 <= value <= high + 1e-12)
               for value, (low, high) in zip(actual_solvents, bounds)):
            continue
        moles = mass / weights
        mole_fraction = moles / moles.sum()
        row = {"presence_pattern": ";".join(name for name, x in zip(names, amount) if x > 0)}
        for j, name in enumerate(names):
            row[f"{name}_mass_g"] = f"{amount[j] * step:.10g}"
            row[f"{name}_mass_fraction"] = f"{mass[j]:.12g}"
            row[f"{name}_mole_fraction"] = f"{mole_fraction[j]:.12g}"
        rows.append(row)
        seen.add(signature)
    if len(rows) < 2:
        raise ValueError("Fewer than two distinct feasible candidates; increase samples or weighing precision")
    snapshot = dict(config)
    snapshot["components"] = [{k: item[k] for k in ("name", "role", "smiles")}
                              for item in components]
    snapshot["design_name"] = str(config.get("name", "Untitled design"))
    (output_dir / "config_snapshot.json").write_text(
        json.dumps(snapshot, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with (output_dir / "components.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["name", "role", "smiles", "molecular_weight"])
        writer.writeheader()
        writer.writerows(components)
    for name in ("feasible_candidates.csv", "pool.csv"):
        with (output_dir / name).open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    return {"generated": len(samples), "feasible": len(rows), "components": len(components)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(generate(json.loads(args.config.read_text(encoding="utf-8")), args.output_dir)))


if __name__ == "__main__":
    main()
