"""Generate a measurable salt/solvent/additive candidate pool for convert_pool.py.

Salt and additive bounds are fractions of the whole electrolyte; solvent bounds are
fractions of the solvent pool, i.e. of the mass left after salt and additives.

Cosolvents (other than the balance solvent) and additives whose lower bound is zero are
optional: each is left out of a candidate with probability ``sampling.absence_probability``
(default 0.25), so the pool also contains the sub-formulations without them. A component
that is present weighs at least its ``minimum_nonzero_mass_g`` (default: one weighing
step). ``absence_probability`` 0 reproduces the earlier fully-populated sampling.
"""
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

ADDITIVE_ROLES = ("functional_additive", "salt_additive")
DEFAULT_ABSENCE_PROBABILITY = 0.25


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
    additives = config.get("additives") or []
    if not isinstance(additives, list) or len(additives) > 8:
        raise ValueError("Additives must be a list of at most eight components")
    if any(not isinstance(item, dict) or item.get("role", "functional_additive") not in ADDITIVE_ROLES
           for item in additives):
        raise ValueError(f"Additive role must be one of: {', '.join(ADDITIVE_ROLES)}")
    components = [salt, *solvents, *additives]
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
    solvent_names = names[1:1 + len(solvents)]
    additive_names = names[1 + len(solvents):]
    bounds = [_bounds(item.get("pool_mass_fraction_bounds"), f"{name} fraction")
              for name, item in zip(solvent_names, solvents)]
    if sum(x[0] for x in bounds) > 1 + 1e-12 or sum(x[1] for x in bounds) < 1 - 1e-12:
        raise ValueError("Solvent bounds cannot sum to one")
    additive_bounds = [_bounds(item.get("final_mass_fraction_bounds"), f"{name} fraction")
                       for name, item in zip(additive_names, additives)]
    if salt_bounds[1] + sum(x[1] for x in additive_bounds) >= 1:
        raise ValueError("Salt and additive upper bounds must leave mass for the solvents")
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
    absence = float(sampling.get("absence_probability", DEFAULT_ABSENCE_PROBABILITY))
    if not (math.isfinite(absence) and 0 <= absence < 1):
        raise ValueError("absence_probability must satisfy 0 ≤ p < 1")
    minimum = []
    for name, item in zip(names, components):
        value = float(item.get("minimum_nonzero_mass_g", step))
        if not (math.isfinite(value) and value > 0):
            raise ValueError(f"minimum_nonzero_mass_g must be positive for {name}")
        minimum.append(max(value, step) / batch)
    temperature = float(config.get("temperature_C", 25))
    if not math.isfinite(temperature):
        raise ValueError("Temperature must be finite")
    ordered = [{"name": names[0], "role": "primary_salt",
                "smiles": Chem.MolToSmiles(mols[0]), "molecular_weight": Descriptors.MolWt(mols[0])}]
    for item, name, mol in zip(solvents, solvent_names, mols[1:]):
        ordered.append({"name": name, "role": "balance_solvent" if item.get("balance") else "cosolvent",
                        "smiles": Chem.MolToSmiles(mol), "molecular_weight": Descriptors.MolWt(mol)})
    for item, name, mol in zip(additives, additive_names, mols[1 + len(solvents):]):
        ordered.append({"name": name, "role": item.get("role", "functional_additive"),
                        "smiles": Chem.MolToSmiles(mol), "molecular_weight": Descriptors.MolWt(mol)})
    return (ordered, salt_bounds, bounds, additive_bounds, batch, step, units, power, seed,
            absence, np.array(minimum))


def presence_masks(samples, column, bounds, additive_bounds, balance, absence, minimum):
    """Decide per candidate which optional solvents/additives are present.

    Uses one extra Sobol column per optional component, starting at ``column``.
    Returns (solvent_present, additive_present, solvent_optional, additive_optional).
    """
    n, n_solvents, n_add = len(samples), len(bounds), len(additive_bounds)
    solvent_present = np.ones((n, n_solvents), dtype=bool)
    additive_present = np.ones((n, n_add), dtype=bool)
    if absence == 0:
        return solvent_present, additive_present, [], []
    solvent_optional = [i for i, (low, _) in enumerate(bounds) if i != balance and low == 0]
    additive_optional = [k for k, (low, _) in enumerate(additive_bounds) if low == 0]
    for i in solvent_optional:
        solvent_present[:, i] = samples[:, column] >= absence
        column += 1
    for k in additive_optional:
        additive_present[:, k] = ((samples[:, column] >= absence)
                                  & (additive_bounds[k][1] >= minimum[1 + n_solvents + k]))
        column += 1
    highs = np.array([high for _, high in bounds])
    for i in solvent_optional:
        short = (solvent_present * highs).sum(axis=1) < 1 - 1e-12
        solvent_present[short, i] = True
    return solvent_present, additive_present, solvent_optional, additive_optional


def generate(config, output_dir):
    (components, salt_bounds, bounds, additive_bounds,
     batch, step, units, power, seed, absence, minimum) = validate(config)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    n_solvents, n_add = len(bounds), len(additive_bounds)
    balance = next(i for i, item in enumerate(components[1:1 + n_solvents])
                   if item["role"] == "balance_solvent")
    n_optional = 0 if absence == 0 else (
        sum(i != balance and low == 0 for i, (low, _) in enumerate(bounds))
        + sum(low == 0 for low, _ in additive_bounds))
    samples = qmc.Sobol(d=n_solvents + n_add + n_optional, scramble=True,
                        seed=seed).random_base2(power)
    n = len(samples)
    solvent_present, additive_present, solvent_optional, additive_optional = presence_masks(
        samples, n_solvents + n_add, bounds, additive_bounds, balance, absence, minimum)
    salt_ratio = salt_bounds[0] + samples[:, 0] * (salt_bounds[1] - salt_bounds[0])
    additive_ratio = np.zeros((n, n_add))
    for k, (low, high) in enumerate(additive_bounds):
        floor = max(low, minimum[1 + n_solvents + k]) if k in additive_optional else low
        value = floor + samples[:, n_solvents + k] * (high - floor)
        additive_ratio[:, k] = np.where(additive_present[:, k], value, 0.0)
    solvent_share = 1 - salt_ratio - additive_ratio.sum(axis=1)
    lows = np.tile([low for low, _ in bounds], (n, 1))
    for i in solvent_optional:
        lows[:, i] = np.maximum(lows[:, i], minimum[1 + i] / solvent_share)
    lows *= solvent_present
    highs = np.array([high for _, high in bounds]) * solvent_present
    last = n_solvents - 1 - np.argmax(solvent_present[:, ::-1], axis=1)
    solvent_ratio = np.zeros((n, n_solvents))
    remaining = np.ones(n)
    for i in range(n_solvents):
        lower = np.maximum(lows[:, i], remaining - highs[:, i + 1:].sum(axis=1))
        upper = np.minimum(highs[:, i], remaining - lows[:, i + 1:].sum(axis=1))
        value = (lower + samples[:, i + 1] * np.maximum(0, upper - lower)
                 if i < n_solvents - 1 else remaining)
        value = np.where(last == i, remaining, value)
        solvent_ratio[:, i] = np.where(solvent_present[:, i], value, 0.0)
        remaining = remaining - solvent_ratio[:, i]
    fractions = np.column_stack([salt_ratio, solvent_share[:, None] * solvent_ratio, additive_ratio])
    raw_units = fractions * units
    integral = np.floor(raw_units).astype(int)
    remainder = units - integral.sum(axis=1)
    order = np.argsort(-np.where(fractions > 0, raw_units - integral, -1.0), axis=1, kind="stable")
    optional = np.zeros(len(components), dtype=bool)
    optional[[1 + i for i in solvent_optional] + [1 + n_solvents + k for k in additive_optional]] = True
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
        if np.any(optional & (amount > 0) & (mass < minimum - 1e-12)):
            continue
        if any(not (low - 1e-12 <= value <= high + 1e-12)
               for value, (low, high) in zip(mass[1 + n_solvents:], additive_bounds)):
            continue
        solvent_total = mass[1:1 + n_solvents].sum()
        if solvent_total <= 0:
            continue
        actual_solvents = mass[1:1 + n_solvents] / solvent_total
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
    return {"generated": n, "feasible": len(rows), "components": len(components),
            "presence_patterns": len({row["presence_pattern"] for row in rows})}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(generate(json.loads(args.config.read_text(encoding="utf-8")), args.output_dir)))


if __name__ == "__main__":
    main()
