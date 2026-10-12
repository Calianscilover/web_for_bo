"""Generate a measurable salt/solvent/additive candidate pool for convert_pool.py.

A design has one or two lithium salts (``salts``; a single ``salt`` object is also
accepted). Salt and additive bounds are fractions of the whole electrolyte; solvent bounds
are fractions of the solvent pool, i.e. of the mass left after salts and additives.

Cosolvents (other than the balance solvent), additives and, in two-salt designs, salts
whose lower bound is zero are optional: each is left out of a candidate with probability
``sampling.absence_probability`` (default 0.25), so the pool also contains the
sub-formulations without them. Both salts are treated alike; when every optional salt is
drawn absent, the one whose presence draw was highest is kept, so each candidate has a
salt. A component that is present weighs at least its ``minimum_nonzero_mass_g``
(default: one weighing step). ``absence_probability`` 0 reproduces the earlier
fully-populated sampling.

Each row also records ``salt_molality_mol_kg``: moles of all lithium salts per kg of
solvents (balance solvent + cosolvents), computed from the weighed masses. It is reported
only; sampling stays in mass fractions.
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

try:
    from pool.convert_pool import MOLALITY_FIELD, format_molality, salt_molality
except ImportError:
    from convert_pool import MOLALITY_FIELD, format_molality, salt_molality

ADDITIVE_ROLES = ("functional_additive", "salt_additive")
DEFAULT_ABSENCE_PROBABILITY = 0.25
MAX_SALTS = 2


def _bounds(value, label):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{label} must have two bounds")
    low, high = map(float, value)
    if not (math.isfinite(low) and math.isfinite(high) and 0 <= low <= high <= 1):
        raise ValueError(f"{label} must satisfy 0 ≤ lower ≤ upper ≤ 1")
    return low, high


def salt_list(config):
    """The design's lithium salts: ``salts`` (one or two) or a single ``salt`` object."""
    if "salts" in config and "salt" in config:
        raise ValueError("Give either salts or salt, not both")
    salts = config["salts"] if "salts" in config else [config.get("salt") or {}]
    if (not isinstance(salts, list) or not 1 <= len(salts) <= MAX_SALTS
            or not all(isinstance(item, dict) for item in salts)):
        raise ValueError(f"Give one to {MAX_SALTS} lithium salts")
    return salts


def validate(config):
    salts = salt_list(config)
    n_salts = len(salts)
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
    components = [*salts, *solvents, *additives]
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
    salt_names = names[:n_salts]
    solvent_names = names[n_salts:n_salts + len(solvents)]
    additive_names = names[n_salts + len(solvents):]
    salt_bounds = [_bounds(item.get("final_mass_fraction_bounds"), f"{name} fraction")
                   for name, item in zip(salt_names, salts)]
    bounds = [_bounds(item.get("pool_mass_fraction_bounds"), f"{name} fraction")
              for name, item in zip(solvent_names, solvents)]
    if sum(x[0] for x in bounds) > 1 + 1e-12 or sum(x[1] for x in bounds) < 1 - 1e-12:
        raise ValueError("Solvent bounds cannot sum to one")
    additive_bounds = [_bounds(item.get("final_mass_fraction_bounds"), f"{name} fraction")
                       for name, item in zip(additive_names, additives)]
    if sum(x[1] for x in salt_bounds) + sum(x[1] for x in additive_bounds) >= 1:
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
    for j in optional_salts(salt_bounds, absence):
        if salt_bounds[j][1] < minimum[j]:
            raise ValueError(f"{salt_names[j]} upper bound is below its minimum nonzero mass")
    temperature = float(config.get("temperature_C", 25))
    if not math.isfinite(temperature):
        raise ValueError("Temperature must be finite")
    ordered = [{"name": name, "role": "primary_salt",
                "smiles": Chem.MolToSmiles(mol), "molecular_weight": Descriptors.MolWt(mol)}
               for name, mol in zip(salt_names, mols)]
    for item, name, mol in zip(solvents, solvent_names, mols[n_salts:]):
        ordered.append({"name": name, "role": "balance_solvent" if item.get("balance") else "cosolvent",
                        "smiles": Chem.MolToSmiles(mol), "molecular_weight": Descriptors.MolWt(mol)})
    for item, name, mol in zip(additives, additive_names, mols[n_salts + len(solvents):]):
        ordered.append({"name": name, "role": item.get("role", "functional_additive"),
                        "smiles": Chem.MolToSmiles(mol), "molecular_weight": Descriptors.MolWt(mol)})
    return (ordered, salt_bounds, bounds, additive_bounds, batch, step, units, power, seed,
            absence, np.array(minimum))


def salt_molality_range(config):
    """Total-salt molality span (mol/kg of solvents) implied by the mass-fraction bounds.

    Molality grows with every salt and additive fraction, so the lowest value uses all
    lower bounds and the highest all upper bounds; every sampled candidate falls inside.
    """
    salts = salt_list(config)
    weights = []
    for item in salts:
        molecule = Chem.MolFromSmiles(str(item.get("smiles", "")))
        if molecule is None:
            raise ValueError("Invalid SMILES for the lithium salt")
        weights.append(Descriptors.MolWt(molecule))
    salt_bounds = [_bounds(item.get("final_mass_fraction_bounds"), "Salt fraction")
                   for item in salts]
    additives = [_bounds(item.get("final_mass_fraction_bounds"), "Additive fraction")
                 for item in config.get("additives") or []]

    def molality(side):
        salt_fractions = [x[side] for x in salt_bounds]
        solvent = 1 - sum(salt_fractions) - sum(x[side] for x in additives)
        if solvent <= 0:
            raise ValueError("Salt and additive upper bounds must leave mass for the solvents")
        return 1000 * sum(x / weight for x, weight in zip(salt_fractions, weights)) / solvent

    return {"molecular_weights": weights, "low": molality(0), "high": molality(1)}


def optional_salts(salt_bounds, absence):
    """Salts a candidate may leave out: zero lower bound, only in multi-salt designs."""
    if absence == 0 or len(salt_bounds) < 2:
        return []
    return [j for j, (low, _) in enumerate(salt_bounds) if low == 0]


def optional_components(salt_bounds, bounds, additive_bounds, balance, absence):
    """Indices of optional components in component order (salts, solvents, additives)."""
    if absence == 0:
        return []
    n_salts, n_solvents = len(salt_bounds), len(bounds)
    return (optional_salts(salt_bounds, absence)
            + [n_salts + i for i, (low, _) in enumerate(bounds) if i != balance and low == 0]
            + [n_salts + n_solvents + k for k, (low, _) in enumerate(additive_bounds) if low == 0])


def presence_masks(samples, column, optional, salt_bounds, bounds, additive_bounds,
                   absence, minimum):
    """Decide per candidate which components are present.

    Each component in ``optional`` uses one Sobol column, in order, starting at ``column``.
    Returns a (candidates x components) boolean array in component order.
    """
    n_salts, n_solvents = len(salt_bounds), len(bounds)
    solvents = slice(n_salts, n_salts + n_solvents)
    highs = np.array([high for _, high in (*salt_bounds, *bounds, *additive_bounds)])
    # Solvent minimums depend on the solvent share, so solvents are filtered after rounding.
    room = highs >= minimum
    room[solvents] = True
    present = np.ones((len(samples), len(highs)), dtype=bool)
    for offset, c in enumerate(optional):
        present[:, c] = (samples[:, column + offset] >= absence) & room[c]
    salts = optional_salts(salt_bounds, absence)
    if salts:
        draws = samples[:, column:column + len(salts)]
        saltless = ~present[:, :n_salts].any(axis=1)
        present[saltless, np.asarray(salts)[draws[saltless].argmax(axis=1)]] = True
    solvent_highs = highs[solvents]
    for c in optional:
        if n_salts <= c < n_salts + n_solvents:
            short = (present[:, solvents] * solvent_highs).sum(axis=1) < 1 - 1e-12
            present[short, c] = True
    return present


def generate(config, output_dir):
    (components, salt_bounds, bounds, additive_bounds,
     batch, step, units, power, seed, absence, minimum) = validate(config)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    n_salts, n_solvents, n_add = len(salt_bounds), len(bounds), len(additive_bounds)
    solvent_end = n_salts + n_solvents
    balance = next(i for i, item in enumerate(components[n_salts:solvent_end])
                   if item["role"] == "balance_solvent")
    optional = optional_components(salt_bounds, bounds, additive_bounds, balance, absence)
    # Sobol columns: salts, all solvents but the last, additives, one per optional component.
    presence_column = solvent_end - 1 + n_add
    samples = qmc.Sobol(d=presence_column + len(optional), scramble=True,
                        seed=seed).random_base2(power)
    n = len(samples)
    present = presence_masks(samples, presence_column, optional, salt_bounds, bounds,
                             additive_bounds, absence, minimum)
    is_optional = np.zeros(len(components), dtype=bool)
    is_optional[optional] = True

    def direct_ratio(component, column, low, high):
        """Whole-electrolyte fraction of a salt or additive; zero where it is absent."""
        floor = max(low, minimum[component]) if is_optional[component] else low
        return np.where(present[:, component], floor + samples[:, column] * (high - floor), 0.0)

    salt_ratio = np.zeros((n, n_salts))
    for j, (low, high) in enumerate(salt_bounds):
        salt_ratio[:, j] = direct_ratio(j, j, low, high)
    additive_ratio = np.zeros((n, n_add))
    for k, (low, high) in enumerate(additive_bounds):
        additive_ratio[:, k] = direct_ratio(solvent_end + k, solvent_end - 1 + k, low, high)
    solvent_share = 1 - salt_ratio.sum(axis=1) - additive_ratio.sum(axis=1)
    solvent_present = present[:, n_salts:solvent_end]
    lows = np.tile([low for low, _ in bounds], (n, 1))
    for i in range(n_solvents):
        if is_optional[n_salts + i]:
            lows[:, i] = np.maximum(lows[:, i], minimum[n_salts + i] / solvent_share)
    lows *= solvent_present
    highs = np.array([high for _, high in bounds]) * solvent_present
    last = n_solvents - 1 - np.argmax(solvent_present[:, ::-1], axis=1)
    solvent_ratio = np.zeros((n, n_solvents))
    remaining = np.ones(n)
    for i in range(n_solvents):
        lower = np.maximum(lows[:, i], remaining - highs[:, i + 1:].sum(axis=1))
        upper = np.minimum(highs[:, i], remaining - lows[:, i + 1:].sum(axis=1))
        value = (lower + samples[:, n_salts + i] * np.maximum(0, upper - lower)
                 if i < n_solvents - 1 else remaining)
        value = np.where(last == i, remaining, value)
        solvent_ratio[:, i] = np.where(solvent_present[:, i], value, 0.0)
        remaining = remaining - solvent_ratio[:, i]
    fractions = np.column_stack([salt_ratio, solvent_share[:, None] * solvent_ratio, additive_ratio])
    raw_units = fractions * units
    integral = np.floor(raw_units).astype(int)
    remainder = units - integral.sum(axis=1)
    order = np.argsort(-np.where(fractions > 0, raw_units - integral, -1.0), axis=1, kind="stable")
    for row_index, count in enumerate(remainder):
        integral[row_index, order[row_index, :count]] += 1
    names = [item["name"] for item in components]
    weights = np.array([item["molecular_weight"] for item in components])
    fields = ["presence_pattern", MOLALITY_FIELD]
    for name in names:
        fields.extend([f"{name}_mass_g", f"{name}_mass_fraction", f"{name}_mole_fraction"])
    rows = []
    seen = set()
    for amount in integral:
        signature = tuple(int(x) for x in amount)
        if signature in seen:
            continue
        mass = amount / units
        if any(not (low - 1e-12 <= value <= high + 1e-12)
               for value, (low, high) in zip(mass[:n_salts], salt_bounds)):
            continue
        if np.any(is_optional & (amount > 0) & (mass < minimum - 1e-12)):
            continue
        if any(not (low - 1e-12 <= value <= high + 1e-12)
               for value, (low, high) in zip(mass[solvent_end:], additive_bounds)):
            continue
        solvent_total = mass[n_salts:solvent_end].sum()
        if solvent_total <= 0:
            continue
        actual_solvents = mass[n_salts:solvent_end] / solvent_total
        if any(not (low - 1e-12 <= value <= high + 1e-12)
               for value, (low, high) in zip(actual_solvents, bounds)):
            continue
        moles = mass / weights
        mole_fraction = moles / moles.sum()
        row = {"presence_pattern": ";".join(name for name, x in zip(names, amount) if x > 0),
               MOLALITY_FIELD: format_molality(salt_molality(amount * step, components))}
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
