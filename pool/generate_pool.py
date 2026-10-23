"""Generate a measurable salt/solvent/additive candidate pool for convert_pool.py.

A design accepts a library of lithium salts (``salts``; a single ``salt`` object is also
accepted). ``salt_selection.sizes`` enumerates fixed single/dual-salt spaces before
Sobol sampling, sharing a total budget of 2**candidate_power across those spaces.
Unselected salts are zero; selected salts obey their conditional bounds and minimum mass.
Without salt_selection the original one/two-salt presence sampling remains unchanged.
Salt and additive bounds are fractions of the whole electrolyte; solvent bounds
are fractions of the solvent pool, i.e. of the mass left after salts and additives.

In legacy mode, cosolvents (other than the balance solvent), additives and salts
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
from itertools import combinations
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
MAX_SALTS = 12


def _bounds(value, label):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{label} must have two bounds")
    low, high = map(float, value)
    if not (math.isfinite(low) and math.isfinite(high) and 0 <= low <= high <= 1):
        raise ValueError(f"{label} must satisfy 0 ≤ lower ≤ upper ≤ 1")
    return low, high


def salt_list(config):
    """The design's salt library, preserving its feature-column order."""
    if "salts" in config and "salt" in config:
        raise ValueError("Give either salts or salt, not both")
    salts = config["salts"] if "salts" in config else [config.get("salt") or {}]
    if (not isinstance(salts, list) or not 1 <= len(salts) <= MAX_SALTS
            or not all(isinstance(item, dict) for item in salts)):
        raise ValueError(f"Give one to {MAX_SALTS} lithium salts")
    return salts


def salt_spaces(config):
    """Fixed salt-index tuples, or None for the legacy presence sampler.

    In this mode a positive salt lower bound applies only when that salt is selected.
    """
    selection = config.get("salt_selection")
    n_salts = len(salt_list(config))
    if selection is None:
        if n_salts > 2:
            raise ValueError("More than two library salts require salt_selection.sizes")
        return None
    sizes = selection.get("sizes") if isinstance(selection, dict) else None
    if (not isinstance(sizes, list) or not sizes
            or any(type(size) is not int or size not in (1, 2) for size in sizes)
            or len(set(sizes)) != len(sizes)):
        raise ValueError("salt_selection.sizes must be [1], [2], or [1, 2]")
    spaces = [space for size in sorted(sizes) for space in combinations(range(n_salts), size)]
    if not spaces:
        raise ValueError("The chosen salt spaces require at least two library salts")
    return spaces


def validate(config):
    salts = salt_list(config)
    spaces = salt_spaces(config)
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
    allowed = spaces if spaces is not None else [tuple(range(n_salts))]
    largest_salt_share = max(sum(salt_bounds[j][1] for j in space) for space in allowed)
    if largest_salt_share + sum(x[1] for x in additive_bounds) >= 1:
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
    if spaces and 2 ** power < len(spaces):
        raise ValueError(f"采样预算不足：{len(spaces)} 个盐组合至少各需一个点，请提高采样指数")
    if type(sampling.get("require_all_solvents", False)) is not bool:
        raise ValueError("sampling.require_all_solvents must be boolean")
    absence = float(sampling.get("absence_probability", DEFAULT_ABSENCE_PROBABILITY))
    if not (math.isfinite(absence) and 0 <= absence < 1):
        raise ValueError("absence_probability must satisfy 0 ≤ p < 1")
    minimum = []
    for name, item in zip(names, components):
        value = float(item.get("minimum_nonzero_mass_g", step))
        if not (math.isfinite(value) and value > 0):
            raise ValueError(f"minimum_nonzero_mass_g must be positive for {name}")
        minimum.append(max(value, step) / batch)
    for j in (range(n_salts) if spaces is not None else optional_salts(salt_bounds, absence)):
        # Check on the actual weighing lattice, not only in continuous fractions.
        floor = max(salt_bounds[j][0], minimum[j])
        if spaces is not None:
            if math.ceil(floor * units - 1e-9) > math.floor(salt_bounds[j][1] * units + 1e-9):
                raise ValueError(f"{salt_names[j]} has no weighable nonzero mass within its bounds")
        elif salt_bounds[j][1] < minimum[j]:
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

    Molality grows with every salt and additive fraction. Bound it separately in each
    allowed salt space, then take the overall span. This is a continuous-bound preview;
    minimum weighing amounts and solvent feasibility can make the actual span narrower.
    """
    salts = salt_list(config)
    spaces = salt_spaces(config)
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

    def molality(side, space):
        salt_fractions = [x[side] if j in space else 0 for j, x in enumerate(salt_bounds)]
        solvent = 1 - sum(salt_fractions) - sum(x[side] for x in additives)
        if solvent <= 0:
            raise ValueError("Salt and additive upper bounds must leave mass for the solvents")
        return 1000 * sum(x / weight for x, weight in zip(salt_fractions, weights)) / solvent

    allowed = spaces if spaces is not None else [tuple(range(len(salts)))]
    return {"molecular_weights": weights,
            "low": min(molality(0, space) for space in allowed),
            "high": max(molality(1, space) for space in allowed)}


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
                   absence, minimum, fixed_salts=None):
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
    if fixed_salts is not None:
        present[:, :n_salts] = fixed_salts
    salts = optional_salts(salt_bounds, absence) if fixed_salts is None else []
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


def stratified_sobol(spaces, n_salts, presence_column, n_optional, power, seed):
    """Independent scrambled Sobol sequences inside preselected salt spaces.

    Divide the TOTAL budget evenly (difference at most one), and embed local draws
    into the full column layout. Inactive salts consume no local Sobol dimensions.
    Non-power-of-two quotas take a prefix of the next complete Sobol block.
    """
    total = 2 ** power
    samples = np.zeros((total, presence_column + n_optional))
    masks = np.zeros((total, n_salts), dtype=bool)
    space_ids = np.empty(total, dtype=int)
    counts = []
    start = 0
    seeds = np.random.SeedSequence(seed).spawn(len(spaces))
    for index, (space, child_seed) in enumerate(zip(spaces, seeds)):
        count = total // len(spaces) + (index < total % len(spaces))
        columns = [*space, *range(n_salts, presence_column + n_optional)]
        draws = qmc.Sobol(d=len(columns), scramble=True, seed=np.random.default_rng(child_seed))
        block = draws.random_base2((count - 1).bit_length())[:count]
        stop = start + count
        samples[np.ix_(range(start, stop), columns)] = block
        masks[np.ix_(range(start, stop), space)] = True
        space_ids[start:stop] = index
        counts.append({"salt_count": len(space), "generated": count, "feasible": 0})
        start = stop
    return samples, masks, space_ids, counts


def generate(config, output_dir):
    (components, salt_bounds, bounds, additive_bounds,
     batch, step, units, power, seed, absence, minimum) = validate(config)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    n_salts, n_solvents, n_add = len(salt_bounds), len(bounds), len(additive_bounds)
    spaces = salt_spaces(config)
    require_solvents = config.get("sampling", {}).get("require_all_solvents", False)
    solvent_end = n_salts + n_solvents
    balance = next(i for i, item in enumerate(components[n_salts:solvent_end])
                   if item["role"] == "balance_solvent")
    optional = optional_components(salt_bounds, bounds, additive_bounds, balance, absence)
    if spaces is not None:
        optional = [j for j in optional if j >= n_salts]
    if require_solvents:
        optional = [j for j in optional if not n_salts <= j < solvent_end]
    # Sobol columns: salts, all solvents but the last, additives, one per optional component.
    presence_column = solvent_end - 1 + n_add
    fixed_salts = space_ids = space_counts = None
    if spaces is None:
        samples = qmc.Sobol(d=presence_column + len(optional), scramble=True,
                            seed=seed).random_base2(power)
    else:
        samples, fixed_salts, space_ids, space_counts = stratified_sobol(
            spaces, n_salts, presence_column, len(optional), power, seed)
    n = len(samples)
    present = presence_masks(samples, presence_column, optional, salt_bounds, bounds,
                             additive_bounds, absence, minimum, fixed_salts)
    is_optional = np.zeros(len(components), dtype=bool)
    is_optional[optional] = True
    enforce_minimum = is_optional.copy()
    if spaces is not None:
        enforce_minimum[:n_salts] = True
    if require_solvents:
        enforce_minimum[n_salts:solvent_end] = True

    def direct_ratio(component, column, low, high):
        """Whole-electrolyte fraction of a salt or additive; zero where it is absent."""
        floor = max(low, minimum[component]) if enforce_minimum[component] else low
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
        if enforce_minimum[n_salts + i]:
            lows[:, i] = np.maximum(lows[:, i], minimum[n_salts + i] / solvent_share)
    lows *= solvent_present
    highs = np.array([high for _, high in bounds]) * solvent_present
    last = n_solvents - 1 - np.argmax(solvent_present[:, ::-1], axis=1)
    solvent_ratio = np.zeros((n, n_solvents))
    remaining = np.ones(n)
    feasible_draw = np.ones(n, dtype=bool)
    for i in range(n_solvents):
        lower = np.maximum(lows[:, i], remaining - highs[:, i + 1:].sum(axis=1))
        upper = np.minimum(highs[:, i], remaining - lows[:, i + 1:].sum(axis=1))
        feasible_draw &= lower <= upper + 1e-12
        value = (lower + samples[:, n_salts + i] * np.maximum(0, upper - lower)
                 if i < n_solvents - 1 else remaining)
        value = np.where(last == i, remaining, value)
        solvent_ratio[:, i] = np.where(solvent_present[:, i], value, 0.0)
        remaining = remaining - solvent_ratio[:, i]
    fractions = np.column_stack([salt_ratio, solvent_share[:, None] * solvent_ratio, additive_ratio])
    fractions = fractions[feasible_draw]
    if fixed_salts is not None:
        fixed_salts = fixed_salts[feasible_draw]
        space_ids = space_ids[feasible_draw]
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
    for row_index, amount in enumerate(integral):
        signature = tuple(int(x) for x in amount)
        if signature in seen:
            continue
        mass = amount / units
        if fixed_salts is not None and not np.array_equal(amount[:n_salts] > 0, fixed_salts[row_index]):
            continue
        if require_solvents and not (amount[n_salts:solvent_end] > 0).all():
            continue
        if any(not (low - 1e-12 <= value <= high + 1e-12)
               and not (spaces is not None and value == 0)
               for value, (low, high) in zip(mass[:n_salts], salt_bounds)):
            continue
        if np.any(enforce_minimum & (amount > 0) & (mass < minimum - 1e-12)):
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
        if space_counts is not None:
            space_counts[space_ids[row_index]]["feasible"] += 1
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
    stats = {"generated": n, "feasible": len(rows), "components": len(components),
             "presence_patterns": len({row["presence_pattern"] for row in rows})}
    if space_counts is not None:
        for space, count in zip(spaces, space_counts):
            count["salts"] = [names[j] for j in space]
        stats["salt_spaces"] = space_counts
        stats["single_salt_candidates"] = sum(x["feasible"] for x in space_counts if x["salt_count"] == 1)
        stats["dual_salt_candidates"] = sum(x["feasible"] for x in space_counts if x["salt_count"] == 2)
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(generate(json.loads(args.config.read_text(encoding="utf-8")), args.output_dir)))


if __name__ == "__main__":
    main()
