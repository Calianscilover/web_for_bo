"""Two formulation validators followed by a separate target-validation gate."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors

from pool.convert_pool import load_components, transform_rows
from pool.import_experiments import (EXAMPLE_PREFIX, RecipeError, pool_row,
                                     read_csv_text, recipe_composition)


LIBRARY = [
    ("LiDFOB", "primary_salt", "F[B-]1(OC(C(O1)=O)=O)F.[Li+]"),
    ("LiPF6", "primary_salt", "[Li+].F[P-](F)(F)(F)(F)F"),
    ("LiFSI", "primary_salt", "[Li+].O=S(=O)(F)[N-]S(=O)(=O)F"),
    ("LiBF4", "primary_salt", "[Li+].F[B-](F)(F)F"),
    ("LiTFSI", "primary_salt", "[Li+].O=S(=O)(C(F)(F)F)[N-]S(=O)(=O)C(F)(F)F"),
    ("TTE", "balance_solvent", "FC(F)C(F)(F)COC(F)(F)C(F)F"),
    ("FEC", "cosolvent", "O=C1OCC(F)O1"),
    ("DME", "cosolvent", "COCCOC"),
    ("NDFA", "balance_solvent", "CN(C)C(=O)C(F)(F)F"),
    ("LiNO3", "salt_additive", "[Li+].[O-][N+](=O)[O-]"),
    ("VC", "functional_additive", "O=C1OC=CO1"),
    ("TMSP", "functional_additive", "C[Si](C)(C)OP(O[Si](C)(C)C)O[Si](C)(C)C"),
    ("EC", "balance_solvent", "O=C1OCCO1"),
    ("DMC", "cosolvent", "COC(=O)OC"),
    ("EMC", "cosolvent", "CCOC(=O)OC"),
    ("MA", "cosolvent", "O=C(OC)C"),
]
MASS_META = {"reference_batch_mass_g", "batch_mass_g", "total_mass_g"}


class DataError(ValueError):
    def __init__(self, message, details=()):
        super().__init__(message)
        self.details = list(details)[:30]


def save_json(path, value):
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def read_json(path):
    return json.loads(path.read_text())


def write_csv(path, fields, rows):
    temp = path.with_suffix(".tmp")
    with temp.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def table(path):
    return read_csv_text(path.read_text(encoding="utf-8-sig"))


def upload(directory, text, name="uploaded.csv"):
    if not isinstance(text, str) or len(text.encode("utf-8")) > 20_000_000:
        raise DataError("CSV 必须为 20 MB 以内的文本")
    try:
        fields, rows = read_csv_text(text)
    except (RecipeError, ValueError) as exc:
        raise DataError(str(exc), getattr(exc, "details", ()))
    if not rows or any(None in row or None in row.values() for row in rows):
        raise DataError("CSV 为空或存在行列不完整的数据")
    directory.mkdir(parents=True)
    (directory / "source.csv").write_text(text, encoding="utf-8")
    recipe = any(x.endswith(("_mass_g", "_mass_fraction")) and x not in MASS_META for x in fields)
    catalog = "compound_0" in fields and any(x in fields for x in ("mass_ratio_0", "mole_ratio_0", "ratio_0"))
    state = {"dataset_id": directory.name, "status": "uploaded", "name": Path(name).name,
             "rows": len(rows), "fields": fields,
             "formats": (["experiment_recipe"] if recipe else []) + (["candidate_catalog"] if catalog else []),
             "source_sha256": hashlib.sha256(text.encode()).hexdigest(),
             "inspection": inspect_table(fields, rows)}
    save_json(directory / "status.json", state)
    return state


def eligible_rows(rows):
    """Retain source indexes/Excel lines while following online template conventions."""
    active, skipped = [], {"examples": 0, "blank": 0}
    for index, row in enumerate(rows):
        if not any(str(value or "").strip() for value in row.values()):
            skipped["blank"] += 1
        elif row.get("experiment_id", "").strip().upper().startswith(EXAMPLE_PREFIX):
            skipped["examples"] += 1
        else:
            active.append((index, row))
    return active, skipped


def property_columns(fields, rows):
    def candidate(field):
        if (field in MASS_META or field in {"temperature[C]", "temperature_C", "source_row", "mass_basis"}
                or field.endswith(("_mass_g", "_mass_fraction", "_mole_fraction", "_id", "_count",
                                   "_pred_std", "_model_devi", "_std", "_pred_values"))
                or re.match(r"(?:mass_ratio|mole_ratio|ratio|compound|smiles|constraint|normalization|source_mass)_", field)
                or "molality" in field or "molarity" in field):
            return False
        values = [row[field].strip() for row in rows if row[field].strip()]
        if not values:
            return True
        for value in values:
            try:
                float(value)
                return True
            except ValueError:
                continue
        return field in ("CRR[%]", "Conductivity", "logCE")
    return sorted((f for f in fields if candidate(f)),
                  key=lambda f: (f.endswith("_pred"), f not in ("CRR[%]", "Conductivity", "logCE"), fields.index(f)))


def inspect_table(fields, rows):
    active, skipped = eligible_rows(rows)
    masses = [f for f in fields if f.endswith(("_mass_g", "_mass_fraction")) and f not in MASS_META]
    slots = [f for f in fields if re.fullmatch(r"compound_\d+", f)]
    species = sorted({row[f].strip() for _, row in active for f in slots} - {""})
    # Inspection reports file contents without consulting the chemistry registry.
    component_counts = {}
    prefix = next((p for p in ("mass_ratio", "mole_ratio", "ratio") if p+"_0" in fields), None)
    for _, row in active:
        present = set()
        for field in slots if prefix else masses:
            name = row[field].strip() if prefix else field.removesuffix("_mass_fraction").removesuffix("_mass_g")
            try:
                value = row.get(prefix+"_"+field[9:]) if prefix else row[field]
                if name and float(value or 0) > 0:
                    present.add(name)
            except ValueError:
                pass
        for name in present:
            component_counts[name] = component_counts.get(name, 0)+1
    return {"profile": "ratio_slots" if "ratio_0" in fields else "catalog" if slots else "recipe",
            "composition_columns": masses or slots, "species": species,
            "eligible_rows": len(active), "skipped": skipped,
            "ratio_columns": [x for x in ("mass_ratio", "mole_ratio", "ratio") if x+"_0" in fields],
            "optional_components": {n: count for n, count in component_counts.items() if count < len(active)},
            "property_columns": property_columns(fields, [row for _, row in active]),
            "property_counts": {f: sum(bool(row[f].strip()) for _, row in active)
                                for f in property_columns(fields, [row for _, row in active])},
            "prediction_columns": [f for f in fields if f.endswith("_pred")],
            "auxiliary_columns": [f for f in fields if f.endswith(("_std", "_model_devi", "_pred_values"))]}


def identity_resolver(known):
    """Resolve molecular aliases by structure and 1:1 lithium salts by their anion.

    Anions are never treated as independent salts; their Li+ balance is checked per row.
    """
    by_structure, anions = {}, {}
    lithium = Chem.MolToSmiles(Chem.MolFromSmiles("[Li+]"))
    for name, item in known.items():
        mol = Chem.MolFromSmiles(item.get("smiles", ""))
        if mol is None or not item.get("role"):
            raise DataError("组分配置需有有效的 role 和 SMILES", [name])
        canonical = Chem.MolToSmiles(mol)
        if canonical in by_structure:
            raise DataError("组分库中同一化学结构对应多个名称，请合并", [name, by_structure[canonical]])
        by_structure[canonical] = name
        fragments = Chem.GetMolFrags(mol, asMols=True)
        if (item["role"] in ("primary_salt", "salt_additive") and len(fragments) == 2
                and lithium in [Chem.MolToSmiles(x) for x in fragments]):
            other = next(x for x in fragments if Chem.MolToSmiles(x) != lithium)
            if Chem.GetFormalCharge(other) == -1:
                anions[Chem.MolToSmiles(other)] = name

    def resolve(name, smiles):
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            raise DataError(f"{name or '空名称'} 的 SMILES 缺失或无效")
        canonical = Chem.MolToSmiles(mol)
        if name in known:
            if canonical != Chem.MolToSmiles(Chem.MolFromSmiles(known[name]["smiles"])):
                raise DataError(f"{name} 的 SMILES 与绑定身份不一致")
            return name, "molecule"
        if canonical == lithium:
            return "Li+", "cation"
        if canonical in anions:
            return anions[canonical], "anion"
        if canonical in by_structure:
            return by_structure[canonical], "molecule"
        raise DataError("无法识别组分，请绑定设计或提供组分身份及角色 JSON", [name, smiles])
    return resolve


def number(value, label, nonnegative=True):
    try:
        result = float(value)
    except (ValueError, TypeError):
        raise DataError(f"{label} 必须为数值")
    if not math.isfinite(result) or (nonnegative and result < 0):
        raise DataError(f"{label} 必须为有限{'非负' if nonnegative else ''}数值")
    return result


def integer(value, label, low, high):
    result = number(value, label)
    if isinstance(value, bool) or int(result) != result or not low <= result <= high:
        raise DataError(f"{label} 必须为 {low}–{high} 的整数")
    return int(result)


def component_config(config):
    if config.get("components"):
        return config["components"]
    if "salts" in config or "salt" in config:
        salts = config.get("salts", [config.get("salt")])
        return ([dict(x, role="primary_salt") for x in salts]
                + [dict(x, role="balance_solvent" if x.get("balance") else "cosolvent")
                   for x in config.get("solvents", [])]
                + [dict(x, role=x.get("role", "functional_additive")) for x in config.get("additives", [])])
    return []


def validate_formulations(directory, payload):
    state = read_json(directory / "status.json")
    if state["status"] == "frozen":
        raise DataError("数据集已冻结，请重新上传以修改配方")
    state = {k: v for k, v in state.items() if k not in ("target", "N", "target_report", "formulation_report", "property_columns", "available_feature_bases")}
    state["status"] = "uploaded"
    save_json(directory / "status.json", state)
    fields, raw = table(directory / "source.csv")
    active, skipped = eligible_rows(raw)
    if not active:
        raise DataError("没有可校验的配方", ["EXAMPLE 示例和全空行不参与校验，请填写正式实验行。"])
    kind = payload.get("input_format")
    if kind not in state["formats"]:
        raise DataError("请选择表格支持的配方格式")
    mode = payload.get("validation_mode", "chemistry" if payload.get("design_id") else "composition")
    if mode not in ("composition", "chemistry"):
        raise DataError("校验方式必须为仅组分数据或按设计校验")
    state["validation_mode"] = mode
    save_json(directory / "status.json", state)
    if mode == "composition":
        from offline_benchmark.composition import validate
        return validate(directory, payload, state)
    config = payload.get("design_config") or {}
    if not isinstance(config, dict):
        raise DataError("设计配置必须为 JSON 对象")
    custom = component_config(config)
    if not isinstance(custom, list) or any(not isinstance(x, dict) or not x.get("name") for x in custom):
        raise DataError("组分配置必须是含 name、role 和 smiles 的对象列表")
    if custom and len({x["name"] for x in custom}) != len(custom):
        raise DataError("绑定配置中的组分名称重复")
    known = {name: {"name": name, "role": role, "smiles": smiles} for name, role, smiles in LIBRARY}
    if custom:
        known = {x["name"]: x for x in custom}
    mass_fields = [x for x in fields if x.endswith("_mass_g") and x not in MASS_META]
    fraction_fields = [x for x in fields if x.endswith("_mass_fraction")]
    issues, warnings, records = [], [], []
    for field in ("temperature[C]", "temperature_C"):
        if field in fields:
            temperatures = {number(row[field], f"第 {i+2} 行 {field}", False) for i, row in active}
            if len(temperatures) != 1:
                raise DataError("数据包含不同温度，请先按实验条件拆分")
            if "temperature_C" in config and not math.isclose(next(iter(temperatures)), float(config["temperature_C"]), abs_tol=1e-8):
                raise DataError("实验温度与绑定设计不一致")
    if kind == "experiment_recipe":
        chosen = mass_fields or fraction_fields
        suffix = "_mass_g" if mass_fields else "_mass_fraction"
        names = [x[:-len(suffix)] for x in chosen]
        unknown = sorted(set(names) - set(known))
        if unknown:
            raise DataError("组分不在绑定的化学身份库中，请提供设计配置 JSON", unknown)
    else:
        slots = sorted(int(x[9:]) for x in fields if re.fullmatch(r"compound_\d+", x))
        if slots != list(range(len(slots))):
            raise DataError("候选池槽位必须从 0 连续编号")
        basis = payload.get("catalog_basis", "mass")
        if basis not in ("mass", "mole", "unielf_mole", "ratio_mass"):
            raise DataError("请明确比例列的物理含义：质量、摩尔或 UniELF 离子摩尔分数")
        prefix = {"mass": "mass_ratio", "mole": "mole_ratio", "unielf_mole": "ratio", "ratio_mass": "ratio"}[basis]
        for j in slots:
            if f"{prefix}_{j}" not in fields or f"smiles_{j}" not in fields:
                raise DataError(f"槽位 {j} 缺少比例或 SMILES")
        resolve = identity_resolver(known)
        names, mapping, resolution_errors = set(), set(), []
        for i, row in active:
            for j in slots:
                name, smiles = row[f"compound_{j}"].strip(), row[f"smiles_{j}"].strip()
                if not name and not smiles:
                    continue
                try:
                    canonical, mode = resolve(name, smiles)
                    if mode != "cation":
                        names.add(canonical)
                    if canonical != name or mode != "molecule":
                        mapping.add(f"{name} → {canonical}" + ("（阴离子配对）" if mode == "anion" else ""))
                except DataError as exc:
                    resolution_errors.append(f"第 {i+2} 行槽位 {j}：{exc} {'；'.join(exc.details)}")
        if resolution_errors:
            raise DataError("候选组分身份校验未通过", resolution_errors)
        names = sorted(names)
        if mapping:
            warnings.append("按 SMILES 对齐化学身份：" + "；".join(sorted(mapping)))
        if basis in ("mole", "unielf_mole"):
            warnings.append(f"{prefix}_i 按摩尔分数解释；1:1 锂盐阴离子与 Li+ 配对后转换为完整盐的质量／摩尔特征。")
    if len(set(names)) != len(names) or not names:
        raise DataError("组分列为空或重名")
    # Preserve the bound library order; zero-only design components remain in X.
    if custom:
        if set(names) - set(known):
            raise DataError("存在设计之外的组分")
        names = list(known)
    else:
        names = [x[0] for x in LIBRARY if x[0] in names]
    snapshot = {**config, "components": [known[name] for name in names]}
    save_json(directory / "config_snapshot.json", snapshot)
    _, components = load_components(directory / "config_snapshot.json")
    weights = {c["name"]: c["molecular_weight"] for c in components}
    roles = {c["name"]: c["role"] for c in components}
    min_salts = integer(payload.get("min_salts", 1), "最少盐种数", 1, 12)
    max_salts = integer(payload.get("max_salts", 2), "最多盐种数", min_salts, 12)
    required = payload.get("require_all_solvents", True)
    if type(required) is not bool:
        raise DataError("溶剂必含设置必须为布尔值")
    # Bound design settings take precedence over the explicit base-rule controls.
    if config:
        sizes = config.get("salt_selection", {}).get("sizes")
        required = config.get("sampling", {}).get("require_all_solvents", False)
    else:
        sizes = list(range(min_salts, max_salts + 1))
    schema_bad = design_bad = 0
    for index, row in active:
        line = index + 2
        try:
            if kind == "experiment_recipe":
                amounts, _ = recipe_composition(row, names, suffix, line)
                values = dict(zip(names, amounts))
            else:
                values = dict.fromkeys(names, 0.0)
                active_names = set()
                original_slots = []
                lithium, anion_total, split_names, input_total = 0.0, 0.0, set(), 0.0
                for j in slots:
                    name, smiles = row[f"compound_{j}"].strip(), row[f"smiles_{j}"].strip()
                    text = row[f"{prefix}_{j}"].strip()
                    value = number(text or (0 if not name and not smiles else None), f"第 {line} 行槽位 {j}")
                    input_total += value
                    if not name:
                        if value > 0 or smiles:
                            raise DataError(f"第 {line} 行槽位 {j}：SMILES 或非零比例对应空组分名")
                        original_slots.append((j, 1.0, 0.0))
                        continue
                    canonical, mode = resolve(name, smiles)
                    original_slots.append((j, Descriptors.MolWt(Chem.MolFromSmiles(smiles)), value))
                    if (canonical, mode) in active_names:
                        raise DataError(f"第 {line} 行组分 {name} 槽位重复")
                    active_names.add((canonical, mode))
                    if mode != "molecule" and basis not in ("mole", "unielf_mole"):
                        raise DataError(f"第 {line} 行：拆分离子需使用 mole_ratio_i 或 UniELF ratio_i 摩尔分数模式")
                    if mode == "cation":
                        lithium += value
                    else:
                        values[canonical] += value
                        if mode == "anion":
                            anion_total += value
                            split_names.add(canonical)
                if abs(input_total - 1) > 1e-8:
                    raise DataError(f"第 {line} 行候选比例之和不为 1")
                if not math.isclose(lithium, anion_total, rel_tol=1e-7, abs_tol=1e-8):
                    raise DataError(f"第 {line} 行：Li+ 与 1:1 盐阴离子摩尔数不平衡")
                if any((n, "molecule") in active_names and (n, "anion") in active_names for n in split_names):
                    raise DataError(f"第 {line} 行：同一种盐同时以完整盐和拆分离子表示")
                # Cross-check in the ORIGINAL species space. After ion recombination,
                # a salt's mole fraction differs from its anion's input mole fraction.
                molecular_basis = basis in ("mole", "unielf_mole")
                raw_masses = [(j, v*mw if molecular_basis else v) for j, mw, v in original_slots]
                raw_moles = [(j, v if molecular_basis else v/mw) for j, mw, v in original_slots]
                for family, quantities in (("mass", raw_masses), ("mole", raw_moles)):
                    quantity_total = sum(v for _, v in quantities)
                    for j, value in quantities:
                        key = f"{family}_ratio_{j}"
                        if key in fields and not math.isclose(number(row[key] or 0, f"第 {line} 行 {key}"), value/quantity_total, rel_tol=0, abs_tol=1e-8):
                            raise DataError(f"第 {line} 行槽位 {j} 质量／摩尔比例不一致")
                if basis in ("mole", "unielf_mole"):
                    values = {name: value * weights[name] for name, value in values.items()}
            total = sum(values.values())
            if total <= 0:
                raise DataError(f"第 {line} 行总质量为 0")
            ratios = {name: values[name] / total for name in names}
            if kind == "experiment_recipe" and mass_fields and fraction_fields:
                for name in names:
                    key = name + "_mass_fraction"
                    if key in fields and not math.isclose(number(row[key], key), ratios[name], abs_tol=1e-8):
                        raise DataError(f"第 {line} 行 {name} 的质量与分数不一致")
        except (DataError, RecipeError, KeyError) as exc:
            schema_bad += 1
            issues.extend(getattr(exc, "details", None) or [str(exc)])
            continue
        salts = [n for n in names if roles[n] == "primary_salt" and ratios[n] > 0]
        solvents = [n for n in names if roles[n] in ("balance_solvent", "cosolvent")]
        problems = []
        if (sizes is not None and len(salts) not in sizes) or (sizes is None and not min_salts <= len(salts) <= max_salts):
            problems.append("盐种数不符合绑定空间")
        if not solvents or sum(ratios[n] for n in solvents) <= 0:
            problems.append("缺少溶剂")
        if required and any(ratios[n] <= 0 for n in solvents):
            problems.append("未包含全部要求的溶剂（"+"、".join(n for n in solvents if ratios[n] <= 0)+"）；若本数据允许溶剂组合变化，请取消全部溶剂必含或绑定对应设计")
        solvent_total = sum(ratios[n] for n in solvents)
        for item in config.get("salts", [config["salt"]] if "salt" in config else []) + config.get("additives", []):
            low, high = item["final_mass_fraction_bounds"]
            value = ratios.get(item["name"], 0)
            conditional = "salt_selection" in config and item["name"] in [n for n in names if roles[n] == "primary_salt"]
            if not (conditional and value == 0) and not low - 1e-8 <= value <= high + 1e-8:
                problems.append(f"{item['name']} 超出质量分数范围")
        for item in config.get("solvents", []):
            low, high = item["pool_mass_fraction_bounds"]
            if solvent_total <= 0 or not low - 1e-8 <= ratios.get(item["name"], 0) / solvent_total <= high + 1e-8:
                problems.append(f"{item['name']} 超出溶剂池比例范围")
        if "salt_molality_bounds" in config and solvent_total > 0:
            molality = 1000 * sum(ratios[n]/weights[n] for n in salts) / solvent_total
            low, high = config["salt_molality_bounds"]
            if not low - 1e-8 <= molality <= high + 1e-8:
                problems.append("总盐 mol/kg 超出绑定范围")
        if problems:
            design_bad += 1
            issues.append(f"第 {line} 行：{'；'.join(problems)}")
            continue
        records.append({"source_index": index, "source_id": row.get("sample_id") or row.get("pool_formulation_id") or row.get("formulation_id") or row.get("experiment_id") or f"row_{line}",
                        "fractions": [ratios[n] for n in names]})
    report = {"input_format": kind, "rows": len(raw), "schema_valid": schema_bad == 0,
              "design_valid": design_bad == 0 and schema_bad == 0,
              "schema_invalid_rows": schema_bad, "design_invalid_rows": design_bad,
              "valid_rows": len(records), "eligible_rows": len(active), "skipped": skipped,
              "catalog_basis": basis if kind == "candidate_catalog" else None,
              "details": issues[:30], "warnings": warnings,
              "constraints": {"source": "bound_design" if config else "explicit_base_rules", "min_salts": min_salts,
                              "max_salts": max_salts, "require_all_solvents": required}}
    save_json(directory / "formulation_report.json", report)
    if issues:
        state["formulation_report"] = report
        save_json(directory / "status.json", state)
        raise DataError("配方校验未通过，请修正表格或绑定匹配的设计配置", issues)
    if mass_fields and kind == "experiment_recipe" and "compound_0" in fields:
        warnings.append("混合表已明确采用组分质量列；旧 compound_i / ratio_i 槽位未参与特征。")
    if "mass_basis" in fields:
        warnings.append("保留原质量来源标记；参考质量不视为原始实测称量。")
    if skipped["examples"] or skipped["blank"]:
        warnings.append(f"沿用在线模板规则：跳过 {skipped['examples']} 行 EXAMPLE 示例、{skipped['blank']} 行全空数据；目标缺失在下一步单独检查。")
    grouped = {}
    source_keys = {}
    for record in records:
        key = tuple(round(x, 10) for x in record["fractions"])
        sid = record["source_id"]
        if sid in source_keys and source_keys[sid] != key:
            raise DataError("同一来源 ID 对应不同配方，无法可靠查表", [sid])
        source_keys[sid] = key
        grouped.setdefault(key, []).append(record)
    canonical = []
    groups = []
    for key in sorted(grouped):
        # Use the same precision and mass→mole conversion as online recipe imports.
        fractions = grouped[key][0]["fractions"]
        pool_fields = ["presence_pattern", *(n+s for n in names for s in ("_mass_g", "_mass_fraction", "_mole_fraction"))]
        row = pool_row(pool_fields, components, {"fractions": fractions, "masses": fractions})
        canonical.append(row)
        groups.append(grouped[key])
    _, _, catalog_fields, catalog, _ = transform_rows(canonical, components, "mass", 1e-12, 1e-8, 10)
    write_csv(directory / "formulation_catalog.csv", catalog_fields, catalog)
    save_json(directory / "groups.json", groups)
    mapping = [{"source_row": item["source_index"]+2, "source_id": item["source_id"],
                "sample_id": recipe["sample_id"], "replicates": len(group)}
               for group, recipe in zip(groups, catalog) for item in group]
    write_csv(directory / "formulation_mapping.csv", ["source_row", "source_id", "sample_id", "replicates"], mapping)
    preview = [{"sample_id": row["sample_id"], **{name: row[f"mass_ratio_{j}"] for j, name in enumerate(names)}} for row in catalog[:5]]
    report.update(unique_formulations=len(catalog), components=names, warnings=warnings,
                  replicates_merged=len(records)-len(catalog), preview=preview,
                  component_details=[{"name": n, "role": known[n]["role"], "smiles": known[n]["smiles"]} for n in names])
    save_json(directory / "formulation_report.json", report)
    candidates = property_columns(fields, [row for _, row in active])
    state.update(status="formulation_validated", formulation_report=report, property_columns=candidates,
                 input_format=kind, N=len(catalog), available_feature_bases=["mass", "mole"])
    save_json(directory / "status.json", state)
    return state


def validate_target(directory, payload):
    state = read_json(directory / "status.json")
    if state["status"] not in ("formulation_validated", "target_validated"):
        raise DataError("请先通过配方校验；冻结数据不可修改")
    state["status"] = "formulation_validated"
    state.pop("target", None)
    state.pop("target_report", None)
    save_json(directory / "status.json", state)
    fields, raw = table(directory / "source.csv")
    column = payload.get("column")
    if column not in state["property_columns"]:
        raise DataError("请选择原表中独立于组分特征的目标性质列")
    direction = payload.get("direction", "max")
    source = payload.get("source", "measured")
    if direction not in ("max", "min") or source not in ("measured", "external_simulation"):
        raise DataError("目标方向或来源无效")
    if column.endswith("_pred") and source != "external_simulation":
        raise DataError("预测列属于外部模拟／计算值，请修改目标来源后再校验")
    groups = read_json(directory / "groups.json")
    catalog_fields, catalog = table(directory / "formulation_catalog.csv")
    valid, truth, missing, errors, replicates = [], [], [], [], []
    for group, recipe in zip(groups, catalog):
        values = []
        for item in group:
            text = raw[item["source_index"]][column].strip()
            if text:
                try:
                    value = number(text, f"第 {item['source_index']+2} 行 {column}", False)
                    values.append(value)
                    replicates.append({"sample_id": recipe["sample_id"], "source_id": item["source_id"], "value": value})
                except DataError as exc:
                    errors.append(str(exc))
        if not values:
            missing.append(recipe["sample_id"])
            continue
        valid.append(recipe)
        truth.append({"sample_id": recipe["sample_id"], "target_value": sum(values) / len(values)})
    report = {"column": column, "valid_formulations": len(valid), "missing_formulations": missing,
              "invalid_values": errors, "replicate_measurements": len(replicates)}
    save_json(directory / "target_report.json", report)
    if errors or (missing and not payload.get("exclude_missing", False)):
        raise DataError("目标校验未通过；请修正无效值或明确排除缺失目标配方", errors + missing)
    if len(valid) < 3:
        raise DataError("离线回放至少需要 3 个有目标值的唯一配方")
    values = [x["target_value"] for x in truth]
    if max(values) == min(values):
        raise DataError("全部目标值相同，没有寻优区分度")
    # Keep only feature/recipe metadata in catalog; even the target column name is absent.
    write_csv(directory / "pool_catalog.csv", catalog_fields, valid)
    write_csv(directory / "truth.csv", ["sample_id", "target_value"], truth)
    write_csv(directory / "replicate_records.csv", ["sample_id", "source_id", "value"], replicates)
    target = {"column": column, "engine_column": "offline_target", "direction": direction,
              "source": source, "unit": str(payload.get("unit", ""))[:30]}
    report.update(min=min(values), max=max(values), best=max(values) if direction == "max" else min(values))
    save_json(directory / "target_report.json", report)
    state.update(status="target_validated", target=target, target_report=report, N=len(valid))
    save_json(directory / "status.json", state)
    return state


def freeze(directory):
    state = read_json(directory / "status.json")
    if state["status"] != "target_validated":
        raise DataError("配方和目标都通过校验后才能冻结")
    files = ("source.csv", "config_snapshot.json", "formulation_report.json", "target_report.json", "pool_catalog.csv", "truth.csv")
    hashes = {name: hashlib.sha256((directory / name).read_bytes()).hexdigest() for name in files}
    manifest = {"files": hashes, "target": state["target"], "N": state["N"],
                "available_feature_bases": state.get("available_feature_bases", ["mass", "mole"])}
    save_json(directory / "manifest.json", manifest)
    state.update(status="frozen", fingerprint=hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest())
    save_json(directory / "status.json", state)
    return state


def verify_frozen(directory):
    state = read_json(directory / "status.json")
    if state["status"] != "frozen":
        raise DataError("数据集尚未完成两阶段校验与冻结")
    manifest = read_json(directory / "manifest.json")
    if "available_feature_bases" in manifest and state.get("available_feature_bases", ["mass", "mole"]) != manifest["available_feature_bases"]:
        raise DataError("冻结数据集的可用特征基准已改变")
    for name, digest in manifest["files"].items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != digest:
            raise DataError("冻结数据集已改变，不能继续回放")
    return state
