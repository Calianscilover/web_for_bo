"""File-defined composition features, independent of a salt/solvent identity registry."""
from __future__ import annotations

import hashlib
import json
import math
import re

from rdkit import Chem
from rdkit.Chem import Descriptors


def validate(directory, payload, state):
    # Imported lazily: dataset dispatches here after resetting the downstream gates.
    from offline_benchmark.dataset import (DataError, eligible_rows, number, property_columns,
                                           save_json, table, write_csv, MASS_META)
    from pool.import_experiments import RecipeError, recipe_composition

    fields, raw = table(directory / "source.csv")
    active, skipped = eligible_rows(raw)
    kind = payload.get("input_format")
    if kind not in state["formats"]:
        raise DataError("请选择表格支持的配方格式")
    if not active:
        raise DataError("没有可校验的正式配方")
    if payload.get("design_id"):
        raise DataError("需要绑定已有设计时，请切换到按设计校验模式")
    config = payload.get("design_config") or {}
    if not isinstance(config, dict):
        raise DataError("组分信息必须为 JSON 对象")
    definitions = config.get("components", [])
    if not isinstance(definitions, list) or any(not isinstance(c, dict) or not c.get("name") for c in definitions):
        raise DataError("可选 components 信息必须包含组分名称；换算摩尔分数时需填写 smiles")
    if len({c["name"] for c in definitions}) != len(definitions):
        raise DataError("可选组分信息中的名称重复")
    definitions = {c["name"]: c for c in definitions}
    temperatures = {}
    for field in ("temperature[C]", "temperature_C"):
        if field in fields:
            temperatures[field] = {number(row[field], f"第 {i+2} 行 {field}", False) for i, row in active}
            if len(temperatures[field]) != 1:
                raise DataError("当前模型只使用组分特征；不同温度的数据请先按条件拆分，避免合并为同一配方")
    if len(temperatures) == 2 and temperatures["temperature[C]"] != temperatures["temperature_C"]:
        raise DataError("两列温度记录不一致")

    identities, name_keys, aliases = {}, {}, {}
    parsed, issues, warnings = [], [], []

    def identify(name, smiles, label, require_smiles):
        if not name:
            raise DataError(f"{label} 缺少组分名称")
        if require_smiles and not smiles:
            raise DataError(f"{label} {name} 缺少 SMILES")
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None or mol.GetNumAtoms() == 0:
                raise DataError(f"{label} {name} 的 SMILES 无效")
            canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
            weight = Descriptors.MolWt(mol)
            if not math.isfinite(weight) or weight <= 0:
                raise DataError(f"{label} {name} 的分子量无效，无法换算比例")
            key = "smiles:"+canonical
        else:
            canonical, weight, key = "", None, "name:"+name
        if name in name_keys and name_keys[name] != key:
            raise DataError(f"{label}：同一组分名称 {name} 对应不同 SMILES，请先明确名称")
        name_keys[name] = key
        identities[key] = {"smiles": canonical, "molecular_weight": weight}
        aliases.setdefault(key, set()).add(name)
        return key

    mass_fields = [f for f in fields if f.endswith("_mass_g") and f not in MASS_META]
    fraction_fields = [f for f in fields if f.endswith("_mass_fraction")]
    if kind == "experiment_recipe":
        suffix = "_mass_g" if mass_fields else "_mass_fraction"
        names = [f[:-len(suffix)] for f in (mass_fields or fraction_fields)]
        if not names:
            raise DataError("没有找到组分质量或质量分数列")
        keys = [identify(n, str(definitions.get(n, {}).get("smiles") or "").strip(), "组分列", False) for n in names]
        basis = "mass"
        for index, row in active:
            try:
                _, fractions = recipe_composition(row, names, suffix, index+2)
                values = dict.fromkeys(identities, 0.0)
                for key, fraction in zip(keys, fractions):
                    values[key] += fraction
                if mass_fields and fraction_fields:
                    for name, fraction in zip(names, fractions):
                        field = name+"_mass_fraction"
                        if field in fields and not math.isclose(number(row[field], field), fraction, abs_tol=1e-8):
                            raise DataError(f"第 {index+2} 行 {name} 质量与质量分数不一致")
                parsed.append((index, values))
            except (DataError, RecipeError) as exc:
                issues.extend(getattr(exc, "details", None) or [str(exc)])
    else:
        basis_option = payload.get("catalog_basis", "mass")
        options = {"mass": ("mass_ratio", "mass"), "mole": ("mole_ratio", "mole"),
                   "unielf_mole": ("ratio", "mole"), "ratio_mass": ("ratio", "mass")}
        if basis_option not in options:
            raise DataError("请明确选择比例列及质量／摩尔基准")
        prefix, basis = options[basis_option]
        slots = sorted(int(f[9:]) for f in fields if re.fullmatch(r"compound_\d+", f))
        if not slots or slots != list(range(len(slots))):
            raise DataError("组分槽位必须从 0 连续编号")
        for j in slots:
            missing = [f for f in (f"smiles_{j}", f"{prefix}_{j}") if f not in fields]
            if missing:
                raise DataError(f"槽位 {j} 缺少列", missing)
        for index, row in active:
            line = index+2
            try:
                values, slot_values = {}, []
                for j in slots:
                    name, smiles = row[f"compound_{j}"].strip(), row[f"smiles_{j}"].strip()
                    text = row[f"{prefix}_{j}"].strip()
                    value = number(text or (0 if not name and not smiles else None), f"第 {line} 行槽位 {j}")
                    if not name and not smiles:
                        if value != 0:
                            raise DataError(f"第 {line} 行槽位 {j} 非零比例没有组分信息")
                        slot_values.append((j, 1.0, 0.0))
                        continue
                    key = identify(name, smiles, f"第 {line} 行槽位 {j}", True)
                    values[key] = values.get(key, 0.0)+value
                    slot_values.append((j, identities[key]["molecular_weight"], value))
                total = sum(values.values())
                if abs(total-1) > 1e-8 or total <= 0:
                    raise DataError(f"第 {line} 行比例之和为 {total:.10g}，应为 1")
                # Validate any second explicit ratio family in its original species space.
                for family in ("mass", "mole"):
                    amounts = [(j, v if family == basis else v*mw if family == "mass" else v/mw)
                               for j, mw, v in slot_values]
                    amount_total = sum(v for _, v in amounts)
                    for j, value in amounts:
                        field = f"{family}_ratio_{j}"
                        if field in fields and not math.isclose(number(row[field] or 0, f"第 {line} 行 {field}"), value/amount_total, rel_tol=0, abs_tol=1e-8):
                            raise DataError(f"第 {line} 行槽位 {j} 质量／摩尔比例不一致")
                parsed.append((index, {k: v/total for k, v in values.items()}))
            except DataError as exc:
                issues.append(str(exc))

    report = {"validation_mode": "composition", "input_format": kind, "rows": len(raw),
              "eligible_rows": len(active), "valid_rows": len(parsed), "skipped": skipped,
              "schema_valid": not issues, "design_valid": not issues,
              "schema_invalid_rows": len(active)-len(parsed), "design_invalid_rows": 0,
              "catalog_basis": payload.get("catalog_basis") if kind == "candidate_catalog" else None,
              "constraints": {"source": "composition_only"}, "warnings": warnings, "details": issues[:30]}
    if issues:
        save_json(directory/"formulation_report.json", report)
        state["formulation_report"] = report
        save_json(directory/"status.json", state)
        raise DataError("组分数据校验未通过", issues)

    ordered = sorted(identities)
    names = [sorted(aliases[k])[0] for k in ordered]
    weights = [identities[k]["molecular_weight"] for k in ordered]
    available = ["mass", "mole"] if all(weights) else [basis]
    warnings.append("组分由本次文件建立，不查询预设身份库；未执行盐种数、溶剂必含和设计比例范围校验。")
    if any(len(v)>1 for v in aliases.values()):
        warnings.append("相同规范 SMILES 的不同名称／槽位已合并，原始名称保留在组分映射中。")
    if kind == "candidate_catalog":
        warnings.append("保留文件中的分子／离子表示；离子不强制配对为完整盐，不从名称推断盐或溶剂角色。")
    if len(available) == 1:
        warnings.append("当前文件只有质量信息：可直接进行质量分数优化；未提供 SMILES 的组分不生成摩尔分数。")
    if skipped["examples"] or skipped["blank"]:
        warnings.append(f"跳过 {skipped['examples']} 行 EXAMPLE 示例、{skipped['blank']} 行空数据。")

    grouped, source_keys = {}, {}
    for index, values in parsed:
        native = [values.get(k, 0.0) for k in ordered]
        vectors = {basis: native}
        if len(available) == 2:
            other = "mole" if basis == "mass" else "mass"
            amounts = [v/w if basis == "mass" else v*w for v, w in zip(native, weights)]
            total = sum(amounts)
            vectors[other] = [v/total for v in amounts]
        group_key = tuple(round(x, 10) for x in vectors["mass"])
        row = raw[index]
        source_id = next((row[f].strip() for f in ("sample_id", "pool_formulation_id", "formulation_id", "experiment_id")
                          if row.get(f, "").strip()), f"row_{index+2}")
        if source_id in source_keys and source_keys[source_id] != group_key:
            raise DataError("同一来源 ID 对应不同配方，无法可靠查表", [source_id])
        source_keys[source_id] = group_key
        grouped.setdefault(group_key, []).append({"source_index": index, "source_id": source_id, "vectors": vectors})

    catalog, groups = [], []
    for group_key, records in sorted(grouped.items()):
        vectors = records[0]["vectors"]
        signature = json.dumps(["file_composition_v1", ordered, group_key], ensure_ascii=False)
        sid = "SMP_"+hashlib.sha256(signature.encode()).hexdigest()[:16]
        present = [key for key, value in zip(ordered, vectors["mass"]) if value > 1e-12]
        group_id = "CG_"+hashlib.sha256(json.dumps(present).encode()).hexdigest()[:12]
        entry = {"sample_id": sid, "chem_group_id": group_id,
                 "source_row": records[0]["source_index"]+2,
                 "presence_pattern": ";".join(n for n, v in zip(names, vectors["mass"]) if v > 1e-12)}
        for j, (name, key) in enumerate(zip(names, ordered)):
            entry.update({f"compound_{j}": name, f"smiles_{j}": identities[key]["smiles"]})
            for family in available:
                entry[f"{family}_ratio_{j}"] = vectors[family][j]
        catalog.append(entry)
        groups.append([{k: v for k, v in item.items() if k != "vectors"} for item in records])
    catalog_fields = list(catalog[0])
    write_csv(directory/"formulation_catalog.csv", catalog_fields, catalog)
    save_json(directory/"groups.json", groups)
    mapping = [{"source_row": item["source_index"]+2, "source_id": item["source_id"],
                "sample_id": entry["sample_id"], "replicates": len(group)}
               for group, entry in zip(groups, catalog) for item in group]
    write_csv(directory/"formulation_mapping.csv", ["source_row", "source_id", "sample_id", "replicates"], mapping)
    details = [{"name": n, "smiles": identities[k]["smiles"], "role": "unclassified",
                "aliases": sorted(aliases[k]), "identity_source": "file"} for n, k in zip(names, ordered)]
    save_json(directory/"config_snapshot.json", {"validation_mode": "composition", "components": details,
                                               "available_feature_bases": available, "temperatures": {k: list(v)[0] for k, v in temperatures.items()}})
    report.update(unique_formulations=len(catalog), components=names, component_details=details,
                  available_feature_bases=available, preview_basis="mass",
                  replicates_merged=len(parsed)-len(catalog),
                  preview=[{"sample_id": row["sample_id"], **{n: row[f"mass_ratio_{j}"] for j, n in enumerate(names)}} for row in catalog[:5]])
    save_json(directory/"formulation_report.json", report)
    state.update(status="formulation_validated", validation_mode="composition", input_format=kind,
                 formulation_report=report, available_feature_bases=available, N=len(catalog),
                 property_columns=property_columns(fields, [r for _, r in active]))
    save_json(directory/"status.json", state)
    return state
