# 步骤 1：生成配方空间

`actions.py::generate_formulations(config, design_dir)` 调用 `pool/generate_pool.py`。前端设计 JSON 包含锂盐列表（`salts`：一至两种，每项为名称、SMILES、`final_mass_fraction_bounds`；单盐也可写成 `salt: {...}`）、溶剂名称和 SMILES、每个溶剂在溶剂池中的分数范围、可选的添加剂列表（`additives`：名称、SMILES、`role` 为 `functional_additive` 或 `salt_additive`、`final_mass_fraction_bounds`）、称量步长、Sobol 采样规模和随机种子。盐与添加剂的范围以电解液总质量为分母；溶剂池分数以扣除盐与添加剂后的质量为分母，池内之和为 1。组分顺序为盐、溶剂、添加剂。

`pool.csv` 每行另有 `salt_molality_mol_kg`：全部主锂盐摩尔数之和 ÷ 溶剂（基准溶剂 + 共溶剂）质量，单位 mol/kg，由称量质量精确计算；撒点仍按质量分数进行。`molality_range(config)` 按锂盐与添加剂的质量分数上下限给出该浓度的范围（下限取盐与添加剂下限，上限取两者上限），返回各盐的 `molecular_weights`，供页面在步骤 1 提示。

双盐设计中下限为 0 的盐、下限为 0 的共溶剂（基准溶剂除外）和添加剂是可选组分，每个候选以 `sampling.absence_probability`（默认 0.25，取值 0 ≤ p < 1）的概率不加入该组分，因此 4 元设计里也会有 3 元、2 元子配方，`presence_pattern` 相应变化。两种盐地位相同；某行两种可选盐都抽到缺席时，补回有/无维度取值较大的那一种，保证每个候选至少含一种盐。单盐设计的盐始终存在。Sobol 维度为：每种盐 1 维、非末位溶剂各 1 维、添加剂各 1 维，再加每个可选组分 1 维的“是否加入”（顺序为盐、共溶剂、添加剂）。加入的可选组分质量不低于其 `minimum_nonzero_mass_g`（可写在组分对象中，默认一个称量步长），取整后低于该值的行会被剔除。如果某行剩余溶剂的上限之和不足 1，缺席的共溶剂会被补回。单盐设计的维度布局与旧版相同；`absence_probability` 为 0 时维度也与旧版一致，同一种子得到与旧版相同的候选池。返回统计中的 `presence_patterns` 为候选池中不同组分组合的数量。

输出为 `config_snapshot.json`、`components.csv`、`feasible_candidates.csv`、`pool.csv`。下一步读取 `config_snapshot.json` 和 `pool.csv`。改造采样策略或组分约束时，先保持这两个输出契约，再替换 `generate_pool.py` 的计算。温度仅记录，不参与现有 BO 特征。
