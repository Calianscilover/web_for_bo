# 步骤 1：生成配方空间

`actions.py::generate_formulations(config, design_dir)` 调用 `pool/generate_pool.py`。前端设计 JSON 包含锂盐/溶剂名称和 SMILES、盐最终质量分数范围、每个溶剂在溶剂池中的分数范围、可选的添加剂列表（`additives`：名称、SMILES、`role` 为 `functional_additive` 或 `salt_additive`、`final_mass_fraction_bounds`）、称量步长、Sobol 采样规模和随机种子。盐与添加剂的范围以电解液总质量为分母；溶剂池分数以扣除盐与添加剂后的质量为分母，池内之和为 1。组分顺序为盐、溶剂、添加剂；没有添加剂时 Sobol 维度与旧版一致，同一种子得到相同候选池。

输出为 `config_snapshot.json`、`components.csv`、`feasible_candidates.csv`、`pool.csv`。下一步读取 `config_snapshot.json` 和 `pool.csv`。改造采样策略或组分约束时，先保持这两个输出契约，再替换 `generate_pool.py` 的计算。温度仅记录，不参与现有 BO 特征。
