# 步骤 3：实验导入与回传

`validate_upload` 要求上传 CSV 的列顺序为候选池所有原始列加所选性能列。逐行检查 `sample_id` 属于候选池、配方字段与候选池一致、目标数值有限且无重复。初始实验可只含已测行；推荐回传必须覆盖本轮推荐的每一行和全部所选目标。

`write_table` 使用同目录临时文件原子替换；`synchronize_feedback` 调用 `bo_utils.sync_observations` 合并初始实验与历轮推荐，更新运行目录的 `observation.csv`。`simulate_feedback` 调用原 `simu_experiment.py`，供演示使用。网页同时保留上传原件于 `feedback_uploads/` 并在运行状态中记录 `real` 或 `simulated` 来源。

## 自有配方（不在候选池）

`recipe_template` / `import_recipes` 包装 `pool/import_experiments.py`。配方 CSV 为 `experiment_id`、每个设计组分一列 `<名称>_mass_g`（或 `<名称>_mass_fraction`，二选一，未用组分留空）和目标列。脚本按 10 位小数的质量分数与 `pool.csv` 比对：相同配方沿用原行；新配方追加到 `pool.csv` 末尾，再用原转换设置重跑 `convert_pool`。`sample_id` 与 `chem_group_id` 都是组成的确定性哈希，因此新旧候选编号规则完全一致，原有行编号不变。重复配方的目标值取平均；超出设计范围的配方保留并在 `experiment_mapping.csv` 的 `within_design_bounds` 标记。含设计外组分的列直接报错，需要先在配方空间加入该组分并重新生成。页面的两种模板列完全相同：`recipe_template(..., source="pool")` 列出全部候选（`experiment_id` 预填 `sample_id`，质量取自 `pool.csv`，逐行回读都对应原候选）；`source="examples"` 带两行 `EXAMPLE-01/02` 示例（一行取自候选池、一行为按 0.01 g 称量的池外配方，未用组分留空）。导入时跳过 `EXAMPLE` 行、未填任何目标值的行和全空行，所以实验人员只需给做过的行填值或在示例下方追加。旧的目录格式文件（`compound_i` + `mass_ratio_i`，如七元 `experiment.csv`）也能导入；组分不属于当前设计时报错并列出这些组分。校验失败抛出 `RecipeError`（`ValueError` 子类），`message` 为中文概述，`details` 为逐条说明：非数字、负质量、填了目标值却没有组分、分数和不为 1、多出逗号等逐行问题会检查完整个文件后一并列出（最多 20 处），并注明行号（表头为第 1 行，与 Excel 行号一致）；有任何问题时不写入任何文件。空文件、Excel 文件、非逗号分隔、缺少目标列、缺少组分列、质量列与分数列混用也各有对应提示。命令行用法：

```bash
python pool/import_experiments.py --design-dir bo_test/designs/<id> --targets Conductivity --template recipe_template.csv
python pool/import_experiments.py --design-dir bo_test/designs/<id> --targets Conductivity --recipes lab.csv
```

后续扩展重点：目标单位/列名映射、上传审核、部分批次回传、重复实验和异常值处理。扩展前需确定 `sample_id` 与原始配方的不可变约束。
