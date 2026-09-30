# 步骤 3：实验导入与回传

`validate_upload` 要求上传 CSV 的列顺序为候选池所有原始列加所选性能列。逐行检查 `sample_id` 属于候选池、配方字段与候选池一致、目标数值有限且无重复。初始实验可只含已测行；推荐回传必须覆盖本轮推荐的每一行和全部所选目标。

`write_table` 使用同目录临时文件原子替换；`synchronize_feedback` 调用 `bo_utils.sync_observations` 合并初始实验与历轮推荐，更新运行目录的 `observation.csv`。`simulate_feedback` 调用原 `simu_experiment.py`，供演示使用。网页同时保留上传原件于 `feedback_uploads/` 并在运行状态中记录 `real` 或 `simulated` 来源。

后续扩展重点：目标单位/列名映射、上传审核、部分批次回传、重复实验和异常值处理。扩展前需确定 `sample_id` 与原始配方的不可变约束。
