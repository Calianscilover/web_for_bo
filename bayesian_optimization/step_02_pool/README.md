# 步骤 2：候选池转换

`actions.py::convert_formulations(design_dir, feature_basis="mass")` 调用原 `pool/convert_pool.py`，将步骤 1 的采样结果转换为 `converted/pool_catalog.csv`、`pool_features.csv` 和 `pool_manifest.json`。目录及列结构是步骤 3、4 的输入契约；`sample_id` 是实验回传与候选配方关联的键。

`import_demo` 把现有七元电解液示例的候选池和 `experiment.csv` 复制进一个独立设计目录，用来复现当前已跑通的工作流，不修改源示例。改变输入维度时，优先修改生成/转换逻辑并维持 `pool_catalog.csv` 与实验文件的列对齐，不在前端写死七元特征。
