# 电解液优化本地网页

## 启动

在仓库根目录、已激活 `botorch` 环境（见根目录 `environment.yml`）下执行：

```bash
python bayesian_optimization/run.py --port 8765
```

浏览器打开 <http://127.0.0.1:8765>。服务默认只监听本机 `127.0.0.1`。所有网页任务输出保存在仓库根目录的 `bo_test/`（已加入 `.gitignore`），不会覆盖 `pool/converted/` 中的示例文件。可用 `--output` 改变输出根目录。

网页文件逐模块说明见 [docs/README.md](docs/README.md)；四步算法动作见 [bayesian_optimization/README.md](../bayesian_optimization/README.md)；服务器迁移与文件/SQLite 选择见 [docs/DEPLOYMENT_AND_STORAGE.md](docs/DEPLOYMENT_AND_STORAGE.md)。旧入口 `web/server.py` 仍可直接运行，调用的也是同一套四步动作。

## 流程

1. 填写锂盐、溶剂 SMILES 和范围，点击“生成候选池”；也可点击“载入七元示例”，用现有候选池和初始实验表快速演示。
2. 下载 `pool_catalog.csv` 检查候选配方。从“实验数据”处下载模板，仅保留已实验的行，填写目标性能值后上传。上传 CSV 必须保留候选池的全部原始列及顺序，最后是所选目标列。
3. 选择单目标或双目标、性能列和优化方向，点击“开始优化”。网页异步调用原有 `qLogNEI.py` 或 `qLogNEHVI.py`，在结果卡片预览并下载 `recommendation_1.csv`。
4. 真实实验完成后，可在推荐表下逐条填写全部目标值并保存，也可在下载的推荐 CSV 中填写后上传；“模拟实验回传”只用于演示。接口逐行核对 `sample_id` 和完整配方，保留回传 CSV 原件，并立即合并更新 `observation.csv`。只有本轮全部推荐完成回传后才允许点击“生成下一轮”；下一轮产生 `recommendation_2.csv`，后续轮次同理。
5. 推荐配方下方会显示代理模型训练样本的预测值与实测值散点，以及基于已测值的优化趋势。单目标趋势为累计最佳性能，双目标趋势为固定参考点下的已测 Pareto 超体积（两个目标时为面积）。双目标还可以展开查看已测 Pareto 前沿：图中显示全部已测配方，前沿点实心相连，最新一轮回传的配方加虚线圈，并注明是否有新配方进入前沿。每次真实或模拟回传后，进展、Pareto、候选预测和历史表都会按最新实测重新计算；拟合图对应已训练的模型轮次，在生成下一轮时更新。历史配方与本轮候选预测分开展示；模拟回传单独标记。

可视化脚本可独立运行：

```bash
python algorithms/visualize_bo.py --run-dir bo_test/runs/<run_id>
```

它输出 `visualization_<round>.json`。拟合散点是训练集后验均值诊断，不代表留出集预测精度；超体积只由已测值计算，参考点沿用 `training_summary.json` 中供 qLogNEHVI 使用的固定值。设计参考：[Bgolearn 文档](https://bgolearn.github.io/docs/) 与 [BoTorch 多目标优化教程](https://botorch.org/docs/v0.17.2/tutorials/multi_objective_bo)。

温度目前仅记录在设计配置中，不进入 BO 特征。新建空间支持一个锂盐、两个及以上溶剂和可选添加剂（功能添加剂或锂盐添加剂，范围以电解液总质量为分母）。实验数据可按候选池模板上传，也可按组分称量上传不在候选池中的自有配方，系统会把新配方并入候选池并生成一致的 `sample_id`/`chem_group_id`。目标名以 CSV 列名为准，示例使用 `Conductivity` 和 `logCE`，其中 `logCE` 不会自动等同于 `LCE`。

## 验证

```bash
python web/test_integration.py
```

测试使用临时目录，覆盖新撒点、候选转换、实验上传、单目标与双目标推荐、真实/模拟回传、第二轮及可视化数据。
