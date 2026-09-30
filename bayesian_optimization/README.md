# Bayesian Optimization 四步工作流

本目录把网页调用的算法动作按数据流拆开。`web/server.py` 负责 HTTP、任务状态和目录隔离；这里的 `actions.py` 负责实际的数据和算法动作。`pool/` 下的 `generate_pool.py`、`convert_pool.py` 与 `algorithms/` 下的 `qLogNEI.py`、`qLogNEHVI.py`、`bo_utils.py`、`simu_experiment.py`、`visualize_bo.py` 仍是计算实现，避免在迁移时出现两份算法。

```text
1 配方空间生成 → 2 候选池转换 → 3 实验数据导入/回传 → 4 模型训练/推荐/可视化
                                         ↑                 │
                                         └── 下一轮回传 ───┘
```

| 步骤 | 文件 | 主要动作 | 输入 | 输出 |
| --- | --- | --- | --- | --- |
| 1 | `step_01_formulation/actions.py` | `generate_formulations` | 设计 JSON、设计目录 | `config_snapshot.json`、`components.csv`、`feasible_candidates.csv`、`pool.csv` |
| 2 | `step_02_pool/actions.py` | `convert_formulations`、`import_demo` | 步骤 1 文件或七元示例文件 | `converted/pool_catalog.csv`、`pool_features.csv`、`pool_manifest.json` |
| 3 | `step_03_experiment/actions.py` | `validate_upload`、`write_table`、`synchronize_feedback`、`simulate_feedback` | 候选池、实验 CSV、回传 CSV | `experiment.csv`、回填的 `recommendation_N.csv`、`observation.csv` |
| 4 | `step_04_optimization/actions.py` | `command_for`、`execute_round`、`refresh_visualization` | 实验和候选池、运行参数 | `training_summary.json`、`model.pt`、预测/推荐 CSV、`visualization_N.json` |

开发启动入口（仓库根目录，已激活 `botorch` 环境）：

```bash
python bayesian_optimization/run.py --port 8765
```

仍可用 `web/server.py` 旧入口。两者调用同一套四步动作。指定独立数据根目录：`--output /path/to/persistent/bo_data`。这只是本地开发服务器入口；生产服务部署及存储选择见 [`../web/docs/DEPLOYMENT_AND_STORAGE.md`](../web/docs/DEPLOYMENT_AND_STORAGE.md)。

每轮需要完整的真实实验回传，步骤 3 才会把推荐合并到累计观测；步骤 4 再运行一次便产生 `recommendation_2.csv` 等。模拟数据只用于演示，状态会标记来源。温度目前只保存于配置，不进入模型特征。
