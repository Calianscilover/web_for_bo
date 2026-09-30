# 步骤 4：训练、推荐与可视化

`command_for` 根据 `mode=single/multi` 组装参数数组，分别调用 `qLogNEI.py` 或 `qLogNEHVI.py`；`execute_round` 在独立运行目录中完成训练、候选打分、推荐，并调用 `visualize_bo.py` 生成 `visualization_N.json`。`refresh_visualization` 在实验回传后重新计算实测进展和历史视图。

输入是设计目录的 `experiment.csv`、`converted/pool_catalog.csv`，以及运行状态中的目标、方向、批量数、采样数、拟合迭代数与随机种子。主要输出是 `training_summary.json`、`model.pt`、`candidate_predictions*.csv`、`recommendation_N.csv`、`observation.csv` 和 `visualization_N.json`。`algorithm.log` 保存算法标准输出与错误信息。

迁服务器时可以直接把 `execute_round` 放进后台任务进程；API 负责入队及查询状态。若任务进程与 API 不在同一台机器，运行目录必须位于共享的持久对象存储或通过文件传输复制，且不能依赖当前进程内的任务锁。
