# 步骤 4：训练、推荐与可视化

`command_for` 根据 `mode=single/multi` 组装参数数组，分别调用 `qLogNEI.py` 或 `qLogNEHVI.py`；`execute_round` 在独立运行目录中完成训练、候选打分、推荐，并调用 `visualize_bo.py` 生成 `visualization_N.json`。`refresh_visualization` 在实验回传后重新计算实测进展和历史视图。

输入是设计目录的 `experiment.csv`、`converted/pool_catalog.csv`，以及运行状态中的目标、方向和模型超参数。`command_for` 把两个脚本公开的全部建模参数都传过去：`batch_size`、`mc_samples`、`fit_maxiter`、`seed`、`pool_batch_size`、`feature_basis`（mass/mole）、`kernel`（default/rbf/matern）、`matern_nu`、`ard`、`lengthscale_init`、`noise_std`（单目标 1 个值，双目标每目标 1 个，缺省为自动估计），双目标另有 `ref_point`（缺省沿用脚本的自动参考点）。旧运行的状态缺少这些键时使用脚本默认值，行为与之前一致。数据路径、目标列与方向由服务层决定，`--feature-columns` 不开放（输入列由候选池的比例列自动确定）。主要输出是 `training_summary.json`、`model.pt`、`candidate_predictions*.csv`、`recommendation_N.csv`、`observation.csv` 和 `visualization_N.json`。`algorithm.log` 保存算法标准输出与错误信息。

迁服务器时可以直接把 `execute_round` 放进后台任务进程；API 负责入队及查询状态。若任务进程与 API 不在同一台机器，运行目录必须位于共享的持久对象存储或通过文件传输复制，且不能依赖当前进程内的任务锁。
