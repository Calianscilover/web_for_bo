# `web/server.py`：API 与任务编排

## 当前职责

Flask 应用提供 `/api/v1` 接口。设计路由创建配方空间、加载七元示例、查询状态/候选池/文件/实验模板；观测路由校验并保存初始 `experiment.csv`，或把按组分称量填写的自有配方并入候选池后生成 `experiment.csv` 与 `experiment_mapping.csv`（已有运行的设计拒绝再次导入）；优化路由创建运行、查询状态、下载推荐/预测/可视化；回传路由保存真实或模拟结果并启动下一轮。具体动作委托给 `bayesian_optimization/step_01...step_04`。

`OUTPUT` 默认为仓库根目录的 `bo_test/`。`designs/<design_id>` 保存输入和候选池；`runs/<run_id>` 保存算法结果和轮次文件。`status.json` 是当前任务状态，`feedback_sources` 记录真实/模拟来源。网页通过轮询查询异步任务。ID 和可下载文件名都使用白名单校验。

## 主要 API

| 操作 | 路由 |
| --- | --- |
| 创建/列出/查询设计 | `POST/GET /api/v1/designs`、`GET /api/v1/designs/{id}` |
| 七元示例 | `POST /api/v1/designs/demo` |
| 候选池、文件、实验模板 | `GET /api/v1/designs/{id}/candidates`、`files/{name}`、`experiment-template` |
| 实验模板（列相同：experiment_id + 各组分 g + 目标） | `GET /api/v1/designs/{id}/experiment-template`（全部候选）、`recipe-template`（EXAMPLE 示例） |
| 实验导入（候选池内外配方统一入口） | `POST /api/v1/designs/{id}/experiment-recipes` |
| 旧版目录格式实验上传（仅 API 兼容，页面不再使用） | `POST /api/v1/designs/{id}/observations` |
| 创建/查询优化 | `POST /api/v1/designs/{id}/optimization-runs`、`GET /api/v1/optimization-runs/{id}` |
| 推荐/可视化/文件 | `GET /api/v1/optimization-runs/{id}/recommendations/{round}`、`visualization`、`files/{name}` |
| 实测/模拟回传与下一轮 | `POST /api/v1/optimization-runs/{id}/feedback`、`simulate`、`next-round` |

校验错误统一返回 400 和 `{code, message}`；异常带 `details` 属性时（如实验 CSV 的 `RecipeError`）一并返回 `details` 列表，页面逐条显示。实验导入、目标列名和启动优化前的检查使用中文提示。

## 修改位置

新增业务动作改四步模块的 `actions.py`；增加接口或改变响应字段改本文件并同步 `static/app.js`。新增可下载产物需要更新 `DESIGN_FILES` 或 `RUN_FILE`。改变状态结构时同时考虑已有运行目录的 `status.json` 向后兼容。

## 服务器迁移注意

当前 `ThreadPoolExecutor` 和 `run_locks` 只在一个 Python 进程中有效，重启会丢失队列状态；开发用 `app.run` 也不应当作为生产服务。先把输出根目录变成持久目录，再把后台任务换成独立 worker 和持久任务状态。详见 [部署与存储](DEPLOYMENT_AND_STORAGE.md)。
