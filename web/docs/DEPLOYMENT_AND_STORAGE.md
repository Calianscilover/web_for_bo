# 从本地网页迁到服务器：部署与存储方案

## 1. 当前运行边界

当前入口是 `bayesian_optimization/run.py` 或 `web/server.py`，运行在本机 Flask 开发服务器。`web/server.py` 管理 HTTP、`status.json` 和进程内线程池；四步算法动作在 `bayesian_optimization/`；实际建模仍调用本地 Python 解释器下的 `qLogNEI.py`、`qLogNEHVI.py`。默认输出在仓库根目录的 `bo_test/`。

```text
浏览器 → Flask API → 线程池 → 四步动作 → 同一台机器上的 Python 算法
                         ↘ bo_test/designs 和 bo_test/runs
```

这条链路适合本地验证。服务重启时未执行的线程任务会消失；多个 Web 进程不会共享 `run_locks`；运行状态和文件目前没有用户权限边界。迁移时需要逐项解决。

## 2. 数据应该放在哪里

**首个服务器版本建议：文件放服务器的持久数据目录，SQLite 只保存任务元数据。** 不需要为了保存 CSV/JSON 先把其全部内容塞进数据库。

| 数据 | 建议主存储 | 原因 |
| --- | --- | --- |
| `config_snapshot.json`、`components.csv`、`pool.csv`、`pool_catalog.csv`、`experiment.csv`、`observation.csv`、`recommendation_N.csv`、预测 CSV | 服务器持久文件目录，后期可换对象存储 | 已有算法按文件路径读写；便于下载、版本留存、独立复算 |
| `model.pt`、可视化 JSON、日志、原始上传 CSV | 同一持久文件目录或对象存储 | 模型二进制和较大文件不适合作为 SQLite 记录字段；上传原件要保留供审计 |
| `design_id`、`run_id`、用户/项目、目标列、方向、轮次、状态、任务时间、文件相对路径、校验摘要 | SQLite（单服务器） | 方便查询、权限、筛选、事务化更新及重启恢复 |
| 少量指标摘要，如样本数、最优值、超体积 | SQLite 可选；JSON 仍可作为算法产物 | 列表页和报表查询快，但不取代可复算文件 |

SQLite 不是必需的第一个迁移步骤。仅一台服务器、少量用户时，可以先沿用 `status.json` 与持久目录，把服务器跑通；在出现登录、多用户、任务检索和审计需求后再加 SQLite。若一开始就需要多人协作，建议同步引入 SQLite 元数据。SQLite 的 WAL 模式允许读写并发，但只有一个写入者，且数据库文件应位于同一主机的本地文件系统；不应将 SQLite 文件放在跨主机网络共享目录。参考 [SQLite 官方 WAL 文档](https://www.sqlite.org/wal.html)。

若 API、后台 worker 分布在多台机器，推荐 PostgreSQL 保存元数据，持久对象存储保存 CSV/JSON/模型，队列保存任务。不要让多台机器直接共享一个 SQLite 数据库文件。这里是后续规模化架构选择，不要求当前就建设。

## 3. 推荐的服务器文件布局

可通过 `--output /srv/electrolyte/data` 指向专用持久盘；实际目录可由运维配置，不能由浏览器提交任意绝对路径。

```text
/srv/electrolyte/data/
  designs/<design_id>/
    config_snapshot.json
    components.csv
    feasible_candidates.csv
    pool.csv
    converted/pool_catalog.csv
    converted/pool_features.csv
    converted/pool_manifest.json
    experiment.csv
    status.json                # 迁数据库后可保留导出快照
  runs/<run_id>/
    status.json
    observation.csv
    recommendation_1.csv
    recommendation_2.csv
    candidate_predictions.csv
    training_summary.json
    visualization_1.json
    model.pt
    algorithm.log
    feedback_uploads/round_1_uploaded.csv
```

下载接口只接受白名单产物名和服务端生成的 ID。文件写入采用临时文件后原子替换；跨文件更新（回传 CSV、`observation.csv`、状态）将来需要任务状态机、幂等键和恢复逻辑。按设计/运行目录进行备份，不要只备份数据库；恢复要连同文件和数据库元数据一起验证。实验原始上传、算法版本、随机种子、依赖版本、输入文件摘要应保留，便于复现。

## 4. SQLite 最小元数据模型（需要时再实施）

```text
designs(id, owner_id, name, kind, status, temperature_c, created_at, updated_at,
        config_path, catalog_path, experiment_path)
runs(id, design_id, mode, targets_json, directions_json, batch_size,
     mc_samples, fit_maxiter, seed, hyperparameters_json, current_round, status, error,
     created_at, updated_at)
rounds(run_id, round_number, recommendation_path, observation_path,
       visualization_path, feedback_source, status, created_at, completed_at)
artifacts(id, design_id, run_id, round_number, kind, relative_path,
          sha256, byte_size, created_at)
```

`relative_path` 只在服务器配置的数据根目录内解析。目标列表可先用 JSON 字符串字段，若以后需要跨项目按目标筛选再拆表。`observation.csv` 是本轮训练输入快照，不应由数据库拼装后悄悄覆盖；如果以后把结构化实测数据作为数据库主数据，需要定义明确的导出版本和双向一致性规则。

## 5. 迁移实施顺序

1. **服务器复现环境。** 按仓库根目录的 `environment.yml` 固定 Python 版本、PyTorch/BoTorch/RDKit/Flask 依赖，确定 CPU/GPU 选择；先把七元示例和当前端到端测试跑通。部署包即本仓库全部内容。
2. **外置数据目录。** 启动时传 `--output`，目录赋予运行用户读写权限；将现有 `bo_test/designs`、`bo_test/runs` 原样复制到持久盘并核对文件与 ID。避免将数据写进临时容器层或 Git 仓库。
3. **拆分 Web 与 worker。** API 只创建任务和返回 ID；后台 worker 执行步骤 1、2、4，步骤 3 的回传合并也在同一 run 的串行锁内完成。队列要能持久化与重试；记录任务提交时间、开始/结束时间和失败原因。`ThreadPoolExecutor` 只作本地开发用。
4. **用生产 WSGI 服务和反向代理。** Flask 官方明确开发服务器不用于生产；使用 Gunicorn/uWSGI/Waitress 等 WSGI 服务，前面配置 HTTPS、上传大小限制、超时和静态资源缓存。参考 [Flask 官方部署文档](https://flask.palletsprojects.com/en/stable/deploying/)。当前 `app` 可被 WSGI 导入，但在引入多 worker 前须先解决进程内队列与锁。
5. **加身份和权限。** 设计和运行都归属用户或项目；文件下载、反馈上传和下一轮启动按所有者授权。日志应记录操作者、操作时间和文件摘要，不记录不必要的敏感实验内容。
6. **按需要引入 SQLite。** 将 `status.json` 中的可查询元数据迁入数据库，文件仍保留为算法输入/输出；执行数据库迁移与文件备份恢复演练。多机部署时升级为 PostgreSQL 与对象存储。

## 6. 迁移后仍需保持的接口契约

- 浏览器请求 `/api/v1`；初始实验 CSV 必须与候选池全部原始列及顺序一致，并在末尾添加性能列。
- 每个 `sample_id` 对应不可变配方；推荐 CSV 的目标列起初为空，完整回传后才允许下一轮。
- `mode=single` 调 `qLogNEI.py`，`mode=multi` 调 `qLogNEHVI.py`；目标名和方向保持与训练摘要、图表一致。
- `recommendation_N.csv`、`observation.csv`、`visualization_N.json` 命名与轮次一致；前端只消费 API 和下载链接，不假设文件在本机。
- 模拟实验来源必须标记；真实实验数据不可被模拟回传覆盖。

最小迁移验收：新服务器上创建候选池、上传初始实验、得到第一轮推荐、回填真实值、生成第二轮；重启 Web 进程后仍能查询设计/运行和下载所有历史产物，并能从失败任务识别错误原因。
