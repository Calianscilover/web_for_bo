# web_for_bo：电解液配方贝叶斯优化工作台

本地网页工作台：定义“一个锂盐 + 多个溶剂”的配方空间并撒点生成候选池，按模板上传实验结果，用 BoTorch 的 qLogNEI（单目标）或 qLogNEHVI（双目标）推荐下一批配方，回填实测值后继续多轮优化。

## 环境

所有代码都在专用的 `botorch` 环境中运行，版本锁定在 `environment.yml`（Python 3.12、PyTorch 2.4.0、BoTorch 0.18.1、GPyTorch 1.15.2、RDKit、Flask 等）。

```bash
conda env create -f environment.yml
conda activate botorch
```

不用 conda 时，可在 Python 3.12 的虚拟环境中执行 `pip install -r requirements.txt`。Linux 服务器只用 CPU 时，可先执行 `pip install torch==2.4.0 --index-url https://download.pytorch.org/whl/cpu`，避免下载 CUDA 版本。

## 启动

在仓库根目录执行：

```bash
python bayesian_optimization/run.py --port 8765
```

浏览器打开 <http://127.0.0.1:8765>。点击“载入七元示例”可以直接用自带的候选池与初始实验数据跑一轮。网页任务的输出写入 `bo_test/`，可用 `--output` 改到其他目录。

## 测试

```bash
python web/test_integration.py                                  # 网页 API 端到端流程
python -m unittest test_acquisition test_qLogNEI test_qLogNEHVI  # BO 脚本
(cd pool && python -m unittest test_convert_pool)               # 候选池转换
```

## 目录

| 路径 | 内容 |
| --- | --- |
| `web/` | Flask 服务 `server.py`、静态页面 `static/`、集成测试和模块文档 `docs/` |
| `bayesian_optimization/` | 网页调用的四步动作：撒点、候选池转换、实验上传与回传、训练推荐与可视化 |
| `pool/` | `generate_pool.py` 撒点、`convert_pool.py` 转换；`converted/` 是七元示例候选池和初始实验数据 |
| `qLogNEI.py`、`qLogNEHVI.py`、`bo_utils.py` | 单目标、双目标 BO 入口和共用的 GP、采集函数与数据回流工具 |
| `simu_experiment.py`、`visualize_bo.py` | 模拟实验回传（仅演示）、可视化数据生成 |
| `assets/` | 品牌 Logo 原图 |

## 文档

- [ELECTROLYTE_UI_API_DESIGN.md](ELECTROLYTE_UI_API_DESIGN.md)：前后端设计方案与 API 契约
- [web/README.md](web/README.md)：网页操作流程
- [bayesian_optimization/README.md](bayesian_optimization/README.md)：四步工作流
- [BO_SCRIPTS.md](BO_SCRIPTS.md)：命令行使用 BO 脚本
- [POOL_TO_EXPERIMENT_GUIDE.md](POOL_TO_EXPERIMENT_GUIDE.md)、[pool/POOL_CONVERSION_DESIGN.md](pool/POOL_CONVERSION_DESIGN.md)：候选池到实验的数据流
- [web/docs/DEPLOYMENT_AND_STORAGE.md](web/docs/DEPLOYMENT_AND_STORAGE.md)：迁移到服务器的部署与存储方案
