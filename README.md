# web_for_bo：电解液配方贝叶斯优化工作台

本地网页工作台：定义“一至两种锂盐 + 多个溶剂 + 可选添加剂”的配方空间并撒点生成候选池，按模板上传实验结果，用 BoTorch 采集函数（单目标默认 qLogNEI，多目标 2–4 个默认 qLogNEHVI，可在网页中换成 qUCB、qLogEHVI、qLogNParEGO 等预设）推荐下一批配方，回填实测值后继续多轮优化。

## 环境

所有代码都在专用的 `botorch` 环境中运行，版本锁定在 `environment.yml`（Python 3.12、PyTorch 2.4.0、BoTorch 0.18.1、GPyTorch 1.15.2、RDKit、Flask 等）。

```bash
conda env create -f environment.yml
conda activate botorch
```

不用 conda 时，可在 Python 3.12 的虚拟环境中执行 `pip install -r requirements.txt`。Linux 服务器只用 CPU 时，可先执行 `pip install torch==2.4.0 --index-url https://download.pytorch.org/whl/cpu`，避免下载 CUDA 版本。

## 启动

在仓库根目录、已激活 `botorch` 环境下执行：

```bash
python bayesian_optimization/run.py --port 8765
```

模型优化会用启动服务的同一个 Python 以子进程运行 `algorithms/` 中的脚本，所以必须用装有 BoTorch 的解释器启动。若用 base 环境的 `python` 启动，页面能生成候选池、载入示例，但“开始优化”会失败，日志为 `ModuleNotFoundError: No module named 'botorch'`。可先用 `python -c "import botorch"` 确认当前环境。

浏览器打开 <http://127.0.0.1:8765>。点击“载入七元示例”可以直接用自带的候选池与初始实验数据跑一轮。网页任务的输出写入 `bo_test/`，可用 `--output` 改到其他目录。

## 配方表示与锂盐浓度

撒点与 BO 特征都使用质量分数（或摩尔分数）：锂盐与添加剂的范围以电解液总质量为分母，溶剂范围以溶剂池为分母。为便于按化学习惯阅读，每条候选另外记录锂盐（双盐时为两种盐之和）的质量摩尔浓度（mol/kg）：

$$
b = \frac{\sum_j m_{\text{锂盐}_j} / M_{\text{锂盐}_j}}{(m_\text{基准溶剂} + m_\text{共溶剂}) / 1000}
$$

分母只计溶剂，不计功能添加剂与锂盐添加剂，锂盐添加剂中的 Li⁺ 也不计入。该值写在 `pool.csv` 与 `pool_catalog.csv` 的 `salt_molality_mol_kg` 中（旧版单盐 `pool.csv` 的列名为 `<盐名>_molality_mol_kg`，导入时仍兼容），随推荐与 `observation.csv` 一起输出；它只用于展示，不进入 BO 特征，也不参与 `sample_id`。网页第一部分会按锂盐与添加剂的质量分数上下限实时显示换算的浓度区间。

双盐设计中两种盐地位相同：下限为 0 的盐是可选组分，按 `sampling.absence_probability` 的概率缺席，候选池因此同时包含双盐与单盐子配方；每个候选至少含一种盐。

## 测试

```bash
python web/test_integration.py                     # 网页 API 端到端流程
python -m unittest discover -s tests               # BO 算法脚本
(cd pool && python -m unittest test_convert_pool test_generate_pool)  # 候选池转换与撒点
```

## 目录

```text
web_for_bo/
├── web/                    网页：Flask 服务 server.py、静态页面 static/、集成测试、模块文档 docs/
├── bayesian_optimization/  网页调用的四步动作：撒点 → 候选池转换 → 实验上传/回传 → 训练推荐/可视化
├── algorithms/             BO 算法脚本（四步动作以子进程调用，也可命令行单独运行）
│   ├── qLogNEI.py          单目标推荐入口
│   ├── qLogNEHVI.py        多目标（2–4 个）推荐入口
│   ├── acquisitions.py     采集函数预设注册表（网页与命令行共用）
│   ├── bo_utils.py         GP 拟合、批量选点、数据读取与回流等共用工具
│   ├── simu_experiment.py  模拟实验回传（仅演示）
│   └── visualize_bo.py     生成 visualization_<round>.json
├── pool/                   generate_pool.py 撒点、convert_pool.py 转换、import_experiments.py 导入池外实验配方；converted/ 为七元示例候选池与初始实验
├── tests/                  algorithms/ 的单元测试
└── docs/                   设计文档、命令行说明、池数据指南、Logo 原图
```

## 文档

- [docs/ELECTROLYTE_UI_API_DESIGN.md](docs/ELECTROLYTE_UI_API_DESIGN.md)：前后端设计方案与 API 契约
- [web/README.md](web/README.md)：网页操作流程
- [bayesian_optimization/README.md](bayesian_optimization/README.md)：四步工作流
- [docs/BO_SCRIPTS.md](docs/BO_SCRIPTS.md)：命令行使用 BO 脚本
- [docs/POOL_TO_EXPERIMENT_GUIDE.md](docs/POOL_TO_EXPERIMENT_GUIDE.md)、[pool/POOL_CONVERSION_DESIGN.md](pool/POOL_CONVERSION_DESIGN.md)：候选池到实验的数据流
- [web/docs/DEPLOYMENT_AND_STORAGE.md](web/docs/DEPLOYMENT_AND_STORAGE.md)：迁移到服务器的部署与存储方案
