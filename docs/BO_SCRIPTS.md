# 电解液贝叶斯优化：训练与候选推荐

当前 `qLogNEI.py` 和 `qLogNEHVI.py` 完成实验数据读取、输入归一化、GP 拟合、候选预测、批量推荐和实验结果回流。以下命令均在仓库根目录、已激活 `botorch` 环境下运行。

## 输入

- 实测结果：`pool/converted/experiment.csv`。单目标默认读取 `Conductivity`；双目标默认读取 `logCE` 和 `Conductivity`。
- 候选池：`pool/converted/pool_catalog.csv`。候选输入默认读取 `mass_ratio_0`、`mass_ratio_1` 等连续编号的数值列；`--feature-basis mole` 改为摩尔比例列。
- 两表通过 `sample_id` 对齐，行顺序可以不同。实验表中的输入特征必须与候选池一致。目标值只从实验表读取；候选池即使含有目标列，也不会作为测量值使用。
- 目标为空的实验行不参与模型训练，仍保留在候选池。双目标要求一行的两个目标都齐全。训练至少需要两行完整实测数据。
- 对整个候选池的所选特征做范围归一化；恒定特征被移除。训练记录会列出实际使用的特征及维度。

输入特征尚未确定时，可用 `--feature-columns` 显式指定两表共有的数值列。例如：

```bash
python algorithms/qLogNEI.py \
  --feature-columns mass_ratio_0 mass_ratio_1 mass_ratio_2
```

## 训练与推荐

```bash
python algorithms/qLogNEI.py
python algorithms/qLogNEHVI.py
```

可通过 `--experiment`、`--pool`、`--output` 指定输入和输出路径。单目标可用 `--target` 修改性能列；双目标可用 `--targets A B` 和 `--directions max min` 修改目标列及优化方向。两者保留原有的 `--kernel`、`--matern-nu`、`--ard`、`--lengthscale-init`、`--noise-std` 和 `--fit-maxiter` GP 参数。`--batch-size` 设置每次推荐配方数（默认 3）；`--mc-samples`、`--pool-batch-size` 和 `--seed` 控制采集函数计算。双目标的 `--ref-point A B` 使用原始目标单位；省略时从首次已测目标自动计算，后续轮次复用该参考点。

单目标使用 qLogNEI，双目标使用 qLogNEHVI。采集函数在未测候选中做离散批量选择，历史已测输入作为 baseline。拟合后的模型还会分块预测所有未测候选，预测均值会换回原始目标方向和单位。

默认输出分别写入 `outputs/single_training/` 与 `outputs/multi_training/`：

- `training_summary.json`：实际特征、输入维度、目标名、已测及剩余候选数量。
- `model.pt`：模型参数状态及训练用的 `train_X`、原始单位的 `train_Y`。
- `candidate_predictions.csv`：首轮全部未测候选的 `sample_id`、预测均值和潜在函数后验标准差；后续轮次保存为 `candidate_predictions_2.csv` 等。
- `recommendation_1.csv`：首批推荐的原始候选池配方字段及空白目标列，供实验回填；后续推荐依次保存为 `recommendation_2.csv` 等。
- `observation.csv`：累计实验结果，首次从 `experiment.csv` 建立，后续按 `sample_id` 合并每批回传值。部分目标已测的行会保留已有数值。

若推荐文件仍有待回填目标值，重复运行入口会暂停。填完 `recommendation_1.csv` 的目标列后，用相同命令再次运行，即生成 `recommendation_2.csv`、`candidate_predictions_2.csv`；后续依次编号。完整候选池耗尽时停止推荐。请保持实验与候选池文件、目标设置及输出目录一致。已有的 `recommendations.csv`、`recommendations_2.csv` 等旧文件仍可继续回流。

模拟回传可用单独脚本。它只填本轮推荐文件中尚未获得的目标值，随机值从相应目标的已有实测范围中抽取，随机数由 `--seed` 与配方 `sample_id` 共同决定：不同配方得到不同的值，同一配方重复模拟结果不变（各轮不会重复同一组数值）；若该配方某目标早已实测，则保留原实测值。

```bash
python algorithms/simu_experiment.py --output outputs/single_training --seed 2026
python algorithms/qLogNEI.py
```

双目标只需将输出目录和入口改为 `multi_training`、`qLogNEHVI.py`。也可以用 `--recommendations` 指定某个推荐 CSV。模拟值仅用于验证数据回流，不代表真实实验结果。

当前示例候选池有 9,656 条配方；`experiment.csv` 含 9 条，其中 8 条两个目标都完整。初次拟合使用这 8 条，其余 9,648 条为待评估候选。模型目标只从实验表或已合并的 `observation.csv` 读取，不从候选池读取。

验证：

```bash
python -m unittest discover -s tests
```
