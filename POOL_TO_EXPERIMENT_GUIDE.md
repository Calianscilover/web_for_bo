# 从候选池到实验推荐：简要流程

## 1. 目标

本流程将预先定义的电解液配方空间转换为可追溯的候选池，并使用
qLogNEI 从候选池中逐轮选择最值得实验的配方。

```text
配方空间与约束
      ↓ 撒点
原始 pool.csv
      ↓ 格式转换与校验
pool_catalog.csv
      ↓ qLogNEI 主动学习
实验配方 → 实验测试 → 回填结果 → 下一轮推荐
```

## 2. 第一步：生成原始候选池

上游撒点程序根据组分范围、配方约束和采样规则生成：

- `pool/pool.csv`：所有合法候选配方，包含各组分的质量分数、摩尔分数及称量质量。
- `pool/config_snapshot.json`：组分名称、角色和 SMILES 等配置快照。

当前候选池包含 **9,656 个配方**，覆盖 **32 种组分组合**。`pool.csv` 是原始数据，
后续转换不会覆盖它。

## 3. 第二步：转换为标准格式

在仓库根目录运行：

```bash
python pool/convert_pool.py \
  --pool pool/pool.csv \
  --config pool/config_snapshot.json \
  --output-dir pool/converted \
  --feature-basis mass
```

转换后生成三个文件：

| 文件 | 用途 |
| --- | --- |
| `pool_catalog.csv` | 完整配方目录，包含 sample ID、化学组 ID、名称、SMILES、质量比和摩尔比 |
| `pool_features.csv` | 仅含数值 ratio，供只接受纯数值输入的模型使用 |
| `pool_manifest.json` | 记录输入文件、转换设置、组分顺序、版本和校验信息 |

转换脚本同时完成 RDKit SMILES 规范化、稳定 ID 生成、比例和校验、重复配方检查。
其中 `sample_id` 对应一个具体配方，`chem_group_id` 对应相同的活性组分集合。

## 4. 第三步：生成第一轮实验配方

单目标推荐直接使用完整的 `pool_catalog.csv`：

```bash
python qLogNEI.py \
  --mode experiment \
  --features pool/converted/pool_catalog.csv \
  --feature-basis mass \
  --output outputs/pool_qlognei
```

首次运行时，程序随机选择 4 个初始配方，并生成：

```text
outputs/pool_qlognei/experiment/
├── observations.csv
├── recommendations_000.csv
└── config.json
```

- `recommendations_000.csv`：第一轮实验清单快照。
- `observations.csv`：整个项目持续更新的实验记录表。
- `config.json`：本次运行配置，用于保证后续运行条件一致。

模型只使用 `mass_ratio_0 ... mass_ratio_6` 作为输入特征；配方名称、SMILES、
摩尔比和各类 ID 会完整保留，方便实验执行和结果追溯。

## 5. 实验闭环

1. 实验端根据 `observations.csv` 中的质量比配制并测试配方。
2. 只在待测行的 `y` 列填写实验结果，不修改其他列。
3. 再次执行完全相同的 qLogNEI 命令。
4. 程序读取全部已测数据，拟合高斯过程模型，并推荐下一批 3 个配方。
5. 重复以上步骤，直到达到实验预算或候选池耗尽。

如果某一批仍有未填写的 `y`，程序会暂停等待，不会提前生成下一批推荐。

## 6. 会议要点

- 所有推荐均来自预先校验过的候选池，算法不会生成池外配方。
- 模型学习使用质量分数，实验清单同时保留质量比和摩尔比。
- `sample_id` 贯穿建模与实验，可避免配方错配。
- 每轮实验数据都会累积进入模型，实现“推荐—实验—更新—再推荐”的主动学习闭环。
- `replay` 模式用于历史数据验证；`experiment` 模式用于真实实验推荐。

