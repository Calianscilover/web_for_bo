# Pool CSV 转换脚本设计

## 1. 目标

`convert_pool.py` 将当前撒点结果 `pool.csv` 和
`config_snapshot.json` 转换为两类职责分离的文件：

1. `pool_features.csv`：只保留模型需要的数值 ratio 列，默认使用摩尔比，
   直接作为 `qLogNEI.py` 或 `qLogNEHVI.py` 的 `--features`。
2. `pool_catalog.csv`：保留样品 ID、RDKit 化学组 ID、溶剂、锂盐、添加剂和
   每个组分的质量比与摩尔比，用于实验调度、追溯和展示。

转换不覆盖原始 `pool.csv`。原始撒点数据作为不可变的 source，转换文件
可以随脚本重新生成。

## 2. 当前数据结构

`config_snapshot.json` 定义了 7 个组分：

| 组分 | role | 输出类别 |
| --- | --- | --- |
| LiDFOB | `primary_salt` | `lithium_salts` |
| NDFA | `balance_solvent` | `solvents` |
| TTE | `cosolvent` | `solvents` |
| FEC | `cosolvent` | `solvents` |
| LiNO3 | `salt_additive` | `lithium_salts` |
| VC | `functional_additive` | `functional_additives` |
| TMSP | `functional_additive` | `functional_additives` |

当前 `pool.csv` 有 9,656 行。七个组分的质量分数之和为 1，最大浮点误差约
`1.1e-16`；七维质量分数组合共有 9,656 个唯一值，目前没有重复候选。

## 3. ratio 定义

每个 sample 同时保存两套 ratio：

```text
compound_0   = LiDFOB
smiles_0     = RDKit 规范 LiDFOB SMILES
mass_ratio_0 = LiDFOB_mass_fraction
mole_ratio_0 = LiDFOB_mole_fraction

compound_1   = NDFA
smiles_1     = RDKit 规范 NDFA SMILES
mass_ratio_1 = NDFA_mass_fraction
mole_ratio_1 = NDFA_mole_fraction

... 直到 compound_6 = TMSP
```

使用方式：

- BO/模拟默认使用 `mole_ratio_0 ... mole_ratio_6`。
- 实验下发时按同一 `sample_id` 从 catalog 读取
  `mass_ratio_0 ... mass_ratio_6`。
- 当前配方按 5 g 总质量设计，质量比可直接与实验称量对应。
- 0 明确表示该组分不存在。

脚本保留 `--feature-basis mole|mass` 参数，只用于选择 `pool_features.csv`
的特征基准，默认为 `mole`。无论选择哪种特征，`pool_catalog.csv` 都保存
两套 ratio。

组分槽位顺序由 `config_snapshot.json` 中 `components` 的顺序唯一决定。
列名不再嵌入 LiDFOB、NDFA 等化学名称，化学身份改为每行的
`compound_i` 和 `smiles_i`。

## 4. 输出文件

### 4.1 `pool_features.csv`

仅包含数值列：

```csv
mole_ratio_0,mole_ratio_1,mole_ratio_2,mole_ratio_3,mole_ratio_4,mole_ratio_5,mole_ratio_6
0.188727508,0.811272492,0,0,0,0,0
```

要求：

- 行数和行顺序与原始 `pool.csv` 完全一致。
- 不包含 `sample_id`、`chem_group_id`、组分名称字符串或目标值。
- 可直接传入 `--features pool/pool_features.csv`。
- 每行 ratio 之和必须在容差内等于 1。

### 4.2 `pool_catalog.csv`

建议列顺序：

```text
sample_id
chem_group_id
solvents
lithium_salts
functional_additives
source_row
presence_pattern
compound_0
smiles_0
mass_ratio_0
mole_ratio_0
...
compound_6
smiles_6
mass_ratio_6
mole_ratio_6
```

样例：

```csv
sample_id,chem_group_id,solvents,lithium_salts,functional_additives,source_row,compound_0,smiles_0,mass_ratio_0,mole_ratio_0,compound_1,smiles_1,mass_ratio_1,mole_ratio_1
SMP_a13f...,CG_82c4...,NDFA,LiDFOB,,1,LiDFOB,F[B-]1(OC(C(O1)=O)=O)F.[Li+],0.19162,0.188727508,NDFA,CN(C)C(=O)C(F)(F)F,0.80838,0.811272492
```

`solvents`、`lithium_salts` 和 `functional_additives` 只列出 ratio 大于容差的活性
组分，多个组分使用分号连接。组分顺序以 `config_snapshot.json` 为准，不依赖
字典或字母表的偶然顺序。

### 4.3 `pool_manifest.json`

记录转换条件：

- 输入文件路径和 SHA-256。
- 脚本版本、RDKit 版本和生成时间。
- 特征 ratio 基准、数值容差和 ID 精度。
- 组分顺序、role 映射、规范 SMILES 和 InChIKey。
- 输出行数、唯一 sample 数和唯一 chemical group 数。

## 5. ID 设计

### 5.1 `chem_group_id`

`chem_group_id` 表示“活性化学物种集合”，不包含 ratio。因此，组分种类相同但
比例不同的配方共用一个 `chem_group_id`。

生成步骤：

1. 从 config 取出每个组分的 SMILES。
2. 用 `Chem.MolFromSmiles` 解析，解析失败立即终止。
3. 用 `Chem.MolToSmiles(..., canonical=True, isomericSmiles=True)` 得到规范 SMILES。
4. 如 RDKit 支持 InChI，同时生成 InChIKey；否则使用规范 SMILES。
5. 对当前行 ratio 大于 `presence_tolerance` 的组分，生成
   `role:InChIKey` 列表并按 config 顺序连接。
6. 对完整签名计算 SHA-256，取前 12 个十六进制字符：
   `CG_<12 hex>`。

不使用 Morgan fingerprint 直接产生 ID。Fingerprint 适合相似度搜索，但需要额外的
阈值和聚类规则，不适合作为精确、可重现的化学集合 ID。

### 5.2 `sample_id`

`sample_id` 表示具体配方，需要包含 ratio。为避免重排 CSV 后 ID 变化，不使用
简单行号作为主 ID。

生成步骤：

1. 按 config 顺序取七个组分的 ratio。
2. 使用固定小数位形成标准字符串，默认 10 位小数。
3. 将组分标识、标准质量比和标准摩尔比一起组成签名。
4. 计算 SHA-256 并取前 16 个十六进制字符：`SMP_<16 hex>`。
5. 检查 ID 冲突；如冲突则报错，不静默追加序号。

`source_row` 另行保留原始 CSV 的 1-based 数据行号，用于人工查找，但不参与
`sample_id` 生成。

## 6. 命令行接口

命令行用法：

```bash
python pool/convert_pool.py \
  --pool pool/pool.csv \
  --config pool/config_snapshot.json \
  --output-dir pool/converted \
  --feature-basis mole
```

参数：

| 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `--pool` | `pool.csv` | 原始撒点文件 |
| `--config` | `config_snapshot.json` | 组分和化学信息 |
| `--output-dir` | `converted/` | 输出目录 |
| `--feature-basis` | `mole` | `pool_features.csv` 使用 `mole` 或 `mass` |
| `--presence-tolerance` | `1e-12` | 判断组分是否存在 |
| `--sum-tolerance` | `1e-8` | ratio 加和校验容差 |
| `--id-decimals` | `10` | sample 签名数值精度 |

脚本使用临时文件写出，所有检查通过后再原子替换最终文件，避免中途失败
留下不完整输出。

## 7. 转换流程

```text
config_snapshot.json
        |
        +--> 解析 component / role / SMILES
        +--> RDKit 规范化与化学标识
        |
pool.csv
        |
        +--> 验证必需列和数值
        +--> 提取七维 ratio
        +--> 校验 ratio 和与唯一性
        +--> 生成 sample_id / chem_group_id / 分类摘要
        |
        +--> pool_features.csv
        +--> pool_catalog.csv
        +--> pool_manifest.json
```

## 8. 强制校验

转换脚本在任何输出落盘前执行以下检查：

1. config 中组分 `name` 唯一，所需字段完整。
2. 所有 SMILES 都能被 RDKit 解析。
3. `pool.csv` 存在每个组分对应的 ratio 源列。
4. ratio 必须为有限非负数，每行之和在容差内为 1。
5. ratio 为 0 时对应 `mass_g` 必须为 0；ratio 大于 0 时质量必须大于 0。
6. `sample_id` 必须唯一。
7. `pool_features.csv` 中的特征向量必须唯一，防止 BO 将两个候选映射到同一点。
8. `pool_features.csv` 与 `pool_catalog.csv` 行数、行顺序一致。
9. 输出 ratio 逐值与输入源列一致，不在转换时重归一化或静默修正。
10. 任何 ID 哈希冲突都导致转换失败。

## 9. 与 qLogNEI / qLogNEHVI 的衔接

单目标模型可以直接读取 catalog，并选择质量比作为 GP 特征：

```bash
python algorithms/qLogNEI.py \
  --mode experiment \
  --features pool/converted/pool_catalog.csv \
  --feature-basis mass \
  --output outputs/pool_qlognei
```

未指定 `--initial-observations` 时，qLogNEI 随机派发第 0 轮。生成的
`observations.csv` 保留整行 catalog 元数据，实验端只需填写 `y`。继续运行时，
程序使用 `mass_ratio_0 ... mass_ratio_6` 拟合 GP，但仍将 `sample_id`、SMILES、
质量比和摩尔比全部写入后续推荐。

如果使用已有测量，`--initial-observations` 的列必须是基础观测列加完整
catalog 列，且 `pool_id` 必须与 catalog 行号对应。

`qLogNEHVI.py` 目前仍只接受纯数值 `pool_features.csv`；若需要多目标推荐也直接
输出 catalog 元数据，需再同步相同的加载逻辑。

## 10. 实现与测试边界

建议将实现限制在 Python 标准库、NumPy 和 RDKit，不引入 pandas 依赖。拆分为：

- `load_config()`：读取并验证组分定义。
- `canonicalize_components()`：生成 RDKit 标识。
- `read_pool()`：以字典行读取原始 CSV。
- `extract_ratios()`：提取并验证数值特征。
- `make_chem_group_id()` 和 `make_sample_id()`：纯函数 ID 生成。
- `build_catalog_row()`：生成可追溯记录。
- `validate_outputs()`：输出前交叉检查。
- `write_outputs()`：原子写入三个文件。

最少测试用例：

1. 用当前 9,656 行数据完整转换。
2. 同一输入重复运行产生完全相同的 ID 和 ratio。
3. 调换输入行顺序后，同一配方的 `sample_id` 不变。
4. ratio 不归一、负数、NaN 和缺列时明确报错。
5. 无效 SMILES 时明确报错。
6. 两个 ratio 完全相同的候选导致失败并报出原始行号。
7. `core_only` 与同时含 TTE/FEC/LiNO3/VC/TMSP 的混合配方能生成正确的分类字段。
8. 程序中途失败时不留下部分输出。
