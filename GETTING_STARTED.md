# GETTING STARTED · sc-audit 使用指南

**适用版本：** 0.8.4　**定位：** 在**投稿前**回答一个问题——**"我这份单细胞病例-对照数据，能不能支撑我打算写的那些组间差异结论？"**

sc-audit 不判断"生物学是否存在"，只判断"**当前数据能否支撑所声称的结论**"。
审计结论为"不可定量 / 不可识别"时，正确表述是"**未定（indeterminate）**"，而不是"这个研究错了"。

---

## 1. 30 秒上手（不需要任何真实数据）

```bash
git clone https://github.com/karvenyao-dot/sc-audit.git && cd sc-audit
python3 -m venv .venv && . .venv/bin/activate
pip install -e .

# 生成一份自带"已知真值"的合成演示数据并跑通全流程
bash examples/run_demo.sh
# → 产物在 ./sc_audit_demo_out/（含 report.html 可读报告）
```

演示数据里 **MYL3 被刻意造成"病例组恰好全 0"**，你会直接看到 `UNQUANTIFIABLE (detection floor)` —— 这正是本工具要拦下的那类结论。

---

## 2. 用你自己的数据：三步

### 第 1 步 · 准备输入（最常见卡点）

输入是 `AnnData`（`.h5ad`）或 10x 矩阵。**`obs` 里必须有分组与供者列**，否则工具无法开始：

| 必需列 | 含义 | 例 |
|---|---|---|
| 分组列（默认 `condition`） | 病例/对照标签 | `control` / `disease` |
| 供者列（默认 `donor`） | 供者/样本 ID | `donor_A1` / `ctrl_d0` |

**审计二（可比性）另需**：≥1 个**供者级设计变量**（机构 / 取材方式 / 分选来源 / 区域 / 年龄 / 性别）。没有这类列时，审计一/三/四照常，**审计二会返回 `INDETERMINATE`**（缺注释不是"可比"，而是"无法判断"）。
**审计四另需**：细胞类型注释列（默认 `lineage_names`）。

> ⚠️ **原始 10x 目录 / raw `.h5` 不携带 `obs` 元数据**。请先用脚本写入 `condition` / `donor` / 设计变量并存成 `.h5ad`；工具在 `obs` 为空时会给出明确提示。
> 若 `var_names` 是 Ensembl ID：**命令行会**自动回退到 `var_symbols` / `gene_symbols` 列；
> ⚠️ **Python API 不会自动回退** —— 需先显式调用 `sc_audit.apply_gene_name_fallback(ad, genes)`，否则目标基因会全部落到 `NOT_IN_DATASET`（工具此时会给出明确提示）。

### 第 2 步 · 先预检（只报告，不采纳）

```bash
sc-audit check --input mydata.h5ad --condition-col condition --donor-col donor \
               --outdir check_out/
```
它告诉你：这份数据能支撑哪些审计、哪些设计变量可以拿来匹配、哪些变量有风险（例如与分组共线的取材列）。
**`check` 不写入任何配置** —— 选哪些变量进 `match_vars` 由你决定并写进 config，这一步是"事前锁定"属性的来源。

### 第 3 步 · 跑审计

```bash
sc-audit run \
  --input mydata.h5ad \
  --target-genes MYL3,PLN,TNNT2,ACTB,GAPDH \
  --condition-col condition --donor-col donor --control-label control \
  --reference reference_table.csv \
  --outdir audit_out/
```

**关键参数**

| 参数 | 作用 | 建议 |
|---|---|---|
| `--reference` | 外部参照检出率表（CSV：`gene, pct_ref`） | **强烈建议提供**；缺它则检出比 R 无法计算，判定回落 `INSUFFICIENT (no reference)` |
| `--preset cardiac` | 用内置心肌细胞 marker 面板自动门控 | 10x 原始数据常用 |
| `--matched-donors` | 手工指定匹配用的对照供者 | 默认自动检测 |
| `--ladder 'none;institution\|source;source'` | **事前锁定**审计二的匹配规则阶梯（`;` 分规则、`\|` 分变量、`none` = 不匹配锚点） | 想固定规则以供他人复核时给出 |
| `--n-perm` | 置换次数 | 默认 10000；快速试跑可用 200 |
| `--ratio-basis case\|overall\|control` | R 的分子口径 | 默认 `case` |

### 参照表（`--reference`）格式

最少两列：

```csv
gene,pct_ref
MYL3,60.0
PLN,85.0
```

> ⚠️ 上表数字只是**格式示例**：务必换成你自己**同组织、同模态、上游门控可比**的
> 数据集实测检出率（照抄示例值会让 R 与判定一起失真）。

可选字段（均由代码接受，见 `sc_audit/reference.py`）：

- **多列检出率**（如 `pct_ref_snRNA,pct_ref_scRNA`）：列名以 `_pct` 结尾，或为
  `pct` / `detection_rate`；聚合方式由 `mode` 决定（`max`（默认，最保守）或 `median`）。
  换聚合方式得到的 R 不同，务必在方法学描述里写明。
- **`modality` 列**（`scRNA` / `snRNA` / …）：声明参照模态。跨模态时 R 的绝对值
  不可解释为生物学量，工具会标 `cross-modal` 并告警。
- **旁挂 manifest**：`<reference>.manifest.json`（或同名去后缀后加 `.manifest.json`）
  可登记上游门控参数。缺它时报告会提示“参照表溯源”无法核验。

> ⚠️ 参照表应来自**同组织、同模态、且上游门控参数可比**的数据集；混用会让 R 出现系统性偏移。

---

## 3. 读结果：产物与判定

```
audit_out/
├── report.html                    ← 可读报告（自包含，含判定配色与逐表）
├── floor_audit.csv                ← 审计一：逐基因（n_pos / R / 二项 P / jackknife / q_value_BH / verdict）
├── floor_jackknife.csv            ← 审计一：逐基因的留一供者 P 范围（事件稀疏基因不出此文件）
├── composition.txt                ← 审计二（v0.2 单规则口径）
├── identifiability.txt            ← 审计二 v0.3 可识别性层（花名册 / 共同支持域 / 阶梯 / 判定）
├── identifiability_roster.csv     ← 变量花名册（哪些变量可用、为什么不可用）
├── identifiability_ladder.csv     ← 逐规则 × 逐基因的衰减比
├── identifiability_ladder_summary.csv
├── robustness.txt                 ← 审计三：深度匹配与分层
├── identity.txt                   ← 审计四：细胞身份与门控共线
└── run_metadata.json              ← 运行参数、输入/参照表 SHA256、配置来源、版本
```

### 审计一（检出地板）

| 判定 | 含义 | 你怎么做 |
|---|---|---|
| `UNQUANTIFIABLE (detection floor)` | 检出率 < 5% | **不要**报告该基因的组间比较 |
| `DROPOUT-LIMITED` | 相对参照严重 dropout（R < 0.3） | 组间差异不可解释为生物学 |
| `QUANTIFIABLE` | 检出率 ≥ 20% 且 R ≥ 0.6 | 可进入组间比较 |
| `INSUFFICIENT` | 无参照 / 证据不足 | 只在补充材料报告 |
| `stratum = event_sparse`（独立轴） | 阳性细胞数 < 20 | 该基因**不出**推断统计量（这是事前锁定规则） |

### 审计二 v0.3（可比性 / 可识别性）

| 判定 | 含义 |
|---|---|
| `SUPPORT-BELOW-FLOOR` | 共同支持域太小且病例组在设计变量上恒定 ⇒ **不存在**可调整的组间估计；衰减比被**抑制**不报 |
| `THRESHOLD-SPANNING` | 衰减比跨预设阈值 ⇒ 它不是单一估计量 |
| `COMPARABLE` / `NOT SEPARABLE` | 所有规则同侧 ⇒ 给方向判定，**必附**"阈值落在观测区间之外、一致性检查无信息"的强制警示 |
| `INDETERMINATE` | 无组别结构 / 病例组跨多个设计剖面 / 无目标基因 / **无可用设计变量**（注释局限） |

> ⚠️ 读 `identifiability.txt` 时**必须连 notes 一起读**：单侧判定与"注释局限"的结论都写在 notes 里。

---

## 4. 事前锁定（为什么值得做）

阈值全部写在 `sc_audit/configs/cardiac.yaml`（与 `AuditConfig` 默认值逐字段一致）。用法：

```bash
cp sc_audit/configs/cardiac.yaml my.yaml     # 改阈值/变量声明
sc-audit run --input mydata.h5ad --config my.yaml ...
```

- 改了什么、用的哪份配置、输入与参照表的 SHA256 —— 全部记进 `run_metadata.json`，他人可复跑。
- **不要**在看到结果之后再改阈值；`--ladder` 显式给出即视为运行前锁定。
- 变量角色（结局 / 治疗 / 忽略）可用 `role_overrides` 显式声明，优先于列名启发式。

---

## 5. 命令速查

```bash
sc-audit check            --input x.h5ad                     # 预检（只报告）
sc-audit run              --input x.h5ad --target-genes A,B  # 四项审计 + v0.3 层 + HTML
sc-audit floor            --input x.h5ad --target-genes A,B  # 仅审计一
sc-audit composition      --input x.h5ad                     # 仅审计二（v0.2 单规则）
sc-audit identifiability  --input x.h5ad --ladder 'none;source'  # 仅审计二 v0.3 层
sc-audit robust           --input x.h5ad --housekeeping ACTB,GAPDH
sc-audit identity         --input x.h5ad --celltype-col lineage_names --target-lineage Cardiomyocytes
sc-audit compare          --input x.h5ad --reference ref.csv # 跨数据集检出比 R
```

---

## 6. Python API

```python
import sc_audit

ad = sc_audit.load_data("mydata.h5ad")
sc_audit.compute_qc_metrics(ad)

floor = sc_audit.audit_detection_floor(
    ad, genes=["MYL3", "PLN"], condition_col="condition",
    donor_col="donor", control_label="control")

ident = sc_audit.audit_identifiability(
    ad, genes=["MYL3", "PLN"], condition_col="condition",
    donor_col="donor", control_label="control")

print(floor.table[["gene", "n_pos", "ratio_R", "verdict"]])
print(ident.verdict)
print(ident.roster)          # 变量花名册：哪些可用、为什么不可用
```

阈值推导也可在代码里复核（预注册容差口径；示例 n_case = 200）：

```python
from sc_audit import derive_event_threshold, event_threshold_table
event_threshold_table(200)        # {0.50: 17, 0.434: 21, 0.40: 24, 0.33: 33, 0.25: 50}
derive_event_threshold(200, 0.434)  # 21 —— x_min 随 n_case 变化
```

---

## 7. 常见问题

| 现象 | 原因与处置 |
|---|---|
| 报错 `obs 中缺少列 'condition'` | 你的输入没带 `obs`（常见于原始 10x）。先写入 `condition`/`donor` 再运行 |
| 审计二一律 `INDETERMINATE` | 没有**非供者级唯一**的设计变量（公开数据集里 institution / source 这类设计字段往往缺失）。这是**注释局限**，不是数据可比 |
| 全部基因 `INSUFFICIENT (no reference)` | 没给 `--reference`。提供同组织参照检出率表即可 |
| `NOT_IN_DATASET` | 目标基因不在 `var_names`。若用 Ensembl ID：**CLI 会自动回退到 `gene_symbols`**；**API 需显式调用 `apply_gene_name_fallback()`**。都做了仍失败 → 该基因确实缺失 |
| **第一条命令卡住 30–60 秒** | 正常：首次调用要 import scanpy/numba 并生成编译缓存（约 30–60 s），**不是死机**。同一次进程里后续命令是 1–2 秒；脚本里批量跑时把多次调用放同一进程最省 |
| 运行 1 分钟以上 | 成本集中在审计一的 **jackknife**，随**供者数**增长（与细胞数关系小）。可分开跑：`identifiability`/`composition` 都是秒级；或 `--n-perm 200` 试跑 |
| 命令直接抛 Python 异常 | 输入类错误应打印 `[error] …` 并 `exit 2`；若看到裸 traceback，说明是未预料的输入形态，请带上完整输出反馈（设 `SC_AUDIT_DEBUG=1` 可看到原始堆栈） |
| 结果里的比率看不懂 | 见 `identifiability.txt`：`ratio` 是**效应量比值**（分子分母均取绝对值，对方向不敏感）；方向信息在 `full_pp` 与 `matched_ctrl_pct` 两列 |
| `--target-genes` 给空值 | `run` / `floor` / `robust` / `identifiability` / `compare` 会**拒绝**（exit 2，必须给至少一个目标基因）；`composition` / `identity` 接受空值（这两项不按基因划分子表） |

---

## 8. 自检与复现

```bash
python3 -m pytest tests/ -q          # 409 项（collected）＝ 405 passed + 4 skipped
                                     #  4 skipped = 2 项外部数据集门控 + 2 项内部文档门控
```

纪律：阈值事前锁定；返回表**列集合稳定**（退化路径也补齐核心列）；文档与代码**逐行一致**（维护侧的
由开发仓库自带的文档–代码一致性检查脚本完成（该脚本不随发布包出厂）。

---

## 9. 引用与许可

MIT。若用于研究，请引用本仓库（见 `CITATION.cff`）与 `docs/methods.md` 中的方法学依据。
