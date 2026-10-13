# sc-audit

**Detection-floor and control-composition audit for case-control single-cell studies.**

一个可安装、可运行的命令行工具，用于在**发表前**检查单细胞病例-对照研究的信号可信度。

> 解决的问题：在人心脏等**样本稀缺**的组织中，低表达基因常因技术原因（检出地板/稀疏）
> 而在病例组"恰好"全部未检出，形成看似显著的组间差异；同时，病例组与对照组在
> 机构、取材方式、分选来源上常与生物学条件完全共线。本工具把这两类问题变成
> **可执行的、事前锁定的判据**。

---

## 安装

**一键安装（推荐）**——发布树附带安装脚本，一条命令建环境 + 装包 + 验证：

```bash
bash install.sh                 # macOS / Linux（联网）
bash install.sh --offline      # 用同目录 wheels/ 离线安装
# macOS 也可在 Finder 里双击 install-macos.command
```
```powershell
.\install.ps1                   # Windows PowerShell
.\install.ps1 -Offline
```
```bash
python3 install.py              # 跨平台通用（纯标准库）
```
默认装到 `~/.sc-audit-venv`（可用 `SC_AUDIT_VENV` 改）；完成后用 `~/.sc-audit-venv/bin/sc-audit`
（Windows 在 `Scripts\sc-audit.exe`）。

**手动安装**

```bash
git clone https://github.com/karvenyao-dot/sc-audit.git          # 或直接使用本地目录
cd sc-audit
pip install -e .          # 可编辑安装
# 或
pip install .             # 常规安装
```

依赖：`numpy pandas scipy anndata scanpy matplotlib pyyaml`（Python ≥ 3.9）

验证安装：
```bash
sc-audit --version
python -c "import sc_audit; print(sc_audit.__version__)"
```

---

## 快速开始

```bash
# 四项审计一次跑完，输出 HTML 报告 + CSV/TXT
sc-audit run \
  --input mydata.h5ad \
  --target-genes MYL3,KCNJ2,KCNQ1,KCNH2,SCN5A \
  --reference reference_table.csv \
  --outdir audit_out/
```

产出：
```
audit_out/
├── report.html        人类可读报告（自包含；覆盖审计一/二/三 + 可识别性层；**不含审计四**——审计四的产物是 identity.txt）
├── floor_audit.csv    审计一：逐基因结果
├── floor_jackknife.csv 审计一：逐基因的留一供者 P 范围
├── composition.txt    审计二：交叉表 + 匹配细节 + 效应量衰减
├── identifiability.txt 审计二 v0.3：花名册 / 共同支持域 / 阶梯 / 判定
├── identifiability_roster.csv          变量花名册
├── identifiability_ladder.csv          逐规则 × 逐基因的衰减比
├── identifiability_ladder_summary.csv  逐规则汇总
├── robustness.txt     审计三：深度匹配 + 分层
├── identity.txt       审计四：细胞身份与捕获完整性
└── run_metadata.json  运行参数（可复现）
```

**单跑某一项：**
```bash
sc-audit floor       --input mydata.h5ad --target-genes MYL3,KCNJ2
sc-audit composition --input mydata.h5ad
sc-audit robust      --input mydata.h5ad --housekeeping ACTB,GAPDH
sc-audit identity    --input mydata.h5ad --celltype-col lineage_names \
                     --target-lineage Cardiomyocytes --modality scRNA
sc-audit compare     --input mydata.h5ad --reference ref.csv
```

**`check` 的退出码与 `--ladder`：**
- `check` 退出码：`0` 全部可行 / `1` 有审计不可行（degraded 或 no）/ `2` 存在阻断项（缺分组或供者列、计数不可用）。它**只报告不采纳**，不会写入配置。
- 事前锁定的匹配阶梯在 `identifiability` 与 `run` 上都用 `--ladder`（`'none;institution|source;source'`，`'none'` = 不匹配的锚点）；`identifiability` **不接受** `--matched-donors`（它的匹配集合由阶梯决定）。

**读取 10x 输出 + 自动门控心肌细胞：**
```bash
sc-audit run --input ./filtered_feature_bc_matrix/ \
             --preset cardiac \
             --condition-col group --donor-col sample \
             --target-genes MYL3,KCNJ2
```

---

## 四项审计

### ① 检出地板审计 `sc-audit floor`

判断目标基因是否**具备定量能力**，独立于正在被声称的组间比较。

| 指标 | 含义 |
|---|---|
| `pct_study` | 目标细胞中该基因的检出率 |
| `n_pos` | **绝对阳性细胞数（事件数）** |
| `pct_ref` / `ratio_R` | 外部参照检出率与归一化检出比 |
| `binom_P_case0` | 病例组观察到 0 个阳性的二项概率 |
| `stratum` | 事件数分层（sparse / limited / sufficient） |
| `jackknife_*` | 留一对照供者后的 P 值范围与效应量范围 |

**事前锁定阈值：**

| 判据 | 阈值 | 处置 |
|---|---|---|
| 检出率 | `< 5%` | 不报告组间比较 |
| 检出比 R | `< 0.3` | dropout-limited，差异不可解释为生物学 |
| 检出比 R | `≥ 0.6` 且检出率 `≥ 20%` | 可定量 |
| **事件数** | **`n_pos < 20`** | **仅描述性统计，不做显著性检验** |
| 事件数 | `20 ≤ n_pos < 50` | 比较需附精确 CI 与 jackknife 范围 |

> **为什么以事件数为准？** 1.2% 的检出率在 800 个细胞中只对应约 10 个阳性细胞。
> 此时任何 P 值都反映抽样噪声。该条件**事前锁定，不因结果而调整**。

### ② 对照构成审计 `sc-audit composition`

判断病例组与对照组是否**可比**。

- 条件 × 机构 / 取材方式 / 分选来源 / 区域 的交叉表
- **供者级 TVD**：条件 × 设计变量的分布距离（0 = 分布相同，1 = 完全分离）；
  判定用供者级（每供者一票），细胞级仅作参考
- **同队列匹配**：自动找出在所有设计变量上与病例组一致的对照供者，
  在每个供者上报告"一致变量数"（透明可核查）
- **效应量衰减**：匹配前后效应量之比（衰减 = 混杂份额的证据）

判定：`COMPARABLE` / `NOT SEPARABLE` / `INDETERMINATE`

### ② b 可识别性层 `sc-audit identifiability`　★v0.3 新增

回答一个比"是否可比"更前置的问题：**这个设计里到底存不存在可调整的组间估计？**

- **变量花名册**：obs 每个非 QC 列必须落桶（未归类就报出，不静默跳过）
  · 病例组内恒定（donor 级熛 ≈ 0）→ `design-defining(scope)`：不可匹配，只能作适用范围声明
  · 病例组内有变异 → `adjustable(confounder)`：真正的可调整混杂
- **共同支持域**（可识别性 / positivity，先于匹配）：病例组设计剖面内的对照供者与细胞数
- **匹配规则阶梯**（事前锁定 R0–R5）：逐规则的 full / matched 效应量与衰减比
- **一致性检查**：单侧（全部规则落在阈值同侧）→ 给方向判定并**强制附警示**
  （阈值在观测区间之外时，一致性检查不携带信息）；两档并存 → `THRESHOLD-SPANNING`
- **低于支持域下限 → 抑制衰减比**：不出任何量化证据
- **小样本警示**：输出 `chance-constancy risk`（某变量在 n 个病例供者中偶然恒定的概率）

判定：`SUPPORT-BELOW-FLOOR` / `THRESHOLD-SPANNING` / `NOT SEPARABLE` / `COMPARABLE` / `INDETERMINATE`

```bash
sc-audit identifiability --input mydata.h5ad \
  --target-genes MYL3,PLN,TNNT2 --condition-col condition --donor-col donor
# 事前锁定自己的阶梯（R0 锚点写 none；规则间用 ';'，变量间用 '|'）：
sc-audit identifiability --input mydata.h5ad --ladder 'none;institution|source;source'
```

产物：`identifiability.txt` + `identifiability_roster.csv` +
`identifiability_ladder.csv` + `identifiability_ladder_summary.csv`。

### ④ 细胞身份与捕获完整性 `sc-audit identity`　★v0.2 新增

在**细胞类型内部**检出构成混杂，并把**门控规则本身当作待审计对象**。
由病例–对照单细胞数据的细胞身份取证实践逼出——前三项审计均无法发现它检出的问题。

八项检查：

| # | 检查 | 判据 |
|---|---|---|
| 1 | 自家门控通过率（按注释谱系） | 目标谱系须远高于其他谱系（默认 ≥5×）|
| 2 | 中位 genes/cell 在同数据集谱系内的排名 | 目标谱系不应最低 |
| 3 | 核滞留转录本（MALAT1）排名 | 偏低 → 核含量低/捕获不完整 |
| 4 | 长/低拷贝基因在**复杂度匹配**下的检出 | 分离"捕获受限"与"ambient 污染" |
| 5 | ★ **亚群不平衡**（来源 × 条件，细胞类型内部） | 某亚群只在单一臂出现 |
| 6 | ★ **门控规则的条件共线性** | 阴性标志检出率随臂变化 |
| 7 | 模态匹配（参照表是否有同模态列） | 无 → 所有 R 标 `cross-modal` |
| 8 | 取材方式标志按臂分层 | 内皮/间质 ambient |

**5 与 6 是承重项。** 前者是"审计二下沉到细胞类型内部"，后者把 QC 规则本身
变成待审计对象。判定：`PASS` / `CAPTURE-LIMITED` / `COMPOSITION-CONFOUNDED (within cell type)`。

```bash
sc-audit identity --input mydata.h5ad \
  --celltype-col lineage_names --source-col source \
  --target-lineage Cardiomyocytes --modality scRNA \
  --reference-modalities snRNA --outdir out
```

`--modality` 与 `--reference-modalities` 用于检查 7；跨模态时 R 的**绝对值不可解释**，
管家基因护栏也**只在同模态内有效**。

### ③ 稳健信号对照 `sc-audit robust`

分离真信号与深度假象，并提供**正面对照**（避免"全盘否定"）。

- **深度匹配**：每细胞基因数 / 总计数 / 线粒体比例的中位数对比（各自独立标度）
- **深度分层检出**：按每细胞基因数分箱，箱内比较两组检出率
- **管家基因内部对照**：若管家基因检出率不低于参照，则排除"细胞类型判错"与
  "文库整体劣化"两个替代解释

判定：`DEPTH-COMPARABLE` / `PARTIAL` / `INSUFFICIENT`

> 注：线粒体比例属**细胞状态**指标，不参与深度可比性判定，仅作独立观察报告。

---

## 设计声明与能力闸门（v0.4.0）

工具**默认假设**输入是「二元、横断面、分析单位 = 供者」的病例–对照设计（支持集 `S0`）。
若你的数据是**纵向 / 治疗臂 / 多状态**，请**显式声明**——程序**不推断**，
声明超出 `S0` 时会**拒绝输出判定**（而不是静默算错）：

```bash
sc-audit run --input x.h5ad --target-genes MYL3,PLN \
  --comparison-type two_arm --unit donor   # 显式声明 S0：与不声明等价
```

| 参数 | 取值 | 默认 |
|---|---|---|
| `--comparison-type` | `two_arm` / `paired` / `time_course` / `multi_arm` / `multi_state` | `two_arm` |
| `--unit` | `donor` / `donor_timepoint` / `sample` / `cell` | `donor` |
| `--contrast-a` / `--contrast-b` | 比较的 A/B 组标签 | 由 condition/control 推 |
| `--time-col` / `--ref-time` | 纵向时间列 / 参考时间点 | 无 |
| `--index-time` / `--index-event` | 索引时间 / 事件（**仅记录**） | 无 |

**行为**：

- 声明 ∈ `S0` → 正常执行；默认路径判定与 v0.3.x **逐字节一致**。
- 声明合法但 ∉ `S0` → `DESIGN-NOT-SUPPORTED`，**退出码 3**，不产出任何判定产物。
- 声明取值写错（不在封闭词表内）→ 错误，**退出码 2**。
- 数据形似更复杂设计（如 condition 水平数 > 2）→ 输出 `[design-hint]` 提示，**不改判定**。

> `DESIGN-NOT-SUPPORTED`（设计超出能力、未做审计）与 `INDETERMINATE`（审计已做但无法给结论）
> 是**两回事**。详见 `docs/methods.md` §8。

---

## 参照表契约：跨 modality 与合并（v0.5.0）

不同 modality（**scRNA** vs **snRNA**）的检出率**不可比**，直接算 R 会系统性偏移。
因此参照表相关有两条**硬约束**：

- **R2 跨 modality 默认禁止**。若被审数据声明 `--modality scRNA`、而参照声明 `snRNA`（或在
  参照表 `modality` 列 / manifest 中声明），工具会**报错（退出码 2）**；确需跨模态时用
  `--allow-cross-modality` 显式放行——放行后受影响基因的 R 标 `cross_modality=True`，
  且其**依赖 R 的判定会被降级**为 `INSUFFICIENT (cross-modal reference)`（v0.7.1 起），
  不再输出 `QUANTIFIABLE`；同时新增 `verdict_contract_ok=False` 列。
- **R1 五元组合并守卫**。把多个参照列并成同一参照时，
  `(modality, intron_counting, gate_spec_hash, upstream_qc_hash, cell_source)`
  必须一致（写在参照表的 `*.manifest.json` 的 `column_contract` 里），否则**拒绝合并**
  并列出差异字段。未提供逐列五元组时，工具**显式登记“R1 未强制”**，不静默跳过。

```bash
# 跨模态会被拒绝：
sc-audit floor --input x.h5ad --reference ref.csv --modality scRNA
# [error] 数据模态 'scRNA' 不在参照模态 ['snRNA'] 中 … 退出码 2

# 确认要跨模态时：
sc-audit floor --input x.h5ad --reference ref.csv --modality scRNA --allow-cross-modality
```

---

## v0.6.0 方法学修订（A/B/C）

基于逐层面文献（`文献_逐层面映射_2026-10-12.md`）与代码逐处对照，v0.6.0 落地：

- **A1 参照侧深度可比性**：参照表若提供 `mean_umi`（或 manifest），会与 case 的中位
  `total_counts` 比对；相对差 > `reference_depth_delta_max_pct`（默认 50%）→ 默认**拒绝**
  （`--allow-depth-mismatch` 放行）；`--depth-matched-sensitivity` 额外输出把 case 下采样到
  参照深度后的 R 与判定列（`case_pct_depth_matched` / `verdict_depth_matched`）。
- **A2 positivity / 协变量重叠**：新增**供者级协变量重叠**诊断；`min(overlap) < positivity_overlap_min`
  → `NOT SEPARABLE`。细胞计数现在只是**容量前提**，不再充当可识别性判据。
- **B1 潜在批次层**（`--with-latent-batch`）：供者级 pseudobulk → HVG → PCA → k-means，
  检查条件与供者表达簇的共线（ARI）。**报警层，默认关闭，不当匹配变量**。
- **B2 管家护栏一致性**：新增各管家 R 的离散度（CV）告警与“高基线组织适用域”声明。
- **B3 分层口径敏感性**：额外用 `total_counts` 分层与 `n_genes` 对照，报告持续性是否口径依赖。
- **C1–C5 透明项**：ambient RNA 局限、双细胞率、注释方法依赖、`binom_P_case0` 的细胞级（伪重复）口径、小面板 FDR 功效——均显式写入输出。

详见 `方法学修订提案_v0.6.0_2026-10-12.md`。

---

## v0.7.1 修复（据外部独立评估，经交叉验证后采信）

- **C1 condition 的“标签副本”列不再被当设计变量**：与 condition **一对一重编码**
  （两臂各自恒定且水平 ≥ 2）的列（如 `condition_original`/`Status`）归入新桶
  `excluded:redundant-with-condition`，从 design-defining / 共同支持域 / 候选匹配中剔除。
  修复前它会把支持域打成 **0**（审计二b 对多个真实数据集恒判 `SUPPORT-BELOW-FLOOR`，
  属**假拒绝**）；修复后若支持域真为 0，必有其因，不再伪装成“数据不可识别”。
- **C2 放行路径现在**限制判定**（不再只打标记）：`--allow-cross-modality` /
  `--allow-depth-mismatch` 放行后，受影响基因的依赖 R 的判定降级为
  `INSUFFICIENT (cross-modal reference)` / `INSUFFICIENT (depth mismatch)`，并新增
  `verdict_contract_ok` 列（仅在受影响时新增，未受影响运行输出零改动）。
- **C12 `check` backed 输入不再崩溃**：`_depth_bin_capacity` 改走 backed-safe 的
  `io.x_matrix()`（原先直接 `sp.csr_matrix(ad.X)`，`backed='r'` 时抛 `dtype object`）。
- **C3 审计四可评估性改为动态清单**：删除硬编码 “checks 4/5/6/7/8 FULLY EVALUATED”，
  改为按**本次运行实际结果**生成 `ASSESSABILITY (this run)`（`not evaluable ≠ passed`）；
  检查 4 的静默跳过改为显式登记。

裁决与 C4/C5/C6/C7/C9–C11（登记为 0.7.2 候选）全记录见项目内审材料（不随公开包发布）。

---

## `design` 子命令（v0.8.0，**解耦版 AI 设计解析**）

**不嵌入任何 LLM、不联网**；AI 在工具**之外**（任一大模型/agent 生成 proposal）。工具只定义契约：

```bash
# 1) 导出 inventory（离线；--geotext 指向 GEO 自由文本）
sc-audit design --action export --input data.h5ad \
  --condition-col condition --donor-col donor --control-label control \
  --label my_dataset --geotext geo.txt --outdir out/

# 2) 把 inventory.json 交给 AI 生成 proposal.json（本项目 agent / 以后换的模型 / 手工均可）

# 3) 读回 → 确定性复验 + 人工确认单 + 待确认的配置声明
sc-audit design --action apply --inventory out/inventory.json \
  --proposal proposal.json --outdir out/
```

产物：`inventory.json`、`design_reverify.csv`（逐列复验）、`design_confirm_sheet.md`（风险分层：① 必须逐条 / ② 批量 / ③ 自动+抽检；opt-out）、`design_config.json`（**待人工确认、不自动写入**）、`design_registry.json`（列注册表，重跑不重复问）。

---

## v0.7.2 修复（据同一轮外部独立评估，经交叉验证后采信·第二批）

- **C4 审计四：门控通过率新增**绝对量**与**按臂差**判据**。旧实现只比“目标/其他谱系中位数”
  的相对 margin，单谱系时该值不可用 → 门控项永不触发；通过率 0.3%、两臂差 3.4× 的门控
  规则仍会判 PASS。新增两个事前锁定阈值：`identity_gate_pass_floor_pct`（通过率绝对
  下限；低于此 → `GATE-NOT-USABLE`）与 `identity_gate_arm_delta_max`（臂间比上限；超过
  → 门控与分组共线，计入 `COMPOSITION-CONFOUNDED`）。
- **C5 审计二：新增 `NOT ASSESSED` 档**。无任何 ≥ 2 水平的设设变量时，可比性**根本没被
  评估**，不再输出褒义的 `COMPARABLE … (TVD not assessed)`。
- **C6 `check`：reference 未传入的措辞**。“参照表只覆盖 0/N” → “覆盖率 0 系**未传入**
  所致（不是参照表缺少这些基因）；请用 --reference 传入”。
- **C7 审计二：小匹配集不再用衰减比下判**。匹配细胞数 < `matched_min_cells_warn` 时，
  该比值由过少供者决定（如 1 供者/31 细胞），从判据中移除并显式说明（只留 TVD）——
  与审计二b 对同一数字的抑制一致，消除“同一量两套认识论权重”。
- **C8 审计四：rank 分母**用注释谱系总数（原无 `n_genes` 时显示 `rank 1/0`）。
- **C9 审计一：表内加 `verdict_basis_group` / `verdict_floor_pct` 两列**，注明 verdict 的
  基准组与阈值（`stratum` 按总阳性数分档，verdict 按 basis 组检出率判定，两者基准不同）。
- **C11 审计一：R2 未强制时显式登记**（与 R1 对称）：参照表未声明模态时发出告警并写入
  `summary['r2_enforcement_note']`。
- **C10 文档门**：`check_doc_code_sync.py` 对**包外**辅助文件（`../tools/`）在副本/发布树中
  缺失时降级为 SKIP，不再令“文档=代码”门在用户侧永远为红。

---

## Python API

```python
import sc_audit

ad = sc_audit.load_data("mydata.h5ad")
sc_audit.compute_qc_metrics(ad)

floor = sc_audit.audit_detection_floor(
    ad, genes=["MYL3", "KCNJ2"],
    condition_col="condition", donor_col="donor", control_label="control",
    reference={"MYL3": 60.0},
)
print(floor.to_text())
print(floor.unquantifiable)        # ['MYL3', 'KCNJ2']

comp = sc_audit.audit_control_composition(
    ad, genes=["TNNT2"], condition_col="condition",
    donor_col="donor", control_label="control")

rob = sc_audit.audit_robust_signal(
    ad, genes=["TNNT2"], condition_col="condition",
    donor_col="donor", control_label="control")

from sc_audit.report import render_report
render_report(floor=floor, composition=comp, robust=rob, outdir="out/")
```

---

## 配置

预设：`configs/cardiac.yaml`（阈值与 `docs/methods.md` 的预注册口径一致）。

自定义：复制该文件并修改，然后 `--config my_config.yaml`。

```python
from sc_audit import load_config, AuditConfig
cfg = load_config("cardiac")      # 或 "path/to/my.yaml"
cfg = AuditConfig(event_sparse=15, n_perm=20000)   # 代码内覆盖
```

---

## 解释指南

| 判定 | 含义 | 建议 |
|---|---|---|
| `UNQUANTIFIABLE` | 检出低于地板或事件数过少 | **不要**报告该基因的组间比较 |
| `DROPOUT-LIMITED` | 相对参照严重 dropout | 差异不可解释为生物学 |
| `NOT SEPARABLE` | 条件与设计变量共线 | 生物学信号与技术来源不可分离 |
| `COMPARABLE` | 在可用设计变量上可比 | 仍需检查样本量与统计方法 |
| `DEPTH-COMPARABLE` | 分层后仍分离 | 是**正面对照**，不等于支持生物学结论 |

**重要**：本工具不判断"生物学是否存在"，只判断"**当前数据能否支持所声称的结论**"。
审计结论为 `UNQUANTIFIABLE` / `NOT SEPARABLE` 时，正确表述是"**未定（indeterminate）**"，
而非"该研究错误"。

---

## 测试

```bash
pip install -e ".[dev]"                # 测试依赖（pytest）在 dev extra 里
python -m pytest tests/ -q            # 409 项（collected）＝ 405 passed + 4 skipped
                                      #  4 skipped = 2 项外部数据集门控 + 2 项内部文档门控
python -m pytest tests/ -q -k gating  # 仅门控
```

测试覆盖：阈值边界值、事件数分层、门控逻辑、交叉表非空、匹配检测、
深度判定与 mito 解耦、报告自包含性、端到端串联。

---

## 引用

若使用本工具，请引用：

> Yao Y. sc-audit: an operational audit framework for low-expression genes in
> case-control single-cell studies (v0.8.4). Zenodo (2026). doi:[DOI]

以及方法学依据（见 `docs/methods.md`）。

---

## 许可

MIT
