# Changelog

## 0.8.4 (2026-10-12)

**一键安装**

- 发布树根新增 `install.sh` / `install.ps1` / `install.py` / `install-macos.command`
  （＋ `make_offline_wheels.sh`）：检测 Python ≥ 3.9 → 建**独立 venv**（默认 `~/.sc-audit-venv`，
  可用 `SC_AUDIT_VENV` 覆盖）→ 装随包 wheel → 自检版本；`--offline` 用同目录 `wheels/` 离线安装。
- 装配脚本接入安装脚本拷贝与可执行位设置。

**测试**：**409** 项（本版无行为变更）。

## 0.8.3 (2026-10-12)

**内置提案调用器（可选）**

- 新增 `design --action propose`（`--backend none|local|http`，**默认 `none`**）：把
  `inventory.json` → `proposal.json` 收成一条命令，并自动写出处块（model / prompt_sha256 /
  inventory_sha256 / 时间戳 / 工具版本）。
- `urllib` / `subprocess` **延迟导入** —— 顶层与默认路径仍**零网络、零子进程**；
  联网/子进程仅发生在显式的 `propose`（且 `--backend` 非 `none`）时。

**测试**：405 → **409**（新增 4 项）。

## 0.8.2 (2026-10-12)

**设计解析契约加固（据独立评估反馈）**

- **数据集名不匹配即报错**：`design apply` 若 proposal 里没有当前数据集名，原先会静默取
  proposal 的第一条 → 现**报错退出**（退出码 2），不再静默误用。
- **`export` 校验分组标签**：case/control 标签缺失即报错（此前不校验）。
- **`apply` 登记幻觉列**：输出 `n_hallucinated` / `hallucinated_columns` 与弃权率。
- **复验判据同源**：`redundant`（condition 副本列）改用与设计变量花名册**同一份**判据；
  `donor_unique` 对混杂角色豁免。
- **常量列不自动接受**：单水平列的 `technical` / `design_defining` 提案一律转人工复核
  （不再自动接受）。
- **出处可溯源**：新增 `design --action attest`（proposal 出处块）；`export` 写
  `inventory_sha256`。工具**仍不内置任何模型**。

**测试**：398 → **405**（新增 7 项）。

## 0.8.1 (2026-10-12)

**设计解析闭环**

- **列注册表**：`apply` 以注册表里 `status="confirmed"` 的人工裁决**优先于**自动提案；
  再次运行同一数据集不再重复提问（人工确认单里单列「已确认」区）。
- **人工确认入口**：新增 `design --action confirm`（`--registry` + 逐列决策）。
- **配置装配**：新增 `apply_design_config_to_cfg()`，把确认后的设计配置装进 `AuditConfig`
  （**不自动调用**；是否采用由使用者决定）。

**测试**：396 → **398**（新增 2 项）。

## 0.8.0 (2026-10-12)

**新增 `design` 子命令（解耦版设计解析）**

- **架构**：工具**不内置语言模型**。工具只定义**契约**，模型在工具之外：
  `design --action export` 产出 `inventory.json`（逐列 dtype / 水平数 / 供者级水平数 /
  是否供者唯一 / 两臂熵 / 样例值 ＋ **确定性基线角色**；`--geotext` 只写文本指针路径、不读文本）
  →（工具外的模型出 `proposal.json`）→ `design --action apply` 读回并做**确定性复验**
  （提案角色 vs 结构判据）＋ **人工确认单**（风险分层：必须逐条 / 批量 / 自动＋抽检）＋
  产出一份**待确认**的 `design_config.json`（**不自动写入**）与 `design_registry.json`。
- `apply` 只依赖 `inventory.json` → 可**完全离线**。
- CLI：`--input` 改为可选（`design apply` 不需要输入数据），并加子命令分发守卫。

**测试**：392 → **396**（新增 4 项）。

## 0.7.3 (2026-10-12)

**缺陷修复**

- 供者级连续协变量（每供者 1 值、且各供者取值互不相同，例如供者数很少时的 `age`）会被
  `is_donor_unique` 误判为「供者唯一」而**当标识符丢弃**。现改为：被判为混杂角色的这类列进入
  **新桶** `excluded:donor-level-all-distinct` 并附显式说明（“不可匹配，但也不是标识符”）——
  判据只用**列名角色**，不看统计。判定行为不变。

**测试**：390 → **392**（新增 2 项）。

## 0.7.2 (2026-10-12)

**独立评估反馈修复（第二批）**

- **细胞身份门控门槛**：新增门控通过率**绝对下限**（`identity_gate_pass_floor_pct=1.0`）与
  **按臂差**上限（`identity_gate_arm_delta_max=2.0`）；触发时新增判定档 `GATE-NOT-USABLE`
  （门控近乎失效或与分组共线）。
- **对照构成**：新增 `NOT ASSESSED` 档（数据不具备可评估条件时，显式区分“未评估”与“不可分离”）。
- 匹配细胞数低于 `matched_min_cells_warn` 时，衰减比不再作判据（标 `attenuation_suppressed`），
  避免小匹配集上的过度解读。
- 谱系排名分母改用 `n_lineages`；检出地板表新增 `verdict_basis_group` / `verdict_floor_pct`
  两列（判定依据可回溯）。
- `check` 在无参照时的措辞；文档–代码同步检查对包外文件缺失降级为跳过；参照契约未强制时
  **显式登记**（不静默）。

**测试**：380 → **390**（新增 10 项）。

## 0.7.1 (2026-10-12)

**独立评估反馈修复（第一批）**

- **冗余副本列**：新增桶 `excluded:redundant-with-condition`——若某列在两臂内各自恒定、
  且水平数 ≥ 2（= condition 的完美副本），从设计变量、共同支持域与候选匹配中**剔除**，
  并在报告中声明。
- **放行路径限制判定**：经 `--allow-cross-modality` / `--allow-depth-mismatch` 放行、
  但其 R 依赖跨模态参照或深度不匹配时，判定**降级为 `INSUFFICIENT`**；受影响时新增
  `verdict_contract_ok` 列。
- `check` 在 backed（磁盘驻留）输入上改走 `io.x_matrix()`，避免全量读入内存。
- 细胞身份审计产出改为**动态可评估性清单**（删去硬编码的 `FULLY EVALUATED`）；检查 4
  无法执行时给显式说明，不再静默跳过。

**测试**：373 → **380**（新增 7 项）。

## 0.7.0 (2026-10-12)

**评审 L 组修复（工程/稳健）**

- **L1** B2 的“高基线”阈值不再硬编码：新增 `housekeeping_high_baseline`（0.8，入 config + yaml）。
- **L2** 潜在批次层选 HVG 改用**归一化离散度**（方差/均值，Fano 型），不再被少数超高表达基因主导。
- **L3** 潜在批次层 k-means 改**多初始化**（5 seed）：报告各 seed ARI 中位数并标记初始化敏感性。
- **L4** HTML 报告（`report.html`）现渲染 v0.6.x 新产物：floor 的深度/口径说明、A2 positivity 重叠行、
  以及 **⑤ Latent-batch scan** 段（`render_report(..., latent_batch=)`）。

**测试**：370 → **373**（新增 L1/L2+L3/L4 回归 3 项）。

## 0.6.2 (2026-10-12)

**评审 M 组修复（口径类）**

- **M1 positivity 纳入可调整混杂**：重叠剖面从 `design_defining` 扩为 `design_defining ∪ adjustable`（忽略 adjustable 会低估 positivity 违例）。输出新增 `profile_vars`。
- **M2 深度口径统一为中位数**：`load_reference_depth` 优先 `median_umi`（与 case 侧 median total_counts 对齐；避免 median/mean 混用）。
- **M3 下采样灵敏度不再静默**：计数非整数 → 不能稀释 → 在 `reference_depth_note` 显式声明（不再默默 NaN）。
- **M4 `check` 接入深度可比性**：预检与正式 `run` 一致，深度不匹配且未放行 → 退出码 2。
- **M5 判定链单一实现**：主链改调 `_verdict_from`，避免与敏感性口径漂移（行为逐字节不变）。

**测试**：365 → **370**（新增 M1–M5 回归 5 项）。

## 0.6.1 (2026-10-12)

**必修项（依 `代码逻辑与漏洞评审_v0.6.0_2026-10-12.md`）**

- **H1 删除死配置字段 `depth_metric_primary`**（定义未用；B3 仅作敏感性，主分层仍为 `n_genes`）。
- **H2 `identity.py` C2 段不再静默吞错**：双细胞列解析失败 → 显式 note（不得 `except: pass`）。
- **H3 `latent_batch` 口径与资源保护**：① X 非整数计数 → `INDETERMINATE` + 显式声明（不得把
  归一化输入当 counts 解读）；② 供者×基因超 `MAX_DENSE_ELEMENTS`（5e7）→ 跳过并声明
  （该扫描会稠密化，防内存耗尽）。

**测试**：361 → **365**（新增 H1/H2/H3 回归 4 项）。

## 0.6.0 (2026-10-12)

**方法学修订（A/B/C；据逐层面文献检索 + 代码逐处对照）**

- **A1 参照侧深度可比性**：参照表若提供 `mean_umi`（列或 manifest），与 case 中位 `total_counts`
  比对；相对差 > `reference_depth_delta_max_pct`（默认 50，事前锁定）→ 默认**拒绝**
  （`--allow-depth-mismatch` 放行）。`--depth-matched-sensitivity` 输出二项稀释到参照深度后的
  R/判定列（`case_pct_depth_matched` / `verdict_depth_matched`）。**未提供深度时不改变输出**。
- **A2 positivity / 协变量重叠**（Petersen 2012；Westreich 2010）：新增**供体级协变量重叠**
  诊断（`overlap_case` / `overlap_ctrl` / `overlap_min`）。**设计决策**：A2 作**并列诊断 + 强制警示**，
  **不覆盖**阶梯判定（阶梯已将估计限定在共同支持域；用 NOT SEPARABLE 覆盖会丢失衰减信息）。
  细胞计数退为**容量前提**。
- **B1 潜在批次层**（`--with-latent-batch`，默认关闭）：供体级 pseudobulk → HVG → PCA(2) →
  k-means(2) → 与 condition 的 ARI；`ARI >= latent_batch_ari_warn`（0.5）→ `LATENT-BATCH-RISK`。
  报警层，**不得**把簇当匹配变量。新增产物 `latent_batch.txt`。
- **B2 管家护栏一致性**：新增各管家 R 的 CV 告警（`housekeeping_cv_warn`）与“高基线组织
  适用域”声明（饱和性 `R_S2=(1-sqrt(1-p))/p`）。
- **B3 分层口径敏感性**：额外用 `total_counts` 分层与 `n_genes` 对照，报告持续性是否口径依赖；
  **不默改主判定**。
- **C1–C5 透明项**：ambient RNA 局限（SoupX/DecontX/Janssen 2023）、双细胞率（Scrublet/DoubletFinder）、
  注释方法依赖（Abdelaal 2019）、`binom_P_case0` 的细胞级（伪重复）口径、小面板 FDR 功效 ——
  均显式写入输出（`identity.txt` / `floor` summary）。

**零影响**：未给 `--reference`（或参照无深度/模态声明）时，审计一的数值产物**逐字节不变**；
审计三/四新增的 note 会出现在 `robustness.txt` / `identity.txt`（[增强]）。

**文档**：`README.md` + `docs/methods.md` 附录 A；修订提案 `方法学修订提案_v0.6.0_2026-10-12.md`。

**测试**：347 → **361**（新增 `tests/test_v060_methodology.py` 14 项：A1×5 / A2×3 / B1×2 / B2×1 / B3×1 / C×2）。

## 0.5.0 (2026-10-11)

**新增接口 + 行为变更（参照契约 / REFERENCE_SCHEMA R1+R2）**

- **R2 跨 modality 默认禁止**：不同 modality（scRNA vs snRNA）的检出率不可比——细胞质转录本在
  whole-cell 与 single-nucleus 制备中的丰度天然不同。当**目标模态与参照模态都已声明且不同**时，
  `floor`（及 `run`/`check`）**报错**（`CrossModalityError`，退出码 2）；仅在向
  `--allow-cross-modality`（config `allow_cross_modality`）显式放行时才继续，且受影响基因的 R 标
  `cross_modality=True`——不得用于 `floor_pct`/`ratio_*` 的可定量判定。
  **行为变更**：旧版仅“告警 + 标 cross-modal”，现默认拒绝。任一侧未声明 → `unknown`（不报错、不静默）。
- **R1 五元组合并守卫**：参照旁 `*.manifest.json` 的 `column_contract` 提供逐列五元组
  `(modality, intron_counting, gate_spec_hash, upstream_qc_hash, cell_source)` 时，
  不一致即**拒绝合并**并列出差异字段（`ReferenceContractError`，退出码 2）；
  未提供逐列元数据时**显式登记“R1 未强制”**（绝不静默）。
- `--modality` / `--reference-modalities` 现也注册给 `floor` 子命令（此前仅 `run`/`identity`）；
  参照模态集合统一取自 `--reference-modalities` ∪ 参照表声明的模态（CSV `modality` 列 / manifest）。
- `run_metadata.json` 新增 `reference_contract` 段。

**零影响**：未给 `--reference`（或参照未声明模态）的默认路径**逐字节不变**（无判定改变）。

**文档**：`README.md` 新增「参照表契约：跨 modality 与合并」；`docs/methods.md` §5.1；
`docs/REFERENCE_SCHEMA_v1.md` 更新实现进度（**R3/R4/R5 仍属设计稿**，未实现）。

**测试**：335 → **347**（新增 12 项：R2 默认禁止 / 放行 / 同模态 / 未声明 + R1 一致 / 不一致 / 未强制 + CLI 层）。

## 0.4.0 (2026-10-11)

**新增接口（默认路径零影响）**

- **设计声明 + 能力闸门（Tier 1）**：工具在 v0.3.x 及以前**默认假设**输入是「二元、
  横断面、分析单位 = 供者」的病例–对照设计。对纵向 / 治疗臂 / 多状态设计若照旧计算，
  会**静默给出无效判定**（例如把复测单位当独立观测、把多臂的对照含义搞错）。现新增：
  - **设计声明**（CLI `--comparison-type` / `--unit` / `--contrast-a|b` / `--time-col` /
    `--ref-time` / `--index-time` / `--index-event`，或 config 的 `design_*` 字段）：
    **只记录、不推断**；全部默认值 = 支持集 S0
    = `{two_arm, donor, time=None, arms=2}`。
  - **能力闸门**：声明 ∈ S0 → 正常执行；声明合法但 ∉ S0 → **拒绝输出任何判定标签**
    （`verdict = DESIGN-NOT-SUPPORTED`，**退出码 3**），不做审计、不产出判定产物；
    声明取值不在封闭词表内（写错）→ 错误（**退出码 2**）。
  - **结构提示**：数据形似更复杂设计（如 condition 水平数 > 2）→ `[design-hint]` 警告，
    **绝不改变任何判定标签**。
  - `run_metadata.json` 新增 `design` 段（`declared` / `supported` / `status` /
    `supported_set` / `structural_hints`）。
- **零影响承诺**：不声明 = 与 v0.3.x **逐字节一致**（默认路径判定、跨数据集结论与 v0.3.x 冻结结果不变）；显式声明 `two_arm` + `donor` 与不声明等价。
  **`DESIGN-NOT-SUPPORTED` 与 `INDETERMINATE` 是两回事**：前者「设计超出能力、
  根本没做审计」，后者「支持的设计下、审计做了但无法给结论」。
- **明确排除**（有意推迟）：纵向 / 治疗臂 / 多状态设计的**实际实现**（Tier 2/3）、
  因果角色（confounder / collider）推断。取舍原则：**把复杂度接口化并推迟** ——
  先保证「不懂时明确拒绝」，而不是静默算错。

**输入契约**

- `cardiac.yaml` 补齐与 `AuditConfig()` 的 9 个 `design_*` 字段（全字段逐值一致测试继续强制）。

**文档**

- `README.md` 新增「设计声明与能力闸门」一节；`docs/methods.md` 新增 §8（原 §8 引用与致谢顺延为 §9）。

**测试**：309 → **335**（新增 26 项：闸门语义 / 声明词表 / 结构提示 / CLI 集成与
「默认 vs 显式 S0」零影响回归）。

## 0.3.3 (2026-10-08)

**行为变更（会改变判定标签，务必注意）**

- **支持域 = 整个对照组即不再报 `COMPARABLE`**：旧守卫只在“没有 design-defining 变量”
  （`vacuous`）时触发。若数据集里存在**在该数据集内恒定的零信息列**（例如全队列同一
  `region`），它会被归入 design-defining 桶 → 守卫失效 → 支持域退化为整个对照组、
  所有规则衰减比恒 = 1，于是落进单侧分支被报成 `COMPARABLE` —— 与本工具自己的文案
  （"do not read this state as 'COMPARABLE'"）自相矛盾。现改为：**支持域等于整个对照组
  且没有任何可评估规则真正缩域** → `INDETERMINATE`（注释局限）+ 强制提示。
  影响面：仅该类数据；真分层列（支持域被压到地板以下）仍走 `SUPPORT-BELOW-FLOOR` 不变。
- **小样本警示的精确枚举下限改正**：提示文案原先写 “最小 P ≈ 1/C(n,k)，可能永远达不到
  显著”。精确枚举下最小可达**双侧** P = 2/C(n,k)（两组等大）/ 1/C(n,k)（不等大，且要求
  极值划分唯一）—— 对等大组且阳性全部分离的情形，5 vs 5 也可达 0.0079。现文案给出两种
  情形并去掉绝对断言。

**缺陷修复（文案 ↔ 行为）**

- `floor` 事件稀疏分支的提示原先写 "jackknife is still reported"，实现里却 `continue` 跳过
  了整段 jackknife（示例产物 `floor_jackknife.csv` 里也无该基因）。现改为
  "no full-sample P … and the gene is excluded from the jackknife under the rule"，
  `summary['perm_p_note']` 同步。
- `docs/methods.md` §6.2 的 `CAPTURE-LIMITED` 条件原写作「复杂度最低 / 核标志最低 / 跨模态」
  （读作“或”，且误含“跨模态”）。实现为复杂度低 **且** 核标志低；现表格与实现一致。
- **外部审查发现的文档错误（同日）**：§2.2 的二项精确区间“跨度超过 30 个百分点”实为
  **1.69 个百分点**（k=10, n=800，CP [0.60%, 2.29%]）；§2.3 说二项检验用“整体检出率”，
  实现用**对照组**检出率；§2.5 口径表两处错位（双侧相对半宽有列 `cp_rel_halfwidth_case`；
  `case_pct_upper95` 是**绝对**上界）；§4.3 写管家“不低于参照”，实现阈值是 0.6。
- **接口层（同日）**：`identifiability --matched-donors` 原先被 argparse 注册却**从未生效**
  （静默 no-op）→ 现显式拒绝（exit 2，并指引用 `--ladder`）；`run` 补 `--ladder`（与
  `identifiability` 同语义）；`run_metadata.json` 补 `perm_mode`（可复现留痕）。
- **健壮性（同日）**：审计四检查 6 在“目标谱系只出现于单臂”时对**空集合求均值** → 静默
  NaN + numpy 警告。现改为显式提示“目标谱系在 {臂} 完全缺席（不是 0）”，与检查 5 的处理一致。
  检查 8（取材标志按臂）同口径修正。
- **`check` 的参数覆盖（同日）**：`_run_check` 原先只覆盖 `n_perm`，于是 `--count-layer` /
  `--ratio-basis` / `--hetero-col` / `--match-vars` 在**预检**里被静默忽略 → 预检比正式运行
  **更严**，会把“计数放在非标准 layer 里”误报成阻断项（而 `floor --count-layer` 能跑通）。
  现预检与 `_prepare` 用**同一套**参数覆盖。
- **检查 5 的归因（同日）**：多谱系且未指定 `identity_target` 时，旧实现把跳过原因写成
  “obs 中无 'source'”——**即使该列确实存在**。现分开归因（“未指定目标谱系” vs “缺列”）。
- **`docs/methods.md` §6.2 真值表（同日）**：`INDETERMINATE` 行删除“只有一个谱系且核心检查
  也无法执行”（与实现不符——单谱系只把检查 1/2/3 降为报告模式，不改变标签），并加注说明
  “报告模式降级 / 检查跳过**不改变**判定标签，但此时 `PASS` 只覆盖**已执行**的检查”。

**输入契约**

- `cardiac.yaml` 补齐与 `AuditConfig()` 的 14 个缺失字段（`fdr_alpha` + 13 个 `identity_*`）；
  对应测试由“抽指定字段”改为**全字段**逐值比对（原文件 52/66 却注释称“逐字段一致”）。

**发布卫生**

- 出厂文件不再点明内审平台名与内审产物目录，不含未发表案例的真实供者 ID 或供者级
  治疗/结局列名，也不以未发表案例的真实基因作演示（演示地板基因改为中性基因）；
  代码注释与测试里的历史评审轮次编号改为描述性说明。
- **有意保留**（已在 `DOCS_INDEX.md` 写明）：`SC_AUDIT_EXTERNAL_CASE_H5AD` 门控的
  外部数据冒烟回归，以及 2 项依赖**开发仓库**审计文档/维护脚本的契约测试 —— 二者
  在发布树内自动跳过，属于刻意设计（它们断言的是维护侧一致性，不是工具行为）。

**测试**：301 → **309**（新增 8 项：零信息列守卫 2、事件稀疏 jackknife 文案 1、
小样本下限文案 1、目标谱系单臂缺席 1、identifiability 拒绝 --matched-donors 与 run --ladder 1、
检查 5 归因 1、check --count-layer 生效 1）

## 0.3.2 (2026-10-07)

**行为变更（会改变判定标签，务必注意）**

- **真空支持域不再报 `COMPARABLE`**：当没有 design-defining 变量（共同支持域按构造 = 整个对照组）
  且没有任何可评估匹配规则**真正缩小**支持域时，衰减比恒 = 1（by construction）——这属于
  "注释局限"，现统一判定为 `INDETERMINATE` 并附强制提示
  （"not evidence of comparability … do not read this state as 'COMPARABLE'"）。
  影响面：仅该类数据；支持域非真空的数据集判定不变。
- **阶梯不再生成空合取规则**：变量列表为空的规则等价于"匹配所有"，既不缩域也不承载一致性信息，
  现不再生成（此前会固定贡献一个 1.000 的"一致性证据"）。

**缺陷修复**

- `check` 的 donor-unique 判定：索引未对齐导致供者内分组为空、判据退化为"水平数 ≥ 供者数"，
  使"供者内变动"的列被误判为 donor-unique（并给出事实错误的原因文案）。现加索引对齐守卫；
  无法分组时返回 `False`（宁漏不误）。
- `gate_cells`：`qc_from='obs'` 但 obs 缺少所需 QC 列时会静默回退为从传入矩阵现算（在基因子集上
  会导致门控全灭且无提示）。现缺失列即告警，并说明该指标已回退。

**新增**

- `gate_cells(..., qc_from="x"|"obs")`：默认从传入矩阵现算（语义与 0.3.1 一致）；在**基因子集**上
  门控时可用 `"obs"` 复用 `compute_qc_metrics` 在全基因空间算得的 QC。口径记入
  `ad.uns["sc_audit_gating"]`（`qc_from` / `qc_source`）。
- 公开提示：Python API 不做基因名自动回退时（`var_names` 为 Ensembl ID 且有符号列），
  工具会显式提示改用 `apply_gene_name_fallback()` 或 CLI；文档与之一致。

**测试**：292 → **301**（新增 9 项回归：门控 QC 口径 2、判定语义与索引错位 5、API 回退 2）

## 0.3.1 (2026-10-06)

四项审计（检出地板 / 对照构成 / 稳健信号 / 细胞身份）+ 审计二 v0.3 可识别性层。

- **审计二 v0.3 可识别性层**：变量花名册 → 病例组内设计剖面分桶 → 共同支持域 →
  事前锁定的匹配规则阶梯 → 一致性检查；判定标签
  `COMPARABLE` / `NOT SEPARABLE` / `SUPPORT-BELOW-FLOOR` /
  `THRESHOLD-SPANNING` / `INDETERMINATE`。
- **单侧判定**：阈值落在观测档位区间之外时直接给方向判定，并强制附
  "一致性检查无信息"的警示（可关：`identifiability_one_sided_verdict=False`）。
- **健壮性**：缺设计列、目标基因 0 命中、空阶梯、空面板均不再崩溃，
  退化为显式判定/带列名空表。
- **一致性**：非设计 QC 列名单收敛到单一真值源；报告配色区分正/负面判定。
- 阈值全部**事前锁定**（`config.py`），运行期不改写。

早期版本（v0.1–v0.3.0）为内部迭代，未单独归档。
