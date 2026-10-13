# sc-audit · 方法学说明(操作定义与阈值)

> 本文档可直接用作投稿稿件的 **Supplementary Methods**,或作为工具文档的方法学依据。
> 所有阈值均为**事前锁定(pre-specified)**,不因结果而调整。

---

## 1. 设计原则

工具不是判断"生物学是否存在",而是判断"**当前数据能否支持所声称的结论**"。因此:

1. **独立于声称**:审计在不知道组间比较结果的前提下计算,避免"看着结果选判据"。
2. **事前锁定**:阈值写在配置文件中,分析前固定。
3. **多重判据并列**:任一"不可定量"判定必须由 ≥3 项指标共同支持,不依赖单一数字。
4. **不回溯修改**:分层后若结论反转或不稳健,写入局限,**不回头改主分析**。

---

## 2. 审计一:检出地板

### 2.1 指标定义

```
pct_study   = 目标细胞群中该基因 count > 0 的比例
pct_ref     = 各外部参照背景中最高的检出率(主口径)
              (敏感性分析:中位数口径)
R           = pct_study / pct_ref                 归一化检出比
n_pos       = 绝对阳性细胞数(事件数)
n_pos_case  = 病例组阳性细胞数
```

### 2.2 阈值(事前锁定)

| 判据 | 阈值 | 处置 |
|---|---|---|
| `pct_study` | `< 5%` | 不报告组间比较 |
| `R` | `< 0.3` | dropout-limited |
| `R` 且 `pct_study` | `≥ 0.6` 且 `≥ 20%` | 可定量 |
| **`n_pos`** | **`< 20`** | **仅描述性统计;不做组间显著性检验** |
| `n_pos` | `20-50` | 比较须附精确二项 CI 与 jackknife 范围 |
| `n_pos` | `≥ 50` | 标准分析 |

**事件数阈值的方法学理由**：比例低不等于不可定量，但**绝对事件数低必然不可定量**。
一个 1.2% 的检出率在 800 个细胞中仅对应约 10 个阳性细胞；此时二项分布的 95% 精确区间
（Clopper–Pearson）为 **[0.60%, 2.29%]——跨度 1.69 个百分点、相对半宽约 ±67%**，
任何 P 值反映的是抽样噪声而非生物学。

### 2.3 三重支持量("不可定量"结论必须同时给出)

1. **二项检验**：病例组观察到 0 个阳性的概率（在**对照组**检出率下，即 `(1 − p_ctrl)^n_case`；仅当 `n_pos_case = 0` 时给出，并显式标注其“细胞独立”假设）。
2. **留一供者敏感性**:逐一移除每个对照供者,报告置换 P 的**范围与倍数波动**,
   以及效应量的范围。报告的是范围而非二值显著性判断--因为二值判断会依赖置换种子。
3. **同地板参照基因**:同数据集中检出率相近的其他基因是否有同类信号。
   若某基因"显著"而其同地板邻居不显著,提示地板效应而非生物学梯度。

### 2.4 不得作为结论的情形

- 事件稀疏基因(`n_pos < 20`)的组间比较
- 显著性依赖单一供者的结果
- 检出比 `R < 0.3` 的基因间差异

---


### 2.5 两种区间口径不得混用(外部审查 B1,2026-10-06 核实)

| 量 | 定义 | 用途 | 工具中的列 |
|---|---|---|---|
| **双侧相对半宽** | (upper − lower) / 2 / p̂，CP 精确**双侧** 95% | 预注册容差表（**±43.4%** → x_min 随 n 变化） | `cp_rel_halfwidth_case`（病例组） |
| **单侧绝对上界** | CP 精确**单侧** 95% 上界（**绝对百分点**，不是相对值） | “以此样本量，真实比例至多约 x%” | `case_pct_upper95` / `overall_pct_upper95` |

复算示例(n_case = 200,本机复算):

| k | p̂ | 双侧相对半宽 |
|---|---|---|
| 19 | 0.0950 | 0.454 |
| **20** | **0.1000** | **0.440** |
| 21 | 0.1050 | 0.428 |

**结论**:锁定阈值 20 对应的是**双侧**相对半宽优于 ±43.4% 这一事先声明的容差;
x_min 随 n_case 变化(`sc_audit.floor.event_threshold_table(n)` 可自行复算,与
`derive_event_threshold(n, tol)` 互补)。而单侧 95% 上界 ≡ 双侧 **90%** 上界,
因此**不能**把工具的单侧列代入双侧容差表--那样得到的 x_min 会偏小。这是**口径**
问题,不是阈值问题;阈值与判定逻辑均未改动。

## 3. 审计二:对照构成

### 3.1 指标定义

```
F_hetero    = 对照中来自「非病例组主机构」的细胞比例
匹配供者     = 在所有可用设计变量上与病例组一致的对照供者
衰减比       = |匹配后效应量| / |匹配前效应量|
```

### 3.2 匹配规则(事前锁定)

**严格规则**:同时匹配 `institution ∧ donation_type ∧ source ∧ region`

**报告**:
- 每个对照供者的**一致变量数 / 可用变量数**(完全透明,便于判断匹配严格度)
- 严格口径与"仅同机构"宽松口径的**对照结果**
- 两套口径下的效应量衰减

### 3.3 判定

| 判定 | 条件 |
|---|---|
| `NOT SEPARABLE` | `F_hetero > 0.5` 或 衰减比中位数 `< 0.5`（后者仅在匹配细胞数 ≥ `matched_min_cells_warn` 时作判据，C7/v0.7.2） |
| `COMPARABLE` | 存在 ≥2 水平的设计变量，且上述均不满足 |
| `NOT ASSESSED` | **没有任何 ≥ 2 水平的设计变量** → 可比性无从评估（C5/v0.7.2；不读作 COMPARABLE） |
| `INDETERMINATE` | 无可用设计变量 |

### 3.4 强制限定语

匹配子集通常很小,**不能**用来论证"生物学效应不存在"。
可解释的结果是**效应量衰减的幅度本身**--它量化了设计无法分离的混杂份额。

---

## 4. 审计三:稳健信号对照

### 4.1 指标

```
深度指标   = 每细胞检出基因数、每细胞总计数(中位数对比)
细胞状态   = 线粒体比例(独立报告,不参与深度判定)
深度分层   = 按每细胞基因数分箱,箱内比较两组检出率
管家基因   = ACTB / GAPDH 的检出率与参照的比值
```

### 4.2 判定

| 判定 | 条件 |
|---|---|
| `DEPTH-COMPARABLE` | 深度指标差异均 < 20%,且分层后组间仍分离 |
| `PARTIAL` | 深度指标差异 ≥ 20% |
| `INSUFFICIENT` | 深度可比性无法建立（`INSUFFICIENT — depth comparability could not be established`）；反之，若深度可比但**无任何分层满足每组最小细胞数**，则判 `DEPTH-COMPARABLE (… separation could not be assessed)`（带“分离无法评估”限定），**不**降为 `INSUFFICIENT` |

**关键区分**:`DEPTH-COMPARABLE` 是**正面对照**,说明分离不是深度假象;
它**不等于**支持生物学结论。线粒体比例差异属细胞状态/取材指标,单独报告。

### 4.3 替代解释排除

管家基因的检出比 *R* ≥ `housekeeping_ratio_ok`（默认 **0.6**）时，视为同时排除：
1. 细胞类型判错(确实是目标细胞)
2. 文库整体劣化(测序质量未全局受损)

此时"低表达基因检出异常"应归因于**检出地板**,而非数据整体不可用。

---

## 5. 跨数据集比较

**禁止直接比较绝对检出率**(CellRanger 版本、核 vs 细胞、建库批次差异可使绝对检出率
相差数倍)。统一使用归一化 `R`。

### 5.1 参照契约(硬约束,v0.5.0)

对应 `docs/REFERENCE_SCHEMA_v1.md` 第 III.2 的机器强制规则:

- **R2 跨 modality 默认禁止**:`modality ∈ {scRNA, snRNA}` 决定“抓哪个 RNA 池”,
  细胞质转录本在 whole-cell 与 single-nucleus 制备中的丰度天然不同 → 绝对检出率不可比。
  当**目标模态**与**参照模态都已声明且不同**时,工具**报错(`CrossModalityError`,CLI 退出码 2)**;
  仅在 `--allow-cross-modality`(config: `allow_cross_modality`)显式放行时才继续,
  且受影响基因的 R 标 `cross_modality=True`——**不得**用于 `floor_pct`/`ratio_*` 的可定量判定。
  任一侧未声明模态 → `ref_modality_flag = unknown`,不报错、也不静默。
- **R1 五元组合并守卫**:把多个参照列并成同一参照时,
  `(modality, intron_counting, gate_spec_hash, upstream_qc_hash, cell_source)`
  必须完全一致(写于参照旁的 `*.manifest.json` 的 `column_contract`)——任一项不同则
  **拒绝合并**并列出差异字段(`ReferenceContractError`,退出码 2)。未提供逐列五元组时
  工具**显式登记“R1 未强制”**,绝不静默跳过。
- **尚未实现**(仍属设计稿):R3(`intron_counting=unknown` 降级)、R4(look-outside-one-tissue /LOSO
  自动剔除,长表形态)、R5(剔除后可用列 < 2 → INDETERMINATE)。现行参照表为宽表格式。

`sensitivity`:参照口径(max vs median)、最小数据集排除、
门控严格度(阳性 marker 数)、细胞级 vs 供者级聚合--四项均作为敏感性分析报告,
**不用于替换主口径**。

---

## 6. 审计四:细胞身份与捕获完整性(v0.2 新增)

由病例-对照单细胞数据的细胞身份取证实践逼出。前三项审计**均无法发现**它检出的
问题,因此作为独立模块加入。它在两个层面上是新的:

1. **把审计二下沉到细胞类型内部**--一个被正确注释的细胞簇,其内部仍可能存在
   只在其中一个臂出现的亚群(例如:某细胞类型内、某分选来源的细胞**全部**集中在
   单一组)。
2. **把 QC 规则本身当作待审计对象**--"阴性标志全 0"这条门控规则可能与取材方式/部位
   共线(例如病例组活检取材 vs 对照供心)。此时一致施加门控会**差异性滤除**一个组,
   在"去掉一个构成差异"的名义下**引入一个新的**。

### 6.1 八项检查与阈值

| # | 字段 | 默认阈值 | 含义 |
|---|---|---|---|
| 1 | `identity_gate_margin` | 5.0 | 目标谱系门控通过率 / 其他谱系中位数 |
| 1b | `identity_gate_pass_floor_pct` | 1.0 | 门控通过率**绝对下限**（%）;低于此 → `GATE-NOT-USABLE`（C4/v0.7.2） |
| 1c | `identity_gate_arm_delta_max` | 2.0 | 门控通过率的**臂间比**上限;超过 → 门控与分组共线（计入 COMPOSITION-CONFOUNDED） |
| 2 | `identity_complexity_col` | `n_genes` | 排名;=1(最低)触发告警(需 ≥4 个谱系)|
| 3 | `identity_nuclear_marker` | `MALAT1` | 排名 =1 触发"核含量低" |
| 4 | `identity_complexity_bins` | 200/400/600/800/1200 | 分箱匹配,每箱双侧 ≥ `identity_bin_min_cells` |
| 5 | `identity_imbalance_min_cells` / `_min_tvd` | 20 / 0.5 | 亚群只在单臂出现，**或**该亚群在供者间的分布偏离目标细胞总体分布（**供者级** TVD > 0.5） |
| 6 | `identity_collinear_min_delta_pp` | 10.0 | 阴性标志两臂差 |
| 7 | `--modality` / `--reference-modalities` | - | 无同模态列 → 全部 R 标 `cross-modal` |
| 8 | `identity_procurement_markers` | PECAM1/VWF/DCN/COL1A1 | 按臂报告 |

### 6.2 判定

| 判定 | 条件 |
|---|---|
| `INDETERMINATE` | 缺注释列 / 缺阳性标志（门控可评估性无从谈起）|
| `COMPOSITION-CONFOUNDED (within cell type)` | 亚群不平衡 **或** 门控共线 |
| `CAPTURE-LIMITED` | 复杂度最低 **且** 核标志最低(模态不匹配只加 flag,不是判定条件) |
| `GATE-NOT-USABLE (own-gate pass rate below floor)` | 目标自家门控通过率 < `identity_gate_pass_floor_pct`（C4/v0.7.2） |
| `PASS` | 无上述任一 |

> **“只有一个谱系”与“核心检查跳过”不影响判定标签**（v0.3.3 澄清）：只有 1 个谱系时，
> 检查 1/2/3（相对排名）降级为**报告模式**（在 notes 里逐项写明），检查 5/6/8 仍然执行；
> 若**未指定 `identity_target`** 或缺**来源/分选列**，则 5/6/8 跳过并在 notes 里**分开归因**
> （“未指定目标谱系” vs “缺列”）。注意：此时 `PASS` 只覆盖**已执行**的检查，
> **不构成**“亚群不平衡与门控共线均不存在”的证据。

### 6.3 关键实现约束

- **backed 兼容**：`backed='r'` 的 h5ad 可直接作为输入——`floor` 与 `identity` 会按需物化所需列
  （审计四用批量 `_detect_many` 一次列切片 `to_memory`）。仅当对象无法物化时，
  `gene_detection_vector` 才返回 `None`（此时按“未检出”处理，并在表中反映）。
  （2026-10-08 复核更正：旧文本称“对**任何**基因都返回 None”，与实际不符——实测
  `backed='r'` 下 `gene_detection_vector` 正常返回计数。）
- **单谱系只降级依赖谱系对照的检查(1/2/3)**,不阻断 5/6/7/8--后两者是核心项,
  不该被一个无关条件连带关闭。
- **缺列/缺基因一律显式 note**,绝不静默跳过。
- **门控的 QC 口径**:`gate_cells` 的 `min_counts`/`min_genes` 默认从**传入矩阵**现算(`qc_from="x"`,保证深度指标与被审矩阵一致);若要在**基因子集**上跑门控(例如只保留面板基因),现算深度只覆盖子集、阈值必然不过、门控会**静默全灭** -- 此时改用 `qc_from="obs"`(复用 `compute_qc_metrics` 在全基因空间算得的 QC),工具在 0 通过且 obs 有 QC 列时会给出该提示(2026-10-07 已知答案验证发现)。口径与来源记录在 `ad.uns["sc_audit_gating"]`。
- **弱信号守卫(weak-signal caveat)**:单个捕获信号(复杂度或核滞留
  读数)**不构成**推断--两个信号必须同时出现才报 `CAPTURE-LIMITED`,
  否则只报"低于双信号阈值,未判"。

---

## 7. 审计二 · v0.3 可识别性层(v0.3.0 并入;v0.3.1 按审计修订;v0.3.2 增真空支持域规则;v0.3.3 放宽真空支持域守卫)

> **与 §3 的关系(统一口径,取代早前的自相矛盾表述)**:
> §3 描述的单规则口径(`composition.py`)**仍然存在且行为未变**,作为**历史/对照**保留;
> 本节的层是**推荐判定**:先问"这个设计里是否存在可调整的组间估计",再问"衰减比
> 跨规则是否稳定"。`run` 会同时输出两者;报告里显式交代关系,冲突时以本节为准。

### 7.1 解决的问题

原审计二把用户手写的 `match_vars` 当作混杂变量定义,导致:
(a) 病例组内**恒定**的设计变量(如 institution/region)被当作混杂去匹配,匹配坍缩到
个位数供者;(b) 病例组内**有变异**的可调整变量(如 age/gender)被排除在判定之外;
(c) 单条规则下的一个中位衰减数被当作"混杂份额"的定量证据,而它对规则高度敏感。

### 7.2 操作定义

1. **变量花名册**:obs 的每个非 QC 列必须落桶(未归类即报出,不得静默跳过):

   | 桶 | 判据 | 处置 |
   |---|---|---|
   | `excluded:donor-unique` | 每供者内恒定且水平数 ≥ 供者数 | 匹配无意义 |
   | `excluded:donor-level-all-distinct` | `confounder` 角色且每供者 1 值、取值互异（如供者数少时的 age） | 不可匹配，但**不是标识符**，需进局限讨论（v0.7.3） |
   | `excluded:all-missing` | 整列缺失 | 不进入任何判定 |
   | `excluded:derived-from-expression` | 由表达量派生的评分列(显式声明 + 默认前缀 `score_`,属**命名启发式**) | 条件化 = collider |
   | `excluded:outcome(collider)` | 角色 = outcome(词根 / `role_overrides`) | 禁止匹配 |
   | `excluded:post-baseline` | 角色 = intervention | 暴露/治疗后协变量 |
   | `excluded:technical(nuisance)` | 角色 = technical | 批次/通道等 |
   | `excluded:redundant-with-condition` | **两臂各自恒定且水平 ≥ 2**（与 condition 一对一重编码） | 条件化在它上面 = 条件化在待估对比上（v0.7.1） |
   | `design-defining(scope)` | **病例组内供者级熵 ≈ 0** | 不可匹配,只作适用范围声明 |
   | `adjustable(confounder)` | 病例组内有变异 | 真正的可调整混杂,必须进混杂集 |

   注:`design-defining` 是**本数据集的属性**,不是群体层面的证明;病例供者数少时,
   某个变量可能因抽样恰好恒定(报告会输出 `chance-constancy risk` = p^n+(1-p)^n)。

2. **共同支持域**(可识别性 / positivity,先于匹配):在 design-defining 变量上取
   病例组设计剖面,统计落在该剖面的对照供者与细胞数。若出现多个病例剖面,**不池化**
   (池化即跨层比较),只报告单一最大剖面并报 `INDETERMINATE`。
3. **匹配规则阶梯**:事前锁定的规则序列 R0(不匹配,锚点)→ R1(现行正文口径)→
   R2(全部候选变量合取)→ R3(design-defining 合取)→ R4(单变量)→ R5(推荐窗内合取)。
   逐规则 × 逐基因计算 `full_pp`、`matched_pp`、`ratio = matched/full`;
   只对 `full_pp ≥ attenuation_min_full_effect_pp`(5 pp)的基因取中位。
   `ratio` 是**效应量比值**,**不是**"混杂份额":匹配子集可能在剩余可调整变量上仍不平衡,
   且 `ratio > 1`(放大幅度)在实测中出现(报告会给出放大档数)。
4. **一致性检查**:中位衰减比在阶梯上是否**跨过** `attenuation_min`。

### 7.3 判定真值表(当前有效,v0.3.3)

| 顺序 | 条件 | verdict | 是否报告衰减比 |
|---|---|---|---|
| 1 | `condition` 不含两组 | `INDETERMINATE` | 否(早退) |
| 2 | 支持域细胞数 < `matched_min_cells_warn`(50) | `SUPPORT-BELOW-FLOOR` | **否**(抑制) |
| 3 | 病例组跨多个设计剖面 | `INDETERMINATE` | 否(不池化) |
| 4 | 无可评估规则(含"无目标基因命中") | `INDETERMINATE` | 否 |
| 5 | **无可用设计/可调整变量**,或**支持域 = 整个对照组**(真空,或存在零信息列)**且无任何可评估规则缩域**(v0.3.3) | `INDETERMINATE`(注释局限) | **否**(抑制) |
| 6 | 中位衰减比**跨**阈值(两档并存) | `THRESHOLD-SPANNING` | 报告区间,不给单一结论 |
| 7 | 单侧且 `identifiability_one_sided_verdict=False`(旧保守读法) | `INDETERMINATE` | 区间 |
| 8a | 全部规则**≥** 阈值(单侧) | `COMPARABLE` + **强制警示** | 是(区间) |
| 8b | 全部规则**<** 阈值(单侧) | `NOT SEPARABLE` + **强制警示** | 是(区间) |

**强制警示(单侧情形必附)**:阈值落在观测档位区间**之外**,故一致性检查在此**不携带
任何信息**(所有规则按构造落在同一侧);该判定仅基于量级比较,不得脱离警示引用。

**历史说明**:v0.3.0 的实现里"两档并存"与"跨规则翻转"是同一谓词,
于是两档并存 ⇒ 必判 `RULE-DEPENDENT`、单侧 ⇒ 被压为 `INDETERMINATE`,
`COMPARABLE` / `NOT SEPARABLE` 在默认配置下**永不可达**(真值表第 6 行是死行)。
v0.3.1 改为上表:单侧给方向 + 强制警示。若需要旧的保守读法
(两端都报 `INDETERMINATE`),把 `identifiability_one_sided_verdict` 置 `False`。

**v0.3.2 追加(已知答案验证 + 独立复核)**:真值表第 5 行被显式化--"无可用变量"与
"支持域真空且无任何可评估规则缩域"合并到同一支(都属**注释局限**,衰减比按构造恒 = 1,
一律 `INDETERMINATE` + 强制提示,不再落进单侧分支被报成 `COMPARABLE`);同时阶梯
**不再生成空合取规则**(变量列表为空 ⇒ 匹配所有,不缩域、不承载一致性信息)。

### 7.4 阈值

| 字段 | 默认 | 说明 |
|---|---|---|
| `identifiability_ladder` | `None` → 内置 R0-R5 | 显式给出即视为运行前锁定(写入 `run_metadata.json`);CLI 要求 ≥2 条规则 |
| `identifiability_min_support_cells` | `None` → 复用 `matched_min_cells_warn` | 单一真值源;`0` = 显式不设下限(不得被 `or` 吞掉) |
| `matched_min_cells_warn` | 50 | 沿用 v0.2.5(自 v0.3.0 起保持) |
| `attenuation_min` | 0.5 | 沿用 v0.2.5 |
| `attenuation_min_full_effect_pp` | 5.0 | 沿用 v0.2.5 |
| `identifiability_one_sided_verdict` | `True` | 单侧给方向 + 强制警示(`False` = 两端 INDETERMINATE) |
| `identifiability_suppress_below_floor` | `True` | 低于下限 → 抑制衰减比 |
| `identifiability_derived_prefixes` | `["score_"]` | 表达量派生评分列(命名启发式) |
| `identifiability_nondesign_cols` | `[]` | 追加的非设计列(QC 列已内置) |

### 7.5 稳健性(阈值敏感性,实测)

- `matched_min_cells_warn` 取 10→∞:判定只在"设计不可识别"与"可给判定"之间切换,
  且切换点唯一(恰等于共同支持域的细胞数)--该阈值不在其它位置造成跳变。
- `attenuation_min` 落在观测档位区间**之内**时任取一值 → `THRESHOLD-SPANNING`;
  区间之外 → 单侧判定 + 强制警示。
- 因此结论陈述不依赖具体阈值:只要各规则的中位衰减比跨度大于任何合理判据阈值,
  该判定就是"阈值位置决定"而非"数据决定"。

### 7.6 边界(诚实清单)

- 匹配目标是**逐变量**的病例组众数(非供者级联合剖面);
- `full_pp` / `matched_pp` 是**细胞级 pooled 检出率差**(伪重复),与供者级估计量
  不同尺度,只在规则之间可比;
- 衰减比数值在汇总前按 2 位小数取整,故报告的跨度与真值之间存在 ≤0.005 的量化误差;
- 例外:`identifiability` 层与 `composition` 层**并列输出**(见本节开头的关系说明)。
- **真空支持域不得读作"可比"**(2026-10-07 已知答案验证 + 独立复核):当**没有** design-defining
  变量时,共同支持域按构造即整个对照组(vacuous),且任何可评估规则都不缩小它 → 所有衰减比恒
  = 1(by construction)。此时判定为 `INDETERMINATE`,并在 notes 明写
  "not evidence of comparability ... do not read this state as 'COMPARABLE'";**不**给单侧方向判定。
- **支持域 = 整个对照组即不得读作"可比"**(v0.3.3 放宽):触发条件从"真空支持域"扩展为
  "支持域**等于整个对照组**(`vacuous`,或匹配没有把任何对照挡在支持域之外)**且**没有任何可
  评估规则真正缩小支持域"。后者覆盖一个重要情形:数据集里存在**在该数据集内恒定的零信息列**
  (例如全队列同一 `region`),它会被归入 design-defining 桶,使旧的 `vacuous` 守卫失效、衰减比
  恒 = 1 却被报成 `COMPARABLE`;该情形现与真空支持域同判(`INDETERMINATE` + 强制提示)。
- **空合取规则不进入一致性统计**:阶梯里变量列表为空的规则(例如 `design_defining` 为空时的
  "同管线合取")等价于"匹配所有",既不缩域也不承载一致性信息,故**不再生成**此类规则。
- **细胞级注释列不是供者级设计变量**:如 `cell_annotation`(细胞类型/谱系注释)这类列,即使
  `check` 列为候选,也应在审计二显式排除(`identifiability_nondesign_cols`)--对结局下游的
  注释做条件化属 overadjustment。

---

## 8. 设计声明与能力闸门(v0.4.0 · Tier 1)

### 8.1 解决的问题

工具在 v0.3.x 及以前**默认假设**输入是「二元、横断面、分析单位 = 供者」的病例–对照设计。
对纵向 / 治疗臂 / 多状态设计若照旧计算,会**静默给出无效判定**(例如把复测单位当独立观测、
把多臂的对照含义搞错)。而「比较的类型」是**实验设计知识**,数据**无法**唯一决定 ——
因此程序不应去猜,而应**索取**。

### 8.2 设计声明(只记录,不推断)

调用方**显式声明**(CLI:`--comparison-type` / `--unit` / `--contrast-a` / `--contrast-b` /
`--time-col` / `--ref-time` / `--index-time` / `--index-event`;或 config 中的 `design_*` 字段)。
全部默认值 = 支持集 S0::

    S0 = { comparison_type: "two_arm", unit: "donor", time: None, arms: 2 }

声明**只记录、留痕**(写入 `run_metadata.json` 的 `design` 段),不改变任何既有判定。

### 8.3 能力闸门

- 声明 ∈ S0 → 正常执行,判定与 v0.3.x **逐字节一致**。
- 声明合法但 ∉ S0 → **拒绝输出任何判定标签**,`verdict = DESIGN-NOT-SUPPORTED`,退出码 **3**;
  只回显声明 + 支持集 + 补齐方式。**不做审计、不产出判定产物**。
- 声明取值不在封闭词表内(写错) → **错误**,退出码 **2**。

**`DESIGN-NOT-SUPPORTED` ≠ `INDETERMINATE`**:前者是「设计超出能力、根本没做审计」;
后者是「支持的设计下、审计做了但无法给结论」。

### 8.4 结构提示(只警告,不改判定)

数据形似更复杂设计时(例如 condition 水平数 > 2)输出 `[design-hint]` 提示,
**绝不改变任何判定标签** —— 程序只索取,不推断。

### 8.5 当前支持集与路线

S0 是**唯一**实现。纵向(`paired` / `time_course`)、治疗臂(`multi_arm`)、多状态(`multi_state`)
与「因果角色(confounder / collider)推断」**均未实现**;声明即拒绝(Tier 2/3)。
这是有意的取舍:**把复杂度接口化并推迟** —— 先保证「不懂时明确拒绝」,而不是静默算错。

---

## 9. 引用与致谢

方法学依据:
- 稀疏基因不可定量性:scRNA-seq dropout 文献(按你的面板自行引用)
- 伪重复与供者级置换:见 §2.3 与 §3
- 多重比较:Benjamini-Hochberg FDR(面板内)

---

*本工具的实现与其测试套件位于 `sc-audit/`;阈值配置文件为 `configs/cardiac.yaml`。
所有报告中引用的数值均可由本工具的 `--outdir` 产物复算。*

---

## 附录 A · v0.6.0 方法学修订（A/B/C）

依逐层面文献（`文献_逐层面映射_2026-10-12.md`）与代码逐处对照：

- **A1 参照侧深度可比性**（检出率是深度的陡函数：Hicks 2018 PMID 29121214、Svensson 2017
  PMID 28263961）。若参照表提供 `mean_umi`（列或 manifest），与 case 中位 `total_counts` 比对；
  相对差 > `reference_depth_delta_max_pct`（默认 50，事前锁定）→ 默认**拒绝**
  （`--allow-depth-mismatch` 放行）。`--depth-matched-sensitivity` 输出二项稀释到参照深度后的
  R/判定列。**未提供深度时不改变任何输出**。
- **A2 positivity / 协变量重叠**（Petersen 2012 PMID 21030422；Westreich 2010 PMID 20139125）。
  在 design-defining 变量上计算**供者级剖面重叠**：`overlap_case` / `overlap_ctrl` / `overlap_min`。
  `overlap_min < positivity_overlap_min`（默认 0.5）且支持域细胞数**达标**时 → `NOT SEPARABLE`。
  细胞计数退为**容量前提**（不达标 → `SUPPORT-BELOW-FLOOR`），与 positivity 判据**并列而非混同**。
- **B1 潜在批次层**（Squair 2021 PMID 34584091；Luecken 2022 PMID 34949812；Korsunsky 2019
  PMID 31740819）。供者级 pseudobulk → HVG → PCA(2) → k-means(2) → 与 condition 的 ARI；
  `ARI >= latent_batch_ari_warn`（默认 0.5）→ `LATENT-BATCH-RISK`。**报警层**（`--with-latent-batch`，
  默认关闭），**不得**把簇当作可匹配变量。
- **B2 管家护栏**：新增各管家 R 的 CV 告警（`housekeeping_cv_warn`）与
  “高基线组织适用域”声明（依据 `R_S2=(1-sqrt(1-p))/p` 的饱和性；geNorm-style M-value 传统）。
- **B3 分层口径敏感性**：额外用 `total_counts` 分层与 `n_genes` 对照（Hafemeister 2019
  PMID 31870423），报告持续性是否口径依赖；**不默改主判定**。
- **C1/C2/C3/C4/C5**：ambient RNA（SoupX/DecontX/Janssen 2023）、双细胞（Scrublet/DoubletFinder）、
  注释方法依赖（Abdelaal 2019）、`binom_P_case0` 细胞级（伪重复）口径、小面板 FDR 功效 —— 均显式写入输出。

所有新阈值均为**事前锁定**值，写于 `config.py` + `configs/cardiac.yaml`（全字段一致测试强制）。


---
