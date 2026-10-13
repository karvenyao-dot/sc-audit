# 审计二诊断指南与方法学措辞

**版本：** 2026-10-08（对应 sc-audit v0.8.4；本页示例输出为 v0.3.3 时点）
**适用：** `sc-audit` v0.8.4 —— 子命令 `check` / `floor` / `composition` / `identifiability` / `robust` / `identity` / `compare` / `run` / `design`
**数字来源：** 本页示例取自随包**合成演示数据**（`examples/make_demo_data.py`）的实测输出；不含任何真实研究数据。

---

## 0. 一句话原则

> 审计二评估的可靠性取决于**设计变量的选择**。当设计变量与分组结构性共线时，
> 匹配对比可能输出 `NOT SEPARABLE`——
> **这通常反映数据的真实混杂结构，而不是工具的失败。**

本工具**只透明化信息，不代替决策**：`check` 报告候选与后果，`match_vars` 由使用者写入 config，
该选择随后被记录进 `run_metadata.json`——这就是"事前锁定"属性的来源。

---

## 1. 变量分档（`check` 输出解读）

`check` 用两条规则给每个候选定档：

```
供者级 TVD 与两档阈值（均为 config 中可事前锁定的字段）比较：
  TVD < match_var_tvd_low (0.1)          → safe_but_uninformative
  0.1 ≤ TVD ≤ hetero_max (0.5)           → primary_match
  TVD > hetero_max (0.5)                 → risky
角色为 outcome / intervention             → never / avoid（无论 TVD）
与供者一一对应（donor-unique）            → not_usable
```

### 演示数据实测输出（合成）

| 变量 | 角色 | 供者级 TVD | 档位 |
|---|---|---|---|
| institution | 取材/来源 | 0.833 | risky（> `hetero_max` 0.5） |
| donation_type | 取材/来源 | 0.833 | risky |
| source | 取材/来源 | 0.833 | risky |

（演示数据里三个来源变量都与分组接近分离：6 个对照供者中只有 1 个与病例同机构。
角色为 outcome / intervention 的列会被直接标 `never` / `avoid`，**与 TVD 无关**——
真实数据里常见结局变量与治疗后变量，均不得用于匹配。）

---

## 2. 三套可选策略（工具不排序，只列权衡）

### 策略 A · conservative —— 只用推荐变量匹配

```yaml
match_vars: [source, gender]
```

| | |
|---|---|
| **收益** | 保留可定量细胞数；匹配对比在 source 上进行，合法且可解释 |
| **代价** | 无法控制被排除变量的混杂，须在方法学描述（Limitations）中写明 |
| **适用** | 你判断被排除的变量**不代表**病例/对照的生物学差异 |

### 策略 B · primary_plus_sensitivity —— 主分析 + 分层敏感性分析

```yaml
match_vars: [source, gender]
strata_cols: [donation_type, institution, region, age_range]
```

| | |
|---|---|
| **收益** | 兼顾混杂控制与生物学完整性；结论对分层稳健性可被检验 |
| **代价** | 需报告两组结果；分层后每层细胞数下降，可能触及 `min_cells_per_stratum` |
| **适用** | 你对被排除变量是否代表生物学差异**不确定**（最常见的情形） |

### 策略 C · no_matching —— 放弃匹配，以供者级 TVD 为可比性依据

```yaml
run_matched_comparison: false
```

| | |
|---|---|
| **收益** | 不做无法自洽的匹配；审计二以 TVD 判定，结论不依赖匹配 |
| **代价** | 必须在 Limitations 声明"未做匹配" |
| **适用** | 混杂是**结构性**的（机构/取材与分组几乎完全共线） |

> ⚠️ **策略 C 不是"退而求其次"。** 当混杂是结构性的，
> 强行匹配只会把真实的病例/对照差异一并匹配掉——此时放弃匹配是最诚实的做法。

---

## 3. 方法学措辞（可直接使用）

### 若采用策略 A

> Control and case groups differed in procurement characteristics
> (donor-level total-variation distance, TVD, up to <TVD_high> for <var_high>).
> To control for the imbalance that could be addressed without removing the
> biological contrast, we matched controls to cases on <var_match>
> (donor-level TVD = <TVD_match>). Residual imbalance in <var_residual>
> (TVD <TVD_residual_range>) could not be matched without also
> matching away the case/control contrast and is reported as a limitation.

（把 `<...>` 处换成**你自己数据**的 `check` 输出值；下文两段同理。）

### 若采用策略 B

> We adopted a two-tier strategy. In the primary analysis, controls were matched
> to cases on <var_match> (donor-level TVD = <TVD_match>). In a prespecified
> sensitivity analysis, the cohort was stratified by <var_strata>, each of which
> was strongly associated with group status (donor-level TVD <TVD_residual_range>)
> and could not be matched without removing the biological contrast. Concordance
> between the primary and stratified analyses is reported in Supplementary Table X.

### 若采用策略 C

> Control and case groups differed systematically in procurement variables
> (donor-level TVD up to <TVD_high>). Because no candidate variable could be matched
> without also removing the biological contrast, we did not perform a matched
> comparison; comparability was instead assessed through the donor-level
> total-variation distance (pre-specified margin 0.5), and the imbalance is
> reported as a limitation.

---

## 4. Limitations 段落（可直接使用）

> **Confounding by procurement characteristics.**
> The composition audit evaluates control-case comparability using annotation
> columns available in the dataset (`institution`, `donation_type`, `source`,
> `region`). Where these variables are strongly collinear with group status
> — i.e. when donor-level TVD approaches complete separation —
> a fully matched comparison is not identifiable: matching on such a variable
> removes the case/control contrast itself. We therefore report the resulting
> `NOT SEPARABLE` verdict as a property of the study design rather than as a
> technical failure, and we recommend one of three responses:
> (i) match only on informative, non-collinear variables;
> (ii) report a prespecified stratified sensitivity analysis; or
> (iii) forgo matching and base comparability on the donor-level TVD.
> Datasets lacking these annotation columns entirely yield an `INDETERMINATE`
> verdict; this reflects the state of public metadata rather than a limitation
> of the method, and we quantify annotation coverage in the readiness report.

---

## 5. 关于变量角色的边界（务必阅读）

`check` 的变量角色由**两条规则**判定，**两条都有已知失效模式**：

```
① 列名词根       失效：列名不含常见词根（如 outcome_X、comp_30d）
② 是否只在单一组内变化   失效：**无法区分**结局变量与分组不平衡的取材变量
                        （对照跨 2 个中心 / 病例全在 1 个中心，
                          与"对照全 no / 病例有 yes/no"在数据上同形）
```

因此：
- **统计信号不覆盖角色判定**，只产出 `condition_dependent` 提示（标为「可用·需甄别」）
- 启发式认不出时，用 config **显式声明**（优先级最高）：

```yaml
role_overrides: {outcome_X: outcome, comp_30d: outcome,
                 on_immunosuppressant: treatment, mystery_col: ignore}
extra_outcome_keywords: [troponin_peak]      # 领域词根扩展
extra_treatment_keywords: [colchicine]
```

- **运行期保护**：把结局/治疗变量写进 `match_vars` 会**中止**运行（防止 collider bias）。
  确认语义无误后可用 `--allow-outcome-matching` 显式放行；开关本身也写在 config
  （`forbid_outcome_matching`），以便事前锁定。

---

## 6. 已知未完成项（诚实清单）

| 项 | 状态 |
|---|---|
| 参照表覆盖 | 当前 6 背景 / 17 基因；**未登记生成参照表所用的门控参数**，故 R 的跨数据集可比性无法核验（已知正确性缺口） |
| 科学验证 | 判定阈值仅经**代码审计**，未经独立科学验证 → 应视为"高度自洽的专家启发式" |
| `check` 的角色启发式 | 有已知假阳性/假阴性模式（见 §5），已标注为低置信度 |
| 深度分层容量 | 若仅少数分箱满足两组各 ≥30 细胞，审计三分层结论需按此设限 |
