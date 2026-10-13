"""审计二 · v0.3 可识别性层（identifiability + matching ladder）。

为什么单独成层（2026-10-04 / 2026-10-06 讨论的直接结论）
------------------------------------------------------
原审计二（`composition.py`）把**用户手写的 `match_vars`** 当作"混杂变量定义"，
由此产生三个问题：

  1. 病例组内**恒定**的设计变量（如 `institution` / `region`）被当作混杂去匹配，
     匹配必然坍缩到 n≈1 个供者，得到的"衰减比"由单个供者决定；
  2. 病例组内**有变异**的变量（如 `age_range` / `gender`）反而被排除在判定之外；
  3. 单条匹配规则下的一个中位衰减数被当作"混杂解释了多少"的**定量证据**，
     而它对匹配规则高度敏感（相同规则集换一条匹配变量即可翻向）。

本层把这三件事分别替换为确定性规则：

  · **花名册完备性**：obs 的每个非 QC 列都必须被归入一个桶，未归类即报错
    （杜绝"静默遗漏一个未知混杂列"）；
  · **变量分桶**：病例组内供者级熵 ≈ 0 → `design-defining(scope)`
    （只能在适用范围声明里出现，不可匹配）；有变异 → `adjustable(confounder)`
    （真正可调整的混杂，必须进混杂集）；
    **两臂各自恒定且水平 ≥ 2**（= condition 的一对一重编码）→
    `excluded:redundant-with-condition`（v0.7.1，见下）；
    **confounder 角色的列每供者 1 值且取值互异**（如供者数少时的 age）→
    `excluded:donor-level-all-distinct`（v0.7.3；不可匹配，但不是标识符，归因不同）。
  · **先判可识别性（共同支持域），再谈匹配**；
  · **匹配规则改成事前锁定的阶梯**（R0–R5），逐级报衰减比并做**一致性检查**；
  · **低于事前锁定下限 → 抑制衰减数字**：不出任何量化证据。

判定标签（`verdict`）
--------------------
  · `SUPPORT-BELOW-FLOOR`  共同支持域低于事前锁定下限 → 设计不可识别，不报衰减比
  · `THRESHOLD-SPANNING`   中位衰减比跨预设阈值（阶梯上两侧并存）→ 不构成单一估计量
  · `NOT SEPARABLE`        阶梯上中位衰减比**一致地低于**衰减阈值（单侧；必附阈值位置警示）
  · `COMPARABLE`           阶梯上中位衰减比**一致地不低于**衰减阈值（单侧；必附同上警示）
  · `INDETERMINATE`        无法评估（无可用规则 / 无目标基因 / 多病例设计剖面）

#### 2026-10-06 审计后的修订（关键）

第一版实现有一个结构性缺陷：一致性检查的"两档并存"与"跨规则翻转"是**同一个谓词**，
于是两档并存 ⇒ 必判 `RULE-DEPENDENT`，单侧 ⇒ 被保护①压为 `INDETERMINATE`，
**`COMPARABLE` / `NOT SEPARABLE` 在默认配置下永不可达**。已修订为：

  · 两档并存（阈值落在观测区间内）→ `THRESHOLD-SPANNING`（无法给单一结论）
  · 单侧（阈值落在观测区间之外）→ **直接给出方向判定**（`NOT SEPARABLE` / `COMPARABLE`），
    但**强制附带警示**：此时一致性检查携带零信息（所有规则按构造落在同一侧），
    该判定仅基于量级比较。
  · 开关：`identifiability_one_sided_verdict`（默认 True = 上面的甲方案；
    False = 两端均报 `INDETERMINATE` 的旧保守行为）

#### 保护①（2026-10-06 首次并入；同日审计后改为上述形式）

旧写法 `consistent = (n_low == 0 or n_high == 0)` 有边界伪影：阈值低于全部观测档时
`n_low = 0` 使检查"通过"，再拿中位衰减比去比阈值 → 输出 `COMPARABLE`，等于取消衰减判据。
**现在不再用"拒绝回答"来堵这个洞，而是用强制警示 + 明写一致性检查无效**：
既保住了单侧情形的方向信息，也不让阈值位置被当成证据。

只读属性：本模块不写入任何配置；阶梯规则一旦由 `cfg.identifiability_ladder`
显式给出，即视为运行前锁定（与 `run_metadata.json` 一起留痕）。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .check import _ROLE_ZH, _name_role
from .composition import _group_distribution_distance
from .config import AuditConfig, DEFAULT_NONDESIGN_QC_COLS
from .donor_utils import (chance_constancy_risk, donor_group, donor_level_values,
                          entropy as _entropy, is_donor_unique)
from .io import assert_donors_nested, gene_detection_vector

# ══════════════════ 常量 ══════════════════

# 细胞级 QC / 分组 / 供者 / 细胞身份注释：不是设计变量，不进花名册
_QC_COLS: Tuple[str, ...] = DEFAULT_NONDESIGN_QC_COLS   # 单一真值源（见 config）

BUCKET_DONOR_UNIQUE = "excluded:donor-unique"
# C-bug（v0.7.3）：供者级**协变量**（role=confounder，如 age/sex）在供者数少时会出现
# "每供者取值互异" → 也满足 donor-unique 判据，但它**不是标识符**。二者对匹配的处置相同
# （均不可匹配），但归因不同：标成 donor-unique 会掩盖"存在真实供者级协变量、只是无法平衡"。
BUCKET_DONOR_SATURATED = "excluded:donor-level-all-distinct"
BUCKET_DERIVED = "excluded:derived-from-expression"
BUCKET_OUTCOME = "excluded:outcome(collider)"
BUCKET_INTERVENTION = "excluded:post-baseline"
BUCKET_TECHNICAL = "excluded:technical(nuisance)"
BUCKET_ALL_MISSING = "excluded:all-missing"
# C1（v0.7.1，据第三方独立评估）：与 condition **一对一重编码**的"标签副本"列
# （病例组、对照组各自恒定，且取值不同 → 完美分离）不得当作设计/范围变量。
BUCKET_REDUNDANT_CONDITION = "excluded:redundant-with-condition"
BUCKET_DESIGN_DEFINING = "design-defining(scope)"
BUCKET_ADJUSTABLE = "adjustable(confounder)"

_EXCLUDED_PREFIX = "excluded:"

# ── 判定标签（2026-10-06 修订）──
# 旧名 → 新名：`NON-OVERLAPPING` → `SUPPORT-BELOW-FLOOR`（支持域"小而非空"，旧名是误名）；
#             `RULE-DEPENDENT` → `THRESHOLD-SPANNING`（描述的是中位衰减比跨阈值，而非"规则错"）。
VERDICT_SUPPORT_BELOW_FLOOR = "SUPPORT-BELOW-FLOOR"
VERDICT_THRESHOLD_SPANNING = "THRESHOLD-SPANNING"
VERDICT_NOT_SEPARABLE = "NOT SEPARABLE"
VERDICT_COMPARABLE = "COMPARABLE"
VERDICT_INDETERMINATE = "INDETERMINATE"


# ══════════════════ 基础工具 ══════════════════

def _donor_mode(s: pd.Series) -> str:
    """（兼容包装）见 `donor_utils.donor_mode`。"""
    from .donor_utils import donor_mode
    return donor_mode(s)


def _donor_level(o: pd.DataFrame, col: str, donor: pd.Series) -> pd.Series:
    """（兼容包装）见 `donor_utils.donor_level_values`。"""
    return donor_level_values(o, col, donor)


def _group_of_donor(is_ctrl: np.ndarray, donor: pd.Series) -> pd.Series:
    return donor_group(is_ctrl, donor)


def _nondesign_columns(ad, condition_col: str, donor_col: str,
                       cfg: AuditConfig) -> set:
    s = set(_QC_COLS)
    s |= set(getattr(cfg, "identifiability_nondesign_cols", []) or [])
    s |= {condition_col, donor_col}
    ident = getattr(cfg, "identity_col", None)
    if ident:
        s.add(str(ident))
    return {c for c in s if c in ad.obs.columns}


def _is_derived(col: str, cfg: AuditConfig) -> bool:
    if col in set(getattr(cfg, "identifiability_derived_cols", []) or []):
        return True
    low = col.lower()
    return any(low.startswith(p.lower())
               for p in (getattr(cfg, "identifiability_derived_prefixes", []) or []))


def redundant_with_condition(ent_case: float, ent_ctrl: float, levels: int,
                             eps: float = 1e-9) -> bool:
    """condition 副本列的**单一判据**（v0.8.2）：两臂各自恒定（供者级熵≈0）且水平≥2。

    `_bucket_of` 与 `design_parse._reverify` **共用本函数**，消除两套实现的判据分叉
    （此前 roster 用整体水平数 `levels`、reverify 用 `donor_level_levels`）。
    """
    return bool(ent_case <= eps and ent_ctrl <= eps and int(levels) >= 2)


def _bucket_of(row: pd.Series, cfg: AuditConfig) -> str:
    """确定性分桶。顺序即优先级（与 README / docs/methods.md 逐条对应）。"""
    if row["donor_unique"]:
        # C-bug（v0.7.3）：confounder 角色的列不归入 donor-unique（那不是标识符，
        # 而是"供者级协变量但每供者 1 值 → 不可平衡"）。判据只用**列名角色**，不靠统计。
        if row["role"] == "confounder":
            return BUCKET_DONOR_SATURATED
        return BUCKET_DONOR_UNIQUE
    if row["all_missing"]:
        return BUCKET_ALL_MISSING
    if _is_derived(row["column"], cfg):
        return BUCKET_DERIVED
    if row["role"] == "outcome":
        return BUCKET_OUTCOME
    if row["role"] == "intervention":
        return BUCKET_INTERVENTION
    if row["role"] == "technical":
        return BUCKET_TECHNICAL
    # C1（v0.7.1）：condition 的"标签副本"——在**两个臂内都恒定**（donor 级熵≈0）
    # 且水平数 ≥ 2（两臂取值不同）→ 该列是 condition 的一对一重编码。把它当
    # design-defining 会在共同支持域里条件化在**待估对比本身**上：病例剖面的取值
    # 在对照里不存在 → 支持域恒为 0 → 审计二b 永判 SUPPORT-BELOW-FLOOR（假拒绝）。
    # 判据用"两臂各自恒定 + 水平≥2"（= 完美分离），与列名无关（不靠命名启发式）。
    if redundant_with_condition(row["ent_case"], row["ent_ctrl"], row["levels"]):
        return BUCKET_REDUNDANT_CONDITION
    if row["ent_case"] <= 1e-9:
        return BUCKET_DESIGN_DEFINING
    return BUCKET_ADJUSTABLE


# ══════════════════ ① 变量花名册 ══════════════════

def build_variable_roster(ad, condition_col: str, donor_col: str,
                          control_label: str, cfg: Optional[AuditConfig] = None
                          ) -> Tuple[pd.DataFrame, List[str]]:
    """obs 每个非 QC 列 → 一行花名册；返回 (roster, unclassified)。

    列名作 `columns` 返回的**完备性**由调用方断言：
    `len(roster) + len(nondesign) + len(unclassified) == n_obs_columns`。
    """
    cfg = cfg or AuditConfig()
    o = ad.obs
    cond = o[condition_col].astype(str)
    is_ctrl = (cond == str(control_label)).values
    donor = o[donor_col].astype(str)
    n_don = donor.nunique()
    nondesign = _nondesign_columns(ad, condition_col, donor_col, cfg)

    rows, unclassified = [], []
    for c in o.columns:
        if c in nondesign:
            continue
        try:
            s = o[c].astype(str)
            tvd_d, _cv, _mode, _tbl = _group_distribution_distance(
                s.values, is_ctrl, donor=donor.values)
            dv = _donor_level(o, c, donor)
            dgrp = _group_of_donor(is_ctrl, donor).reindex(dv.index)
            per_donor = s.groupby(donor, observed=True).nunique()
            donor_unique = is_donor_unique(s, donor, n_don)
            rows.append(dict(
                column=c, role=_name_role(c, cfg),
                role_zh=_ROLE_ZH.get(_name_role(c, cfg), _name_role(c, cfg)),
                tvd_donor=round(float(tvd_d), 3),
                ent_case=round(_entropy(dv[dgrp == "case"]), 3),
                ent_ctrl=round(_entropy(dv[dgrp == "ctrl"]), 3),
                donor_unique=donor_unique,
                all_missing=bool(o[c].isna().all()),
                levels=int(s.nunique())))
        except Exception as exc:  # noqa: BLE001 —— 未归类必须显式浮出，不得静默跳过
            unclassified.append(f"{c} ({type(exc).__name__}: {exc})")
    R = pd.DataFrame(rows)
    if len(R):
        R["bucket"] = R.apply(lambda r: _bucket_of(r, cfg), axis=1)
    return R, unclassified


# ══════════════════ ② 共同支持域 ══════════════════

def common_support(ad, design_defining: Sequence[str], condition_col: str,
                   donor_col: str, control_label: str, cfg: Optional[AuditConfig] = None
                   ) -> Dict[str, object]:
    """在 design-defining 变量上取病例组剖面，返回落在该剖面的对照供者与细胞数。

    病例组在这些变量上恒定（这是"design-defining"的定义），故剖面通常只有一个；
    若病例组出现多个剖面（阈值边界情形），则取**全部**病例剖面的并集。
    """
    cfg = cfg or AuditConfig()
    o = ad.obs
    cond = o[condition_col].astype(str)
    is_ctrl = (cond == str(control_label)).values
    donor = o[donor_col].astype(str)
    cell = o.groupby(donor, observed=True).size()

    if not design_defining:
        ctrl_donors = sorted(set(donor[is_ctrl]))
        return dict(case_profiles=["<none: no design-defining variable>"],
                    support_donors=ctrl_donors,
                    support_cells=int(cell.reindex(ctrl_donors).sum()),
                    n_control_donors=len(ctrl_donors), vacuous=True)

    per = pd.DataFrame({c: _donor_level(o, c, donor) for c in design_defining})
    per["grp"] = _group_of_donor(is_ctrl, donor).reindex(per.index)
    key = per[list(design_defining)].agg("|".join, axis=1)
    case_prof = sorted(set(key[per.grp == "case"]))
    by_profile = {p: sorted(per.index[(per.grp == "ctrl") & (key == p)]) for p in case_prof}
    # A1-F6（2026-10-06 审计）：病例组跨多个设计剖面时，**不得**把各剖面的对照并成一个
    # "共同支持域"（那是跨层比较，不是标准的 positivity 检查）。此处取**单一最大剖面**
    # 作为报告值；判定层会直接报 INDETERMINATE，不给出池化估计。
    if len(case_prof) > 1:
        prof = max(case_prof, key=lambda p: cell.reindex(by_profile[p]).sum())
        support = by_profile[prof]
    elif case_prof:
        prof = case_prof[0]
        support = by_profile[prof]
    else:
        prof, support = None, []
    return dict(case_profiles=case_prof,
                support_donors=sorted(support),
                support_cells=int(cell.reindex(support).sum()),
                support_by_profile={k: v for k, v in by_profile.items()},
                reported_profile=prof,
                n_control_donors=int(len(set(donor[is_ctrl]))), vacuous=False)


def positivity_overlap(ad, design_defining: Sequence[str], condition_col: str,
                       donor_col: str, control_label: str,
                       adjustable: Sequence[str] = ()) -> Dict[str, object]:
    """A2（v0.6.0）：供体级协变量重叠／positivity 诊断。

    在 design-defining 变量构成的**供体级剖面**上：
      overlap_case = 病例供体中，其剖面也出现在对照组的比例；
      overlap_ctrl = 对照供体中，其剖面也出现在病例组的比例；
      overlap_min  = min(两者)。空 design_defining → None（不可判定）。

    动机（Petersen 2012 PMID 21030422；Westreich 2010 PMID 20139125）：positivity 定义在
    **协变量空间**，不是细胞计数。细胞足够多 ≠ 协变量重叠好。
    """
    o = ad.obs
    cond = o[condition_col].astype(str)
    is_ctrl = (cond == str(control_label)).values
    donor = o[donor_col].astype(str)
    # M1（v0.6.2）：剖面 = design_defining ∪ adjustable（可调整混杂也在协变量空间里）
    prof_cols = [c for c in list(design_defining) + list(adjustable)
                 if c in o.columns and c != str(donor_col)]
    if not prof_cols:
        return dict(overlap_case=None, overlap_ctrl=None, overlap_min=None,
                    n_case_donors=0, n_ctrl_donors=0, profile_vars=[])
    per = pd.DataFrame({c: _donor_level(o, c, donor) for c in prof_cols})
    per["grp"] = _group_of_donor(is_ctrl, donor).reindex(per.index)
    key = per[prof_cols].agg("|".join, axis=1)
    ck = key[per.grp == "case"]
    tk = key[per.grp == "ctrl"]
    cs, ts = set(ck), set(tk)
    ov_case = (float(ck.isin(ts).mean()) if len(ck) else None)
    ov_ctrl = (float(tk.isin(cs).mean()) if len(tk) else None)
    both = [v for v in (ov_case, ov_ctrl) if v is not None]
    return dict(overlap_case=(None if ov_case is None else round(ov_case, 3)),
                overlap_ctrl=(None if ov_ctrl is None else round(ov_ctrl, 3)),
                overlap_min=(None if not both else round(float(min(both)), 3)),
                n_case_donors=int(len(ck)), n_ctrl_donors=int(len(tk)),
                profile_vars=list(prof_cols))


def adjustable_balance(ad, adjustable: Sequence[str], support_donors: Sequence[str],
                       condition_col: str, donor_col: str, control_label: str
                       ) -> Dict[str, pd.DataFrame]:
    """支持域内每个 adjustable 变量按层平衡表：某层病例有、支持域对照为 0 → 不可平衡。"""
    o = ad.obs
    cond = o[condition_col].astype(str)
    is_ctrl = (cond == str(control_label)).values
    donor = o[donor_col].astype(str)
    grp = _group_of_donor(is_ctrl, donor)
    out: Dict[str, pd.DataFrame] = {}
    for c in adjustable:
        dv = _donor_level(o, c, donor)
        t = pd.DataFrame({"case": dv[grp == "case"].value_counts()})
        t["ctrl_support"] = dv.reindex(list(support_donors)).value_counts()
        t["ctrl_support"] = t["ctrl_support"].fillna(0).astype(int)
        t["case"] = t["case"].fillna(0).astype(int)
        t["case_pct"] = (t["case"] / max(1, t["case"].sum()) * 100).round(1)
        t["unbalanceable"] = (t["case"] > 0) & (t["ctrl_support"] == 0)
        out[c] = t
    return out


# ══════════════════ ③ 匹配规则阶梯 ══════════════════

def default_ladder(design_defining: Sequence[str], adjustable: Sequence[str],
                   in_window: Sequence[str]) -> List[Tuple[str, List[str]]]:
    """内置事前锁定阶梯（与 2026-10-04 原型逐条一致）。

    R0 不匹配（锚点）→ R1 现行正文口径 → R2 全部候选变量合取
    → R3 design-defining 合取 → R4 单变量（design-defining + adjustable）
    → R5 推荐窗内变量合取。
    """
    # C.4（2026-10-07）：**空合取规则**（变量列表为空 = 匹配所有）不是"不匹配的锚点"，
    # 它既不缩域也不承载一致性信息，却会被计入 medians 并稳定贡献一个 1.000。
    # 因此当 design_defining 为空时，不再生成 R3；R2 也随之退化为"仅 adjustable 合取"。
    ladder: List[Tuple[str, List[str]]] = [
        ("R0 不匹配（锚点）", []),
        ("R1 现行口径 institution∧donation_type∧source∧region",
         ["institution", "donation_type", "source", "region"]),
    ]
    if list(design_defining) or list(adjustable):
        ladder.append(("R2 全部候选变量合取", list(design_defining) + list(adjustable)))
    if list(design_defining):
        ladder.append(("R3 同管线 design-defining 合取 (%s)" % "∧".join(design_defining),
                       list(design_defining)))
    for c in list(design_defining) + list(adjustable):
        ladder.append((f"R4 单变量 {c}", [c]))
    if in_window:
        ladder.append((f"R5 推荐窗内变量合取 ({'∧'.join(in_window)})", list(in_window)))
    return ladder


def _match_donors(o: pd.DataFrame, is_ctrl: np.ndarray, donor: pd.Series,
                  vars_: Sequence[str]) -> List[str]:
    """保留"主取值等于病例组众数（逐变量）"的对照供者。

    语义说明（已知近似）：目标值是**逐变量**的病例组众数，不是供者级联合剖面。
    因此 R2 在共线数据上可能直接得到 0 个供者（这本身即可识别性失败的另一种表述）。

    P0-1（2026-10-06 审计）：`vars_` 里的列**可能不存在于 obs**（内置阶梯的 R1 引用
    institution/donation_type/source/region，而真实数据集常缺 region/source）——
    此处先过滤到实有列，避免 KeyError（调用方另有 notes 记录被丢弃的规则变量）。
    """
    avail = [c for c in vars_ if c in o.columns]
    dropped = [c for c in vars_ if c not in o.columns]
    if vars_ and not avail:
        # 规则声明的变量在本数据集里一个都不存在 → 该规则不可评估（0 供者），
        # 而不是静默退化成"不匹配"（后者会把 ratio 写成 1.0，误导读者）
        return []
    if dropped:
        # 部分缺失：用实有列继续（并在调用方的 notes 里留痕）
        vars_ = avail
    if not vars_:
        return sorted(set(donor[is_ctrl]))
    case = ~is_ctrl
    targets: Dict[str, str] = {}
    for c in vars_:
        md = o.loc[case, c].astype(str).mode()
        if len(md):
            targets[c] = str(md.iloc[0])
    kept: List[str] = []
    for d in sorted(set(donor[is_ctrl])):
        sub = o.loc[(donor == d).values & is_ctrl]
        if not len(sub):
            continue
        ok = True
        for c, t in targets.items():
            md = sub[c].astype(str).mode()
            if len(md) and str(md.iloc[0]) != t:
                ok = False
                break
        if ok:
            kept.append(str(d))
    return kept


def build_matching_ladder(ad, genes: Sequence[str], ladder: Sequence[Tuple[str, List[str]]],
                          condition_col: str, donor_col: str, control_label: str,
                          cfg: Optional[AuditConfig] = None
                          ) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """逐规则 × 逐基因算 full / matched 效应量与衰减比。

    返回 (long 表, 规则级汇总表, 参与汇总的衰减比向量)。

    ⚠️ ratio 的口径（外部审查 B6）：分子分母都取**绝对值**
    （`ratio = |matched_pp| / |full_pp|`），因此它是**效应量比值**、对**方向不敏感**：
    若匹配后组间方向反号（病例原本更高、匹配后更低），ratio 仍是正数（如 0.4），
    不会出现负数。方向信息保存在 `full_pp` 与 `matched_ctrl_pct` 两列中（可自行比对），
    需要方向时可据此还原；不要把 ratio 当作带符号的效应量变化。
    """
    cfg = cfg or AuditConfig()
    o = ad.obs
    cond = o[condition_col].astype(str)
    is_ctrl = (cond == str(control_label)).values
    donor = o[donor_col].astype(str)
    cell = o.groupby(donor, observed=True).size()

    panel = [g for g in dict.fromkeys(genes)]          # P1-5：去重保序
    vec = {}
    for g in panel:
        v = gene_detection_vector(ad, g)
        if v is not None:
            vec[g] = v
    panel = [g for g in panel if g in vec]
    # P0-2：目标基因全部不在 var_names / 阶梯为空时，不得构造"无列"的空 DataFrame
    #（否则下游 `Lad.rule` 抛 AttributeError，并连带 run 崩溃、report 不产出）
    LONG_COLS = ["rule", "n_donors", "n_cells", "gene", "full_pp", "matched_pp",
                 "matched_ctrl_pct", "ratio"]
    SUM_COLS = ["rule", "n_donors", "n_cells", "n_usable", "n_amplified",
                "median_ratio"]
    if not panel or not list(ladder):
        return (pd.DataFrame(columns=LONG_COLS), pd.DataFrame(columns=SUM_COLS),
                np.array([], dtype=float))
    case_pct = {g: float((vec[g][~is_ctrl] > 0).mean() * 100) for g in panel}
    base_ctrl = {g: float((vec[g][is_ctrl] > 0).mean() * 100) for g in panel}

    rows = []
    for label, vars_ in ladder:
        matched = _match_donors(o, is_ctrl, donor, vars_)
        m_cells = int(cell.reindex(matched).sum()) if matched else 0
        if vars_:
            m_mask = is_ctrl & donor.isin(matched).values
        else:
            m_mask = is_ctrl.copy()
        for g in panel:
            ctrl_pct = (float((vec[g][m_mask] > 0).mean() * 100)
                        if int(m_mask.sum()) else np.nan)
            full = abs(case_pct[g] - base_ctrl[g])
            matched_pp = abs(case_pct[g] - ctrl_pct) if ctrl_pct == ctrl_pct else np.nan
            rows.append(dict(rule=label, n_donors=len(matched), n_cells=m_cells,
                             gene=g, full_pp=round(full, 2),
                             matched_pp=(round(matched_pp, 2) if matched_pp == matched_pp
                                         else np.nan),
                             matched_ctrl_pct=(round(ctrl_pct, 1)
                                               if ctrl_pct == ctrl_pct else np.nan),
                             ratio=(round(matched_pp / full, 2)
                                    if (full > 0 and matched_pp == matched_pp) else np.nan)))
    Lad = pd.DataFrame(rows)

    summ = []
    for label, _ in ladder:
        s = Lad[Lad.rule == label]
        usable = s.loc[s.full_pp >= cfg.attenuation_min_full_effect_pp, "ratio"].dropna()
        summ.append(dict(rule=label,
                         n_donors=int(s.n_donors.iloc[0]) if len(s) else 0,
                         n_cells=int(s.n_cells.iloc[0]) if len(s) else 0,
                         n_usable=int(len(usable)),
                         n_amplified=int((usable > 1).sum()),
                         median_ratio=float(usable.median()) if len(usable) else np.nan))
    S = pd.DataFrame(summ)
    meds = S.loc[~S.rule.str.startswith("R0") & (S.median_ratio == S.median_ratio),
                 "median_ratio"].to_numpy(dtype=float)
    return Lad, S, meds


# ══════════════════ 结果容器 ══════════════════

@dataclass
class IdentifiabilityResult:
    # ① 花名册
    roster: pd.DataFrame = field(default_factory=pd.DataFrame)
    unclassified: List[str] = field(default_factory=list)
    n_obs_columns: int = 0
    n_nondesign_columns: int = 0
    design_defining: List[str] = field(default_factory=list)
    adjustable: List[str] = field(default_factory=list)
    excluded: Dict[str, List[str]] = field(default_factory=dict)
    # ② 支持域
    case_profiles: List[str] = field(default_factory=list)
    support_donors: List[str] = field(default_factory=list)
    support_cells: int = 0
    n_control_donors: int = 0
    balance: Dict[str, pd.DataFrame] = field(default_factory=dict)
    unbalanceable: List[str] = field(default_factory=list)
    # ③ 阶梯
    ladder: pd.DataFrame = field(default_factory=pd.DataFrame)
    ladder_summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    ladder_rules: List[str] = field(default_factory=list)
    # ④ 判定
    cell_floor: int = 50
    excluded_from_consistency: List[str] = field(default_factory=list)
    verdict: str = ""
    details: Dict[str, object] = field(default_factory=dict)
    attenuated_suppressed: bool = False
    notes: List[str] = field(default_factory=list)

    def to_text(self) -> str:
        L = ["== Audit 2 (v0.3) · Identifiability & matching ladder ==", ""]
        dv = self.details
        L.append(f"cells {dv.get('n_cells')} · donors {dv.get('n_donors')} · "
                 f"control {dv.get('n_control_cells')}/{dv.get('n_control_donors')} · "
                 f"case {dv.get('n_case_cells')}/{dv.get('n_case_donors')}")
        L.append("")
        L.append("-- ① variable roster (every non-QC obs column must be bucketed) --")
        if len(self.roster):
            disp = self.roster.sort_values(["bucket", "tvd_donor"],
                                           ascending=[True, False])[
                ["column", "role_zh", "tvd_donor", "ent_case", "ent_ctrl", "bucket"]]
            L.append(disp.to_string(index=False))
        else:
            L.append("    (no design-variable candidate in obs)")
        L.append(f"    completeness: {self.n_obs_columns} columns = "
                 f"non-design {self.n_nondesign_columns} + bucketed {len(self.roster)} "
                 f"+ unclassified {len(self.unclassified)}")
        if self.unclassified:
            L.append(f"    UNCLASSIFIED (must be resolved before use): {self.unclassified}")
        L.append(f"    design-defining (constant within cases → scope statement, "
                 f"not matchable): {self.design_defining}")
        L.append(f"    adjustable (varies within cases → true confounders, must be in "
                 f"the adjustment set): {self.adjustable}")
        for k, v in self.excluded.items():
            if v:
                L.append(f"    excluded · {k}: {v}")
        L.append("")
        L.append("-- ② common support (identifiability is checked BEFORE matching) --")
        L.append(f"    case design profile(s): {self.case_profiles}")
        L.append(f"    control donors inside that profile: {self.support_donors} → "
                 f"{self.support_cells} cells (of {self.n_control_donors} control donors)")
        _po = (self.details or {}).get("positivity_overlap")
        if _po:
            L.append(f"    donor-level positivity overlap: case {_po.get('overlap_case')} · "
                     f"ctrl {_po.get('overlap_ctrl')} · min {_po.get('overlap_min')} "
                     f"(A2；低于 positivity_overlap_min 阈值 → 强制警示，不覆盖阶梯判定)")
        for c, t in self.balance.items():
            L.append(f"    · {c}: {int(t['unbalanceable'].sum())} case level(s) have no "
                     f"support-domain control donor → "
                     f"{'cannot be balanced' if bool(t['unbalanceable'].any()) else 'balanceable'}")
            for lvl, row in t.iterrows():
                L.append(f"        {str(lvl):<26} case {int(row['case']):>2} "
                         f"({row['case_pct']:>5.1f}%)  ctrl_support {int(row['ctrl_support'])}")
        L.append("")
        L.append("-- ③ prespecified matching ladder --")
        if len(self.ladder_summary):
            L.append(self.ladder_summary.to_string(index=False))
        if len(self.ladder):
            piv = self.ladder.pivot_table(index="rule", columns="gene",
                                          values="ratio", aggfunc="first")
            L.append("")
            L.append("attenuation ratio (matched / full effect), rules × genes "
                     "(R0 = anchor, ratio ≡ 1 by definition):")
            L.append(piv.to_string())
        L.append("")
        L.append("-- ④ verdict --")
        d = self.details
        L.append(f"    pre-specified floor: support cells >= {self.cell_floor} "
                 f"(from matched_min_cells_warn)")
        L.append(f"    observed support cells = {self.support_cells}")
        if self.excluded_from_consistency:
            L.append(f"    rules not evaluable (no matched control cells): "
                     f"{self.excluded_from_consistency}")
        meds = d.get("medians", np.array([]))
        if len(meds):
            L.append(f"    ladder medians: {', '.join(f'{m:.3f}' for m in meds)} "
                     f"(median {float(np.median(meds)):.3f}, spread "
                     f"{float(meds.max() - meds.min()):.3f})")
            L.append(f"    below attenuation_min={d.get('attenuation_min')}: "
                     f"{d.get('n_low')} rule(s); at/above: {d.get('n_high')} rule(s)")
        L.append(f"    VERDICT: {self.verdict}")
        if self.attenuated_suppressed:
            L.append("    * attenuation ratios are SUPPRESSED: any matched subset is "
                     "decided by too few donors to constitute evidence of the "
                     "confounding share.")
        for n in self.notes:
            L.append(f"  · {n}")
        return "\n".join(L)


# ══════════════════ 主函数 ══════════════════

def audit_identifiability(ad, genes: Sequence[str], condition_col: str, donor_col: str,
                          control_label: str, cfg: Optional[AuditConfig] = None,
                          ladder: Optional[Sequence[Sequence[str]]] = None
                          ) -> IdentifiabilityResult:
    """审计二 · v0.3 可识别性层。

    参数
    ----
    ladder : 事前锁定的匹配规则阶梯。`None` → 使用内置阶梯（R0–R5）。
             显式给出时为"变量名列表的列表"，空列表代表"不匹配（锚点）"。
             一旦显式给出，即视为运行前锁定，应记入 `run_metadata.json`。
    """
    from .io import warn_if_gene_names_unresolved as _warn_genes
    _warn_genes(ad, genes, context="audit_identifiability")
    cfg = cfg or AuditConfig()
    o = ad.obs
    cond = o[condition_col].astype(str).values
    donor = o[donor_col].astype(str)
    is_ctrl = cond == str(control_label)
    assert_donors_nested(donor.values, cond, donor_col, condition_col)

    res = IdentifiabilityResult()
    # P2-10（审计 A2-F9）：condition 只有一个水平时，病例/对照之一为空 →
    # 后续所有比较都无意义（且会产生 NaN 传播与错误归因），此处显式早退。
    if not is_ctrl.any() or is_ctrl.all():
        n_ctrl_ = int(is_ctrl.sum())
        res.verdict = (f"{VERDICT_INDETERMINATE} — the condition column does not contain "
                       f"both groups (control cells = {n_ctrl_}, case cells = "
                       f"{int((~is_ctrl).sum())}); no between-group structure to audit")
        res.notes.append(
            f"check --control-label / --condition-col: control label '{control_label}' "
            f"matched {n_ctrl_} of {o.shape[0]} cells.")
        res.n_obs_columns = int(o.shape[1])
        return res
    res.n_obs_columns = int(o.shape[1])
    res.n_nondesign_columns = int(len(_nondesign_columns(ad, condition_col, donor_col, cfg)))
    res.details.update(dict(
        n_cells=int(o.shape[0]), n_donors=int(donor.nunique()),
        n_control_cells=int(is_ctrl.sum()), n_case_cells=int((~is_ctrl).sum()),
        n_control_donors=int(len(set(donor[is_ctrl]))),
        n_case_donors=int(len(set(donor[~is_ctrl]))),
        attenuation_min=cfg.attenuation_min,
        attenuation_min_full_effect_pp=cfg.attenuation_min_full_effect_pp))

    # ── ① 花名册 ──
    R, unclassified = build_variable_roster(ad, condition_col, donor_col,
                                            control_label, cfg)
    res.roster, res.unclassified = R, unclassified
    if len(R):
        res.design_defining = R.loc[R.bucket == BUCKET_DESIGN_DEFINING, "column"].tolist()
        res.adjustable = R.loc[R.bucket == BUCKET_ADJUSTABLE, "column"].tolist()
        for b in sorted({x for x in R.bucket if x.startswith(_EXCLUDED_PREFIX)}):
            res.excluded[b.replace(_EXCLUDED_PREFIX, "")] = \
                R.loc[R.bucket == b, "column"].tolist()
        # C1（v0.7.1）：显式声明被剔除的 condition 副本列——否则"支持域被打到 0"的
        # 真实原因不可见，读者会把它当成"数据不可识别"的证据。
        _redundant = R.loc[R.bucket == BUCKET_REDUNDANT_CONDITION, "column"].tolist()
        if _redundant:
            res.notes.append(
                f"redundant-with-condition column(s) {_redundant} were bucketed as "
                f"`{BUCKET_REDUNDANT_CONDITION}` and excluded from design-defining and "
                "from the common-support domain: each is constant within BOTH arms with "
                "a different value per arm, i.e. a one-to-one re-encoding of the "
                "condition contrast. Leaving it in design-defining would force the case "
                "profile to have zero control support (SUPPORT-BELOW-FLOOR by artefact). "
                "If a column listed here is NOT a label duplicate, rename it or declare "
                "it explicitly via identifiability_nondesign_cols.")
        # C-bug（v0.7.3）：供者级协变量（每供者 1 值、取值互异）不是标识符 —— 显式声明。
        _saturated = R.loc[R.bucket == BUCKET_DONOR_SATURATED, "column"].tolist()
        if _saturated:
            res.notes.append(
                f"donor-level covariate(s) {_saturated} have one distinct value per donor "
                f"and were bucketed `{BUCKET_DONOR_SATURATED}` (NOT an identifier): they "
                "cannot be balanced by matching (every donor is its own level), but they "
                "DO exist and should be carried into the limitation/estimand discussion "
                "rather than read as an id column.")
    if res.unclassified:
        res.notes.append(
            f"{len(res.unclassified)} obs column(s) could not be bucketed and were NOT "
            f"silently skipped: {res.unclassified} — resolve (rename / role_overrides / "
            f"identifiability_nondesign_cols) before interpreting this audit.")

    # ── ② 共同支持域 ──
    sup = common_support(ad, res.design_defining, condition_col, donor_col,
                         control_label, cfg)
    res.case_profiles = list(sup["case_profiles"])
    res.support_donors = list(sup["support_donors"])
    res.support_cells = int(sup["support_cells"])
    res.n_control_donors = int(sup["n_control_donors"])
    res.balance = adjustable_balance(ad, res.adjustable, res.support_donors,
                                     condition_col, donor_col, control_label)
    res.unbalanceable = [c for c, t in res.balance.items() if bool(t["unbalanceable"].any())]
    # ── A2（v0.6.0）：positivity / 协变量重叠（与细胞计数下限并列）──
    _pos = positivity_overlap(ad, res.design_defining, condition_col, donor_col,
                              control_label, adjustable=res.adjustable)
    res.details["positivity_overlap"] = _pos
    if sup.get("vacuous"):
        res.notes.append(
            "no design-defining variable was found: the common support is the whole "
            "control group by construction (vacuous), so the identifiability check "
            "cannot restrict anything — interpret the ladder accordingly.")

    # ── ③ 阶梯 ──
    in_window = R.loc[(R.tvd_donor >= cfg.match_var_tvd_low) &
                      (R.tvd_donor <= cfg.hetero_max) &
                      (~R.bucket.str.startswith(_EXCLUDED_PREFIX)),
                      "column"].tolist() if len(R) else []
    if ladder is not None:
        locked = [(f"R{i} " + ("不匹配（锚点）" if not v else "∧".join(v)), list(v))
                  for i, v in enumerate(ladder)]
        res.notes.append(
            "matching ladder supplied explicitly and treated as pre-specified: "
            + " | ".join(k for k, _ in locked))
    else:
        locked = default_ladder(res.design_defining, res.adjustable, in_window)
    # P0-1（审计）：规则里引用了 obs 不存在的列时，不得静默忽略 —— 显式告知
    _missing_rule_vars = sorted({c for _, v in locked for c in v
                                 if c not in o.columns})
    if _missing_rule_vars:
        res.notes.append(
            f"rule variable(s) {_missing_rule_vars} are absent from obs: any rule whose "
            f"variables are all absent is reported as not evaluable (0 matched cells) "
            f"rather than silently degraded to an unmatched anchor.")
    res.ladder_rules = [k for k, _ in locked]
    # 可适用性修正（2026-10-06 跨数据集实跑发现）：若阶梯里**没有任何一条规则引用
    # 花名册中真正可用的变量**，则 R2/R3/R4 都退化为"不匹配锚点"，衰减比恒 = 1.000
    # → 单侧 → 会输出 COMPARABLE。这等于把"没有元数据"当成"两组可比"。
    # ⚠️ 判据必须基于**可用变量**（design-defining + adjustable），而不是"列是否存在"：
    # 有的数据集里 `source` 列存在但属供者级唯一（不可用），拿它当"用了变量"会漏判。
    _usable_vars = set(res.design_defining) | set(res.adjustable)
    _uses_real_var = any((c in _usable_vars) for _, v in locked for c in v)
    res.details["ladder_uses_real_variable"] = bool(_uses_real_var)
    res.details["usable_variables"] = sorted(_usable_vars)
    Lad, S, meds = build_matching_ladder(ad, genes, locked, condition_col,
                                         donor_col, control_label, cfg)
    res.ladder, res.ladder_summary = Lad, S

    # ── ④ 判定 ──
    # P2-8（审计）：`min_support_cells = 0` 是合法值（“不设下限”），不得被 `or` 吞掉
    _minsup = getattr(cfg, "identifiability_min_support_cells", None)
    floor = int(cfg.matched_min_cells_warn if _minsup is None else _minsup)
    res.cell_floor = floor
    below_floor = res.support_cells < floor
    _real_profiles = [p for p in res.case_profiles if not str(p).startswith("<none")]
    multi_profile = len(_real_profiles) > 1
    n_low = int((meds < cfg.attenuation_min).sum()) if len(meds) else 0
    n_high = int((meds >= cfg.attenuation_min).sum()) if len(meds) else 0
    spanning = bool(n_low > 0 and n_high > 0)
    med_of_med = float(np.median(meds)) if len(meds) else np.nan
    _cells_by_rule = dict(zip(S["rule"], S["n_cells"])) if len(S) else {}
    res.excluded_from_consistency = [k for k, _ in locked
                                     if int(_cells_by_rule.get(k, 0)) == 0]
    _n_genes = int(Lad["gene"].nunique()) if len(Lad) else 0
    res.details.update(dict(n_low=n_low, n_high=n_high, threshold_spanning=spanning,
                            median_of_medians=med_of_med, medians=meds,
                            below_floor=bool(below_floor),
                            multi_case_profile=bool(multi_profile),
                            n_panel_genes_found=_n_genes,
                            support_cells=res.support_cells))

    _ANNOTATION_LIMITATION_NOTE = (
        "this is an ANNOTATION limitation, not evidence of comparability: with no "
        "design-defining variable (or with matching that does not restrict the support "
        "domain) the ratio is 1 by construction. The v0.2 single-rule layer returns "
        "INDETERMINATE on the same data; do not read this state as 'COMPARABLE'. "
        "Supply donor-level design metadata (institution, donation type, sorting "
        "source, region, ...) or declare candidate variables via "
        "identifiability_nondesign_cols / role_overrides.")
    # C.4：支持域真空 + 无任何可评估规则真正缩域 → 一致性检查无信息（见下方判定链）
    # 只统计**可评估**（n_cells > 0）的非锚点规则：n_cells == 0 的规则根本不进入一致性
    # 统计（它们"不匹配任何东西"），不能算作"缩了支持域"。
    _evaluable_cells = [int(v) for k, v in _cells_by_rule.items()
                        if int(v) > 0 and not str(k).startswith("R0")]
    # 触发条件（v0.3.3 放宽）：支持域**等于整个对照组**——要么 `vacuous`（没有任何
    # design-defining 变量），要么匹配没把任何对照挡在支持域之外——**且**没有任何可评估
    # 规则真正缩小支持域。只测 `vacuous` 会漏掉一个重要情形：数据集里存在**在该数据集内
    # 恒定的零信息列**（如全队列同一 `region`），它会被归入 design-defining 桶，于是
    # `vacuous=False`，守卫失效 → 衰减比恒 = 1 却被报成 `COMPARABLE`，与本工具自己的
    # 文案（"do not read this state as 'COMPARABLE'"）自相矛盾。
    _ctrl_cells_total = int(res.details.get("n_control_cells", 0) or 0)
    _support_is_full_control = (_ctrl_cells_total > 0
                               and int(res.support_cells) >= _ctrl_cells_total)
    _vacuous_no_restriction = (bool(sup.get("vacuous")) or _support_is_full_control) and all(
        c >= int(res.support_cells) for c in _evaluable_cells)
    one_sided_ok = bool(getattr(cfg, "identifiability_one_sided_verdict", True))
    _rng = (f"[{meds.min():.3f}, {meds.max():.3f}]" if len(meds) else "[n/a]")

    _pos_ov_min = _pos.get("overlap_min")
    _pos_bad = bool(_pos_ov_min is not None
                    and float(_pos_ov_min) < float(getattr(cfg, "positivity_overlap_min",
                                                          0.5)))
    if below_floor:
        res.attenuated_suppressed = bool(getattr(cfg, "identifiability_suppress_below_floor",
                                                 True))
        res.verdict = (f"{VERDICT_SUPPORT_BELOW_FLOOR} — design not identifiable "
                       f"(common support {res.support_cells} cells < pre-specified floor "
                       f"{floor}, and the case group is constant on "
                       f"{'/'.join(res.design_defining) or 'no design variable'})")
        res.notes.append(
            "no adjusted between-group estimate exists: matching on a variable that is "
            "constant within cases conditions on the very contrast being estimated. "
            "The deliverable is the design table above plus a scope statement.")
        if res.attenuated_suppressed and len(meds):
            _comp = S.loc[~S.rule.str.startswith("R0") & (S.n_cells > 0), "n_donors"]
            _maxd = int(_comp.max()) if len(_comp) else 0
            res.notes.append(
                f"attenuation ratios are NOT reported as evidence: they range "
                f"{meds.min():.3f}–{meds.max():.3f} across the ladder and every matched "
                f"subset is decided by ≤ {_maxd} donor(s).")
    elif multi_profile:
        # A1-F6：多个病例设计剖面 → 并入一个支持域会跨层比较，不是标准共同支持域
        res.attenuated_suppressed = True
        res.verdict = (f"{VERDICT_INDETERMINATE} — the case group spans multiple design "
                       f"profiles ({_real_profiles}); a pooled common support would "
                       f"compare across strata")
        res.notes.append(
            "no pooled support was formed: with >1 case profile the between-group "
            "estimate would mix strata. Run the audit separately within each profile "
            "(e.g. subset the data), or declare the profile as part of the estimand.")
    elif not len(meds):
        if not _uses_real_var:
            # 没有任何可用设计/可调整变量 → 阶梯里不存在可评估规则。此处仍用"注释局限"的
            # 措辞（比"no evaluable matching rule"更能说明**为什么**），并抑制衰减比。
            res.attenuated_suppressed = True
            res.verdict = (f"{VERDICT_INDETERMINATE} — no usable design or adjustable variable "
                           f"(none of the ladder rules references a column that is both present "
                           f"and non-degenerate); no evaluable matching rule exists, so the "
                           f"attenuation branch cannot be evaluated")
            res.notes.append(_ANNOTATION_LIMITATION_NOTE)
        else:
            res.verdict = (f"{VERDICT_INDETERMINATE} — no evaluable matching rule in the "
                           f"ladder" + (f" (and no target gene found among the {_n_genes} "
                                        f"resolved)" if _n_genes == 0 else ""))
        if _n_genes == 0:
            res.notes.append(
                "no target gene was found in the dataset — check the gene symbols / "
                "the gene-name fallback (ENSG → symbol) before interpreting anything.")
        if not res.ladder_rules:
            res.notes.append(
                "the ladder is empty (no rule supplied and none could be built): "
                "the attenuation branch cannot be evaluated at all.")
    elif (not _uses_real_var) or _vacuous_no_restriction:
        # 两种情形都不得读作"可比"（C.4，2026-10-07 已知答案验证 + 独立复核）：
        #   ① 无任何可用设计/可调整变量 → 衰减分支根本不可评估；
        #   ② 支持域**真空**（无 design-defining 变量 → 共同支持域按构造 = 整个对照组）
        #      且**没有任何可评估规则真正缩小支持域** → 所有衰减比恒 = 1（by construction），
        #      一致性检查无信息。原实现只测"规则是否**引用**了可用列"，于是情形②会落进
        #      单侧分支被报成 COMPARABLE，与本工具自己的文案（"do not read this state as
        #      'COMPARABLE'"）自相矛盾。
        res.attenuated_suppressed = True
        if not _uses_real_var:
            res.verdict = (f"{VERDICT_INDETERMINATE} — no usable design or adjustable variable "
                           f"(none of the ladder rules references a column that is both present "
                           f"and non-degenerate); the attenuation branch cannot be evaluated")
        else:
            res.verdict = (f"{VERDICT_INDETERMINATE} — the common support equals the whole "
                           f"control group ({res.support_cells} cells) and no evaluable "
                           f"rule restricts it (every matched subset is the whole control "
                           f"group); attenuation ratios are 1 by construction")
        res.notes.append(_ANNOTATION_LIMITATION_NOTE)
    elif spanning:
        res.verdict = (f"{VERDICT_THRESHOLD_SPANNING} — the across-ladder median moves "
                       f"across the prespecified threshold ({n_low} rule(s) below, "
                       f"{n_high} at/above {cfg.attenuation_min}); the attenuation "
                       f"ratio is not a single estimand across rules")
        res.notes.append(
            f"per-rule medians: "
            + ", ".join(f"{k}={v:.3f}" for k, v in
                        zip(S.loc[~S.rule.str.startswith('R0'), 'rule'],
                            S.loc[~S.rule.str.startswith('R0'), 'median_ratio'])))
    elif not one_sided_ok:
        # 乙方案（保守；“两个方向都不给结论”）：默认关闭
        res.verdict = (f"{VERDICT_INDETERMINATE} — the prespecified attenuation threshold "
                       f"({cfg.attenuation_min}) lies outside the observed ladder range "
                       f"{_rng}, so all rules fall on one side")
        res.notes.append(
            "identifiability_one_sided_verdict=False: one-sided outcomes are reported "
            "as INDETERMINATE by configuration (the more conservative reading).")
    else:
        # 甲方案（默认）：单侧直接给方向判定 + **强制警示**
        _side = "below" if n_high == 0 else "at or above"
        _lab = VERDICT_NOT_SEPARABLE if n_high == 0 else VERDICT_COMPARABLE
        res.verdict = (f"{_lab} — under every evaluable rule the across-ladder median "
                       f"attenuation ratio is {_side} {cfg.attenuation_min} "
                       f"(n={len(meds)} rule(s); range {_rng})")
        res.notes.append(
            f"ONE-SIDED CAVEAT (mandatory): the threshold lies outside the observed "
            f"ladder range {_rng}, so the consistency check carries NO information here "
            f"(all rules fall on the same side by construction). This verdict rests on "
            f"the magnitude comparison alone and cannot be reported without this "
            f"caveat.")
        if n_low == 0:
            res.notes.append(
                "COMPARABLE here means 'matching did not change the effect size under "
                "any prespecified rule' — it does NOT establish comparability of the "
                "groups: matched subsets too small (or too imbalanced on the remaining "
                "adjustable variables) also fail to change the effect size.")

    # ── A2（v0.6.0）：positivity / 协变量重叠（**不强覆阶梯判定**，强制警示）──
    # 设计决策：阶梯已经把估计限定在共同支持域上；用 NOT SEPARABLE 覆盖会丢失衰减信息。
    # 故 A2 作**并列诊断**：报告供体级重叠，并在不足时给强制警示。
    if _pos_bad:
        res.notes.append(
            f"POSITIVITY-OVERLAP CAVEAT (mandatory): donor-level covariate overlap is "
            f"insufficient (case {_pos.get('overlap_case')} / control "
            f"{_pos.get('overlap_ctrl')}; min < positivity_overlap_min="
            f"{getattr(cfg, 'positivity_overlap_min', 0.5)}). The ladder estimate is "
            f"defined only on the (possibly small) common-support population; the "
            f"overall case-vs-control contrast violates positivity and is NOT identified "
            f"for the full population (Petersen 2012, PMID 21030422). Cell counts are a "
            f"CAPACITY prerequisite, not the identifiability criterion.")

    # ── 附加报告（审计 A1-F3 / P2-1 / P2-2 / P2-5 / P2-6）──
    if len(meds) and not below_floor:
        _amp = S.loc[~S.rule.str.startswith("R0"), "n_amplified"]
        if int(_amp.sum()) > 0:
            res.notes.append(
                f"effect-size ratio > 1 in {int(_amp.sum())} rule×gene cell(s) "
                f"(max {int(_amp.max())} in a single rule): matching made the contrast "
                f"LARGER. The ratio is therefore an effect-size ratio, not a monotone "
                f"measure of the confounding share (it can exceed 1).")
        res.notes.append(
            "scale note: full/matched effects are cell-pooled detection-rate differences "
            "(pseudo-replicated), not donor-level estimates; they are comparable "
            "across rules but are not the donor-level estimand used in the manuscript.")
    _derived = [c for c in (res.excluded.get("derived-from-expression") or [])]
    if _derived:
        res.notes.append(
            f"columns excluded as expression-derived ({_derived}) were selected by "
            f"column-NAME prefix ({getattr(cfg, 'identifiability_derived_prefixes', [])}) "
            f"plus explicit declaration — this is a naming heuristic, not a causal "
            f"property; declare such columns explicitly if the prefix is misleading.")
    # 小样本风险不要重叠到"没有变量可用"的情形上
    if res.design_defining and res.details.get("n_case_donors"):
        # 修正（2026-10-06 晚，据外部审查 B5 的**症状**）：原实现把"病例组内的观测主占比"
        # 当作参照占比 p —— 而 design-defining 变量在病例组内的主占比恒 = 1，
        # 于是 P ≡ 1^n + 0^n = 1.00（恒定退化，恒定触发警告且无信息量）。
        # 现改为用**参照人群**：① 整个供者池 ② 对照组臂；取两者中较大（更保守）。
        # 两者都为 0/1（参照人群里也无变异）→ 不可评估，不报数字。
        _d = int(res.details["n_case_donors"])
        _rows = []
        for _c in res.design_defining:
            _dv = _donor_level(o, _c, donor)
            _grp = _group_of_donor(is_ctrl, donor).reindex(_dv.index)
            _case_vals = _dv[_grp == "case"].astype(str)
            if not len(_case_vals):
                continue
            _lvl = _case_vals.value_counts().idxmax()
            _all_vals = _dv.astype(str)
            _ctrl_vals = _all_vals[_grp == "ctrl"]
            _p_pool = float((_all_vals == _lvl).mean())
            _p_ctrl = float((_ctrl_vals == _lvl).mean()) if len(_ctrl_vals) else float("nan")
            _r_pool = chance_constancy_risk(_p_pool, _d)
            _r_ctrl = chance_constancy_risk(_p_ctrl, _d)
            _rows.append((_c, _p_pool, _r_pool, _p_ctrl, _r_ctrl))
        _assessable = [r for r in _rows
                       if not (np.isnan(r[2]) and np.isnan(r[4]))]
        if _assessable:
            _worst = max(max(r[2], r[4]) if not np.isnan(r[4]) else r[2]
                         for r in _assessable)
            if _worst > 0.05:
                _detail = "; ".join(
                    f"{r[0]}: P={max(r[2], r[4]) if not np.isnan(r[4]) else r[2]:.3f} "
                    f"(pooled share {r[1]:.2f}; control-arm share "
                    f"{'n/a' if np.isnan(r[3]) else f'{r[3]:.2f}'})"
                    for r in _assessable)
                res.notes.append(
                    f"chance-constancy risk: with only {_d} case donor(s), a design "
                    f"variable can be constant within cases by sampling alone "
                    f"({_detail}); the reference prevalences are taken from the pooled "
                    f"donor pool and from the control arm, not from the case group "
                    f"(whose share is 1 by construction). 'design-defining' is a "
                    f"property of this dataset, not proof of population-level "
                    f"constancy.")
        else:
            res.notes.append(
                f"chance-constancy risk is not assessable for "
                f"{res.design_defining}: the variable has no variation in the pooled "
                f"donor pool either, so a sampling-based risk cannot be quantified.")

    if res.case_profiles and res.details.get("multi_case_profile"):
        res.notes.append(
            "support reported for the largest single case profile only: pooling the "
            "other profile(s) would compare across strata.")
    if res.unbalanceable:
        res.notes.append(
            f"adjustable variable(s) {res.unbalanceable} cannot be balanced inside the "
            f"common support (a case level has no support-domain control donor).")
    if len(meds) and spanning:
        res.notes.append(
            f"even disregarding identifiability, the attenuation ratio moves between "
            f"{meds.min():.3f} and {meds.max():.3f} across rules (spread "
            f"{meds.max() - meds.min():.3f}), which is larger than any plausible "
            f"threshold — the conclusion therefore does not rest on its exact value.")
    if res.ladder_rules and S.loc[S.rule.str.startswith("R1"), "median_ratio"].notna().any():
        r1 = S.loc[S.rule.str.startswith("R1")].iloc[0]
        res.notes.append(
            f"for comparison, the original single-rule convention (R1) uses "
            f"{int(r1.n_donors)} donor(s) / {int(r1.n_cells)} cell(s) and yields a median "
            f"attenuation of {r1.median_ratio:.3f} — v0.2 would have reported that number "
            f"as the confounding share.")
    return res
