"""审计二：对照构成审计。

判断病例组与对照组是否可比，或生物学条件是否与机构/取材/分选等设计变量完全共线。

修订记录（2026-10-01）：
  · S1  F_hetero 的 mode 逻辑替换为**组间分布距离**（TVD + Cramér's V）。
        原定义在"三中心完全平衡"设计下会误判 NOT SEPARABLE（假阳性），
        因为它衡量的是"对照与病例的众数是否一致"，而非两组分布是否相同。
        f_hetero 字段保留（语义改为 TVD）以兼容旧 API。
  · S2  衰减比只在 full effect >= 阈值（默认 5pp）的基因上计算，
        避免检出地板基因的抽样噪声主导中位数；并报告 ratio > 1 的情形。
  · B2  匹配变量与展示变量分离（cfg.match_vars vs cfg.strata_cols），
        并自动剔除样本级唯一的分层列（batch/lane 等），否则严格匹配恒失败。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import AuditConfig
from .donor_utils import donor_mode, is_donor_unique
from .io import gene_detection_vector


# ══════════════════ 分布距离 ══════════════════

def _group_distribution_distance(values: np.ndarray, is_ctrl: np.ndarray,
                                 donor: Optional[np.ndarray] = None):
    """返回 (TVD, Cramér's V, case_mode, 每水平比例表)。

    TVD = 0.5 * Σ|p_case(level) − p_ctrl(level)|，取值 0（分布相同）～1（完全分离）。
    平衡的多中心设计 TVD≈0；病例全在 A、对照全在 B 时 TVD=1。

    donor 非 None 时按**供者**加权（L1）：设计变量是供者属性，
    细胞级加权会让一个大供者主导 TVD，把"细胞数不均"误读为"设计不均"。
    """
    if donor is not None:
        # 每供者取其在该变量上的主取值，然后每供者一票
        d = pd.Series(donor).astype(str)
        v = pd.Series(values).astype(str)
        per_donor = v.groupby(d).agg(donor_mode)          # 单一实现：donor_utils.donor_mode
        dmask = pd.Series(np.where(is_ctrl, "ctrl", "case"),
                          index=pd.Series(donor).index).groupby(d).first()
        values = per_donor.values
        is_ctrl = (dmask.reindex(per_donor.index).values == "ctrl")

    s = pd.Series(values).astype(str)
    p_case = s[~is_ctrl].value_counts(normalize=True)
    p_ctrl = s[is_ctrl].value_counts(normalize=True)
    levels = sorted(set(p_case.index) | set(p_ctrl.index))
    pc = p_case.reindex(levels, fill_value=0.0)
    pk = p_ctrl.reindex(levels, fill_value=0.0)
    tvd = float(0.5 * (pc - pk).abs().sum())

    ct = pd.crosstab(pd.Series(np.where(is_ctrl, "control", "case")), s)
    n = int(ct.values.sum())
    cv = np.nan
    if ct.shape[0] >= 2 and ct.shape[1] >= 2 and n > 0:
        exp = np.outer(ct.sum(1), ct.sum(0)) / n
        chi2 = float(((ct.values - exp) ** 2 / np.where(exp > 0, exp, 1)).sum())
        cv = float(np.sqrt(chi2 / (n * (min(ct.shape) - 1))))
    case_mode = pc.idxmax() if len(pc) else None
    table = pd.DataFrame({"p_case": pc.round(4), "p_ctrl": pk.round(4),
                          "abs_diff": (pc - pk).abs().round(4)})
    return tvd, cv, case_mode, table


# ══════════════════ 结果容器 ══════════════════

@dataclass
class CompositionResult:
    crosstabs: Dict[str, pd.DataFrame] = field(default_factory=dict)
    f_hetero: Optional[float] = None          # 语义已改为 TVD（保持字段名兼容）
    f_hetero_cell: Optional[float] = None     # 细胞级 TVD（仅参考，L1）
    case_major_institution: Optional[str] = None
    hetero_table: Optional[pd.DataFrame] = None
    design_imbalance: Dict[str, Dict] = field(default_factory=dict)
    matched_donors: List[str] = field(default_factory=list)
    matched_cells: int = 0
    match_detail: Dict[str, Dict] = field(default_factory=dict)
    match_vars_used: List[str] = field(default_factory=list)
    match_vars_dropped: List[str] = field(default_factory=list)
    # P2-a：区分"尝试匹配但无人匹配"与"根本没有可用变量、未尝试匹配"
    match_attempted: bool = False
    match_not_possible_on: List[str] = field(default_factory=list)
    attenuation: Optional[pd.DataFrame] = None
    attenuation_suppressed: bool = False    # C7（v0.7.2）：匹配集过小 → 衰减比不当判据
    verdict: str = ""
    notes: List[str] = field(default_factory=list)

    def to_text(self) -> str:
        L = ["== Audit 2 · Control composition ==", ""]
        for k, v in self.crosstabs.items():
            L.append(f"-- condition × {k} (cells)")
            L.append(v.to_string())
            L.append("")
        if self.design_imbalance:
            L.append("-- condition-wise distribution distance (TVD / Cramér's V) --")
            for k, d in self.design_imbalance.items():
                cv = d.get("cramers_v")
                cv_s = f"{cv:.3f}" if isinstance(cv, (int, float)) and cv == cv else "n/a"
                tvd_d = d.get("tvd_donor", d.get("tvd"))
                L.append(f"    {k:<18} TVD(donor)={tvd_d:.3f}  TVD(cell)={d['tvd']:.3f}   V={cv_s}")
            L.append("")
        if self.match_vars_used:
            L.append(f"matching variables: {self.match_vars_used}")
        if self.match_vars_dropped:
            L.append(f"dropped (sample-unique / non-design): {self.match_vars_dropped}")
        if self.match_detail:
            L.append("  donor agreement (agree/total vars, cells):")
            for d, v in sorted(self.match_detail.items(),
                               key=lambda x: (-x[1]["agree"], x[0])):
                L.append(f"    {d:<14} {v['agree']}/{v['total']}   {v['cells']} cells")
        L.append("")
        L.append(f"-- matched comparison: donors={self.matched_donors} "
                 f"({self.matched_cells} cells) --")
        if self.attenuation is not None and len(self.attenuation):
            L.append(self.attenuation.to_string(index=False))
        L.append("")
        L.append(f"VERDICT: {self.verdict}")
        for n in self.notes:
            L.append(f"  · {n}")
        return "\n".join(L)


# ══════════════════ 匹配 ══════════════════

def _pick_match_vars(ad, cfg: AuditConfig, donor_col: str) -> Tuple[List[str], List[str]]:
    """确定匹配变量：优先 cfg.match_vars 中存在的列；自动剔除**样本/供者级唯一**的列（B2/N2）。

    N2 修复：原判据 `nunique >= 0.5 * n_cells` 只识别"细胞级唯一"，
    而真实数据里 batch/lane/library 是**供者级唯一**（10 供者 / 2000 细胞 → nunique=10），
    因此不会被剔除，严格匹配依然恒失败。
    新判据：该列在每个供者内恒定（per-donor nunique <= 1）**且**
            水平数 >= 供者数 → 该列与供者一一对应，不可能用于匹配。
    """
    o = ad.obs
    cand = [c for c in cfg.match_vars if c in o.columns]
    dropped = []
    if cfg.auto_filter_stratum_cols and cand:
        donor = o[donor_col].astype(str)
        n_donors = donor.nunique()
        keep = []
        for c in cand:
            s = o[c].astype(str)
            n_levels = s.nunique()
            # 去重（2026-10-06 审计 P2-12）：单一实现见 donor_utils.is_donor_unique
            donor_unique = is_donor_unique(s, donor, n_donors)
            cell_unique = n_levels >= max(2, int(0.5 * len(o)))
            if donor_unique or cell_unique:
                dropped.append(c)
            else:
                keep.append(c)
        cand = keep
    return cand, dropped


def _find_matched_donors(ad, condition_col, donor_col, control_label, match_vars):
    """找出与病例组设计变量一致的对照供者。

    返回 (matched, detail)：detail = {donor: {agree, total, cells}}，供人工判断。
    """
    o = ad.obs
    cond = o[condition_col].astype(str)
    donor = o[donor_col].astype(str)
    is_case = cond != str(control_label)

    if not match_vars:
        return [], {}

    targets = {}
    for c in match_vars:
        md = o.loc[is_case, c].astype(str).mode()
        if len(md):
            targets[c] = md.iloc[0]

    matched, detail = [], {}
    for d in sorted(set(donor[~is_case])):
        sub = o.loc[(donor == d) & (~is_case)]
        if not len(sub):
            continue
        agree = 0
        for c, t in targets.items():
            md = sub[c].astype(str).mode()
            if len(md) and md.iloc[0] == t:
                agree += 1
        detail[str(d)] = {"agree": agree, "total": len(targets), "cells": int(len(sub))}
        if agree == len(targets):
            matched.append(str(d))
    return matched, detail


# ══════════════════ 主函数 ══════════════════

def audit_control_composition(ad, genes: Sequence[str], condition_col: str,
                              donor_col: str, control_label: str,
                              matched_donors: Optional[Sequence[str]] = None,
                              cfg: Optional[AuditConfig] = None) -> CompositionResult:
    from .io import warn_if_gene_names_unresolved as _warn_genes
    _warn_genes(ad, genes, context="audit_control_composition")
    cfg = cfg or AuditConfig()
    o = ad.obs
    cond = o[condition_col].astype(str).values
    donor = o[donor_col].astype(str).values
    is_ctrl = cond == str(control_label)
    # P2-3：API 路径也须校验供者 ⊂ 组
    from .io import assert_donors_nested
    assert_donors_nested(donor, cond, donor_col, condition_col)
    res = CompositionResult()

    # ── 交叉表（展示用；两者取自同一 obs 保证索引对齐）──
    for c in cfg.strata_cols:
        if c in o.columns:
            res.crosstabs[c] = pd.crosstab(o[condition_col].astype(str),
                                           o[c].astype(str))

    # ── 组间分布距离（S1 / L1：判定用供者级，细胞级作参考）──
    for c in cfg.strata_cols:
        if c in o.columns:
            tvd_c, cv, _, _ = _group_distribution_distance(o[c].values, is_ctrl)
            tvd_d, _, _, _ = _group_distribution_distance(o[c].values, is_ctrl,
                                                          donor=donor)
            res.design_imbalance[c] = {
                "tvd": round(tvd_c, 3),           # 细胞级（参考）
                "tvd_donor": round(tvd_d, 3),     # 供者级（判定用）
                "cramers_v": round(cv, 3) if cv == cv else np.nan}

    hcol = cfg.hetero_col
    if hcol in o.columns:
        # 判定用供者级 TVD（L1）；细胞级作为参考一并输出
        tvd_d, cv_d, case_mode, _ = _group_distribution_distance(
            o[hcol].values, is_ctrl, donor=donor)
        tvd_c, cv_c, _, tbl = _group_distribution_distance(o[hcol].values, is_ctrl)
        res.f_hetero = round(tvd_d, 3)        # 字段名保留，语义已改为供者级 TVD
        res.f_hetero_cell = round(tvd_c, 3)
        res.case_major_institution = case_mode
        res.hetero_table = tbl
        res.notes.append(
            f"condition × {hcol}: donor-level TVD = {tvd_d:.3f}"
            + (f" (Cramér's V = {cv_d:.3f})" if cv_d == cv_d else "")
            + f"; cell-level TVD = {tvd_c:.3f} (reference only)"
            + "  — 0 = identical distributions, 1 = fully collinear")

    # ── 匹配变量（B2）──
    match_vars, dropped = _pick_match_vars(ad, cfg, donor_col)
    res.match_vars_used = match_vars
    res.match_vars_dropped = dropped
    if dropped:
        res.notes.append(
            f"columns excluded from matching because they are sample-level unique or "
            f"not design variables: {dropped}")

    # ── 匹配供者（P2-a：三分支，避免给出错误的原因）──
    if not getattr(cfg, "run_matched_comparison", True):
        # 策略 3：使用者显式关闭匹配（check 的 no_matching 方案）
        res.match_attempted = False
        res.notes.append(
            "Matched comparison was DISABLED by configuration "
            "(run_matched_comparison: false). Comparability is therefore "
            "assessed on the donor-level TVD alone; this choice must be "
            "stated in the limitations — it is the honest option when "
            "confounding with group is structural rather than incidental.")
    elif matched_donors:
        res.matched_donors = [str(d) for d in matched_donors]
        res.match_attempted = True
    elif not match_vars:
        # 分支 A：没有任何可用匹配变量 → 未尝试匹配
        res.match_attempted = False
        res.match_not_possible_on = list(dropped)
        res.notes.append(
            f"NO matching was ATTEMPTED: none of the declared match variables "
            f"({list(cfg.match_vars)}) is usable in this dataset "
            f"(all were dropped as sample/donor-unique: {dropped}). "
            f"Any 'matched' comparison would be trivially satisfied — use the "
            f"donor-level TVD as the primary evidence instead.")
    else:
        # 分支 B：有变量、真的尝试了
        res.match_attempted = True
        # P2-a：**部分剔除**才是更常见的真实场景——
        # 若某个设计变量因供者级唯一被剔除，剩下的可能只是单水平列（如 region=heart），
        # 于是"所有对照都匹配"，衰减比≈1 而 verdict 仍 COMPARABLE，
        # 名不副实。此处显式记录被剔除的设计变量并降级说明。
        design_dropped_b = [c for c in dropped if c in cfg.match_vars]
        if design_dropped_b:
            res.match_not_possible_on = list(design_dropped_b)
        trivially = [c for c in match_vars if o[c].astype(str).nunique() == 1]
        res.matched_donors, res.match_detail = _find_matched_donors(
            ad, condition_col, donor_col, control_label, match_vars)
        if res.matched_donors:
            res.notes.append(
                f"matched donors auto-detected on {len(match_vars)} design variable(s): "
                f"{', '.join(match_vars)}")
            if design_dropped_b:
                res.notes.append(
                    f"CAVEAT: matching was performed on {match_vars} only; "
                    f"{design_dropped_b} could not be matched (donor-unique). "
                    f"The matched comparison is therefore WEAKER than a same-design "
                    f"comparison — weigh the donor-level TVD alongside it.")
            if trivially and len(trivially) == len(match_vars):
                res.notes.append(
                    f"CAVEAT: every remaining match variable is single-level "
                    f"({trivially}) — the match is TRIVIALLY satisfied by all controls "
                    f"and carries no design information; rely on the donor-level TVD.")
        else:
            res.notes.append(
                f"matching WAS attempted on {match_vars} but NO control donor matched "
                f"the case group on all of them — no matched comparison can be formed "
                f"(see donor agreement table)")

    m_mask = (is_ctrl & np.isin(donor, res.matched_donors)) if res.matched_donors \
        else np.zeros(len(cond), dtype=bool)
    res.matched_cells = int(m_mask.sum())

    # ── 效应量衰减（S2：只在 full effect 达标的基因上计算）──
    if res.matched_cells > 0:
        rows = []
        for g in genes:
            v = gene_detection_vector(ad, g)
            if v is None:
                continue
            full_c = v[is_ctrl].mean() * 100
            full_s = v[~is_ctrl].mean() * 100
            mat_c = v[m_mask].mean() * 100
            eff_full = abs(full_c - full_s)
            eff_mat = abs(mat_c - full_s)
            rows.append({
                "gene": g,
                "full_ctrl_pct": round(full_c, 1),
                "full_case_pct": round(full_s, 1),
                "full_effect_pp": round(eff_full, 1),
                "matched_ctrl_pct": round(mat_c, 1),
                "matched_effect_pp": round(eff_mat, 1),
                "attenuation_ratio": (round(eff_mat / eff_full, 2) if eff_full > 0 else np.nan),
                "direction_preserved": bool((mat_c - full_s) * (full_c - full_s) >= 0),
            })
        res.attenuation = pd.DataFrame(rows)

    # ── 判定 ──
    if res.f_hetero is None and res.matched_cells == 0:
        res.verdict = "INDETERMINATE — no design variables available to assess comparability"
    else:
        flags = []
        # P2-d/P3：供者级 TVD 的粒度由【每组供者数的较小者】决定
        # （2 ctrl + 8 case 时最小非零步长是 1/2，而不是 1/10）
        _cond_arr = np.asarray(donor)
        _n_ctrl_d = len({d for d in set(_cond_arr) if is_ctrl[_cond_arr == d][0]})
        _n_case_d = len(set(_cond_arr)) - _n_ctrl_d
        _n_min = min(_n_ctrl_d, _n_case_d) if _n_ctrl_d and _n_case_d else len(set(_cond_arr))
        if _n_min < 5:
            res.notes.append(
                f"TVD resolution note: {_n_ctrl_d} control + {_n_case_d} case donors → "
                f"the smallest non-zero donor-level TVD is ~1/{_n_min}; values near the "
                f"threshold should not be over-interpreted")
        # ⚠️ 判定须用**全部** strata_cols 的供者级 TVD，而非仅 hetero_col。
        # 原实现只看 hetero_col（默认 institution），使 source / donation_type
        # 等列的极端不平衡完全不进入判定 —— 而分选来源列恰恰可能是核心混杂。
        # 该缺口此前只能靠审计四检查 5 兜底。
        tvd_items = []
        _donor_s = o[donor_col].astype(str)
        _n_don = _donor_s.nunique()
        # C5（v0.7.2）：一个列只有在**供者级有 ≥ 2 个水平**时才算“可评估的设计变量”；
        # 单水平列（如全队列同一 region）的 TVD 恒为 0，不携带可比性信息。
        def _lev(c):
            return int(o[c].astype(str).groupby(_donor_s).first().nunique())
        _assessable_cols = []
        if res.f_hetero is not None:
            tvd_items.append((hcol, res.f_hetero))
            if _lev(hcol) >= 2:
                _assessable_cols.append(hcol)
        # N2 规则同样适用于**判定**：供者级唯一（或细胞级唯一）的列不是设计
        # 变量 —— 每个供者自成一类，其 TVD 必然 = 1，与条件无关。原实现只在
        # 匹配时排除它们（_pick_match_vars），判定时未排除，会把 batch/lane/
        # library 这类列误报为"与条件共线"。
        for _c, _d in (res.design_imbalance or {}).items():
            _v = _d.get("tvd_donor")
            if _c == hcol or _v is None or _v != _v:
                continue
            _s = o[_c].astype(str)
            if is_donor_unique(_s, _donor_s, _n_don):
                continue                      # 样本/供者级唯一 → 非设计变量（单一实现）
            tvd_items.append((_c, _v))
            if _lev(_c) >= 2:
                _assessable_cols.append(_c)
        for _c, _v in tvd_items:
            if _v > cfg.hetero_max:
                flags.append(f"TVD({_c})={_v:.2f} (>{cfg.hetero_max})")
        _tvd_assessed = bool(_assessable_cols)   # C5：是否真的评估过（而非全为单水平列）
        if res.f_hetero is None:
            res.notes.append(
                f"TVD not assessed: column '{hcol}' is absent; "
                "comparability rests on the matched comparison alone")

        if res.attenuation is not None and len(res.attenuation):
            usable = res.attenuation.loc[
                res.attenuation["full_effect_pp"] >= cfg.attenuation_min_full_effect_pp,
                "attenuation_ratio"].dropna()
            n_excl = int(len(res.attenuation) - len(usable))
            if n_excl:
                res.notes.append(
                    f"{n_excl} gene(s) excluded from the attenuation summary "
                    f"(full effect < {cfg.attenuation_min_full_effect_pp} pp)")
            if len(usable):
                med = float(usable.median())
                res.notes.append(
                    f"median attenuation ratio over {len(usable)} gene(s) with "
                    f"full effect >= {cfg.attenuation_min_full_effect_pp} pp = {med:.2f}")
                if med < cfg.attenuation_min:
                    # C7（v0.7.2）：匹配集细胞数 < matched_min_cells_warn 时，该数字由
                    # 过少供者决定（如心脏案例 1 供者/31 细胞），**不得**当作判定依据——
                    # 与审计二b 对同一数字的抑制保持一致，避免“同一量两套认识论权重”。
                    if res.matched_cells and res.matched_cells < cfg.matched_min_cells_warn:
                        res.attenuation_suppressed = True
                        res.notes.append(
                            f"attenuation ratio {med:.2f} (<{cfg.attenuation_min}) is NOT "
                            f"used as a verdict criterion: the matched subset is only "
                            f"{res.matched_cells} cells < matched_min_cells_warn="
                            f"{cfg.matched_min_cells_warn} (decided by too few donors). "
                            "Comparability rests on the donor-level TVD alone — consistent "
                            "with audit-2b, which suppresses the same number.")
                    else:
                        flags.append(f"median attenuation ratio={med:.2f} (<{cfg.attenuation_min})")
                amp = usable[usable > 1.0]
                if len(amp):
                    res.notes.append(
                        f"{len(amp)} gene(s) show effect-size AMPLIFICATION after matching "
                        "(ratio > 1) — expected by chance in small matched subsets")

        if flags:
            res.verdict = ("NOT SEPARABLE — condition is collinear with design variables ("
                           + "; ".join(flags) + ")")
        elif not _tvd_assessed:
            # C5（v0.7.2）：没有任何 ≥ 2 水平的设计变量 → 可比性**根本没被评估**，
            # 不得输出褒义的 COMPARABLE（旧实现会给“COMPARABLE … (TVD not assessed)”）。
            res.verdict = ("NOT ASSESSED — no design variable with >=2 levels was "
                           "available, so comparability was not assessed; this is NOT "
                           "evidence of comparability")
            res.notes.append(
                "no design variable had >=2 donor-level values, so the donor-level TVD "
                "could not assess comparability. This differs from COMPARABLE: nothing "
                "was assessed. Supply design metadata (institution / donation type / "
                "region / sorting source) or declare candidate variables.")
        else:
            res.verdict = ("COMPARABLE on the available design variables"
                           + ("" if res.f_hetero is not None else " (TVD not assessed)"))

    if res.matched_cells and res.matched_cells < cfg.matched_min_cells_warn:
        res.notes.append(
            f"matched comparison uses only {res.matched_cells} cells — low power; "
            "the interpretable result is the attenuation, not the matched estimate itself")
    return res
