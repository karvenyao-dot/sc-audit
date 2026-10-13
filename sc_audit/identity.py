"""审计四：细胞身份与捕获完整性。

由病例–对照单细胞数据的细胞身份取证实践逼出。前三项审计均无法
发现它检出的问题，因此作为独立模块加入。

它在两个层面上是新的：
  · 把审计二**下沉到细胞类型内部**——一个被正确注释的细胞簇，其内部仍可能存在
    只在其中一个臂出现的亚群（例如：某细胞类型内、某分选来源的细胞全部集中在
    单一组）。
  · 把 **QC 规则本身当作待审计对象**——"阴性标志全 0"这条门控规则与取材方式/部位
    共线（例如病例组活检取材 vs 对照供心），一致施加会差异性过滤。

八项检查：
  1 自家门控通过率（按注释谱系）
  2 中位 genes/cell 在同数据集谱系内的排名
  3 核滞留转录本（MALAT1）相对排名
  4 长/低拷贝基因在**复杂度匹配**下的检出
  5 ★ 亚群不平衡检测（来源 × 条件，供者级 TVD）
  6 ★ 门控规则的条件共线性检验
  7 模态匹配检查（参照表是否存在同模态列）
  8 取材方式标志按臂分层

设计原则：**测量与判定分离**——所有量先算出来放进 dataclass，判定层只按 config
阈值给结论；任何一项缺列或缺基因都**显式降级并记录原因**，绝不静默跳过。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import AuditConfig
from .io import gene_detection_vector


# ══════════════════ 结果容器 ══════════════════

@dataclass
class IdentityResult:
    # 1 门控通过率
    gate_by_lineage: Optional[pd.DataFrame] = None
    target_gate_pass: Optional[float] = None
    gate_margin: Optional[float] = None          # 目标 / 其他谱系中位数
    # 2 复杂度排名
    complexity_by_lineage: Optional[pd.DataFrame] = None
    complexity_rank: Optional[int] = None        # 目标谱系在谱系中的名次（1=最低）
    complexity_n_lineages: int = 0
    n_lineages: int = 0                          # 注释谱系总数（C8：rank 分母，独立于复杂度列）
    # 3 核标志
    nuclear_marker_pct: Optional[float] = None
    nuclear_marker_rank: Optional[int] = None
    # 4 复杂度匹配
    matched_complexity: Optional[pd.DataFrame] = None
    long_gene_deficit: Optional[Dict[str, float]] = None
    # 5 亚群不平衡
    subset_by_arm: Optional[pd.DataFrame] = None
    subset_tvd: Dict[str, float] = field(default_factory=dict)
    imbalanced_subsets: List[str] = field(default_factory=list)
    # 6 门控共线
    negative_marker_by_arm: Optional[pd.DataFrame] = None
    collinear_markers: List[str] = field(default_factory=list)
    gate_pass_by_arm: Optional[pd.DataFrame] = None
    # 7 模态
    modality: Optional[str] = None
    reference_modalities: List[str] = field(default_factory=list)
    modality_matched: Optional[bool] = None
    # 8 取材标志
    procurement_by_arm: Optional[pd.DataFrame] = None
    # 判定
    verdict: str = ""
    notes: List[str] = field(default_factory=list)

    def to_text(self) -> str:
        L = ["== Audit 4 · Cell identity and capture completeness ==", ""]
        if self.gate_by_lineage is not None:
            L += ["-- 1. own-gate pass rate, by annotated lineage --",
                  self.gate_by_lineage.to_string(), ""]
        if self.complexity_by_lineage is not None:
            L += ["-- 2. median genes/cell, by lineage --",
                  self.complexity_by_lineage.to_string(), ""]
        if self.nuclear_marker_pct is not None:
            _nl = self.n_lineages or self.complexity_n_lineages
            L += [f"-- 3. {self.nuclear_marker_pct:.1f}% nuclear-marker detection"
                  f" (rank {self.nuclear_marker_rank}/{_nl})", ""]
        if self.matched_complexity is not None:
            L += ["-- 4. long/low-copy genes at matched complexity --",
                  self.matched_complexity.to_string(), ""]
        if self.subset_by_arm is not None:
            L += ["-- 5. subset (source) x arm --",
                  self.subset_by_arm.to_string(), ""]
            if self.subset_tvd:
                L.append(f"   subset TVD: {self.subset_tvd}")
            L.append("")
        if self.negative_marker_by_arm is not None:
            L += ["-- 6. negative markers by arm (gate collinearity) --",
                  self.negative_marker_by_arm.to_string(), ""]
        if self.gate_pass_by_arm is not None:
            L += ["-- 6b. gate pass rate by arm --",
                  self.gate_pass_by_arm.to_string(), ""]
        if self.procurement_by_arm is not None:
            L += ["-- 8. procurement markers by arm --",
                  self.procurement_by_arm.to_string(), ""]
        L += [f"VERDICT: {self.verdict}"]
        for n in self.notes:
            L.append(f"  note: {n}")
        return "\n".join(L)


# ══════════════════ 工具 ══════════════════

def _detect_many(ad, genes: Sequence[str]) -> Dict[str, np.ndarray]:
    """一次取出多个基因的 0/1 检出向量。

    与 gene_detection_vector 的区别：本函数**支持 backed 模式**（`.h5ad` 以
    `backed='r'` 打开时 ad.X 是 h5py dataset，`sp.issparse` 判定会失败，
    gene_detection_vector 因此在 backed 数据上对任何基因都返回 None —— 这会
    让审计四静默地"找不到标志基因"）。这里改为一次列切片 to_memory，
    既兼容 backed 又避免逐基因重复读盘。

    基因缺失时不返回该键（调用方须处理缺键，不得静默填零）。
    """
    # 单一实现：io.detection_vectors（backed-safe；单/多基因同路径）。
    # 原实现在单基因时回退到整矩阵 ad.X 并返回第 0 列 —— 静默给出错误基因。
    from .io import detection_vectors
    return detection_vectors(ad, genes)


def _detect(ad, gene: str):
    """单个基因的检出向量；缺失返回 None。标量便捷入口（内部仍走批量）。"""
    d = _detect_many(ad, [gene])
    return d.get(gene)


def _tvd(a: pd.Series, b: pd.Series) -> float:
    """两分布的总变差距离，定义在对齐后的取值集合上。"""
    idx = sorted(set(a.index) | set(b.index))
    pa = a.reindex(idx, fill_value=0).astype(float)
    pb = b.reindex(idx, fill_value=0).astype(float)
    pa = pa / pa.sum() if pa.sum() else pa
    pb = pb / pb.sum() if pb.sum() else pb
    return float(0.5 * np.abs(pa - pb).sum())


# ══════════════════ 主审计 ══════════════════

def audit_cell_identity(ad, condition_col: str, donor_col: str, control_label: str,
                        cfg: Optional[AuditConfig] = None,
                        reference_modalities: Optional[Sequence[str]] = None,
                        modality: Optional[str] = None) -> IdentityResult:
    """执行审计四的八项检查。

    reference_modalities: 参照表中已有的模态列表（用于检查 7）。
    modality:            被审计数据集的模态（用于检查 7）。
    """
    cfg = cfg or AuditConfig()
    res = IdentityResult(modality=modality)
    o = ad.obs
    cond = o[condition_col].astype(str).values
    is_ctrl = cond == str(control_label)

    # ── 前置：注释列 ──
    if cfg.identity_col not in o.columns:
        res.verdict = "INDETERMINATE (no cell-type annotation column)"
        res.notes.append(f"obs 中无 '{cfg.identity_col}'；审计四无法执行。"
                         f"设 cfg.identity_col 指向注释列，或改用可门控方案（见 gating.py）。")
        return res
    lin = o[cfg.identity_col].astype(str).values
    lineages = sorted(pd.unique(lin))
    res.n_lineages = len(lineages)
    # 单谱系只降级【依赖谱系对照】的检查（1/2/3 的"相对排名"），
    # 不阻断 5/6/7/8 —— 亚群不平衡与门控共线性与谱系数无关，
    # 且它们是本审计的核心项，不该被一个无关条件连带关闭。
    single_lineage = len(lineages) < 2
    if single_lineage:
        # 显式列出**哪些可评估、哪些仅部分可评估** —— 否则读者会把"不阻断 5/6/8"
        # 误读为"完整审计通过"。（不阻断是对的，但须防止误读。）
        res.notes.append(
            "ASSESSABILITY — only 1 annotated lineage: checks 1 (gate margin), "
            "2 (complexity rank) and 3 (nuclear-marker rank) are DEGRADED to "
            "reporting-only, because a rank of 1 is vacuous without a comparison "
            "set. Cross-lineage statements cannot be made, so this is NOT a "
            "complete identity/capture assessment. Which of checks 4-8 were "
            "ACTUALLY evaluable in this run is listed in the dynamic "
            "'ASSESSABILITY (this run)' note at the end — do not assume they all ran.")

    target = cfg.identity_target or None
    if target is None and single_lineage:
        target = lineages[0]        # 唯一谱系即目标，使 5/6/8 可执行
    if target is not None and target not in set(lineages):
        res.notes.append(f"cfg.identity_target='{target}' 不在注释谱系中；"
                         f"检查 1/2/3 降级为报告模式。")
        target = None

    # ── 阳性/阴性标志 ──
    pos_m = list(cfg.pos_markers)
    neg_m = list(cfg.neg_markers)
    _alldet = _detect_many(ad, list(pos_m) + list(neg_m)
                           + list(cfg.identity_long_genes)
                           + [cfg.identity_nuclear_marker]
                           + list(cfg.identity_procurement_markers))
    pos_v = {g: _alldet.get(g) for g in pos_m}
    neg_v = {g: _alldet.get(g) for g in neg_m}
    pos_ok = [g for g in pos_m if pos_v[g] is not None]
    neg_ok = [g for g in neg_m if neg_v[g] is not None]
    if not pos_ok:
        res.verdict = "INDETERMINATE (no positive markers present)"
        res.notes.append(f"阳性标志 {pos_m} 在数据中一个都不存在。")
        return res
    if len(neg_ok) < len(neg_m):
        res.notes.append(f"阴性标志缺失 {sorted(set(neg_m) - set(neg_ok))}；"
                         f"门控按现有 {len(neg_ok)} 个执行，通过率会被高估。")

    npos = np.sum([pos_v[g] for g in pos_ok], axis=0)
    nneg = np.sum([neg_v[g] for g in neg_ok], axis=0)
    gate_pass = (npos >= cfg.min_pos_markers) & (nneg == 0)

    # ══════════ 1. 自家门控通过率（按谱系）══════════
    g1 = (pd.DataFrame({"lin": lin, "pass": gate_pass})
          .groupby("lin")["pass"]
          .agg(n="size", n_pass="sum")
          .assign(pct_pass=lambda d: 100 * d.n_pass / d.n)
          .sort_values("pct_pass", ascending=False))
    res.gate_by_lineage = g1
    if target is not None:
        res.target_gate_pass = float(g1.loc[target, "pct_pass"])
        others = g1.drop(index=target)["pct_pass"]
        if len(others) == 0:
            res.gate_margin = None      # 单谱系无法定义 margin（判定层会跳过该项）
        else:
            med = float(others.median())
            res.gate_margin = float(res.target_gate_pass / med) if med > 0 else np.inf

    # ══════════ 2. 复杂度排名 ══════════
    if cfg.identity_complexity_col in o.columns:
        g2 = (o.groupby(cfg.identity_col, observed=True)[cfg.identity_complexity_col]
              .median().sort_values()
              .to_frame("median_genes_per_cell"))
        g2["rank_low_to_high"] = range(1, len(g2) + 1)
        res.complexity_by_lineage = g2
        res.complexity_n_lineages = len(g2)
        if target is not None:
            res.complexity_rank = int(g2.loc[target, "rank_low_to_high"])
    else:
        res.notes.append(f"obs 中无 '{cfg.identity_complexity_col}'；检查 2/4 跳过。")

    # ══════════ 3. 核滞留标志 ══════════
    nm = cfg.identity_nuclear_marker
    vnm = _alldet.get(nm)
    if vnm is not None and target is not None:
        tab = (pd.DataFrame({"lin": lin, "d": vnm})
               .groupby("lin")["d"].mean().mul(100).sort_values())
        res.nuclear_marker_pct = float(tab.loc[target])
        res.nuclear_marker_rank = int(list(tab.index).index(target) + 1)
    elif vnm is None:
        res.notes.append(f"核标志 '{nm}' 不存在；检查 3 跳过。")

    # ══════════ 4. 复杂度匹配下的长基因检出 ══════════
    if (target is not None and cfg.identity_complexity_col in o.columns):
        ng = pd.Series(o[cfg.identity_complexity_col].values)
        edges = list(cfg.identity_complexity_bins)
        bins = pd.cut(ng, bins=[-np.inf] + edges + [np.inf], right=False)
        bin_arr = bins.values
        rows = []
        deficit = {}
        for b in bins.cat.categories:
            m = (bin_arr == b)
            mt, mo = m & (lin == target), m & (lin != target)
            if mt.sum() < cfg.identity_bin_min_cells or mo.sum() < cfg.identity_bin_min_cells:
                continue
            row = {"bin": str(b), "n_target": int(mt.sum()), "n_other": int(mo.sum())}
            for g in cfg.identity_long_genes:
                v = _alldet.get(g)
                if v is None:
                    continue
                pt, po = 100 * v[mt].mean(), 100 * v[mo].mean()
                row[f"{g}_target"] = round(pt, 2)
                row[f"{g}_other"] = round(po, 2)
                deficit.setdefault(g, []).append(pt - po)
            rows.append(row)
        if rows:
            res.matched_complexity = pd.DataFrame(rows)
            res.long_gene_deficit = {g: float(np.mean(d)) for g, d in deficit.items()}
        else:
            # C3（v0.7.1）：单谱系时 mo 恒空 → 所有复杂度分箱被跳过 → 检查 4 静默不落
            # 任何产物。本模块自述原则是"绝不静默跳过"，故显式登记原因。
            res.notes.append(
                "check 4 (long/low-copy genes at matched complexity) NOT evaluated in "
                "this run: no complexity bin contained both target and non-target "
                "cells (>= identity_bin_min_cells each), so no matched-complexity "
                "comparison exists (this is the expected single-lineage case).")

    # ══════════ 5. ★ 亚群不平衡（细胞类型内部）══════════
    src_col = cfg.identity_source_col
    if target is not None and src_col in o.columns:
        sub = o[[src_col, condition_col, donor_col]].copy()
        sub["_in"] = (lin == target)
        cells = sub[sub["_in"]]
        res.subset_by_arm = pd.crosstab(cells[src_col].astype(str),
                                        cells[condition_col].astype(str),
                                        margins=True)
        # 每亚群：跨臂计数 + 跨供者分布，与整体供者分布求 TVD
        src = cells[src_col].astype(str)
        overall_by_donor = cells.groupby(donor_col, observed=True).size()
        ct_arm = pd.crosstab(src, cells[condition_col].astype(str))
        for sv in ct_arm.index:
            in_s = (src == sv)
            by_donor = cells[in_s].groupby(donor_col, observed=True).size()
            res.subset_tvd[sv] = float(_tvd(by_donor, overall_by_donor))
            n_s = int(in_s.sum())
            if n_s < cfg.identity_imbalance_min_cells:
                continue
            # ★ 结构性检查：该亚群在某臂**完全缺席**（本例的核心发现）
            zero_arms = [a for a in ct_arm.columns if ct_arm.loc[sv, a] == 0]
            if zero_arms:
                res.imbalanced_subsets.append(
                    f"{sv}: {n_s} cells inside the annotated '{target}' cluster but "
                    f"absent from arm(s) {zero_arms}")
            # 运算符与审计二统一为**严格 >**（P2-d 已定）；原为 >=，
            # 使 2 对照供者时 TVD 恰为 0.5 的场景在审计四报、审计二不报。
            elif res.subset_tvd[sv] > cfg.identity_imbalance_min_tvd:
                res.imbalanced_subsets.append(
                    f"{sv}: {n_s} cells, donor-level TVD "
                    f"{res.subset_tvd[sv]:.3f} > {cfg.identity_imbalance_min_tvd}")
    elif target is None:
        # v0.3.3：多谱系 + 未指定目标谱系时，旧实现也走到下面的 else 分支，
        # 于是把原因写成“obs 中无 'source'”——即使该列**确实存在**（归因错误）。
        # 现分开归因：这是“未指定目标”，不是“缺列”。
        res.notes.append(
            f"检查 5/6/8 跳过：注释列有 {len(lineages)} 个谱系但**未指定目标谱系**"
            f"（配置 `identity_target` 或 CLI `--target-lineage`）；"
            f"“亚群是否局限于单臂”是本审计的核心问题，未指定目标时无法回答——"
            f"该状态**不构成**“身份与捕获合格”的证据。")
    else:
        res.notes.append(f"obs 中无 '{src_col}'；检查 5 跳过（亚群不平衡是本审计的核心项，"
                         f"建议提供来源/分选列）。")

    # ══════════ 6. ★ 门控规则的条件共线性 ══════════
    if target is not None and neg_ok:
        rows = []
        _absent_arms = []
        for g in neg_ok:
            v = neg_v[g]
            m_c = is_ctrl & (lin == target)
            m_k = (~is_ctrl) & (lin == target)
            if not m_c.any() or not m_k.any():
                # v0.3.3：目标谱系在某一臂完全缺席 → 不报“每臂检出率”（那会对空集合求
                # 均值，得到 NaN 与 numpy 警告，读者反而看不出“缺席”的事实）。
                # 与检查 5 的处理一致：显式说“absent from arm”。
                for _arm, _m in (("control", m_c), ("case", m_k)):
                    if not _m.any():
                        _absent_arms.append((g, _arm))
                continue
            pc = 100 * v[m_c].mean()
            pk = 100 * v[m_k].mean()
            rows.append({"marker": g, "control_pct": round(pc, 2),
                         "case_pct": round(pk, 2), "delta_pp": round(pk - pc, 2)})
        if _absent_arms:
            _g = sorted({g for g, _ in _absent_arms})
            _arms = sorted({a for _, a in _absent_arms})
            res.notes.append(
                f"检查 6：目标谱系 '{target}' 在 {_arms} 臂完全缺席——这些阴性标志的"
                f"逐臂检出率无定义（不是 0），已从表中排除而非报 NaN：{_g}。"
                f"“亚群局限于单一臂”本身就是审计四要报的结构性问题（见检查 5）。")
        nm_df = pd.DataFrame(rows).set_index("marker") if rows else pd.DataFrame(
            columns=["control_pct", "case_pct", "delta_pp"])
        res.negative_marker_by_arm = nm_df
        res.collinear_markers = sorted(
            nm_df.index[nm_df["delta_pp"].abs() >= cfg.identity_collinear_min_delta_pp]) \
            if len(nm_df) else []
        # 通过率按臂
        arm = (pd.DataFrame({"arm": np.where(is_ctrl, "control", "case"),
                             "pass": gate_pass, "lin": lin})
               .query("lin == @target")
               .groupby("arm")["pass"].agg(n="size", n_pass="sum"))
        arm["pct_pass"] = 100 * arm.n_pass / arm.n
        res.gate_pass_by_arm = arm

    # ══════════ 7. 模态匹配 ══════════
    if reference_modalities is not None:
        res.reference_modalities = sorted({str(m) for m in reference_modalities})
        if modality is not None:
            res.modality_matched = str(modality) in res.reference_modalities
        else:
            res.notes.append("未提供被审计数据集的 modality；检查 7 无法判定。")
    else:
        res.notes.append("未提供参照表模态列；检查 7 跳过。"
                         "强烈建议提供——跨模态比较会使 R 的绝对值不可解释。")

    # ══════════ 8. 取材方式标志按臂 ══════════
    if target is not None:
        rows = []
        _absent8 = []
        for g in cfg.identity_procurement_markers:
            v = _alldet.get(g)
            if v is None:
                continue
            m_c = is_ctrl & (lin == target)
            m_k = (~is_ctrl) & (lin == target)
            if not m_c.any() or not m_k.any():
                # 与检查 5/6 同口径（v0.3.3）：目标谱系在某臂缺席时，逐臂检出率无定义；
                # 不得对空集合求均值（会得到 NaN + numpy 警告）。
                _absent8.append(g)
                continue
            pc = 100 * v[m_c].mean()
            pk = 100 * v[m_k].mean()
            rows.append({"marker": g, "control_pct": round(pc, 2),
                         "case_pct": round(pk, 2), "delta_pp": round(pk - pc, 2)})
        if _absent8:
            res.notes.append(
                f"检查 8：目标谱系 '{target}' 在某一臂缺席，以下取材标志的逐臂检出率无定义"
                f"（已排除而非报 NaN）：{sorted(_absent8)}。")
        if rows:
            res.procurement_by_arm = pd.DataFrame(rows).set_index("marker")

    # ── C3（v0.6.0）：注释方法依赖声明 ──
    res.notes.append(
        "C3（annotation dependence）：本层的身份结论**以输入注释列为前提**；作者注释/门控本身的误差"
        "不由本工具排除（自动注释方法的准确率差异见 Abdelaal 2019, PMID 31500660）。")
    # ── C1（v0.6.0）：ambient RNA 局限声明 ──
    res.notes.append(
        "C1（ambient RNA）：未做去污染时，环境 RNA 会让错簇出现'标记信号'，污染身份判定与检出率"
        "（SoupX PMID 33367645 / DecontX PMID 32138770 / Janssen 2023 PMID 37337297）。"
        "本工具不校正 ambient RNA，该局限不由本层排除。")
    # ── C2（v0.6.0）：双细胞率报告（已采集则用，否则不编造）──
    _obs = ad.obs
    for _dc in ("predicted_doublet", "doublet_score"):
        if _dc in _obs.columns:
            try:
                if _dc == "predicted_doublet":
                    _v = _obs[_dc].astype(str).str.lower().isin(["true", "1", "yes"]).values
                    _pc = float(_v[is_ctrl].mean() * 100) if is_ctrl.sum() else float("nan")
                    _pk = float(_v[~is_ctrl].mean() * 100) if (~is_ctrl).sum() else float("nan")
                    res.notes.append(
                        f"C2（doublets）：预测双细胞率 control {_pc:.1f}% / case {_pk:.1f}%"
                        "（Scrublet PMID 30954476 / DoubletFinder PMID 30954475）；双臂差异大时可能污染身份判定。")
                else:
                    _m = pd.to_numeric(_obs[_dc], errors="coerce")
                    _pc = float(_m[is_ctrl].median()) if is_ctrl.sum() else float("nan")
                    _pk = float(_m[~is_ctrl].median()) if (~is_ctrl).sum() else float("nan")
                    res.notes.append(
                        f"C2（doublets）：`{_dc}` 中位 control {_pc:.3f} / case {_pk:.3f}"
                        "（Scrublet / DoubletFinder）；差异大时需检查双细胞混杂。")
            except Exception as _e:  # noqa: BLE001 —— 显式 note，绝不静默跳过
                res.notes.append(
                    f"C2（doublets）：列 `{_dc}` 解析失败"
                    f"（{type(_e).__name__}: {_e}）——双细胞率不可得（未静默跳过）。")
            break

    # ══════════════════ 判定层 ══════════════════
    flags: List[str] = []

    # 捕获信号须**分别计数**：单靠"复杂度排名最低"是弱信号——同模态健康心脏中
    # 同模态健康心脏中心肌细胞也可能排在低位。故 CAPTURE-LIMITED
    # 要求低复杂度与低核滞留信号**同时**出现：单一读数不构成推断，组合才承载结论。
    complexity_low = bool(res.complexity_rank is not None
                          and res.complexity_rank == 1
                          and res.complexity_n_lineages >= 4)
    # ⚠️ 须要求谱系数 >= 2：只有一个注释谱系时 nuclear_marker_rank 恒为 1，
    # "核标志最低"是**空洞的真**（没有可比较的对象）→ 会误报 nuclear-poor capture。
    # 用 >=2 而非 >=4：排名只要求在 >=2 个谱系间才有意义；而整体 CAPTURE-LIMITED
    # 判定仍需 complexity_low（其自身要求 >=4 个谱系），故不会在少谱系数据上误判。
    nuclear_low = bool(res.nuclear_marker_rank == 1
                       and res.nuclear_marker_pct is not None
                       and res.complexity_n_lineages >= 2)
    # 核标志缺失时 nuclear_low 恒为 False → CAPTURE-LIMITED 在**本轮不可达**。
    # 若只在 else 分支说 "single capture signal present"，会误导读者以为
    # 第二个信号被评估过但未达标。故须显式声明"不可评估"。
    nuclear_marker_present = res.nuclear_marker_pct is not None

    if res.gate_margin is not None and res.gate_margin < cfg.identity_gate_margin:
        flags.append(f"gate is not lineage-specific (margin {res.gate_margin:.1f}x "
                     f"< {cfg.identity_gate_margin}x)")
    if complexity_low:
        flags.append(f"annotated '{target}' has the LOWEST median genes/cell of "
                     f"{res.complexity_n_lineages} lineages -> capture signal (weak alone)")
    if nuclear_low:
        flags.append(f"nuclear marker '{nm}' lowest in '{target}' -> nuclear-poor capture")
    elif res.nuclear_marker_rank == 1 and res.nuclear_marker_pct is not None:
        # 排名为 1 但谱系数 < 2 → 空洞的真，不构成信号
        res.notes.append(
            f"nuclear marker '{nm}' ranks 1 but only "
            f"{res.n_lineages} annotated lineage(s) -> rank is vacuous, "
            "not reported as a capture signal")
    if res.imbalanced_subsets:
        flags.append(f"subset imbalance within '{target}': {res.imbalanced_subsets}")
    if res.collinear_markers:
        flags.append(f"gate rule is condition-collinear on {res.collinear_markers}")
    if res.modality_matched is False:
        flags.append(f"no same-modality reference ({modality} vs "
                     f"{res.reference_modalities}) -> R values are cross-modal")

    # C4（v0.7.2，据第三方独立评估）：门控通过率的**绝对量**与**按臂差**判据。
    # 旧实现只比较“目标/其他谱系中位数”的相对 margin，单谱系时 margin=None → 门控项
    # 永不触发；于是通过率 0.3%、两臂差 3.4× 的门控规则仍判 PASS（漏报）。
    _gate_floor_bad = bool(
        res.target_gate_pass is not None
        and res.target_gate_pass < cfg.identity_gate_pass_floor_pct)
    _gate_arm_bad = False
    _gate_arm_ratio = None
    if res.gate_pass_by_arm is not None and len(res.gate_pass_by_arm):
        _pcts = res.gate_pass_by_arm["pct_pass"].astype(float)
        _lo, _hi = float(_pcts.min()), float(_pcts.max())
        if _lo > 0:
            _gate_arm_ratio = _hi / _lo
            _gate_arm_bad = _gate_arm_ratio > cfg.identity_gate_arm_delta_max
        elif _hi > 0:
            _gate_arm_ratio = float("inf")
            _gate_arm_bad = True
    if _gate_floor_bad:
        flags.append(
            f"own-gate pass rate {res.target_gate_pass:.3f}% < floor "
            f"{cfg.identity_gate_pass_floor_pct}% -> the gate is effectively "
            "non-functional in this dataset")
    if _gate_arm_bad:
        _a = res.gate_pass_by_arm
        _fmt = lambda arm: (f"{_a.loc[arm, 'pct_pass']:.3f}%"
                            if arm in _a.index else "n/a")
        flags.append(
            f"gate pass rate is condition-collinear: differs {_gate_arm_ratio:.1f}x "
            f"between arms (case {_fmt('case')} / control {_fmt('control')})")

    if not flags:
        res.verdict = "PASS (identity and capture adequate)"
    elif any("subset imbalance" in f or "condition-collinear" in f for f in flags):
        res.verdict = "COMPOSITION-CONFOUNDED (within cell type)"
    elif complexity_low and nuclear_low:
        res.verdict = "CAPTURE-LIMITED"
    elif _gate_floor_bad:
        # 门控在本数据集近乎失效（通过率 < 绝对下限）→ 不足以判 PASS
        res.verdict = "GATE-NOT-USABLE (own-gate pass rate below floor)"
    else:
        # 只有单一捕获信号（或仅有门控/模态类 note）→ 不足以判 CAPTURE-LIMITED
        if not nuclear_marker_present:
            res.notes.append(
                f"CAPTURE-LIMITED **not assessable** in this run: nuclear marker "
                f"'{cfg.identity_nuclear_marker}' is absent, so only one of the two "
                "required capture signals can ever be evaluated")
        else:
            res.notes.append(
                "single capture signal present (below the 2-signal threshold) -> not "
                "reported as CAPTURE-LIMITED; see docs/methods.md §6.3 on the "
                "weak-signal caveat")
        res.verdict = "PASS (identity and capture adequate)"
    res.notes.extend(flags)
    # ── C3（v0.7.1，据第三方独立评估）：**动态**报告本次运行实际可评估的检查 ──
    # 旧实现硬编码 "Checks 4/5/6/7/8 remain FULLY EVALUATED"（在单谱系分支），
    # 而检查 4/5/7 在缺列/单谱系时会跳过 → 文案与产物矛盾。此处改成按**实际结果**
    # 生成，并存原则：`not evaluable` ≠ `passed`。
    _chk = [
        ("1 gate pass-rate by lineage", res.gate_by_lineage is not None),
        ("2 complexity rank", res.complexity_rank is not None),
        ("3 nuclear-marker rank", res.nuclear_marker_rank is not None),
        ("4 matched-complexity long genes", res.matched_complexity is not None),
        ("5 subset imbalance", res.subset_by_arm is not None),
        ("6 gate collinearity", res.negative_marker_by_arm is not None),
        ("7 modality match", res.modality_matched is not None),
        ("8 procurement by arm", res.procurement_by_arm is not None),
    ]
    _evald = [n for n, ok in _chk if ok]
    _skipped = [n for n, ok in _chk if not ok]
    res.notes.append(
        "ASSESSABILITY (this run): evaluable = "
        + (", ".join(_evald) if _evald else "none")
        + "; NOT evaluable = "
        + (", ".join(_skipped) if _skipped else "none")
        + " (reasons in the notes above; a NOT-evaluable check is not a passed check).")
    return res
