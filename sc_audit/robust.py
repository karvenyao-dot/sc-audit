"""审计三：稳健信号对照。

分离真信号与深度/技术假象，并给出正面对照（避免"全盘否定"）。

修订记录（2026-10-01）：
  · S3a 判定不再把"深度可比"与"分离持续"混为一谈；k=0 时给出 warn/bad 而非 ok
  · S3b "分离持续"要求**方向一致**（原代码只看 |Δ|>5，符号翻转也算持续）
  · S3c 分层要求每箱两组各 >= cfg.min_cells_per_stratum（默认 30），
        并用 Wilson CI 不重叠作为"分离"判据（原 5 细胞/箱、5pp 阈值无统计意义）
  · S3d 移除 try/except pass，深度指标计算失败时显式记录并降级判定
  · 4.6 mito% 差异作为独立档 "CELL-STATE DIFFERS" 报告，且用 pp 而非相对 %
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import AuditConfig
from .io import gene_detection_vector, compute_qc_metrics


def _wilson(k: int, n: int, level: float = 0.95) -> Tuple[float, float]:
    if n <= 0:
        return (np.nan, np.nan)
    from scipy.stats import norm
    import math
    z = norm.ppf(1 - (1 - level) / 2)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return (max(0.0, centre - half) * 100, min(1.0, centre + half) * 100)


@dataclass
class RobustSignalResult:
    depth_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    depth_delta_pct: Dict[str, float] = field(default_factory=dict)
    depth_delta_pp: Dict[str, float] = field(default_factory=dict)
    stratified: Dict[str, pd.DataFrame] = field(default_factory=dict)
    genes_persisting: List[str] = field(default_factory=list)
    strata_ok: Dict[str, str] = field(default_factory=dict)   # {gene: "k/n strata"}
    housekeeping: Optional[pd.DataFrame] = None
    verdict: str = ""
    verdict_class: str = "warn"        # ok | warn | bad（供报告着色）
    depth_verdict: str = ""
    cell_state_verdict: str = ""
    housekeeping_cv: float = float("nan")      # B2：各管家 R 的离散度（CV）
    depth_metric_sensitivity: Dict[str, object] = field(default_factory=dict)  # B3
    errors: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def to_text(self) -> str:
        L = ["== Audit 3 · Robust-signal control ==", "",
             "-- depth comparability --", self.depth_table.to_string(), ""]
        L.append(f"depth verdict   : {self.depth_verdict}")
        if self.cell_state_verdict:
            L.append(f"cell-state      : {self.cell_state_verdict}")
        L.append("")
        for g, t in self.stratified.items():
            L.append(f"-- depth-stratified detection: {g} --")
            L.append(t.to_string(index=False))
            L.append("")
        if self.housekeeping is not None and len(self.housekeeping):
            L.append("-- housekeeping internal control --")
            L.append(self.housekeeping.to_string(index=False))
            L.append("")
        L.append(f"VERDICT: {self.verdict}")
        for e in self.errors:
            L.append(f"  ! {e}")
        for n in self.notes:
            L.append(f"  · {n}")
        return "\n".join(L)


def _hk_msg(modality_matched: Optional[bool]) -> str:
    """管家护栏的解读文案 —— 由模态匹配结论决定，不再让读者自己判断。"""
    base = ("housekeeping genes detected at comparable or higher rates than the "
            "reference(s)")
    if modality_matched is True:
        return (base + ". The comparison is **within a modality**, so this argues "
                "against cell-type misassignment and global library degradation")
    if modality_matched is False:
        # ⚠️ 措辞限定在**本比较框架**内 —— 不得写成
        # "跨模态一定不能解释质量" 这类普适规律，只能写"在该框架下更符合
        # 模态差异的预期，因而不能直接作为质量证据"。
        return (base + ". The comparison is **cross-modal**. Within this comparison "
                "framework, a ratio >= 1 is more consistent with the expected "
                "modality difference (cytoplasmic transcripts are intrinsically more "
                "abundant per cell in whole-cell than in nuclear preparations) than "
                "with variation in library quality, and therefore **cannot be used "
                "directly as evidence about library quality**. This is a statement "
                "about the present framework, not a general rule.")
    return (base + ". **Within a modality** this argues against cell-type "
            "misassignment and global library degradation; **across modalities** it "
            "does not — cytoplasmic transcripts are expected to be enriched in "
            "whole-cell relative to single-nucleus preparations, so a ratio >= 1 is "
            "then the expected modality signature, not evidence about library "
            "quality (see docs/methods.md §4)")


def audit_robust_signal(ad, genes: Sequence[str], condition_col: str,
                        donor_col: str, control_label: str,
                        housekeeping: Optional[Sequence[str]] = None,
                        reference: Optional[Dict[str, float]] = None,
                        cfg: Optional[AuditConfig] = None,
                        modality_matched: Optional[bool] = None,
                        ) -> RobustSignalResult:
    """modality_matched: 由审计四算出的模态匹配结论。

    True  -> 同模态：管家 R >= 1 支持"无文库劣化"。
    False -> 跨模态：管家 R >= 1 是模态差的预期表现，**不是**质量证据。
    None  -> 未知：两种解读并列输出（原行为）。
    """
    cfg = cfg or AuditConfig()
    # 给了参照却未声明模态 → 会静默落入 None 分支（两解并列），管家护栏的
    # 结论强度无法确定。显式告警，避免"看起来跑了但没有定论"。
    if reference and modality_matched is None:
        import warnings as _w
        _w.warn(
            "reference was supplied but modality_matched was not declared: the "
            "housekeeping guard will report both readings side by side instead "
            "of a conclusion. Pass modality_matched=True if the reference is "
            "same-modality, False otherwise.", stacklevel=2)
    ad = compute_qc_metrics(ad)
    o = ad.obs
    cond = o[condition_col].astype(str).values
    donor = o[donor_col].astype(str).values
    # P2-3：API 路径也须校验供者 ⊂ 组
    from .io import assert_donors_nested
    assert_donors_nested(donor, cond, donor_col, condition_col)
    is_ctrl = cond == str(control_label)
    res = RobustSignalResult()

    have_mito = "percent_mito" in o.columns and o["percent_mito"].notna().any()

    # ── 深度匹配（显式计算，不用 try/except 静默吞错）──
    if is_ctrl.sum() == 0 or (~is_ctrl).sum() == 0:
        res.errors.append("one of the two groups is empty — depth metrics not computable")
    else:
        aggs = {"n": ("n_genes", "size"),
                "n_genes_median": ("n_genes", "median"),
                "total_counts_median": ("total_counts", "median")}
        if have_mito:
            aggs["percent_mito_median"] = ("percent_mito", "median")
        try:
            dt = (o.assign(_g=np.where(is_ctrl, "control", "case"))
                    .groupby("_g").agg(**aggs).T)
            res.depth_table = dt
            for metric, key in [("n_genes_median", "genes_per_cell"),
                                ("total_counts_median", "counts_per_cell")]:
                c, k = float(dt.loc[metric, "control"]), float(dt.loc[metric, "case"])
                res.depth_delta_pct[key] = round((k - c) / c * 100, 1) if c else np.nan
                res.depth_delta_pp[key] = round(k - c, 2)
            if have_mito:
                c = float(dt.loc["percent_mito_median", "control"])
                k = float(dt.loc["percent_mito_median", "case"])
                res.depth_delta_pct["percent_mito"] = round((k - c) / c * 100, 1) if c else np.nan
                res.depth_delta_pp["percent_mito"] = round(k - c, 2)
        except Exception as exc:  # 显式记录，不静默
            res.errors.append(f"depth metric computation failed: {type(exc).__name__}: {exc}")

    # ── 深度判定（S3a：不默认通过）──
    depth_metrics = {k: v for k, v in res.depth_delta_pct.items()
                     if k != "percent_mito" and v == v and not np.isnan(v)}
    if not depth_metrics:
        res.depth_verdict = "INSUFFICIENT (depth metrics unavailable)"
    elif any(abs(v) >= cfg.depth_delta_max_pct for v in depth_metrics.values()):
        res.depth_verdict = ("PARTIAL (depth differs: "
                             + ", ".join(f"{k}={v:+.1f}%" for k, v in depth_metrics.items()) + ")")
    else:
        res.depth_verdict = ("DEPTH-COMPARABLE ("
                             + ", ".join(f"{k}={v:+.1f}%" for k, v in depth_metrics.items()) + ")")

    # ── 细胞状态（mito，独立档；用 pp）──
    if have_mito:
        dv_pp = res.depth_delta_pp.get("percent_mito", np.nan)
        if dv_pp == dv_pp and abs(dv_pp) >= cfg.mito_delta_warn_pp:
            res.cell_state_verdict = (
                f"CELL-STATE DIFFERS (mitochondrial fraction {dv_pp:+.1f} pp) — "
                "may dilute non-mitochondrial UMI detection; interpret gene-level "
                "detection rates with care")

    # ── 深度分层（S3b/S3c：方向一致 + Wilson CI 不重叠）──
    ng = o["n_genes"].values
    edges = list(cfg.depth_edges)
    for g in genes:
        v = gene_detection_vector(ad, g)
        if v is None:
            continue
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (ng > lo) & (ng <= hi)
            mc, mk = m & is_ctrl, m & (~is_ctrl)
            if mc.sum() >= cfg.min_cells_per_stratum and mk.sum() >= cfg.min_cells_per_stratum:
                kc, kk = int(v[mc].sum()), int(v[mk].sum())
                nc, nk = int(mc.sum()), int(mk.sum())
                cl, ch = _wilson(kc, nc, cfg.strat_ci_level)
                kl, kh = _wilson(kk, nk, cfg.strat_ci_level)
                disjoint = (ch < kl) or (kh < cl)
                rows.append({
                    "stratum": f"({lo},{hi}]",
                    "ctrl_pct": round(kc / nc * 100, 1), "n_ctrl": nc,
                    "case_pct": round(kk / nk * 100, 1), "n_case": nk,
                    "delta_pp": round((kc / nc - kk / nk) * 100, 1),
                    "ci_disjoint": bool(disjoint),
                })
        if rows:
            res.stratified[g] = pd.DataFrame(rows)

    # 阳性对照基因集合（管家基因）——先定义，供下文两处使用
    hk = list(housekeeping) if housekeeping is not None else list(cfg.housekeeping_genes)
    hk_upper = {str(g).upper() for g in hk}

    # 分离持续性：方向一致 + 幅度达标 + CI 不重叠比例达标（L2：阈值入 config）
    for g, t in res.stratified.items():
        # P3：管家基因是阳性对照，不是被检验的假设 —— 不计入“分离持续”
        # （其组间差异反映两组细胞状态/深度，而非生物学信号）
        if str(g).upper() in hk_upper:
            res.notes.append(
                f"{g} excluded from the persistence count (housekeeping / positive control)")
            continue
        d = t["delta_pp"]
        same_dir = bool((d > 0).all() or (d < 0).all())
        big = bool((d.abs() >= cfg.strat_min_delta_pp).all())
        frac = float(t["ci_disjoint"].mean())
        ci_ok = bool(frac >= cfg.strat_ci_frac_ok)
        if same_dir and big and ci_ok:
            res.genes_persisting.append(g)
            res.strata_ok[g] = f"{int(t['ci_disjoint'].sum())}/{len(t)}"

    # ── B3（v0.6.0）：分层口径敏感性（n_genes ↔ total_counts）──
    # 历史口径用 n_genes 分层；但 n_genes 与目标测量并非正交（Hafemeister 2019
    # PMID 31870423）。此处**不改主判定**，只并列报告另一口径下的持续性是否一致。
    _tc = o["total_counts"].values if "total_counts" in o.columns else None
    if _tc is not None and len(res.stratified):
        _alt: List[str] = []
        for g in res.stratified:
            if str(g).upper() in hk_upper:
                continue
            v = gene_detection_vector(ad, g)
            if v is None:
                continue
            rows2 = []
            for lo, hi in zip(edges[:-1], edges[1:]):
                m = (_tc > lo) & (_tc <= hi)
                mc, mk = m & is_ctrl, m & (~is_ctrl)
                if (mc.sum() >= cfg.min_cells_per_stratum
                        and mk.sum() >= cfg.min_cells_per_stratum):
                    kc, kk = int(v[mc].sum()), int(v[mk].sum())
                    nc, nk = int(mc.sum()), int(mk.sum())
                    cl, ch = _wilson(kc, nc, cfg.strat_ci_level)
                    kl, kh = _wilson(kk, nk, cfg.strat_ci_level)
                    rows2.append(((kc / nc - kk / nk) * 100, (ch < kl) or (kh < cl)))
            if not rows2:
                continue
            dd = np.array([r[0] for r in rows2])
            ddj = np.array([r[1] for r in rows2])
            _ok = (bool((dd > 0).all() or (dd < 0).all())
                   and bool((np.abs(dd) >= cfg.strat_min_delta_pp).all())
                   and bool(ddj.mean() >= cfg.strat_ci_frac_ok))
            if _ok:
                _alt.append(g)
        res.depth_metric_sensitivity = {
            "primary_metric": "n_genes",
            "alt_metric": "total_counts",
            "primary_persisting": list(res.genes_persisting),
            "alt_persisting": _alt,
            "consistent": set(_alt) == set(res.genes_persisting),
        }
        if set(_alt) != set(res.genes_persisting):
            res.notes.append(
                "depth-stratification metric sensitivity: the persistence verdict DIFFERS "
                "between n_genes and total_counts strata — treat 'separation persists' as "
                "metric-dependent (Hafemeister 2019, PMID 31870423).")
        else:
            res.notes.append(
                "depth-stratification metric sensitivity: n_genes and total_counts strata "
                "give the SAME persistence verdict.")

    # ── 管家基因内部对照 ──
    if reference is not None and hk:
        # P1-c：管家基因的 R 必须与目标基因用**同一分子口径**（cfg.ratio_basis）；
        # 原实现固定用整体检出率，与目标基因（默认 case）不可比。
        _cond = o[condition_col].astype(str).values
        _ctrl = _cond == str(control_label)
        rows = []
        for g in hk:
            v = gene_detection_vector(ad, g)
            if v is None:
                continue
            basis_pct = {
                "case": (v[~_ctrl].mean() * 100 if (~_ctrl).sum() else np.nan),
                "overall": v.mean() * 100,
                "control": (v[_ctrl].mean() * 100 if _ctrl.sum() else np.nan),
            }.get(cfg.ratio_basis, v.mean() * 100)
            pct = basis_pct
            ref = reference.get(g, np.nan)
            rows.append({"gene": g,
                         "basis": cfg.ratio_basis,
                         "pct_in_dataset": round(pct, 2) if pct == pct else np.nan,
                         "best_reference_pct": round(ref, 2) if not np.isnan(ref) else np.nan,
                         "ratio_R": round(pct / ref, 3) if (ref and ref > 0) else np.nan})
        if rows:
            res.housekeeping = pd.DataFrame(rows)

    # ── 综合判定（S3a：k=0 不得为 ok）──
    n_g = len(res.stratified)
    n_p = len(res.genes_persisting)
    if res.depth_verdict.startswith("DEPTH-COMPARABLE") and n_g > 0 and n_p > 0:
        detail = ", ".join(f"{g}: {res.strata_ok.get(g, '?')} strata CI-disjoint"
                           for g in res.genes_persisting)
        res.verdict = (f"DEPTH-COMPARABLE — separation persists (direction-consistent, "
                       f"magnitude >= {cfg.strat_min_delta_pp} pp, "
                       f"CI-disjoint in >= {cfg.strat_ci_frac_ok:.0%} of strata) "
                       f"for {n_p}/{n_g} gene(s) tested [{detail}]")
        res.verdict_class = "ok"
    elif res.depth_verdict.startswith("DEPTH-COMPARABLE") and n_g > 0 and n_p == 0:
        res.verdict = (f"DEPTH-COMPARABLE but separation does NOT persist consistently "
                       f"(0/{n_g} genes met direction + magnitude + CI criteria) — "
                       "the separation may be depth-dependent")
        res.verdict_class = "warn"
    elif res.depth_verdict.startswith("DEPTH-COMPARABLE"):
        res.verdict = ("DEPTH-COMPARABLE but no depth stratum met the minimum cell count "
                       f"({cfg.min_cells_per_stratum}/group) — separation could not be assessed")
        res.verdict_class = "warn"
    elif res.depth_verdict.startswith("PARTIAL"):
        res.verdict = ("PARTIAL — depth differs between groups; interpret stratified "
                       "results with care")
        res.verdict_class = "warn"
    else:
        res.verdict = "INSUFFICIENT — depth comparability could not be established"
        res.verdict_class = "bad"

    if res.housekeeping is not None and len(res.housekeeping):
        hk_r = res.housekeeping["ratio_R"].dropna()
        if len(hk_r) and (hk_r >= cfg.housekeeping_ratio_ok).all():
            res.notes.append(
                _hk_msg(modality_matched))
        elif len(hk_r):
            res.notes.append(
                f"housekeeping ratio(s) below {cfg.housekeeping_ratio_ok}: "
                f"{hk_r.round(3).tolist()} — cell-type assignment or library quality "
                "should be re-examined before attributing low detection to a floor effect")
        # ── B2（v0.6.0）：多管家一致性 + 适用域声明 ──
        if len(hk_r) >= 2 and float(hk_r.mean()) != 0:
            _cv = float(hk_r.std(ddof=0) / hk_r.mean())
            res.housekeeping_cv = round(_cv, 3) if _cv == _cv else float("nan")
            if _cv > float(getattr(cfg, "housekeeping_cv_warn", 0.5)):
                res.notes.append(
                    f"housekeeping genes are mutually inconsistent (CV={_cv:.3f} > "
                    f"{cfg.housekeeping_cv_warn}): the single-threshold guardrail reading "
                    "is unreliable here — inspect individual housekeeping genes "
                    "(housekeeping stability is tissue/depth dependent, geNorm-style M-value).")
        if res.housekeeping is not None and "pct_in_dataset" in res.housekeeping.columns:
            _base = res.housekeeping["pct_in_dataset"].dropna() / 100.0
            if len(_base) and float(_base.min()) >= float(
                    getattr(cfg, "housekeeping_high_baseline", 0.8)):
                res.notes.append(
                    "housekeeping-applicability: all housekeeping genes have high baseline "
                    "detection (>= 80%); at these baselines a global count reduction (the "
                    "S2 scenario) tends to SATURATE and may not trigger the guardrail "
                    "(R_S2 = (1 - sqrt(1 - p)) / p -> 1 as p -> 1). The guardrail is less "
                    "sensitive in this tissue — declare this applicability limit.")
    if res.cell_state_verdict:
        res.notes.append(res.cell_state_verdict)
    return res
