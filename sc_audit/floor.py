"""审计一：检出地板审计。

判断目标基因在本数据集中是否【具备定量能力】，独立于正在被声称的组间比较。

修订记录（2026-10-01）：
  · S6/N4  verdict 与 R 使用**同一分子口径**（唯一字段 cfg.ratio_basis；
        原 verdict_basis 已删除），无参照时不再默认 QUANTIFIABLE，
        改为 INSUFFICIENT (no reference)
  · S7  unquantifiable property 与 summary['n_unquantifiable'] 统一按 verdict 定义
  · S4  二项检验改为仅在 case_pct == 0 时报告，并显式标注其细胞独立假设
  · S5  新增供者加权效应量（置换统计量的观测值），与 perm_P 配套
  · 4.9 BH-FDR 多重比较校正（原工具缺失）
  · 4.3/L6  组合数 <= cfg.exact_max_combinations 且单位数 <= cfg.exact_enum_max_units
        时**精确枚举**（消除种子依赖；原参数名 exact_max_donors 已改名）
  · 4.3 双侧 leave-one-donor（对照 + 病例）
  · 4.10 Wilson/Clopper-Pearson CI（事件有限层）
  · 4.10 事件稀疏阈值三档敏感性（10/20/50）
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import warnings

import numpy as np
import pandas as pd

from .reference import CrossModalityError, DepthMismatchError

from .config import AuditConfig
from .io import gene_detection_vector


# ══════════════════ 统计工具 ══════════════════

def _wilson_ci(k: int, n: int, level: float = 0.95) -> Tuple[float, float]:
    """Wilson score interval（比例）。n=0 时返回 (nan, nan)。"""
    if n <= 0:
        return (np.nan, np.nan)
    from scipy.stats import norm
    z = norm.ppf(1 - (1 - level) / 2)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return (max(0.0, centre - half) * 100, min(1.0, centre + half) * 100)


def _cp_upper(k: int, n: int, level: float = 0.95) -> float:
    """Clopper–Pearson **单侧**上界（L3 修复：统一单侧口径）。

    ⚠️ 两种区间口径不得混用（外部审查 B1）：容差表（±43.4% → x_min 随 n 变化）
    用的是**双侧**相对半宽 (upper−lower)/2/p̂；本函数给的是**单侧 95% 上界**
    （≡ 双侧 90% 区间的上界）。两者是不同的量，不可互相代入。
    复算范式（n_case = 200）：k = 19/20/21 的双侧相对半宽 = 0.454 / 0.440 / 0.428；
    若误用单侧列会得到偏小的 x_min。详见 docs/methods.md §2.5。


    原实现对 k=0 用单侧、k>0 用双侧 97.5% 分位，同一列混两种口径；
    现统一为单侧 beta.ppf(level, k+1, n-k)（k=0 时退化为 1-(1-level)^(1/n)）。
    该列语义为"以此样本量，真实比例至多约 xx%"。
    """
    if n <= 0:
        return np.nan
    from scipy.stats import beta
    if k <= 0:
        return float((1 - (1 - level) ** (1 / n)) * 100)
    if k >= n:
        # 全部阳性：Clopper–Pearson 单侧上界恰为 1（无法排除 p=1）。
        # 原实现 b = n - k = 0 → beta.ppf(..., 0) = NaN，该类基因上界列全空。
        return 100.0
    return float(beta.ppf(level, k + 1, n - k)) * 100


def _cp_bounds(k: int, n: int, level: float = 0.95) -> Tuple[float, float]:
    """Clopper–Pearson 精确**双侧**区间（返回比例 0–1）。

    这是预注册容差表所用的口径（相对半宽 = (upper−lower)/2/p̂）。
    与 `_cp_upper`（**单侧**上界，≡ 双侧 90% 上界）是**不同的量**，不得混用。
    """
    if n <= 0:
        return (np.nan, np.nan)
    from scipy.stats import beta
    a = 1 - level
    lo = 0.0 if k <= 0 else float(beta.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k >= n else float(beta.ppf(1 - a / 2, k + 1, n - k))
    return lo, hi


def relative_halfwidth_cp(k: int, n: int, level: float = 0.95) -> float:
    """双侧相对半宽 (upper − lower) / 2 / p̂ —— 预注册容差口径。

    p̂ = k/n。k = 0 或 n ≤ 0 时不可定义 → NaN。
    """
    if n <= 0 or k <= 0:
        return np.nan
    lo, hi = _cp_bounds(k, n, level)
    p = k / n
    return float((hi - lo) / 2 / p)


def derive_event_threshold(n: int, tolerance: float, level: float = 0.95) -> Optional[int]:
    """达到给定**双侧相对精度**容差所需的最小阳性细胞数 x_min（预注册容差口径）。

    本函数的目的是把"锁定 20"从散文变成**可执行、可复算**的推导，例如
      derive_event_threshold(200, 0.434) == 21
    ⚠️ 它**不**参与任何自动调参——`event_sparse` 仍是事前锁定值，本函数只用于
    文档校核与测试断言（防"数字与口径错配"，例如把 ±50% 列的 x_min 当成 ±43.4% 的）。
    """
    if n <= 0:
        return None
    for x in range(1, n + 1):
        r = relative_halfwidth_cp(x, n, level)
        if r == r and r <= tolerance:
            return x
    return None


def event_threshold_table(n: int, level: float = 0.95,
                          tolerances: Sequence[float] = (0.50, 0.434, 0.40, 0.33, 0.25)
                          ) -> Dict[float, Optional[int]]:
    """`容差 → x_min` 表（代码侧单一实现；x_min 随 n 变化）。"""
    return {float(t): derive_event_threshold(n, float(t), level) for t in tolerances}


def _bh_fdr(pvals: Sequence[float]) -> List[float]:
    """Benjamini–Hochberg q 值。NaN 不参与 m 的计数。"""
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    q = np.full(len(p), np.nan)
    if m == 0:
        return q.tolist()
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    qv = ranked * m / np.arange(1, m + 1)
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    q[order] = np.clip(qv, 0, 1)
    return q.tolist()


def _donor_means(v: np.ndarray, donor: np.ndarray) -> Tuple[np.ndarray, list]:
    """返回 (每供者检出均值, 排序后的供者列表)。"""
    donors = sorted(set(donor))
    D = np.array([v[donor == d].mean() for d in donors])
    return D, donors


def perm_test(v: np.ndarray, donor: np.ndarray, is_ctrl: np.ndarray,
              n_perm: int = 10000, seed: int = 42,
              exact_max_units: int = 40,
              exact_max_combinations: int = 5_000_000) -> Dict:
    """供者级置换检验。

    组合数 <= exact_max_combinations 且单位数 <= exact_max_units 时**精确枚举**，
    否则随机置换。（L6 修复：原 exact_max_donors=12 被 max(...,40) 覆盖，
    实际门只有组合数，参数名误导；现参数名与语义一致。）
    返回 dict：p / observed_effect（供者级）/ n_units / mode / n_combinations
    """
    D, donors = _donor_means(v, donor)
    ctrl = np.array([bool(is_ctrl[donor == d][0]) for d in donors])
    n, k = len(donors), int(ctrl.sum())
    if k == 0 or k == n or n < 2:
        return {"p": np.nan, "observed_effect": np.nan, "n_units": n,
                "mode": "undefined", "n_combinations": 0}

    obs = abs(D[~ctrl].mean() - D[ctrl].mean())

    n_comb = math.comb(n, k)
    if n_comb <= exact_max_combinations and n <= exact_max_units:
        cnt = 0
        for idx in itertools.combinations(range(n), k):
            m = np.zeros(n, dtype=bool)
            m[list(idx)] = True
            if abs(D[~m].mean() - D[m].mean()) >= obs - 1e-12:
                cnt += 1
        return {"p": cnt / n_comb, "observed_effect": float(obs), "n_units": n,
                "mode": "exact", "n_combinations": n_comb}

    rng = np.random.default_rng(seed)
    cnt = 0
    for _ in range(n_perm):
        i = rng.permutation(n)
        if abs(D[i][~ctrl].mean() - D[i][ctrl].mean()) >= obs - 1e-12:
            cnt += 1
    return {"p": (1 + cnt) / (1 + n_perm), "observed_effect": float(obs),
            "n_units": n, "mode": "random", "n_combinations": n_perm}


def leave_one_donor(v: np.ndarray, donor: np.ndarray, is_ctrl: np.ndarray,
                    n_perm: int = 10000, seed: int = 42,
                    exact_max_units: int = 40,
                    exact_max_combinations: int = 5_000_000,
                    sides: str = "both") -> Dict:
    """留一供者敏感性（**整供者级** leave-one-out，非细胞级）。

    · 单位是**供者**：每次剔除一个完整供者（其全部细胞一并移出），重新做供者级置换；
      这不等于"随机剔除同等数量的细胞"，因此其 P 值范围反映的是**供者层**的敏感性。
    · sides='both'（默认）同时扫对照臂与病例臂；'control' / 'case' 只扫单侧。
    · 报告的是**范围**（min–max）而非二值结论：任何单个 p 都依赖该供者是否在场。

    sides: 'control' | 'case' | 'both'
    """
    donors = sorted(set(donor))
    ctrl_donors = [d for d in donors if is_ctrl[donor == d][0]]
    case_donors = [d for d in donors if not is_ctrl[donor == d][0]]
    targets = []
    if sides in ("control", "both"):
        targets += [("control", d) for d in ctrl_donors]
    if sides in ("case", "both"):
        targets += [("case", d) for d in case_donors]

    rows = []
    for side, d in targets:
        keep = donor != d
        r = perm_test(v[keep], donor[keep], is_ctrl[keep],
                      n_perm=n_perm, seed=seed,
                      exact_max_units=exact_max_units,
                      exact_max_combinations=exact_max_combinations)
        # P2-c：jackknife 与主检验必须用同一 exact/random 门限，否则两者的
        # perm_mode 可能不一致（用户下调组合数阈值时尤其明显）
        rows.append({"side": side, "holdout_donor": d,
                     "n_pos_in_holdout": int(v[donor == d].sum()),
                     "cells_in_holdout": int((donor == d).sum()),
                     "p": r["p"], "effect": r["observed_effect"], "mode": r["mode"]})
    df = pd.DataFrame(rows)
    if not len(df):
        return {"table": df}
    ps = df["p"].dropna()
    es = df["effect"].dropna()
    return {
        "table": df,
        "n_ctrl_donors": len(ctrl_donors),
        "n_case_donors": len(case_donors),
        "p_min": round(float(ps.min()), 4) if len(ps) else np.nan,
        "p_max": round(float(ps.max()), 4) if len(ps) else np.nan,
        "effect_min": round(float(es.min()), 4) if len(es) else np.nan,
        "effect_max": round(float(es.max()), 4) if len(es) else np.nan,
        "n_still_significant": int((ps < 0.05).sum()) if len(ps) else np.nan,
        "n_holdouts": len(ps),
        "sides": sides,
    }


# ══════════════════ 结果容器 ══════════════════

@dataclass
class DetectionFloorResult:
    table: pd.DataFrame
    summary: Dict
    n_cells: int
    n_control: int
    n_case: int
    n_donors: int
    jackknife: Optional[Dict] = None
    config: Optional[AuditConfig] = None

    @property
    def unquantifiable(self) -> List[str]:
        """判定为 UNQUANTIFIABLE 的基因（与 summary['n_unquantifiable'] 同定义，S7 修复）。"""
        return self.table.loc[
            self.table["verdict"].astype(str).str.startswith("UNQUANTIFIABLE"),
            "gene"].tolist()

    @property
    def event_sparse(self) -> List[str]:
        return self.table.loc[self.table["stratum"] == "event_sparse", "gene"].tolist()

    def to_text(self) -> str:
        L = ["== Audit 1 · Detection floor =="]
        L.append(f"cells={self.n_cells} (control={self.n_control}, case={self.n_case})"
                 f"  donors={self.n_donors}")
        s = self.summary
        L.append(f"basis={s.get('ratio_basis')}  reference={'yes' if s.get('has_reference') else 'NO'}"
                 f"  permutation={s.get('perm_mode')}  fdr={s.get('fdr_method')}")
        L.append("")
        L.append(self.table.to_string(index=False))
        L.append("")
        L.append(f"UNQUANTIFIABLE: {s['n_unquantifiable']}/{s['n_genes']}   "
                 f"event-sparse: {s['n_event_sparse']}/{s['n_genes']}   "
                 f"(thresholds {s['event_sparse_threshold']})")
        if s.get("event_threshold_sensitivity"):
            L.append("threshold sensitivity (n event-sparse): " +
                     ", ".join(f"<{k}: {v}" for k, v in s["event_threshold_sensitivity"].items()))
        return "\n".join(L)


# ══════════════════ 主函数 ══════════════════

def _verdict_from(basis_pct, R, cfg, has_ref: bool) -> str:
    """检出地板判定链（S6：统一分子；无参照不默认 QUANTIFIABLE）。"""
    if np.isnan(basis_pct):
        return ("INDETERMINATE (no control cells)" if cfg.ratio_basis == "control"
                else "INDETERMINATE (no case cells)")
    if basis_pct < cfg.floor_pct:
        return "UNQUANTIFIABLE (detection floor)"
    if (not np.isnan(R)) and R < cfg.ratio_min:
        return "DROPOUT-LIMITED"
    if np.isnan(R):
        return "INSUFFICIENT (no reference)"
    if basis_pct >= cfg.pct_ok and R >= cfg.ratio_ok:
        return "QUANTIFIABLE"
    return "INSUFFICIENT"


def _depth_matched_pct(ad, gene: str, p: float, rng) -> float:
    """A1：把某基因的原始计数按 p 做二项稀释，返回稀释后的检出率（%）。

    仅当计数为非负整数时有效；否则返回 NaN（不编造）。
    """
    try:
        if gene not in set(map(str, ad.var_names)):
            return np.nan
        col = ad[:, [str(gene)]].X
        import scipy.sparse as _sp
        arr = np.asarray(col.todense()).ravel() if _sp.issparse(col) else np.asarray(col).ravel()
        if not np.allclose(arr, np.round(arr)):
            return np.nan
        thin = rng.binomial(np.round(arr).astype(np.int64), float(p))
        return float((thin > 0).mean() * 100.0)
    except Exception:  # noqa: BLE001
        return np.nan


def audit_detection_floor(ad, genes: Sequence[str], condition_col: str,
                          donor_col: str, control_label: str,
                          cfg: Optional[AuditConfig] = None,
                          reference: Optional[Dict[str, float]] = None,
                          run_jackknife: Optional[bool] = None,
                          modality: Optional[str] = None,
                          reference_modalities: Optional[Sequence[str]] = None,
                          reference_provenance: Optional[Dict[str, str]] = None,
                          allow_cross_modality: bool = False,
                          reference_depth: Optional[float] = None,
                          case_depth: Optional[float] = None,
                          allow_depth_mismatch: bool = False,
                          depth_matched_sensitivity: bool = False,
                          ) -> DetectionFloorResult:
    """modality / reference_modalities / reference_provenance（参照溯源）：

    三者用于在输出表中显式标注 **R 的分母来自哪一列、哪个模态**。
    原实现只给 pct_ref 数值 —— 使用者看到 R 偏低却无从知道分母是核参照列，
    跨模态比较被静默掩盖。分别在表中加两列：
      pct_ref_source      该基因的参照值取自哪一列（mode="max" 时即 argmax 列）
      ref_modality_flag   same-modality / cross-modal / unknown

    v0.5.0 硬约束（REFERENCE_SCHEMA R2）：
      · 目标 modality 与参照 modality **都已声明且不同** → 默认**报错**
        （`CrossModalityError`，CLI 转退出码 2）；
      · 仅当 `allow_cross_modality=True` 时放行，并**逐基因**标 `cross_modality=True`；
      · 任一侧未声明 → `unknown`（无法判定，不报错，也不静默）。
    """
    # 重复基因会产出重复行，任何下游 mean/len 都会双计该基因 —— 且**静默**。
    # 与 io 的重复 var_names 守卫同一纪律：宁可报错，不可静默给错数。
    _dup = [g for g in set(genes) if list(genes).count(g) > 1]
    if _dup:
        raise ValueError(
            f"genes contains duplicates {sorted(_dup)}: this would produce duplicate "
            f"rows and silently double-count them in any downstream mean or count")

    """执行检出地板审计。

    reference : {gene: pct_ref}，外部参照检出率（用于 R）。可选。
    run_jackknife : None=自动（事件数达标者做），True/False 强制。
    """
    from .io import warn_if_gene_names_unresolved as _warn_genes
    _warn_genes(ad, genes, context="audit_detection_floor")
    cfg = cfg or AuditConfig()
    o = ad.obs
    cond = o[condition_col].astype(str).values
    donor = o[donor_col].astype(str).values
    is_ctrl = cond == str(control_label)
    # P2-3：Python API 直接调用时不经 ensure_obs_fields，必须在此显式校验供者 ⊂ 组
    from .io import assert_donors_nested
    assert_donors_nested(donor, cond, donor_col, condition_col)
    n_case = int((~is_ctrl).sum())
    has_ref = bool(reference)

    # 空面板（genes=[]）：rows 为空 → df 无列 → 下游 df["gene"] 抛 KeyError
    # （CLI `--target-genes ""` 曾因此崩）。返回**带列名**的空表，与
    # identifiability 层的空面板处理保持一致（公开就绪度复核 2026-10-07）。
    if not list(genes):
        _cols = ["gene", "n_pos", "n_pos_case", "pct_study", "ctrl_pct", "case_pct",
                 "basis_pct", "pct_ref", "ratio_R", "effect_size_pp", "binom_P_case0",
                 "case_pct_upper95", "cp_rel_halfwidth_case", "overall_pct_upper95",
                 "wilson_lo", "wilson_hi", "stratum", "perm_P", "perm_mode",
                 "q_value_BH", "verdict"]
        return DetectionFloorResult(
            table=pd.DataFrame(columns=_cols),
            summary={"n_genes": 0, "has_reference": has_ref, "empty_panel": True,
                     "note": ("no target gene was supplied: the floor table is empty. "
                              "Pass --target-genes with at least one gene symbol.")},
            n_cells=int(ad.n_obs), n_control=int(is_ctrl.sum()), n_case=n_case,
            n_donors=int(len(set(donor))), config=cfg)

    # 统一分子口径（S6 + N4）：verdict 与 R 共用 cfg.ratio_basis 这一个字段
    basis_name = cfg.ratio_basis
    rows = []
    for g in genes:
        v = gene_detection_vector(ad, g)
        if v is None:
            rows.append({"gene": g, "verdict": "NOT_IN_DATASET",
                         "stratum": "n/a", "note": "gene not found in var_names"})
            continue
        vs = v.astype(float)
        pct_all = vs.mean() * 100
        n_pos = int(v.sum())
        pct_ctrl = vs[is_ctrl].mean() * 100 if is_ctrl.sum() else np.nan
        pct_case = vs[~is_ctrl].mean() * 100 if n_case else np.nan
        n_pos_case = int(v[~is_ctrl].sum())

        # 统一分子（S6）
        basis_pct = {"case": pct_case, "overall": pct_all, "control": pct_ctrl}.get(
            basis_name, pct_case)

        ref = (reference or {}).get(g, np.nan)
        R = (basis_pct / ref) if (ref and not np.isnan(ref) and ref > 0
                                  and not np.isnan(basis_pct)) else np.nan

        # 二项检验（S4 修复：仅在 case_pct == 0 时报告，并标注细胞独立假设）
        binom_p = np.nan
        if n_case and not np.isnan(pct_ctrl) and not np.isnan(pct_case) and pct_case == 0:
            p_ctrl = pct_ctrl / 100
            binom_p = float((1 - p_ctrl) ** n_case)
        cp_upper = _cp_upper(n_pos_case, n_case) if n_case else np.nan
        # 预注册容差口径：病例组检出率的**双侧**相对半宽（供读者逐基因核验阈值推导，见 docs/methods.md §2.2）
        cp_rel_hw = (relative_halfwidth_cp(n_pos_case, n_case)
                     if n_case else np.nan)
        cp_upper_all = _cp_upper(n_pos, len(vs))

        # Wilson CI（事件有限层需要）
        w_lo, w_hi = _wilson_ci(n_pos, len(vs), cfg.strat_ci_level)

        # 事件数分层（事前锁定）
        if n_pos < cfg.event_sparse:
            stratum = "event_sparse"
        elif n_pos < cfg.event_limited:
            stratum = "event_limited"
        else:
            stratum = "event_sufficient"

        # ── 判定（S6：统一分子；无参照不再默认 QUANTIFIABLE）──
        # M5（v0.6.2）：判定链单一实现（`_verdict_from`），避免与敏感性口径漂移。
        verdict = _verdict_from(basis_pct, R, cfg, has_ref)
        rows.append({
            "gene": g, "n_pos": n_pos, "n_pos_case": n_pos_case,
            "pct_study": round(pct_all, 2),
            "ctrl_pct": round(pct_ctrl, 2) if not np.isnan(pct_ctrl) else np.nan,
            "case_pct": round(pct_case, 2) if not np.isnan(pct_case) else np.nan,
            "basis_pct": round(basis_pct, 2) if not np.isnan(basis_pct) else np.nan,
            "pct_ref": round(ref, 2) if (ref and not np.isnan(ref)) else np.nan,
            "ratio_R": round(R, 3) if not np.isnan(R) else np.nan,
            "effect_size_pp": (round(abs(pct_ctrl - pct_case), 2)
                               if not (np.isnan(pct_ctrl) or np.isnan(pct_case)) else np.nan),
            "binom_P_case0": round(binom_p, 3) if not np.isnan(binom_p) else np.nan,
            # P1-d：两列上界语义与命名对齐
            #   case_pct_upper95     ← 病例组检出率的单侧上界
            #   overall_pct_upper95  ← 全体细胞检出率的单侧上界
            #   cp_rel_halfwidth_case ← 病例组检出率的双侧相对半宽（预注册容差口径）
            "case_pct_upper95": round(cp_upper, 2) if not np.isnan(cp_upper) else np.nan,
            "cp_rel_halfwidth_case": (round(cp_rel_hw, 4)
                                      if cp_rel_hw == cp_rel_hw else np.nan),
            "overall_pct_upper95": round(cp_upper_all, 2) if not np.isnan(cp_upper_all) else np.nan,
            "wilson_lo": round(w_lo, 1) if not np.isnan(w_lo) else np.nan,
            "wilson_hi": round(w_hi, 1) if not np.isnan(w_hi) else np.nan,
            # C9（v0.7.2）：显式注明 verdict 的**基准组与阈值** —— stratum 按总阳性细胞数
            # 分档，verdict 按 basis（默认 case）组检出率判定；两者基准不同，并排呈现
            # 会诱发“221 个阳性细胞为何不可定量”的误读。
            "verdict_basis_group": cfg.ratio_basis,
            "verdict_floor_pct": cfg.floor_pct,
            "stratum": stratum,
            "verdict": verdict,
        })

    df = pd.DataFrame(rows)

    # ── 事件稀疏阈值三档敏感性（4.2）──
    thr_sens = {}
    if "n_pos" in df.columns:
        for t in cfg.event_sparse_grid:
            # NOT_IN_DATASET 行不计入事件稀疏数（其 n_pos 为 NaN，不是 0 个阳性细胞）
            _np = df.loc[df["n_pos"].notna(), "n_pos"]
            thr_sens[t] = int((_np < t).sum())

    # ── 置换检验 + jackknife（仅对事件数达标者）──
    perm_ps, perm_modes = [], []
    jk_all = {}
    do_jk = run_jackknife if run_jackknife is not None else True
    for i, row in df.iterrows():
        g = row["gene"]
        if row.get("verdict") == "NOT_IN_DATASET":
            perm_ps.append(np.nan); perm_modes.append("n/a"); continue
        v = gene_detection_vector(ad, g)
        if v is None or int(row["n_pos"]) < cfg.event_sparse:
            perm_ps.append(np.nan); perm_modes.append("skipped (event-sparse: below the pre-specified event "
                              "threshold, so no full-sample P is computed and the gene is "
                              "excluded from the jackknife under the rule)")
            continue
        vs = v.astype(float)
        r = perm_test(vs, donor, is_ctrl, n_perm=cfg.n_perm, seed=cfg.seeds[0],
                      exact_max_units=cfg.exact_enum_max_units,
                      exact_max_combinations=cfg.exact_max_combinations)
        perm_ps.append(r["p"]); perm_modes.append(r["mode"])
        df.loc[i, "donor_weighted_effect"] = round(r["observed_effect"] * 100, 2) if \
            not np.isnan(r["observed_effect"]) else np.nan
        df.loc[i, "perm_units"] = r["n_units"]
        if do_jk:
            jk_all[g] = leave_one_donor(vs, donor, is_ctrl, n_perm=cfg.n_perm,
                                        seed=cfg.seeds[0],
                                        exact_max_units=cfg.exact_enum_max_units,
                                        exact_max_combinations=cfg.exact_max_combinations)
            j = jk_all[g]
            df.loc[i, "jackknife_p_range"] = (f"{j['p_min']:.4f}–{j['p_max']:.4f}"
                                              if not np.isnan(j.get("p_min", np.nan)) else "")
            # 语义注明：本列用**未校正** P<0.05（历史列，保留供对照）
            df.loc[i, "jackknife_n_sig"] = j.get("n_still_significant", np.nan)
            df.loc[i, "jackknife_n_sig_note"] = "uncorrected p<0.05"

    df["perm_P"] = [round(p, 4) if (p == p and not np.isnan(p)) else np.nan for p in perm_ps]
    df["perm_mode"] = perm_modes

    # ── FDR（4.9 + L4）──
    fdr_set, fdr_m = [], 0
    crit_map = {}
    if cfg.fdr_method == "bh" and "perm_P" in df.columns:
        # L4 修复：阳性对照基因（管家基因）不是假设，不参与 m
        exclude = set(str(g).upper() for g in (cfg.fdr_exclude_from_m or []))
        p_for_fdr = [np.nan if str(g).upper() in exclude else p
                     for g, p in zip(df["gene"], df["perm_P"])]
        qs = _bh_fdr(p_for_fdr)
        df["q_value_BH"] = [round(q, 4) if q == q else np.nan for q in qs]
        fdr_set = [g for g, p in zip(df["gene"], p_for_fdr) if p == p]
        fdr_m = len(fdr_set)

        # ── L5：jackknife 的 BH 临界值统计（采纳方案）──
        # 方法学要点：留一拟合是**同一个假设**的重复估计，不应各自做 BH。
        # 正确做法是用**主分析**的 BH 临界值 p_crit = 0.05 · rank_g / m
        # 去数"有多少个留一拟合仍低于该临界值"。
        # 被 fdr_exclude_from_m 排除的基因没有 p_crit → 该列记 NaN。
        # ⚠️ 必须与 q 值用**同一检验族**：q 值由 _bh_fdr(p_for_fdr) 计算，其 m
        # 与 rank 都**排除** fdr_exclude_from_m（管家基因）。原实现用未排除的
        # df["perm_P"] 建 order，m 与 rank 均被管家基因抬高（管家基因通常 P 最小，
        # 会把所有目标基因的 rank 推后），临界值系统性偏松，与 q 值口径不一致。
        # p_for_fdr 已按 df 行序对齐、被排除者为 NaN → 直接用它建 order。
        alpha = float(getattr(cfg, "fdr_alpha", 0.05))
        order = sorted(((i, p) for i, p in enumerate(p_for_fdr) if p == p),
                       key=lambda t: t[1])
        m_tests = len(order)
        for rank, (idx, _p) in enumerate(order, start=1):
            crit_map[df.loc[idx, "gene"]] = alpha * rank / m_tests if m_tests else np.nan
    if jk_all:
        df["jackknife_n_below_bh_crit"] = np.nan
        for g, j in jk_all.items():
            crit = crit_map.get(g, np.nan)
            tbl = j.get("table")
            if crit != crit or tbl is None or not len(tbl) or "p" not in tbl.columns:
                continue
            df.loc[df["gene"] == g, "jackknife_n_below_bh_crit"] = int(
                (tbl["p"] < crit).sum())

    # ── 参照溯源：R 的分母来自哪一列、哪个模态（参照溯源）──
    # 原实现只输出 pct_ref 数值，跨模态比较被静默掩盖：使用者看到 R 偏低，
    # 却无从知道分母是核参照列。此处显式标注，使一份报告即可自洽解读。
    if reference_provenance:
        df["pct_ref_source"] = [reference_provenance.get(str(g), np.nan)
                                for g in df["gene"]]
    else:
        df["pct_ref_source"] = np.nan
    r2_note = None      # C11（v0.7.2）：R2 未强制登记（未声明模态时非空）
    if modality and reference_modalities:
        _mods = {str(m) for m in reference_modalities if str(m).strip()}
        _match = str(modality) in _mods
        if not _match:
            # v0.5.0：跨 modality 默认禁止（R2 硬约束）
            if not allow_cross_modality:
                raise CrossModalityError(
                    f"数据模态 {modality!r} 不在参照模态 {sorted(_mods)} 中："
                    "跨 modality 比较默认禁止（不同 modality 的检出率不可比，R 的绝对值不可解释为生物学量）。"
                    "如确认要跨模态，请用 --allow-cross-modality 显式放行；放行后受影响基因的 R 会标 "
                    "cross_modality=True，且不得用于可定量判定。")
            df["ref_modality_flag"] = "cross-modal"
            df["cross_modality"] = True      # 逐基因标记（仅供放行的跨模态运行）
            warnings.warn(
                f"数据模态 {modality!r} 不在参照模态 {sorted(_mods)} 中："
                "R 为跨模态比较，绝对值不可解释为生物学量（表中已标 cross-modal / cross_modality=True）",
                stacklevel=2)
        else:
            df["ref_modality_flag"] = "same-modality"
    else:
        df["ref_modality_flag"] = "unknown (modality not declared)"
        # C11（v0.7.2）：R2 与 R1 对称——“未强制”时必须显式登记，不得静默。
        r2_note = ("R2 (cross-modality guard) NOT ENFORCED: the modality was not declared "
                   f"on both sides (data modality={modality!r}, reference modalities="
                   f"{list(reference_modalities or [])!r}); per-gene column set to "
                   "'unknown (modality not declared)'. Declare --modality and the "
                   "reference modality (column/manifest) to enforce R2.")
        warnings.warn(r2_note, stacklevel=2)

    # ── A1（v0.6.0）：参照侧深度可比性 ──
    # 仅当参照表提供深度时才加列（保证未提供深度的运行输出不变）。
    depth_note = None
    if reference_depth is not None and reference_depth > 0:
        df["ref_depth"] = round(float(reference_depth), 1)
        if case_depth is not None and case_depth > 0:
            ddp = (float(case_depth) - float(reference_depth)) / float(reference_depth) * 100.0
            df["case_depth"] = round(float(case_depth), 1)
            df["depth_delta_pct"] = round(ddp, 1)
            _lim = float(getattr(cfg, "reference_depth_delta_max_pct", 50.0))
            if abs(ddp) > _lim:
                if not allow_depth_mismatch:
                    raise DepthMismatchError(
                        f"case 深度（中位 total_counts={case_depth:.0f}）与参照深度"
                        f"（median/mean_umi={reference_depth:.0f}）相差 {ddp:+.1f}%，超过事前锁定阈值 "
                        f"{_lim}%。检出率是深度的陡函数，跨深度算 R 会被系统性污染"
                        "（Hicks 2018 / Svensson 2017）。建议换深度可比的参照表，或把 case "
                        "下采样到参照深度（--depth-matched-sensitivity）；若已确认，"
                        "用 --allow-depth-mismatch 显式放行。")
                depth_note = (f"REFERENCE-DEPTH-MISMATCH: case 与参照深度相差 {ddp:+.1f}%"
                              f"（已放行，R 不得用于可定量判定）")
                warnings.warn(depth_note, stacklevel=2)
            if depth_matched_sensitivity:
                _p = min(1.0, float(reference_depth) / float(case_depth))
                rng = np.random.default_rng(int((cfg.seeds or [42])[0]))
                _pcts, _verds = {}, {}
                for g in df["gene"].astype(str):
                    pc = _depth_matched_pct(ad, g, _p, rng)
                    _pcts[g] = pc
                    rg = (reference or {}).get(g, np.nan)
                    R_dm = (pc / rg) if (pc == pc and rg and not np.isnan(rg) and rg > 0) else np.nan
                    _verds[g] = _verdict_from(pc, R_dm, cfg, has_ref)
                df["case_pct_depth_matched"] = [
                    (round(_pcts.get(g), 2) if _pcts.get(g) == _pcts.get(g) else np.nan)
                    for g in df["gene"].astype(str)]
                df["verdict_depth_matched"] = [
                    _verds.get(g) for g in df["gene"].astype(str)]
                # M3（v0.6.2）：计数非整数时稀释不可行 → 全部 NaN，必须显式声明（不得静默）
                if not bool(df["case_pct_depth_matched"].notna().any()):
                    _m3 = ("depth-matched sensitivity SKIPPED: counts are not integers "
                           "(binomial thinning requires integer counts)")
                    depth_note = ((depth_note + " | ") if depth_note else "") + _m3
        else:
            df["case_depth"] = np.nan
            df["depth_delta_pct"] = np.nan
            depth_note = "参照表提供了 mean_umi，但 case 深度不可得：未做深度可比性检查"
    # ── C2（v0.7.1，据第三方独立评估）：放行路径必须**限制判定**，而非只打标记 ──
    # 文档/代码注释承诺"放行后 R 不得用于可定量判定"，但旧实现只在表里写
    # `cross_modality=True`，`_verdict_from()` 根本不看它 → 下游仍读到 QUANTIFIABLE。
    # 此处对**受影响的基因**（跨模态 或 深度差超事前阈值）在放行后：把依赖 R 的判定
    # 降级为 INSUFFICIENT，并加 `verdict_contract_ok` 列（仅在有受影响基因时新增，
    # 保持未受影响运行的输出**逐字节不变**）。
    _R_DEPENDENT = {"QUANTIFIABLE", "DROPOUT-LIMITED", "INSUFFICIENT"}
    _cmod_affected = ("cross_modality" in df.columns and bool(
        df["cross_modality"].fillna(False).astype(bool).any()))
    _lim_d = float(getattr(cfg, "reference_depth_delta_max_pct", 50.0))
    _depth_affected = ("depth_delta_pct" in df.columns and bool(
        df["depth_delta_pct"].notna().any()
        and (df["depth_delta_pct"].abs() > _lim_d).any()))
    if len(df) and "verdict" in df.columns and (_cmod_affected or _depth_affected):
        _contract_ok = pd.Series(True, index=df.index)
        if _cmod_affected:
            _m = df["cross_modality"].fillna(False).astype(bool)
            _contract_ok &= ~_m
            df.loc[_m & df["verdict"].isin(_R_DEPENDENT), "verdict"] = \
                "INSUFFICIENT (cross-modal reference)"
        if _depth_affected:
            _m = df["depth_delta_pct"].notna() & (df["depth_delta_pct"].abs() > _lim_d)
            _contract_ok &= ~_m
            df.loc[_m & df["verdict"].isin(_R_DEPENDENT), "verdict"] = \
                "INSUFFICIENT (depth mismatch)"
        df["verdict_contract_ok"] = _contract_ok
    # 事件稀疏基因一个也没有时，下列列不会被创建，而报告层/白名单会引用它们
    # → KeyError。此处统一补齐为 NaN。
    # 2026-10-06（跨数据集实跑发现）：还有一条退化路径会丢列 —— **目标基因全部不在
    # var_names**（例：数据集只有 Ensembl ID 且未做基因名回退）时，rows 只含
    # {gene, verdict, stratum, note}，于是 table 缺少 n_pos / ratio_R 等核心列。
    # 任何下游消费者（`df["n_pos"]`）都会 KeyError。此处把**核心列**也一并补齐，
    # 保证返回表的**模式稳定**（值可为 NaN，但列不会时有时无）。
    for _col in ("n_pos", "n_pos_case", "pct_study", "ctrl_pct", "case_pct",
                 "basis_pct", "pct_ref", "ratio_R", "effect_size_pp", "binom_P_case0",
                 "case_pct_upper95", "cp_rel_halfwidth_case", "overall_pct_upper95",
                 "wilson_lo", "wilson_hi",
                 "donor_weighted_effect", "perm_units", "jackknife_p_range",
                 "jackknife_n_sig", "jackknife_n_sig_note",
                 "jackknife_n_below_bh_crit"):
        if _col not in df.columns:
            df[_col] = np.nan

    summary = {
        # 置换 P 的性质显式标注：使用者用它复现原作者分析时，该 P 不是本研究
        # 判据下的推断统计量。不标注会被误读为"本研究的检验结果"。
        # jackknife 的 P 波动反映样本量限制，不是结果不稳定；初学者常混淆二者
        # （2026-10-02 建议显式说明）。
        "jackknife_note": ("the leave-one-donor-out P range reflects the number "
                           "of positive cells, not instability of the estimate: "
                           "with few positive cells the range is wide by construction"),
        "perm_p_note": ("permutation P reproduces the authors' own test; it is not "
                        "an inferential statistic under this study's pre-specified "
                        "rule (event-sparse genes: no full-sample P and no jackknife "
                        "under the rule)"),
        "n_genes": int(len(df)),
        "n_unquantifiable": int(df["verdict"].astype(str)
                                .str.startswith("UNQUANTIFIABLE").sum()),
        "n_event_sparse": int((df["stratum"] == "event_sparse").sum()),
        "event_sparse_threshold": cfg.event_sparse,
        "event_threshold_sensitivity": thr_sens,
        "floor_pct_threshold": cfg.floor_pct,
        "ratio_min_threshold": cfg.ratio_min,
        "ratio_basis": cfg.ratio_basis,
        "basis_note": "R 与判定共用 ratio_basis（N4 修复：原 verdict_basis 为死参数）",
        "fdr_gene_set": fdr_set,
        "fdr_n_tests": fdr_m,
        "fdr_excluded_from_m": list(cfg.fdr_exclude_from_m or []),
        "has_reference": has_ref,
        "n_perm": cfg.n_perm,
        "perm_mode": ("exact enumeration" if any(m == "exact" for m in perm_modes)
                      else "random permutation"),
        "fdr_method": cfg.fdr_method,
        "jackknife_sides": "both (control + case)" if do_jk else "disabled",
        # C4（v0.6.0）：binom_P_case0 的口径显式标注（细胞级 = 伪重复，非供体级推断）
        "binom_p_note": ("binom_P_case0 is a CELL-LEVEL statistic (cells treated as "
                          "independent); cells within a donor are pseudoreplicated "
                          "(Squair 2021, PMID 34584091), so it is not a donor-level "
                          "inferential statistic. Donor-level evidence: the jackknife "
                          "range and perm_units."),
        # C5（v0.6.0）：小面板 FDR 功效提示
        "fdr_small_panel_note": ("BH-FDR on a small target panel (m={}) has limited "
                                 "power; q values are reported for completeness, not "
                                 "as a discovery claim.".format(fdr_m)),
        # A1（v0.6.0）：深度可比性（仅在提供了深度时非空）
        "reference_depth_note": depth_note,
        # C11（v0.7.2）：R2 未强制登记（未声明模态时非空）
        "r2_enforcement_note": r2_note,
    }
    return DetectionFloorResult(table=df, summary=summary, n_cells=int(ad.n_obs),
                                n_control=int(is_ctrl.sum()), n_case=n_case,
                                n_donors=int(len(set(donor))),
                                jackknife=jk_all or None, config=cfg)
