"""B1（v0.6.0）：数据驱动的潜在批次 / 未登记混杂检测（只报警，不替代审计二）。

**动机**：审计二完全依赖**声明的元数据**；而批次效应常常**不在元数据里**
（Luecken 2022 PMID 34949812；Squair 2021 PMID 34584091；Korsunsky 2019 PMID 31740819）。
本层用**供者级 pseudobulk** 的表达结构检测"供者技术簇"，并检查它与 condition 的共线程度。

**边界（必须随结果引用）**：
  · 本层是本工具里**证据最弱**的一层，是**报警**不是判定；默认关闭（`--with-latent-batch`）；
  · 技术簇**可能含生物学**（如细胞组成差异），因此**绝不**把它当作可匹配变量；
  · 结果依赖 k-means 初始化（种子写入 metadata），属探索性证据。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.cluster.vq import kmeans2

# H3（v0.6.1）：供者×基因稠密矩阵的上限（与 floor 同量级的病态保护）
MAX_DENSE_ELEMENTS = 50_000_000


@dataclass
class LatentBatchResult:
    verdict: str                      # LATENT-BATCH-RISK / NO-LATENT-BATCH-SIGNAL / INDETERMINATE
    ari: float
    n_donors: int
    n_genes_used: int
    n_hvg: int
    seed: int
    ari_by_seed: List[float] = field(default_factory=list)   # L3：多初始化稳定性
    cluster_labels: Dict[str, int] = field(default_factory=dict)
    condition_by_donor: Dict[str, str] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    def to_text(self) -> str:
        L = ["-- latent-batch scan (B1, v0.6.0; exploratory, alarm-only) --",
             f"    verdict        = {self.verdict}",
             f"    donors         = {self.n_donors}  (genes {self.n_genes_used}, HVG {self.n_hvg})",
             f"    ARI(cluster, condition) = {self.ari:.3f}" if self.ari == self.ari else
             "    ARI(cluster, condition) = n/a",
             f"    k-means seed   = {self.seed}"]
        if self.cluster_labels:
            L.append("    cluster by donor:")
            for d, c in sorted(self.cluster_labels.items()):
                L.append(f"        {d:<24} cluster {c}  ({self.condition_by_donor.get(d, '?')})")
        for n in self.notes:
            L.append(f"  · {n}")
        return "\n".join(L)


def _ari(labels_a: np.ndarray, labels_b: np.ndarray) -> float:
    """调整后兰德指数（不依赖 sklearn）。"""
    a = pd.factorize(labels_a)[0]
    b = pd.factorize(labels_b)[0]
    n = len(a)
    if n == 0:
        return float("nan")
    cont = pd.crosstab(a, b).values.astype(float)
    comb2 = lambda x: x * (x - 1) / 2.0
    sum_ij = comb2(cont).sum()
    sum_i = comb2(cont.sum(axis=1)).sum()
    sum_j = comb2(cont.sum(axis=0)).sum()
    total = comb2(n)
    exp = sum_i * sum_j / total if total else 0.0
    mx = (sum_i + sum_j) / 2.0
    if mx == exp:
        return 1.0 if sum_ij == mx else 0.0
    return float((sum_ij - exp) / (mx - exp))


def audit_latent_batch(ad, condition_col: str, donor_col: str, control_label: str,
                       cfg=None, n_hvg: int = 2000, seed: Optional[int] = None
                       ) -> LatentBatchResult:
    """供者级 pseudobulk → HVG → PCA(2) → k-means(2) → 与 condition 的 ARI。"""
    seed = int(seed if seed is not None else (getattr(cfg, "seeds", [42]) or [42])[0])
    warn = float(getattr(cfg, "latent_batch_ari_warn", 0.5))
    o = ad.obs
    cond = o[condition_col].astype(str)
    is_ctrl = (cond == str(control_label)).values
    donor = o[donor_col].astype(str).values
    donors = sorted(set(donor))
    n_d = len(donors)
    res = LatentBatchResult(verdict="INDETERMINATE", ari=float("nan"), n_donors=n_d,
                            n_genes_used=0, n_hvg=0, seed=seed)

    n_ctrl_d = len({d for d in donors if is_ctrl[donor == d][0]})
    n_case_d = n_d - n_ctrl_d
    if n_ctrl_d < 3 or n_case_d < 3:
        res.notes.append(
            f"too few donors per arm (control {n_ctrl_d}, case {n_case_d}) for an "
            f"exploratory latent-batch scan; it is NOT performed.")
        return res

    # ── H3（v0.6.1）：资源上限 + 计数口径校验 ──
    n_genes_total = int(ad.n_vars)
    if n_d * n_genes_total > MAX_DENSE_ELEMENTS:
        res.notes.append(
            f"latent-batch scan SKIPPED: donor×gene matrix would be {n_d}×{n_genes_total} "
            f"(> {MAX_DENSE_ELEMENTS} elements); this scan densifies it and could exhaust "
            f"memory. Subset genes/donors before enabling --with-latent-batch.")
        res.n_genes_used = n_genes_total
        return res
    _X0 = ad.X
    try:
        _data = _X0.data if sp.issparse(_X0) else np.asarray(_X0).ravel()
        _head = _data[: min(len(_data), 2_000_000)]
        _int_ok = bool(np.allclose(_head, np.round(_head)))
    except Exception:  # noqa: BLE001
        _int_ok = False
    if not _int_ok:
        res.verdict = "INDETERMINATE"
        res.notes.append(
            "X does not look like integer counts (normalised / subset / logged input): "
            "the latent-batch scan assumes raw counts, so the donor profiles and the "
            "ARI would NOT be interpretable as technical structure. Declare raw counts "
            "(a counts layer) before enabling --with-latent-batch.")
        return res

    # 供者级 pseudobulk（sum / cells → log1p）
    X = ad.X if sp.issparse(ad.X) else sp.csr_matrix(ad.X)
    idx = {d: i for i, d in enumerate(donors)}
    rows = np.array([idx[d] for d in donor])
    D = sp.csr_matrix((np.ones(len(rows)), (rows, np.arange(len(rows)))),
                      shape=(n_d, len(rows)))
    sums = np.asarray((D @ X).todense())
    counts = D.sum(axis=1).A1
    with np.errstate(invalid="ignore", divide="ignore"):
        means = sums / np.maximum(counts, 1)[:, None]
    pb = np.log1p(means)                       # n_donors × n_genes
    res.n_genes_used = int(pb.shape[1])
    if res.n_genes_used < 50:
        res.notes.append("too few genes to characterise donor profiles.")
        return res

    # HVG（跨供者方差最大者）
    # HVG：L2（v0.7.0）改用**归一化离散度**（Fano 型，方差/均值）而非原始方差，
    # 避免被少数超高表达基因主导（Hafemeister 2019 PMID 31870423 精神）。
    # 评审 R1：必须先**过滤近零表达基因**——否则 mean≈0 的噪声基因（只在 1–2 个供者出现）
    # 会得到极大离散度而被选为“HVG”，使扫描被噪声驱动。
    _detected = (means > 0).sum(axis=0)
    _min_det = max(2, int(np.ceil(0.2 * n_d)))
    _eligible = _detected >= _min_det
    if int(_eligible.sum()) < max(50, int(0.02 * res.n_genes_used)):
        res.verdict = "INDETERMINATE"
        res.notes.append(
            f"too few genes pass the expression filter (detected in >= {_min_det} donors): "
            f"only {int(_eligible.sum())} gene(s); donor profiles cannot be characterised.")
        return res
    _disp = np.where(_eligible,
                     means.var(axis=0) / (np.abs(means.mean(axis=0)) + 1e-9), -np.inf)
    k = min(int(n_hvg), int(_eligible.sum()))
    top = np.argsort(_disp)[::-1][:k]
    res.n_hvg = k
    M = pb[:, top]
    M = (M - M.mean(axis=0)) / (M.std(axis=0) + 1e-9)
    # SVD → 2 PCs
    U, s, _ = np.linalg.svd(M - M.mean(axis=0), full_matrices=False)
    Z = U[:, :2] * s[:2]
    cond_by_donor = np.array([0 if is_ctrl[donor == d][0] else 1 for d in donors])
    # L3（v0.7.0）：k-means 多初始值稳定性（单次初始化是弱证据）。
    _aris, _best, _best_score = [], None, np.inf
    for _s in (seed, seed + 1, seed + 2, seed + 3, seed + 4):
        try:
            _c, _lab = kmeans2(Z, 2, minit="++", seed=int(_s))
        except Exception:  # noqa: BLE001
            continue
        _score = float(np.sum((Z - _c[_lab]) ** 2))
        _aris.append(float(_ari(np.asarray(_lab), cond_by_donor)))
        if _score < _best_score:
            _best_score, _best = _score, _lab
    if _best is None:
        res.notes.append("k-means failed on all seeds; scan not completed.")
        return res
    lab = _best
    ari = float(np.median(_aris))
    res.ari_by_seed = _aris
    if len(_aris) >= 2 and (max(_aris) - min(_aris)) > 0.2:
        res.notes.append(
            f"k-means initialisation is sensitive: ARI across {len(_aris)} seeds spans "
            f"{min(_aris):.3f}–{max(_aris):.3f}; treat the scan as exploratory only.")
    res.ari = float(ari)
    res.cluster_labels = {d: int(lab[i]) for i, d in enumerate(donors)}
    res.condition_by_donor = {d: ("control" if is_ctrl[donor == d][0] else "case")
                              for d in donors}
    if ari == ari and ari >= warn:
        res.verdict = "LATENT-BATCH-RISK"
        res.notes.append(
            f"the case/control split is strongly collinear with the donor expression "
            f"clusters (ARI={ari:.3f} >= {warn}). This is compatible with a latent, "
            f"unregistered batch/technical structure — verify whether a technical axis "
            f"tracks the condition (Luecken 2022; Squair 2021). ALARM ONLY: do NOT use "
            f"the clusters as matching variables (they may encode biology).")
    else:
        res.verdict = "NO-LATENT-BATCH-SIGNAL"
        res.notes.append(
            "no strong collinearity between condition and donor expression clusters. "
            "This does NOT certify absence of confounding: it only rules out one "
            "specific (expression-cluster) pattern, and the scan depends on k-means "
            "initialisation.")
    return res
