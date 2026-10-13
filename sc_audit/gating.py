"""细胞类型门控（可配置 marker 面板）。

修订记录（2026-10-01）：
  · 修复"两套真值"问题：CARDIAC_PRESETS 与 AuditConfig 原各自维护阈值，
    CLI 用的是前者、yaml 里的 min_counts/pos_markers 从未生效。
    现统一为 preset_from_config()：阈值一律来自 AuditConfig，
    preset 只决定 marker 面板。
"""
from __future__ import annotations

from typing import Dict, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .io import _to_csr, x_matrix


# 预设 → marker 面板（仅面板不同；QC/阈值一律取自 AuditConfig）
PRESET_PANELS: Dict[str, Dict[str, Sequence[str]]] = {
    "cardiac": {},                       # 直接用 cfg.pos_markers / cfg.neg_markers
    "cardiomyocyte_atrial": {"extra_pos": ["NPPA"]},
    "t_cell": {
        "pos": ["CD3D", "CD3E", "CD2", "TRAC", "IL7R"],
        "neg": ["PECAM1", "COL1A1", "DCN", "MYH7", "TTN"],
    },
}

# 兼容旧 API（只读视图；实际阈值请用 preset_from_config）
CARDIAC_PRESETS = {
# API 使用者请勿直接使用本字典：它只是内置默认值的快照，会绕过 config 的
# 覆盖（YAML / CLI 参数）。永远经 preset_from_config(cfg, name) 取门控参数，
# 那是单一真值源。（2026-10-02 指出该绕过路径。）
    "cardiac": {
        "pos_markers": ["TTN", "MYH7", "MYBPC3", "ACTN2", "TNNI3", "MYL2"],
        "neg_markers": ["PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"],
        "min_pos_markers": 2, "min_counts": 500, "min_genes": 200,
    },
    "cardiomyocyte_atrial": {
        "pos_markers": ["TTN", "MYH7", "MYBPC3", "ACTN2", "TNNI3", "MYL2", "NPPA"],
        "neg_markers": ["PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"],
        "min_pos_markers": 2, "min_counts": 500, "min_genes": 200,
    },
    "t_cell": {
        "pos_markers": ["CD3D", "CD3E", "CD2", "TRAC", "IL7R"],
        "neg_markers": ["PECAM1", "COL1A1", "DCN", "MYH7", "TTN"],
        "min_pos_markers": 2, "min_counts": 500, "min_genes": 200,
    },
}


def preset_from_config(cfg, name: str = "cardiac") -> dict:
    """单一真值源：阈值取自 AuditConfig，preset 只决定 marker 面板。"""
    if name not in PRESET_PANELS:
        raise KeyError(f"未知预设 {name!r}；可选：{sorted(PRESET_PANELS)}")
    panel = PRESET_PANELS[name]
    if "pos" in panel:
        pos, neg = list(panel["pos"]), list(panel["neg"])
    else:
        pos = list(cfg.pos_markers) + list(panel.get("extra_pos", []))
        neg = list(cfg.neg_markers)
    return {
        "pos_markers": pos,
        "neg_markers": neg,
        "min_pos_markers": cfg.min_pos_markers,
        "min_counts": cfg.min_counts,
        "min_genes": cfg.min_genes,
    }


def gate_cells(ad, pos_markers: Sequence[str], neg_markers: Sequence[str],
               min_pos_markers: int = 2, min_counts: int = 0,
               min_genes: int = 0, preset: str = None, qc_from: str = "x"):
    """返回布尔 mask（True = 通过门控）。

    规则:
      QC:        n_counts >= min_counts  AND  n_genes >= min_genes
      identity:  阳性标志中 >= min_pos_markers 个表达(count>0)
                 且 阴性标志全部不表达

    缺失的 marker 基因按"未检出"处理，并在 ad.uns 中报告。
    

    ⚠️ QC 阈值（min_counts / min_genes）默认从**传入矩阵**现算（`qc_from="x"`）。
    若你是在**基因子集**（例如只保留面板基因）上做门控，现算深度只覆盖子集 → 阈值必然不过、
    门控会**静默全灭**；此时请传 `qc_from="obs"`（复用 obs 的 total_counts/n_genes，
    通常由 `compute_qc_metrics` 在**全基因空间**算得），或干脆先别预裁基因。
    本次口径与 `qc_from` 均记录在 `ad.uns["sc_audit_gating"]` 中。
    """
    if preset:
        import warnings as _w
        _w.warn("gate_cells(preset=...) 使用模块内置的硬编码阈值，绕过 AuditConfig。"
                "推荐改用 preset_from_config(cfg, name) + 显式传参，使阈值保持单一真值源。",
                DeprecationWarning, stacklevel=2)
        if preset not in CARDIAC_PRESETS:
            raise KeyError(f"未知预设 {preset!r}；可选：{list(CARDIAC_PRESETS)}")
        p = CARDIAC_PRESETS[preset]
        pos_markers = p["pos_markers"]
        neg_markers = p["neg_markers"]
        min_pos_markers = p.get("min_pos_markers", min_pos_markers)
        min_counts = p.get("min_counts", min_counts)
        min_genes = p.get("min_genes", min_genes)

    X = x_matrix(ad)
    genes = pd.Index(ad.var_names)
    n = ad.n_obs

    def col(g):
        if g not in genes:
            return np.zeros(n, dtype=bool)
        return np.asarray(X[:, genes.get_loc(g)].todense()).ravel() > 0

    mask = np.ones(n, dtype=bool)
    # QC 阈值：**优先复用 obs 中已登记的 QC 列**（2026-10-07「已知答案验证」发现：
    # 原实现一律用**传入矩阵**现算 X.sum / (X>0).sum，调用方若为省内存先把基因裁到面板，
    # 行和只剩 ~7 → min_counts=500 必然不过 → 门控**静默全灭**，
    # 下游只报 INDETERMINATE (no case cells)。复用 obs（通常由 compute_qc_metrics
    # 在**全基因空间**算得）可消除该陷阱。）
    if qc_from not in ("x", "obs"):
        raise ValueError(f"qc_from 必须是 'x' 或 'obs'，收到 {qc_from!r}")
    qc_source = []
    _fell_back = []          # qc_from="obs" 但所需列缺失 → 静默回退到 X（A-5：必须提示）
    _x_tc = _x_ng = None
    if min_counts > 0:
        if qc_from == "obs" and "total_counts" in ad.obs.columns:
            _tc = np.asarray(ad.obs["total_counts"].values, dtype=float)
            qc_source.append("obs.total_counts")
        else:
            if qc_from == "obs":
                _fell_back.append("total_counts")
            _x_tc = np.asarray(X.sum(axis=1)).ravel()
            _tc = _x_tc
            qc_source.append("computed_from_X.total_counts")
        mask &= _tc >= min_counts
    if min_genes > 0:
        if qc_from == "obs" and "n_genes" in ad.obs.columns:
            _ng = np.asarray(ad.obs["n_genes"].values, dtype=float)
            qc_source.append("obs.n_genes")
        else:
            if qc_from == "obs":
                _fell_back.append("n_genes")
            _x_ng = np.asarray((X > 0).sum(axis=1)).ravel()
            _ng = _x_ng
            qc_source.append("computed_from_X.n_genes")
        mask &= _ng >= min_genes
    if _fell_back:
        import warnings as _w
        _w.warn(
            f"gate_cells(qc_from='obs')：obs 中缺少 QC 列 {_fell_back}，该指标已**回退**为"
            "从传入矩阵现算 —— 若传的是**基因子集**，该指标仍会失真并可能导致门控全灭。"
            "请先用 compute_qc_metrics(ad) 在**全基因空间**补齐 obs 的 QC 列。", stacklevel=2)

    pos_hits = np.sum([col(g) for g in pos_markers], axis=0) if pos_markers else np.zeros(n, int)
    neg_hits = np.sum([col(g) for g in neg_markers], axis=0) if neg_markers else np.zeros(n, int)
    mask &= (pos_hits >= min_pos_markers) & (neg_hits == 0)

    # 诊断（2026-10-07 已知答案验证发现）：在**基因子集**上做门控时，X 现算的深度只覆盖
    # 该子集（行和可能只有个位数）→ 阈值必然不过、门控**静默全灭**。此处给出可执行的提示。
    if int(mask.sum()) == 0 and (
            any(x.startswith("computed_from_X") for x in qc_source) and (
                qc_from == "x" or _fell_back)) and (
            "total_counts" in ad.obs.columns or "n_genes" in ad.obs.columns):
        import warnings as _w
        _w.warn(
            "gate_cells: 0 个细胞通过门控，而 obs 中**部分或全部**已登记 QC 列（total_counts/n_genes）。"
            "若你是在**基因子集**（如只保留面板基因）上做门控，X 现算的深度只覆盖该子集，"
            "阈值必然不过 —— 请改用 qc_from='obs'（复用全基因空间的 QC），"
            "或先在全基因空间调用 compute_qc_metrics(ad)。", stacklevel=2)

    ad.uns["sc_audit_gating"] = {
        "pos_markers": list(pos_markers),
        "neg_markers": list(neg_markers),
        "min_pos_markers": int(min_pos_markers),
        "min_counts": int(min_counts),
        "min_genes": int(min_genes),
        "qc_from": qc_from,             # "x"(默认，从传入矩阵现算) | "obs"(复用 obs QC)
        "qc_source": qc_source,          # 诊断"门控全灭"是 QC 口径还是 marker 所致
        "missing_pos_markers": [g for g in pos_markers if g not in genes],
        "missing_neg_markers": [g for g in neg_markers if g not in genes],
        "n_pass": int(mask.sum()),
        "n_total": int(n),
        "pct_pass": round(float(mask.mean() * 100), 2),
    }
    return mask

