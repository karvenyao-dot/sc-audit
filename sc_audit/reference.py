"""跨数据集检出比 R（normalised detection ratio）。

绝对检出率不可跨数据集直接比较（CellRanger 版本、核 vs 细胞、建库批次不同），
因此统一使用 R = pct_study / pct_ref。

本模块同时实现 `docs/REFERENCE_SCHEMA_v1.md` 的两条**硬约束**（v0.5.0）：
  · R1 五元组合并守卫：把多个参照列并成同一参照时，
    `(modality, intron_counting, gate_spec_hash, upstream_qc_hash, cell_source)`
    必须完全一致，否则**拒绝合并**并列出差异字段。
  · R2 跨 modality 默认禁止：需 `--allow-cross-modality` 显式放行；放行后受影响
    基因的 R 标 `cross_modality=True`，且不得用于可定量判定。
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd


# 参照列五元组（R1 合并判据；字段名即 long-schema 的列名）
REFERENCE_FIVE_TUPLE = ("modality", "intron_counting", "gate_spec_hash",
                        "upstream_qc_hash", "cell_source")

_CROSS_HINT = ("跨 modality 比较默认**禁止**：不同 modality 的检出率不可比，"
               "R 的绝对值不可解释为生物学量。如确认要跨模态，请用 "
               "`--allow-cross-modality` 显式放行；放行后受影响基因的 R 会标 "
               "`cross_modality=True`，且**不得**用于 `floor_pct`/`ratio_*` 的可定量判定。")


class CrossModalityError(ValueError):
    """R2：跨 modality 比较被禁止（可用 `--allow-cross-modality` 显式放行）。"""


class ReferenceContractError(ValueError):
    """R1：参照列合并违反五元组一致性。"""


class DepthMismatchError(ValueError):
    """A1（v0.6.0）：case 与参照的测序深度相差过大，R 的分子分母不可比。

    检出率是深度的陡函数（Hicks 2018 PMID 29121214；Svensson 2017 PMID 28263961），
    故跨深度算 R 会被系统性污染。可用 `--allow-depth-mismatch` 显式放行。
    """


@dataclass
class ReferenceContract:
    """参照表契约检查结果（R1 + R2）。"""
    path: Optional[str] = None
    ref_modalities: List[str] = field(default_factory=list)
    target_modality: Optional[str] = None
    same_modality: Optional[bool] = None      # None = 无法判定（缺声明）
    allow_cross_modality: bool = False
    r1_enforced: bool = False                  # 是否真的检查了五元组
    r1_violations: List[Dict[str, Any]] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def load_reference_manifest(path: str):
    """读取参照表旁的 manifest（`.manifest.json`）；返回 (dict|None, 路径|None)。"""
    for cand in (path + ".manifest.json", os.path.splitext(path)[0] + ".manifest.json"):
        if os.path.exists(cand):
            try:
                with open(cand, encoding="utf-8") as fh:
                    return json.load(fh), cand
            except Exception:  # noqa: BLE001
                return None, cand
    return None, None


def load_reference_depth(path: str) -> Optional[float]:
    """读取参照表的深度信息（用于 A1 深度可比性）。

    优先级：manifest 的 `mean_umi` → CSV 的 `mean_umi` / `median_umi` / `depth` 列的中位数。
    无任何深度信息 → None（不报错；调用方按“未声明”处理）。
    """
    man, _ = load_reference_manifest(path)
    if isinstance(man, dict):
        # M2（v0.6.2）：优先**中位数**口径（与 case 侧的 median total_counts 对齐，避免 median/mean 混用）
        for k in ("median_umi", "mean_umi", "depth"):
            v = man.get(k)
            if v is not None:
                try:
                    return float(v)
                except (TypeError, ValueError):
                    pass
    try:
        df = pd.read_csv(path)
    except Exception:  # noqa: BLE001
        return None
    for c in ("median_umi", "mean_umi", "depth"):
        if c in df.columns:
            s = pd.to_numeric(df[c], errors="coerce").dropna()
            if len(s):
                return float(s.median())
    return None


def reference_modality_set(path: str) -> List[str]:
    """参照表**声明的模态**（CSV 的 `modality` 列 ∪ manifest 的 modality / column_contract）。"""
    mods = set(reference_modalities(path))
    man, _ = load_reference_manifest(path)
    if isinstance(man, dict):
        for key in ("modality", "reference_modality"):
            v = man.get(key)
            if v:
                mods.add(str(v))
        cc = man.get("column_contract")
        if isinstance(cc, dict):
            for v in cc.values():
                if isinstance(v, dict) and v.get("modality"):
                    mods.add(str(v["modality"]))
    return sorted(m for m in mods if str(m).strip() and str(m).lower() != "nan")


def validate_reference_contract(path: str, target_modality: Optional[str] = None,
                                allow_cross_modality: bool = False) -> ReferenceContract:
    """检查参照表的 R1（五元组合并）与 R2（跨 modality）。**不抛异常**，返回结果。

    · R1：仅当 manifest 提供 `column_contract`（逐列五元组）时才能强制；
      若未提供，在 notes 里**显式登记**“未强制”（绝不静默）。
    · R2：仅当“目标模态”与“参照模态”**都已声明**且确实不同时，才算跨模态。
    """
    res = ReferenceContract(path=path, target_modality=target_modality,
                            allow_cross_modality=bool(allow_cross_modality))
    res.ref_modalities = reference_modality_set(path)

    # ── R1：五元组一致性（宽表下视为把 column_contract 里的各列并成同一参照）──
    man, _ = load_reference_manifest(path)
    cc = man.get("column_contract") if isinstance(man, dict) else None
    if isinstance(cc, dict) and cc:
        res.r1_enforced = True
        tup = {}
        for col, meta in cc.items():
            if isinstance(meta, dict):
                tup[col] = {k: meta.get(k) for k in REFERENCE_FIVE_TUPLE}
        distinct = {tuple(sorted(v.items())) for v in tup.values()}
        if len(distinct) > 1:
            # 列出“与首列不一致”的字段
            base_col = next(iter(tup))
            for col, v in tup.items():
                diff = {k: {"vs_%s" % base_col: tup[base_col][k], col: v[k]}
                        for k in REFERENCE_FIVE_TUPLE if v[k] != tup[base_col][k]}
                if diff:
                    res.r1_violations.append({"column": col, "differs": diff})
    else:
        res.notes.append(
            "R1 未强制：参照表未提供逐列五元组（manifest.column_contract）；"
            "仅登记，不拒绝合并。")

    # ── R2：跨 modality ──
    if target_modality and res.ref_modalities:
        res.same_modality = str(target_modality) in set(res.ref_modalities)
    else:
        res.same_modality = None
    return res


def enforce_reference_contract(contract: ReferenceContract) -> None:
    """按结果抛错：R1 违规 → ReferenceContractError；跨模态且未放行 → CrossModalityError。

    调用方（CLI）负责把它转成 `[error]` + 退出码 2。
    """
    if contract.r1_violations:
        raise ReferenceContractError(
            "R1 五元组不一致：参照列不可合并到同一参照。\n"
            + "\n".join(
                f"  - {v['column']}: " + "; ".join(
                    f"{k}({d}) " for k, d in v["differs"].items())
                for v in contract.r1_violations))
    if contract.same_modality is False and not contract.allow_cross_modality:
        raise CrossModalityError(
            f"数据模态 {contract.target_modality!r} 不在参照模态 "
            f"{contract.ref_modalities} 中。" + _CROSS_HINT)


def reference_modalities(path: str) -> list:
    """读取参考表的 `modality` 列（若存在）。

    用途：`check` 子命令据此在参照表混合模态时主动告警 —— 无需等 schema
    字段全部填充完成。返回去重后的模态列表；无该列时返回空列表。
    """
    df = pd.read_csv(path, nrows=5000)
    col = next((c for c in df.columns if str(c).lower() == "modality"), None)
    if col is None:
        return []
    return sorted({str(x) for x in df[col].dropna().unique()
                   if str(x).strip() and str(x).lower() != "nan"})


def load_reference_with_provenance(path: str, mode: str = "max"):
    """同 load_reference_table，但**并返回每个基因的参照值来自哪一列**。

    动机（2026-10-02）：`mode="max"` 在混合模态参照表上取最大值时，
    核富集基因的最大值必然来自 snRNA 列 —— 而原实现只返回数值，floor 层
    因此**无法知道 R 的分母是哪一列、哪个模态**，跨模态比较被静默掩盖。

    返回 (values, provenance)：
      values[gene]      = 聚合后的参照检出率
      provenance[gene]  = 贡献该值的列名（mode="max" 时即 argmax 列）
    """
    vals = load_reference_table(path, mode=mode)
    df = pd.read_csv(path)
    if "gene" not in df.columns:
        return vals, {}
    cand = [c for c in df.columns if c == "pct_ref"]
    if not cand:
        cand = [c for c in df.columns
                if c.lower() in ("pct", "detection_rate") or c.endswith("_pct")]
    if not cand:
        cand = [c for c in df.select_dtypes("number").columns if c != "gene"]
    if not cand:
        return vals, {}
    num = df[cand].apply(pd.to_numeric, errors="coerce")
    prov = {}
    for i, g in enumerate(df["gene"].astype(str)):
        row = num.iloc[i]
        if row.notna().any():
            if mode == "max":
                prov[g] = str(row.idxmax())
            else:
                prov[g] = "+".join(sorted({str(c) for c in cand
                                           if row[c] == row[c]}))
    return vals, prov


def load_reference_table(path: str, mode: str = "max") -> Dict[str, float]:
    """读取外部参考表（CSV，需含 gene 列与检出率列）。

    接受的列名：pct_ref / pct / detection_rate / 任一以 _pct 结尾的列
    mode='max'    多列时取每基因最大值（主口径，最保守的参照）
    mode='median' 多列时取中位数（敏感性分析）
    """
    df = pd.read_csv(path)
    if "gene" not in df.columns:
        raise ValueError(f"参考表缺少 'gene' 列：{list(df.columns)}")

    cand_cols = [c for c in df.columns if c == "pct_ref"]
    if not cand_cols:
        cand_cols = [c for c in df.columns
                     if c.lower() in ("pct", "detection_rate") or c.endswith("_pct")]
    if not cand_cols:
        num = df.select_dtypes("number").columns.tolist()
        if "gene" in num:
            num.remove("gene")
        cand_cols = num
    if not cand_cols:
        raise ValueError(f"参考表中未找到检出率列：{list(df.columns)}")

    sub = df[["gene"] + cand_cols].copy()
    for c in cand_cols:
        sub[c] = pd.to_numeric(sub[c], errors="coerce")
    sub = sub.set_index("gene")[cand_cols]
    if mode == "median":
        ref = sub.median(axis=1, skipna=True)
    elif mode == "max":
        ref = sub.max(axis=1, skipna=True)
    else:
        raise ValueError(f"mode 必须是 'max' 或 'median'，收到 {mode!r}")
    return {str(k): float(v) for k, v in ref.dropna().items()}


def compute_detection_ratio(study: Dict[str, float],
                            reference: Dict[str, float]) -> pd.DataFrame:
    """study: {gene: pct}；reference: {gene: pct}。

    注意：聚合方式（max / median）在**构建参照时**就已决定
    （见 load_reference_table / build_reference_from_datasets 的 mode 参数）；
    本函数只做逐基因除法。原先此处还有一个 reference_mode 形参，函数体从不使用
    （P3 死参数，已删除）。
    """
    rows = []
    for g, pct in study.items():
        ref = reference.get(g, np.nan)
        R = (pct / ref) if (ref and not np.isnan(ref) and ref > 0) else np.nan
        rows.append({
            "gene": g,
            "pct_study": round(float(pct), 2),
            "pct_ref": round(float(ref), 2) if not np.isnan(ref) else np.nan,
            "ratio_R": round(float(R), 3) if not np.isnan(R) else np.nan,
        })
    out = pd.DataFrame(rows).sort_values("ratio_R", na_position="last")
    return out


def build_reference_from_datasets(dataset_pct: Dict[str, Dict[str, float]],
                                  reference_mode: str = "max") -> Dict[str, float]:
    """从多个数据集的 {gene: pct} 字典构建参照。

    dataset_pct: {'atlas_A': {'MYL3': 38.2, ...}, 'atlas_B': {...}}
    """
    genes = sorted({g for d in dataset_pct.values() for g in d})
    ref = {}
    for g in genes:
        vals = [d[g] for d in dataset_pct.values() if g in d and not np.isnan(d[g])]
        if not vals:
            continue
        ref[g] = float(np.max(vals) if reference_mode == "max" else np.median(vals))
    return ref
