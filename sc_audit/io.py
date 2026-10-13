"""数据读取与基本校验。

设计目标：让工具能读【别人的】数据，而不是只认本项目的数据集。

修订记录（2026-10-01）：
  · B1  用户显式指定的 control_label 不再被自动推断覆盖；指定值必须存在
  · B3  多水平 condition 时给出警告，并明确 case_label = "<all non-control>"
  · B4  供者跨组（违反"供者 ⊂ 组"）在入口处报错，而非在置换检验中静默取首细胞
  · B5  校验 X 是否为原始计数；否者优先切到 counts layer，再否则报错
  · P2  基因名回退：var_names 非基因符号（如 Ensembl ID）时尝试 gene_symbols/feature_name
"""
from __future__ import annotations

import os
import warnings
from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp


def x_matrix(ad):
    """返回 ad.X 的 CSR 形态；**backed 时先 to_memory**。

    供**需要全矩阵**的调用使用（门控、QC 指标）——这些无法像
    `detection_vectors` 那样按基因列切片。原实现在 backed 下对 h5py dataset
    调 sp.csr_matrix() 会抛 `ValueError: scipy.sparse does not support dtype
    object`，使 API 用户在 backed 打开数据后调门控/QC 直接崩溃（#4 的边界）。

    ⚠️ 该函数会把矩阵 materialize 到内存；对超大文件请改用 detection_vectors。
    ⚠️ backed 时 `to_memory()` 返回**新对象**：返回的矩阵与原 `ad` **不共享内存**。
    调用方若需写回，必须显式赋值（`ad = ad.to_memory()`），否则改动会落在副本上。
    """
    if getattr(ad, "isbacked", False):
        ad = ad.to_memory()
    return _to_csr(ad.X)


def _to_csr(X):
    return X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)


def _looks_like_counts(X) -> bool:
    """判断矩阵是否像原始计数：非负且近似整数。"""
    try:
        data = X.data if sp.issparse(X) else np.asarray(X).ravel()
    except Exception:
        return False
    data = np.asarray(data)
    if data.size == 0:
        return True
    if np.any(data < 0):
        return False
    return bool(np.allclose(data, np.round(data), atol=1e-6))


def _any_noninteger(X) -> bool:
    return not _looks_like_counts(X)


def load_data(path: str, backed: Optional[str] = None):
    """读取单细胞数据，自动识别常见格式。

    支持：
      · .h5ad                     → AnnData
      · 目录（含 matrix.mtx(.gz)）→ 10x CellRanger 输出
      · *_filtered_feature_bc_matrix.h5 / *.h5 → 10x HDF5
    """
    import scanpy as sc
    sc.settings.verbosity = 0

    if not os.path.exists(path):
        raise FileNotFoundError(f"输入不存在: {path}")

    if path.endswith(".h5ad"):
        ad = sc.read_h5ad(path) if backed is None else sc.read_h5ad(path, backed=backed)
    elif os.path.isdir(path):
        if not any(os.path.exists(os.path.join(path, f))
                   for f in ("matrix.mtx", "matrix.mtx.gz")):
            raise ValueError(f"目录中未找到 10x 矩阵文件: {path}")
        ad = sc.read_10x_mtx(path, var_names="gene_symbols", make_unique=True)
    elif path.endswith(".h5"):
        ad = sc.read_10x_h5(path)
        ad.var_names_make_unique()
    else:
        raise ValueError(
            f"无法识别的格式: {path}（支持 .h5ad / 10x 目录 / 10x .h5）")

    ad.var_names_make_unique()
    return ad


def resolve_gene_names(ad, target_genes: Sequence[str]) -> dict:
    """**只检测**基因名回退可能性，不修改 ad（纯查询）。

    返回 {'resolved': {requested: actual_var_name}, 'missing': [...],
          'fallback_column': str|None}。
    真正执行替换请用 apply_gene_name_fallback()（CLI 与 Python API 共用同一条路径）。
    """
    genes = set(map(str, ad.var_names))
    requested = [str(g) for g in target_genes]
    missing = [g for g in requested if g not in genes]
    if not missing:
        return {"resolved": {g: g for g in requested}, "missing": []}

    for col in ("gene_symbols", "feature_name", "gene_name", "symbol", "gene"):
        if col in getattr(ad.var, "columns", []):
            syms = ad.var[col].astype(str).values
            mapping = {}
            for g in requested:
                hit = np.where(syms == g)[0]
                if len(hit):
                    mapping[g] = str(ad.var_names[hit[0]])
            if mapping:
                still = [g for g in requested if g not in mapping]
                return {"resolved": mapping, "missing": still,
                        "fallback_column": col}
    # 无可用的符号列：resolved 退回"逐字命中"的恒等映射，missing 只报**真正不在
    # var_names 里**的基因。原实现返回 `missing=requested`（全量），会把确实存在的
    # 基因也报成"不在数据集"（2026-10-07 发现：demo 实含 9/14 个默认目标基因，
    # 预检却打印"缺失 14 个"）。
    return {"resolved": {g: g for g in requested if g in genes},
            "missing": [g for g in requested if g not in genes],
            "fallback_column": None}


def warn_if_gene_names_unresolved(ad, target_genes: Sequence[str],
                                  context: str = "audit") -> None:
    """目标基因在 var_names 中**全部未命中**、但 var 中存在符号列时，显式提示。

    为什么需要（2026-10-07「已知答案验证」双方独立发现）：
      CLI 路径会自动做基因名回退（`cli.py` 调 `apply_gene_name_fallback`），
      但**纯 Python API 不会**；而 `GETTING_STARTED.md` 曾写"工具会自动回退（无需手工处理）"，
      于是 API 使用者会直接拿到一整列 `NOT_IN_DATASET`，且不知道原因。
    本函数**不改变任何行为**，只把该状态说出来并给出修复动作（不改 var_names）。
    """
    requested = [str(g) for g in target_genes]
    if not requested:
        return
    names = set(map(str, ad.var_names))
    if any(g in names for g in requested):
        return
    info = resolve_gene_names(ad, requested)
    col = info.get("fallback_column")
    if col:
        import warnings as _w
        _w.warn(
            f"[{context}] 目标基因在 var_names 中一个都未命中，但 var[{col!r}] 提供基因符号。"
            f"CLI 会自动回退；**Python API 需显式调用** "
            f"`sc_audit.apply_gene_name_fallback(ad, genes)` 后再跑（或改用 CLI）。"
            f"本次未回退，结果可能全部为 NOT_IN_DATASET。", stacklevel=2)


def apply_gene_name_fallback(ad, target_genes: Sequence[str]) -> dict:
    """把 var_names 替换为基因符号（若 var 中有符号列），**原地修改 ad**。

    为什么必须是独立函数：
      N1 第一版把替换写死成一句赋值，导致替换逻辑只在 CLI 里存在；
      Python API 直接调用 audit_detection_floor() 时，ENSG 索引的数据集
      仍会全部 NOT_IN_DATASET。现两条路径共用本函数。

    必须在 compute_qc_metrics() **之前**调用：后者靠 var_names 的前缀
    （默认 "MT-"）识别线粒体基因，若此时 var_names 还是 ENSG…，percent_mito 会算成 NaN
    （P1-a 回归）。

    返回 {'applied': bool, 'fallback_column': str|None, 'missing': [...],
          'n_replaced': int}
    """
    info = resolve_gene_names(ad, target_genes)
    col = info.get("fallback_column")
    if not col:
        info["applied"] = False
        info["n_replaced"] = 0
        return info
    ad.var["_orig_var_names"] = ad.var_names.astype(str)
    ad.var_names = ad.var[col].astype(str).values
    ad.var_names_make_unique()
    genes = set(map(str, ad.var_names))
    info["missing"] = [str(g) for g in target_genes if str(g) not in genes]
    info["applied"] = True
    info["n_replaced"] = len(ad.var_names)
    return info


def ensure_obs_fields(ad, condition_col: str, donor_col: str,
                      control_label: Optional[str] = None,
                      case_label: Optional[str] = None) -> dict:
    """校验/推断 obs 中的条件列与供者列。

    B1: 用户显式给出的 control_label 一律不被覆盖；不存在则报错。
    B3: 多水平 condition 时，case 组包含所有非对照取值，并给出警告。
    B4: 供者不得同时出现在多个 condition 组（否则置换/jackknife 的独立单位假设被破坏）。
    """
    for col in (condition_col, donor_col):
        if col not in ad.obs.columns:
            _hint = ""
            if ad.obs.shape[1] == 0:
                _hint = ("\n        提示：该输入的 obs 完全是空的 —— 这是**原始 10x 矩阵**"
                         "（mtx 目录 / raw .h5）的常见情形，它们不携带样本元数据。\n"
                         "        请先用脚本写入 condition / donor（以及 institution / "
                         "donation_type / source / region 等设计变量）并保存为 .h5ad。")
            raise KeyError(
                f"obs 中缺少列 {col!r}；可用列：{list(ad.obs.columns)[:20]}" + _hint)

    vals = ad.obs[condition_col].astype(str)
    uniq = sorted(vals.unique())
    # 空数据（门控/子集把细胞全过滤掉）时给出可诊断的错误，
    # 而不是 "control_label 不在 [] 中" 这种误导信息
    if not uniq:
        raise ValueError(
            f"{condition_col} 在剩余细胞中没有任何取值 —— 数据已被过滤为空。"
            f"若使用了 --preset，请检查 marker 阈值是否与本数据集匹配；"
            f"若使用了 --cell-subset，请检查表达式。")

    # ── ① 对照组标签 ──────────────────────────────────
    if control_label is None:
        ctrl_cands = ["control", "ctrl", "normal", "non_disease", "nondisease",
                      "healthy", "donor", "wt", "nf", "non_failing", "nonfailing"]
        found = next((v for v in uniq if v.lower() in ctrl_cands), None)
        if found is None:
            found = next((v for v in uniq
                          if any(k in v.lower() for k in
                                 ("ctrl", "control", "normal", "healthy", "non_fail"))),
                         None)
        if found is None:
            raise ValueError(
                f"无法自动识别对照组标签（{condition_col} 取值: {uniq}）。"
                f"请用 --control-label 指定。")
        control_label = found
    else:
        control_label = str(control_label)
        if control_label not in uniq:
            raise ValueError(
                f"control_label={control_label!r} 不在 {condition_col} 的取值 {uniq} 中")

    # ── ② 病例组标签 ──────────────────────────────────
    rest = [v for v in uniq if v != control_label]
    if not rest:
        raise ValueError(
            f"除对照组 {control_label!r} 外没有其他取值（{condition_col} 取值: {uniq}）")
    if case_label is None:
        case_label = rest[0] if len(rest) == 1 else "<all non-control>"
    else:
        case_label = str(case_label)
        # P1-2：哨兵 "<all non-control>" 表示“所有非对照水平合并”，
        # 它**不是**数据中的一个真实取值。CLI 在门控/子集后会把这个值原样传回，
        # 若按真实标签校验会误报 ValueError（多水平 + --preset/--cell-subset 崩溃）。
        if case_label == "<all non-control>":
            case_label = rest[0] if len(rest) == 1 else "<all non-control>"
        elif case_label not in uniq:
            raise ValueError(
                f"case_label={case_label!r} 不在 {condition_col} 的取值 {uniq} 中")

    if len(rest) > 1:
        warnings.warn(
            f"{condition_col} 有 {len(rest)} 个非对照取值 {rest}；"
            f"所有审计把它们合并为 case 组（case_label={case_label!r}）",
            UserWarning, stacklevel=2)

    # ── ③ 供者不得跨组（B4）──────────────────────────
    ct = pd.crosstab(ad.obs[donor_col].astype(str), vals)
    cross = ct.index[(ct > 0).sum(axis=1) > 1].tolist()
    if cross:
        raise ValueError(
            f"以下供者同时出现在多个 {condition_col} 组中，违反“供者 ⊂ 组”假设："
            f"{cross[:10]}{' ...' if len(cross) > 10 else ''}。"
            f"请用 --cell-subset 限定或重新编码 {donor_col} 列。")

    mask_ctrl = vals == control_label
    n_ctrl, n_case = int(mask_ctrl.sum()), int((~mask_ctrl).sum())
    if n_ctrl == 0 or n_case == 0:
        raise ValueError(
            f"分组后为空：control={n_ctrl}（label={control_label!r}）, "
            f"case={n_case}（label={case_label!r}）；"
            f"{condition_col} 全部取值 {uniq}")

    return {
        "control_label": control_label,
        "case_label": case_label,
        "case_scope": ("all non-control" if len(rest) > 1 else "single level"),
        "n_condition_levels": len(uniq),
        "condition_levels": uniq,
        "n_control_cells": n_ctrl,
        "n_case_cells": n_case,
        "n_donors": int(pd.Series(ad.obs[donor_col]).nunique()),
    }


def resolve_counts(ad, cfg=None) -> Tuple[object, str]:
    """确保 ad.X 承载**原始计数**（B5）。

    顺序：X（若像计数）→ 指定 layer → 'counts'/'raw_counts' layer → ad.raw → 报错。
    返回 (ad, source_description)。
    """
    layer_pref = []
    if cfg is not None and getattr(cfg, "count_layer", None):
        layer_pref.append(cfg.count_layer)
    layer_pref += ["counts", "raw_counts", "umi_counts"]

    if _looks_like_counts(ad.X):
        return ad, "X"

    for lname in layer_pref:
        if lname and lname in getattr(ad, "layers", {}):
            M = ad.layers[lname]
            if _looks_like_counts(M):
                ad.X = M.copy()
                return ad, f"layers[{lname}]"

    if getattr(ad, "raw", None) is not None and _looks_like_counts(ad.raw.X):
        import anndata as adata
        ad2 = adata.AnnData(X=ad.raw.X.copy(), obs=ad.obs.copy(),
                            var=ad.raw.var.copy())
        return ad2, "raw"

    # 无法得到计数
    data = ad.X.data if sp.issparse(ad.X) else np.asarray(ad.X).ravel()
    raise ValueError(
        "ad.X 不是原始计数（含负值或非整数值），且未找到可用的 counts layer / raw。"
        f"X 的取值示例: {np.asarray(data)[:5].tolist()}…\n"
        "sc-audit 的检出率、门控与深度指标都要求原始计数。"
        "请用 --count-layer 指定存放计数的 layer 名，或先执行 ad.X = ad.layers['counts']。")


def compute_qc_metrics(ad, mito_prefix: str = "MT-", force: bool = False):
    """计算 n_genes / total_counts / percent_mito。

    force=False（默认）：仅在列缺失时计算。
    force=True：强制重算（N6 修复：换源到 counts layer / raw 后，
    上游在 log-normalized X 上算过的 n_genes 会与新 X 不一致，
    导致深度分层使用与检出向量不匹配的深度指标）。
    """
    o = ad.obs
    X = x_matrix(ad)

    if force or "n_genes" not in o.columns:
        o["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel().astype(int)
    if force or "total_counts" not in o.columns:
        o["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
    if force or "percent_mito" not in o.columns:
        genes = list(ad.var_names)
        mt = [i for i, g in enumerate(genes) if str(g).upper().startswith(mito_prefix.upper())]
        if mt and o["total_counts"].sum() > 0:
            mt_counts = np.asarray(X[:, mt].sum(axis=1)).ravel()
            o["percent_mito"] = np.where(o["total_counts"] > 0,
                                         mt_counts / np.maximum(o["total_counts"], 1) * 100,
                                         0.0)
        else:
            o["percent_mito"] = np.nan
    return ad


def detection_vectors(ad, genes) -> dict:
    """一次取出多个基因的 0/1 检出向量（**backed-safe**）。

    这是本包内**唯一**的检出向量读取实现——`gene_detection_vector` 与
    审计四的 `_detect_many` 都经由它，以保证：
      · backed 模式（`.h5ad` 以 `backed='r'` 打开时 `ad.X` 是 h5py
        dataset）下行为一致。直接对 h5py dataset 调 `sp.csr_matrix()` 会抛
        `ValueError: scipy.sparse does not support dtype object`，原实现使
        审计一/三在 backed 数据上**整体失效**且不报错（#4）。
      · 单基因与多基因走**同一条**列切片路径。原实现在单基因时回退到整矩阵
        `ad.X` 并返回第 0 列，会静默给出**错误的基因**，且把 611 MB 矩阵
        整体稠密化（OOM）（#1）。

    基因缺失时**不返回该键**；调用方须处理缺键，不得静默填零。
    """
    # var_names 重名时 get_loc(g) 返回 slice/布尔掩码而非整数位置，
    # `ad[:, idx]` 会失败（或在部分版本上给出错误列）。CLI 路径在 load_data 中
    # 调用了 var_names_make_unique()，但 API 路径没有 —— 故在此显式拒绝。
    if not ad.var_names.is_unique:
        dup = ad.var_names[ad.var_names.duplicated()].unique()[:5].tolist()
        raise ValueError(
            "var_names 不唯一，基因定位不明确（例如 " + ", ".join(map(str, dup)) +
            "）；请先调用 ad.var_names_make_unique() 或提供唯一的基因标识")
    want = [g for g in dict.fromkeys(genes) if g in ad.var_names]
    if not want:
        return {}
    idx = [ad.var_names.get_loc(g) for g in want]
    sub = ad[:, idx]
    if getattr(ad, "isbacked", False):
        sub = sub.to_memory()
    X = sub.X
    M = X.toarray() if sp.issparse(X) else np.asarray(X)
    if M.ndim == 1:
        M = M.reshape(-1, 1)
    return {g: (M[:, i] > 0).ravel() for i, g in enumerate(want)}


def gene_detection_vector(ad, gene: str) -> Optional[np.ndarray]:
    """返回该基因的布尔检出向量；基因不存在时返回 None。

    backed-safe：经 `detection_vectors` 的列切片路径（#4 修复）。
    """
    if gene not in ad.var_names:
        return None
    return detection_vectors(ad, [gene]).get(gene)


def assert_donors_nested(donor, cond, donor_col: str = "donor",
                         condition_col: str = "condition") -> None:
    """B4 的可复用断言：供者必须 ⊂ 组（否则置换/jackknife 的独立单位假设被破坏）。

    P2-3：CLI 走 ensure_obs_fields 时会检查，但 **Python API 直接调用
    audit_detection_floor / audit_control_composition / audit_robust_signal
    会绕过全部入口校验**（B4 / B5 / N1 同时失效）。故三个 audit_* 内部各自调用本函数。
    """
    ct = pd.crosstab(pd.Series(donor, name=donor_col),
                     pd.Series(cond, name=condition_col))
    cross = ct.index[(ct > 0).sum(axis=1) > 1].tolist()
    if cross:
        raise ValueError(
            f"以下供者同时出现在多个 {condition_col} 组中，违反“供者 ⊂ 组”假设："
            f"{cross[:10]}{' ...' if len(cross) > 10 else ''}。"
            f"请先重新编码 {donor_col} 列，或只保留每个供者的一个组。")
