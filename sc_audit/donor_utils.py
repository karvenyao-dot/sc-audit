"""供者级统计的**单一实现**（2026-10-06 审计 A3/P2-12 去重）。

为什么需要这个模块
------------------
审计发现同一个算法被复制到 3–4 处（供者级众数、Shannon 熵、供者级唯一列判据、
供者分组标签），且其中部分副本的**边界行为不一致**——这是"改了 A 忘改 B"的
典型温床。此处集中定义，其余模块从这里导入。

约定
----
· "供者级取值" = 该供者在某列上的**主取值（众数）**；供者内无众数（全 NaN）时退化为首个取值。
· "供者级唯一" = 该列在每个供者内恒定，且水平数 ≥ 供者数（即该列与供者一一对应）。
  该判据**不涉及分组**，因此可用于判定层（防止把 batch/lane 误报为与条件共线）。
· 熵为 log2 熵，0 表示单一取值。
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def donor_mode(s: pd.Series) -> str:
    """供者内主取值（众数）；无众数时取首个取值。"""
    m = s.mode()
    return str(m.iloc[0]) if len(m) else str(s.iloc[0])


def donor_level_values(o: pd.DataFrame, col: str, donor: pd.Series) -> pd.Series:
    """每供者一票：返回以 donor 为索引的供者级取值（str）。"""
    return o[col].astype(str).groupby(donor, observed=True).agg(donor_mode)


def donor_group(is_ctrl: np.ndarray, donor: pd.Series) -> pd.Series:
    """每个供者属于 'case' 还是 'ctrl'。"""
    return (pd.Series(np.where(is_ctrl, "ctrl", "case"), index=donor.index)
            .groupby(donor, observed=True).first())


def entropy(values: pd.Series) -> float:
    """log2 熵；0 = 单一取值。"""
    p = pd.Series(values).astype(str).value_counts(normalize=True).values
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum()) if len(p) else 0.0


def is_donor_unique(s: pd.Series, donor: pd.Series, n_donors: int) -> bool:
    """该列是否与供者一一对应（供者内恒定且水平数 ≥ 供者数）。

    真 → 不是设计变量：每个供者自成一类，其供者级 TVD 必然 = 1，与条件无关。
    """
    # 索引对齐守卫（2026-10-07 发现）：若 donor 与 s 的索引不对齐（例如调用方传了
    # 无索引的 numpy/Series），groupby 会得到**空分组** → `(per_donor <= 1).all()`
    # 平凡为真 → 误判为 donor-unique（并给出事实上错误的原因）。此处显式对齐；
    # 无法对齐时返回 False（宁可漏判，不可误判）。
    try:
        donor = pd.Series(donor)
        if not donor.index.equals(s.index):
            if len(donor) != len(s):
                return False
            donor = pd.Series(donor.to_numpy(), index=s.index)
        per_donor = s.groupby(donor, observed=True).nunique()
    except Exception:  # pragma: no cover - 保守兜底
        return False
    if len(per_donor) == 0:
        return False
    return bool((per_donor <= 1).all()) and int(s.nunique()) >= int(n_donors)


def chance_constancy_risk(p_level: float, n_donors: int) -> float:
    """给定**参照人群**中主取值的占比 p 与供者数 n，求"n 个供者恰好全部同值"的概率。

    P = p^n + (1−p)^n（允许全部等于任一侧）。

    ⚠️ p **必须取自参照人群**（整体供者池或对照组），**不能**取病例组内的观测占比：
    病例组内 design-defining 变量的主占比恒 = 1，代入后 P ≡ 1（恒定退化）。
    若参照占比为 0 或 1（该变量在参照人群里无变异），概率不可评估 → 返回 NaN。
    """
    p = float(p_level)
    if p != p or not (0.0 < p < 1.0):
        return float('nan')          # 参照人群无变异 → 无法评估偶然性
    return float(p ** n_donors + (1.0 - p) ** n_donors)
