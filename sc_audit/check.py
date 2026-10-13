"""预检（readiness check）：在跑审计【之前】告诉使用者这份数据能支撑哪些审计。

设计原则（与工具的“信息流而非决策流”约定一致）：
  本命令是【信息流】，不是决策流——
    · 它只**报告**候选设计变量及其可匹配性，**绝不**写入 cfg.match_vars
    · 决策权在使用者：看完输出后由人把选定变量写进 config
    · 因此它不破坏工具的"事前锁定"属性（审计二说服力的来源）
  与"自动探测并采纳设计变量"的区别正在于此：后者是循环论证
  （用同一份数据既选变量又判变量），本命令不做这件事。

输出：
    · 控制台文本（人类可读）
    · check_report.json（机器可读，便于流水线/补充材料）
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .composition import _group_distribution_distance
from .config import AuditConfig, DEFAULT_NONDESIGN_QC_COLS
from .donor_utils import is_donor_unique
from .io import ensure_obs_fields, resolve_counts, resolve_gene_names

# obs 中常见的非设计变量（细胞级 QC 列），自动扫描时排除
_QC_COLS = set(DEFAULT_NONDESIGN_QC_COLS)   # 单一真值源（见 config）


def _is_categorical_like(s: pd.Series) -> bool:
    """字符串 / 类别 / 布尔列 —— 可作候选设计变量的类型。

    ⚠️ pandas 3.0 起字符串列的 dtype 不再是 `object`，而是专用 `str` 扩展 dtype
    （PDEP-14）。原先的 `s.dtype == object` 判定在 pandas 3 下对字符串列返回 False，
    导致 `check` 在内存构造的 AnnData 上**静默**给出空候选列表
    （2026-10-07 在干净环境 pandas 3.0.6 下复现：自带测试 4 项失败）。
    此处改为按"语义类型"判定，兼容 object / str / string / category / bool。
    """
    dt = s.dtype
    try:
        import pandas.api.types as ptypes
        if (ptypes.is_object_dtype(dt) or ptypes.is_string_dtype(dt)
                or ptypes.is_bool_dtype(dt)):
            return True
    except Exception:  # pragma: no cover - pandas API 异常时退回 dtype 字面判定
        pass
    if isinstance(dt, pd.CategoricalDtype) or dt == object or dt == bool:
        return True
    return str(dt).startswith(("category", "string", "str"))


def _entropy_bits(series: pd.Series) -> float:
    """类别的 Shannon 熵（bit）。0 = 常量列，越高信息量越大。"""
    p = series.astype(str).value_counts(normalize=True).values
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum()) if len(p) else 0.0


# ── 变量角色启发式（关键词 + 统计信号双判据）────────────────────────────
# 重要：这是**低置信度的提示**，不是判定。分类错误的风险由使用者承担，
#       因此分类结果只用于【分组展示与排序建议】，绝不用于自动剔除或采纳。
_ROLE_PATTERNS = (
    ("outcome", ("fatal", "death", "died", "deceased", "mortality", "surviv",
                 "outcome", "response", "responder", "relapse", "recurrence",
                 "endpoint", "severity", "complication", "readmission",
                 "ventilat", "icu", "mace", "hospitaliz")),
    ("intervention", ("steroid", "treat", "therap", "drug", "medication",
                      "immunosuppress", "dose", "dosing", "tocilizumab", "ivig",
                      "infliximab", "abatacept", "regimen", "chemo", "vaccin")),
    ("technical", ("batch", "lane", "library", "channel", "pool", "plate", "well",
                   "chemistry", "version", "barcode", "sample_id", "chip",
                   "sequenc", "10x", "run_id")),
    ("design", ("institution", "center", "centre", "site", "hospital", "donation",
                "source", "region", "biopsy", "autopsy", "procure", "collection",
                "tissue", "bank", "country", "city", "specimen")),
    ("confounder", ("age", "sex", "gender", "bmi", "ethnic", "race", "smok",
                    "hba1c", "comorbid", "diabet", "hypertens")),
)


_ROLE_ALIASES = {"treatment": "intervention", "therapy": "intervention",
                 "outcome_var": "outcome", "tech": "technical",
                 "procurement": "design"}


def _name_role(col: str, cfg=None) -> str:
    """判断变量角色。优先级：**显式声明 > 领域词根 > 内置词根 > unknown**。

    显式声明（cfg.role_overrides）优先级最高——启发式认不出的结局变量
    （如 outcome_X、complication_30d）由使用者显式标注，避免被误选进 match_vars。
    """
    if cfg is not None:
        ov = getattr(cfg, "role_overrides", None) or {}
        if col in ov:
            v = str(ov[col]).strip().lower()
            return _ROLE_ALIASES.get(v, v)
    low = col.lower()
    if cfg is not None:
        for k in (getattr(cfg, "extra_outcome_keywords", None) or []):
            if k and str(k).lower() in low:
                return "outcome"
        for k in (getattr(cfg, "extra_treatment_keywords", None) or []):
            if k and str(k).lower() in low:
                return "intervention"
    for r, pats in _ROLE_PATTERNS:
        if any(p in low for p in pats):
            return r
    return "unknown"


def classify_variable(col: str, s: pd.Series, is_ctrl: np.ndarray, cfg=None):
    """启发式判断变量角色。返回 (role, reasons, condition_dependent)。

    ⚠️ 角色只由【列名词根】决定；统计信号**不覆盖角色**，只作为附加提示。

    为什么不让统计信号覆盖角色（第一版的错误）：
      "只在单一组内变化"这个信号**同时**匹配两类完全不同的变量——
        · 结局变量：对照全为 no（不会因没有的病而死），病例有 yes/no
        · 混杂的设计变量：对照跨 2 个中心、病例全在 1 个中心
      二者在数据上无法区分。若据此把 institution 判为"疑似结局"，
      就会**把最该被匹配的变量藏起来**——假阳性比假阴性更糟。
      因此统计信号只产出 condition_dependent 提示，交由使用者甄别。
    """
    role = _name_role(col, cfg)
    reasons = []
    hit = None
    low = col.lower()
    _ov = getattr(cfg, "role_overrides", None) or {} if cfg is not None else {}
    if col in _ov:
        reasons.append(f"由配置显式声明（role_overrides）→ 角色：{_ROLE_ZH.get(role, role)}")
    else:
        _extra = []
        if cfg is not None:
            _extra = list(getattr(cfg, "extra_outcome_keywords", None) or []) + \
                     list(getattr(cfg, "extra_treatment_keywords", None) or [])
        if _extra:
            hit = next((k for k in _extra if k and str(k).lower() in low), None)
        if not hit:
            for _r, pats in _ROLE_PATTERNS:
                hit = next((p for p in pats if p in low), None)
                if hit:
                    break
        if hit:
            reasons.append(f"列名含「{hit}」词根 → 角色：{_ROLE_ZH.get(role, role)}")

    cv = int(pd.Series(s.values[is_ctrl]).nunique(dropna=True))
    kv = int(pd.Series(s.values[~is_ctrl]).nunique(dropna=True))
    cond_dep = (cv <= 1 < kv) or (kv <= 1 < cv)
    if cond_dep:
        reasons.append(
            f"仅在单一组内变化（对照 {cv} 个水平 / 病例 {kv} 个水平）→ "
            f"既可能是结局/条件依赖变量，也可能只是分组不平衡的取材变量，"
            f"数据上无法区分，需人工甄别")
    if not reasons:
        reasons.append("列名与分布均无明确信号 → 角色未知，请自行判断")
    return role, reasons, cond_dep


def recommend_role(role: str, tvd_donor: float, donor_unique: bool,
                   cfg: "AuditConfig", condition_dependent: bool = False):
    """给出**建议**（不写入配置、不参与判定）。返回 (recommendation, reasoning)。"""
    if donor_unique:
        return ("not_usable",
                "与供者一一对应，匹配无意义（每个供者一个取值）")
    if role in ("outcome", "outcome_suspect"):
        return ("never",
                "疑似结局/条件依赖变量——纳入匹配等于条件化在结局上，"
                "会引入 collider bias，使暴露-结局关系无法解释")
    if role == "intervention":
        return ("avoid", "治疗/干预变量，属暴露或治疗后协变量，不是取材差异")
    if role == "ignore":
        return ("avoid", "由配置显式声明为忽略（role_overrides: ignore）")
    if role == "technical":
        return ("avoid", "技术特征（批次/通道等），不代表取材来源差异")
    if tvd_donor > cfg.hetero_max:
        return ("risky",
                f"供者级 TVD={tvd_donor:.3f} > hetero_max={cfg.hetero_max}："
                f"与分组接近分离 → 若纳入并命中，判定将直接指向 NOT SEPARABLE。"
                f"匹配它会在控制来源差异的同时抹掉真实的病例/对照差异"
                f"（把生物学一起匹配走）")
    if condition_dependent:
        return ("caution_condition_dependent",
                f"只在单一组内变化：可能是结局/条件依赖变量（不可匹配），"
                f"也可能只是分组不平衡的取材变量（正是要匹配的对象）。"
                f"工具无法自动区分 → 请人工确认该变量在两组中语义是否相同")
    if tvd_donor < cfg.match_var_tvd_low:
        return ("safe_but_uninformative",
                f"供者级 TVD={tvd_donor:.3f} < {cfg.match_var_tvd_low}："
                f"与分组近乎独立，匹配它几乎不控制任何东西（安全但无信息量）")
    return ("primary_match",
            f"供者级 TVD={tvd_donor:.3f} 落在 "
            f"[{cfg.match_var_tvd_low}, {cfg.hetero_max}] 之间："
            f"既与分组相关、又未接近完全分离——匹配它能控制部分来源差异"
            f"而不至于把生物学一并抹去。代价：匹配会丢弃细胞，TVD 越低收益越小"
            f"（接近 {cfg.match_var_tvd_low} 时已接近'匹配也无所谓'）")


_STATUS_BY_REC = {
    "primary_match": "可用·推荐",
    "safe_but_uninformative": "可用·低信息",
    "caution_condition_dependent": "可用·需甄别",
    "risky": "可用·需注意",
    "not_usable": "剔除",
    "avoid": "不建议",
    "never": "禁用",
}
_ROLE_ZH = {"design": "取材/来源", "confounder": "人口学协变量",
            "technical": "技术特征", "intervention": "治疗/干预",
            "outcome": "结局", "outcome_suspect": "疑似结局/条件依赖",
            "unknown": "未识别"}


@dataclass
class ReadinessResult:
    # ── 输入概况 ──
    n_cells: int = 0
    n_genes: int = 0
    counts_source: str = ""
    counts_ok: bool = False
    gene_lookup: Dict = field(default_factory=dict)
    gating: Optional[Dict] = None

    # ── 分组与供者 ──
    obs_ok: bool = False
    condition_levels: List[str] = field(default_factory=list)
    case_scope: str = ""
    control_label: Optional[str] = None
    case_label: Optional[str] = None
    n_control_cells: int = 0
    n_case_cells: int = 0
    n_donors_ctrl: int = 0
    n_donors_case: int = 0

    # ── 设计变量评估（仅报告，不采纳）──
    design_vars: List[Dict] = field(default_factory=list)
    # 疑似结局/干预变量单独列出（仅提示，仍会完整展示，不静默剔除）
    suspect_vars: List[Dict] = field(default_factory=list)

    # ── 可行性 ──
    feasibility: Dict = field(default_factory=dict)
    audit_2_recommendation: Dict = field(default_factory=dict)
    depth_bins_ok: Dict = field(default_factory=dict)
    reference_coverage: Dict = field(default_factory=dict)
    reference_provenance: Dict = field(default_factory=dict)

    # ── 诊断 ──
    blockers: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    recommendations: List[str] = field(default_factory=list)
    reference_genes: int = 0

    def to_json(self) -> Dict:
        return {
            "schema": "sc-audit/readiness/1",
            "input": {"n_cells": self.n_cells, "n_genes": self.n_genes,
                      "counts_source": self.counts_source,
                      "counts_ok": self.counts_ok,
                      "gene_lookup": self.gene_lookup,
                      "gating": self.gating},
            "groups": {"obs_ok": self.obs_ok,
                       "condition_levels": self.condition_levels,
                       "case_scope": self.case_scope,
                       "control_label": self.control_label,
                       "case_label": self.case_label,
                       "n_control_cells": self.n_control_cells,
                       "n_case_cells": self.n_case_cells,
                       "n_donors_ctrl": self.n_donors_ctrl,
                       "n_donors_case": self.n_donors_case},
            "design_variables": self.design_vars,
            "suspect_design_variables": self.suspect_vars,
            "feasibility": self.feasibility,
            "audit_2_recommendation": self.audit_2_recommendation,
            "depth_bins": self.depth_bins_ok,
            "reference_coverage": self.reference_coverage,
            "reference_provenance": self.reference_provenance,
            "blockers": self.blockers,
            "warnings": self.warnings,
            "recommendations": self.recommendations,
            "reference_genes": self.reference_genes,
            "note": ("This command only REPORTS. It never writes to the config. "
                     "Choose your match variables yourself and write them into the "
                     "config yaml before running composition/run — that is what keeps "
                     "the audit pre-registered."),
        }

    def to_text(self) -> str:
        L: List[str] = []
        L.append("=" * 74)
        L.append("sc-audit readiness check（预检 · 只报告，不采纳）")
        L.append("=" * 74)

        L.append("\n【输入】")
        L.append(f"  细胞 {self.n_cells:,} · 基因 {self.n_genes:,}")
        L.append(f"  计数来源: {self.counts_source or '未检查'}"
                 f"{'' if self.counts_ok else '   ⚠️ 非原始计数'}")
        if self.gating:
            L.append(f"  门控({self.gating.get('preset', 'preset')}): "
                     f"{self.gating.get('n_pass')}/{self.gating.get('n_total')} 通过 "
                     f"({self.gating.get('pct_pass')}%)  "
                     f"缺失 marker: pos={self.gating.get('missing_pos_markers')}, "
                     f"neg={self.gating.get('missing_neg_markers')}")
        gl = self.gene_lookup or {}
        if gl:
            L.append(f"  基因名: {'符号（直接可用）' if not gl.get('applied') else '已回退 → ' + str(gl.get('fallback_column'))}"
                     f"   缺失 {len(gl.get('missing') or [])} 个")

        L.append("\n【分组与供者】")
        if not self.obs_ok:
            L.append("  ❌ 无法解析分组/供者（见下方 blocker）")
        else:
            L.append(f"  条件水平: {self.condition_levels}")
            L.append(f"  case 口径: {self.case_scope}")
            L.append(f"  对照 {self.control_label!r}: {self.n_control_cells:,} 细胞 / "
                     f"{self.n_donors_ctrl} 供者")
            L.append(f"  病例 {self.case_label!r}: {self.n_case_cells:,} 细胞 / "
                     f"{self.n_donors_case} 供者")

        L.append("\n【候选设计变量（按建议优先级排序；本命令不会自动采用）】")
        if not self.design_vars:
            L.append("  未找到可评估的候选列")
        else:
            L.append(f"  {'变量':<18}{'角色':<12}{'水平':>4}{'供者级TVD':>10}"
                     f"{'V':>7}{'熵':>6}   建议")
            for d in self.design_vars:
                cv_s = f"{d['cramers_v']:.3f}" if isinstance(d.get("cramers_v"),
                                                             (int, float)) else "  n/a"
                L.append(f"  {d['column']:<18}{d.get('role_zh', ''):<12}"
                         f"{d['n_levels']:>4}{d['tvd_donor']:>10.3f}"
                         f"{cv_s:>7}{d['entropy_bits']:>6.2f}   {d['status']}")
            L.append("")
            for d in self.design_vars:
                if d.get("reasoning"):
                    L.append(f"  · [{d['status']}] {d['column']}：{d['reasoning']}")

        if self.suspect_vars:
            L.append("\n【疑似结局 / 治疗变量 — 不建议用于匹配（仅提示，未隐藏）】")
            L.append("  这些列的角色由「列名词根」与「是否只在单一组内变化」双判据推测，"
                     "属低置信度提示：")
            for d in self.suspect_vars:
                L.append(f"  ✗ {d['column']}（{d.get('role_zh', '')}）"
                         f"  供者级TVD={d['tvd_donor']:.3f}  → {d['status']}")
                for r in d.get("role_reasons", []):
                    L.append(f"        判据: {r}")
                if d.get("reasoning"):
                    L.append(f"        后果: {d['reasoning']}")

        L.append("\n【三项审计可行性】")
        for k in ("audit_1_floor", "audit_2_composition", "audit_3_robust"):
            f_ = self.feasibility.get(k, {})
            mark = {"ok": "✅", "degraded": "⚠️", "no": "❌"}.get(f_.get("level"), "?")
            L.append(f"  {mark} {k:<20} {f_.get('status', '')}")
            for n in f_.get("notes", []):
                L.append(f"       - {n}")

        a2r = self.audit_2_recommendation or {}
        if a2r:
            L.append("\n【审计二 · 总体建议】")
            L.append(f"  建议做法: {a2r.get('method')}   可行: {a2r.get('feasible')}")
            if a2r.get("recommended_variables"):
                L.append(f"  推荐变量: {a2r['recommended_variables']}")
            if a2r.get("excluded_by_role"):
                L.append(f"  按角色排除（不建议用于匹配）: {a2r['excluded_by_role']}")
            if a2r.get("excluded_as_unusable"):
                L.append(f"  不可用于匹配: {a2r['excluded_as_unusable']}")
            if a2r.get("statement_for_paper"):
                L.append("  可直接用于方法学描述的措辞：")
                L.append(f"    “{a2r['statement_for_paper']}”")
            if a2r.get("options"):
                L.append("  可选策略（工具不排序，只列权衡）：")
                for o in a2r["options"]:
                    L.append(f"    · [{o['id']}] {o['label']}")
                    L.append(f"        配置提示: {o['config_hint']}")
                    L.append(f"        代价/取舍: {o['trade_off']}")
                if a2r.get("options_note"):
                    L.append(f"    {a2r['options_note']}")

        rc = self.reference_coverage or {}
        if rc and rc.get("n_target_genes"):
            L.append("\n【参照表覆盖】")
            L.append(f"  参照表共 {rc['n_reference_genes']} 个基因 · 覆盖目标面板 "
                     f"{len(rc.get('target_covered', []))}/{rc['n_target_genes']}"
                     f"（{rc.get('pct_target_covered')}%）")
            if rc.get("target_missing"):
                L.append(f"  未覆盖 → 将报 INSUFFICIENT (no reference): "
                         f"{rc['target_missing']}")

        rp = self.reference_provenance or {}
        if rp.get("manifest_found"):
            L.append("\n【参照表溯源（门控参数登记）】")
            L.append(f"  背景数: {rp['n_backgrounds']}")
            L.append(f"  pct_ref 规则: {rp.get('pct_ref_rule','')}")
            if rp.get("gate_unknown_backgrounds"):
                L.append(f"  ⚠️ 门控未登记的背景: {rp['gate_unknown_backgrounds']}"
                         f" → R 分母含口径未知分量")
            for g in (rp.get("known_gaps") or []):
                L.append(f"  · {g}")
        elif rp:
            L.append("\n【参照表溯源】⚠️ 未找到 manifest → R 可比性无法核验")

        L.append("\n【深度分层容量】")
        db = self.depth_bins_ok or {}
        if db:
            L.append(f"  分箱数 {db.get('n_bins')} · "
                     f"两组均 >= {db.get('min_cells_per_stratum')} 细胞的箱: "
                     f"{db.get('n_bins_ok')} 个")
            for b in db.get("bins", []):
                flag = "✓" if b["ok"] else " "
                L.append(f"    [{flag}] {b['range']:<18} 对照 {b['n_ctrl']:>6} "
                         f"病例 {b['n_case']:>6}")

        if self.blockers:
            L.append("\n【阻断项 — 必须先解决】")
            for b in self.blockers:
                L.append(f"  ❌ {b}")
        if self.warnings:
            L.append("\n【警告】")
            for w in self.warnings:
                L.append(f"  ⚠️ {w}")
        if self.recommendations:
            L.append("\n【建议】")
            for r in self.recommendations:
                L.append(f"  → {r}")

        L.append("\n" + "-" * 74)
        L.append("提示：本命令不写入任何配置。请自行把选定的设计变量写进 config yaml，")
        L.append("      再运行 composition / run —— 这一步是审计二「事前锁定」属性的来源。")
        L.append("=" * 74)
        return "\n".join(L)


def _scan_candidate_columns(ad, condition_col: str, donor_col: str,
                            explicit: Optional[Sequence[str]]) -> List[str]:
    """确定要评估的候选设计变量列。

    优先使用使用者显式给出的列表；否则自动扫描 obs 中的类别型列（仅用于报告）。
    """
    o = ad.obs
    if explicit:
        return [c for c in explicit if c in o.columns]

    out: List[str] = []
    for c in o.columns:
        if c in (condition_col, donor_col) or c in _QC_COLS:
            continue
        s = o[c]
        if _is_categorical_like(s):
            if s.nunique(dropna=True) < 2:      # 常量列不进候选
                continue
            if s.nunique(dropna=True) >= ad.n_obs:   # 细胞级唯一 → 无信息
                continue
            out.append(c)
    # 让常见的设计变量排在前面，便于阅读
    pref = ["institution", "center", "site", "hospital", "donation_type",
            "source", "region", "sex", "gender", "age_group", "batch", "lane",
            "library", "sample"]
    out.sort(key=lambda c: (pref.index(c) if c in pref else len(pref), c))
    return out


def _depth_bin_capacity(ad, is_ctrl: np.ndarray, cfg: AuditConfig) -> Dict:
    """按 cfg.depth_edges 统计各箱两组的细胞数（审计三的容量前提）。"""
    o = ad.obs
    if "n_genes" not in o.columns:
        # C12（v0.7.1，据第三方独立评估）：原实现直接 `sp.csr_matrix(ad.X)`，
        # 在 backed='r' 打开时 ad.X 是 h5py dataset → 抛
        # `ValueError: scipy.sparse does not support dtype object`。改为走
        # `io.x_matrix()`（backed-safe，与本包其他 need-full-matrix 调用一致）。
        from .io import x_matrix as _x_matrix
        X = _x_matrix(ad)
        ng = np.asarray((X > 0).sum(axis=1)).ravel()
    else:
        ng = o["n_genes"].values
    edges = list(cfg.depth_edges)
    bins = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (ng >= lo) & (ng < hi)
        n_c = int((m & is_ctrl).sum())
        n_k = int((m & ~is_ctrl).sum())
        ok = (n_c >= cfg.min_cells_per_stratum) and (n_k >= cfg.min_cells_per_stratum)
        bins.append({"range": f"({lo},{hi}]" if hi < 100000 else f"({lo},inf)",
                     "n_ctrl": n_c, "n_case": n_k, "ok": bool(ok)})
    return {"n_bins": len(bins),
            "min_cells_per_stratum": cfg.min_cells_per_stratum,
            "n_bins_ok": sum(1 for b in bins if b["ok"]),
            "bins": bins}


def check_readiness(ad, condition_col: str = "condition", donor_col: str = "donor",
                    control_label: Optional[str] = None,
                    cfg: Optional[AuditConfig] = None,
                    reference: Optional[Dict[str, float]] = None,
                    candidate_vars: Optional[Sequence[str]] = None,
                    target_genes: Optional[Sequence[str]] = None,
                    reference_manifest: Optional[Dict] = None,
                    check_counts: bool = True) -> ReadinessResult:
    """生成预检报告（不修改任何配置）。"""
    cfg = cfg or AuditConfig()
    res = ReadinessResult()
    res.n_cells = int(ad.n_obs)
    res.n_genes = int(ad.n_vars)
    res.reference_genes = len(reference or {})

    # ── ① 原始计数 ──
    # check_counts=False 时跳过（对应 CLI 的 --no-count-check，保持与其他子命令一致）
    if check_counts:
        try:
            _, src = resolve_counts(ad, cfg)
            res.counts_source = src
            res.counts_ok = (src in ("X", "X (assumed)")
                             or "layer" in src or "raw" in src)
        except Exception as exc:                   # noqa: BLE001
            res.counts_source = "unavailable"
            res.counts_ok = False
            res.blockers.append(f"原始计数不可用：{exc}")
    else:
        res.counts_source = "skipped (--no-count-check)"
        res.counts_ok = True

    # ── ② 基因名（只检测，不改动）──
    if target_genes:
        try:
            res.gene_lookup = resolve_gene_names(ad, list(target_genes))
        except Exception as exc:                   # noqa: BLE001
            res.gene_lookup = {"missing": list(target_genes), "error": str(exc)}

    # ── ③ 分组与供者 ──
    try:
        info = ensure_obs_fields(ad, condition_col, donor_col, control_label, None)
        res.obs_ok = True
        res.condition_levels = info.get("condition_levels", [])
        res.case_scope = info.get("case_scope", "")
        res.control_label = info.get("control_label")
        res.case_label = info.get("case_label")
        res.n_control_cells = int(info.get("n_control_cells", 0))
        res.n_case_cells = int(info.get("n_case_cells", 0))
    except Exception as exc:                       # noqa: BLE001
        res.blockers.append(f"分组/供者解析失败：{exc}")
        res.feasibility = {
            "audit_1_floor": {"level": "no", "status": "不可运行",
                              "notes": ["需先修复分组/供者列"]},
            "audit_2_composition": {"level": "no", "status": "不可运行", "notes": []},
            "audit_3_robust": {"level": "no", "status": "不可运行", "notes": []},
        }
        return res

    cond = ad.obs[condition_col].astype(str).values
    donor = ad.obs[donor_col].astype(str).values
    is_ctrl = cond == str(res.control_label)
    res.n_donors_ctrl = len({d for d in set(donor) if is_ctrl[donor == d][0]})
    res.n_donors_case = len(set(donor)) - res.n_donors_ctrl

    if res.case_scope == "all non-control":
        res.warnings.append(
            f"condition 有 {len(res.condition_levels)} 个水平；所有非对照将被合并为 case")
    for side, n in (("对照", res.n_donors_ctrl), ("病例", res.n_donors_case)):
        if n < 5:
            res.warnings.append(
                f"{side}组只有 {n} 个供者 → 置换检验功效不足"
                f"（精确枚举下的最小可达双侧 P 为 2/C(n,k)（两组等大）或 1/C(n,k)"
                f"（两组不等大，且要求极值划分唯一）；若该下限 ≥ 0.05，则即使两组"
                f"完全分离也给不出显著结论）")

    # ── ④ 设计变量评估（只报告）──
    cands = _scan_candidate_columns(ad, condition_col, donor_col, candidate_vars)
    # 使用者显式声明 ignore 的列不进入候选
    _ign = {c for c in cands if _name_role(c, cfg) == "ignore"}
    if _ign:
        cands = [c for c in cands if c not in _ign]
    n_donors_total = len(set(donor))
    for c in cands:
        s = ad.obs[c].astype(str)
        # 去重（2026-10-06 审计 P2-12）：单一实现见 donor_utils.is_donor_unique
        donor_unique = is_donor_unique(s, pd.Series(donor, index=s.index), n_donors_total)
        cell_unique = s.nunique() >= max(2, int(0.5 * ad.n_obs))
        tvd_d, cv, _, _ = _group_distribution_distance(s.values, is_ctrl, donor=donor)
        tvd_c, _, _, _ = _group_distribution_distance(s.values, is_ctrl)
        ent = _entropy_bits(s)

        role, role_reasons, cond_dep = classify_variable(c, s, is_ctrl, cfg)
        rec, rec_reason = recommend_role(role, float(tvd_d), donor_unique, cfg,
                                         condition_dependent=cond_dep)
        if cell_unique and not donor_unique:
            rec, rec_reason = "not_usable", "细胞级唯一（每个细胞一个取值），无匹配意义"
        if s.nunique() < 2:
            rec, rec_reason = "not_usable", "常量列"

        row = {
            "column": c, "n_levels": int(s.nunique()),
            "role": role, "role_zh": _ROLE_ZH.get(role, role),
            "role_reasons": role_reasons,
            "condition_dependent": bool(cond_dep),
            "recommendation": rec, "reasoning": rec_reason,
            "tvd_donor": round(float(tvd_d), 3), "tvd_cell": round(float(tvd_c), 3),
            "cramers_v": round(float(cv), 3) if cv == cv else None,
            "entropy_bits": round(ent, 2), "donor_unique": donor_unique,
            "status": _STATUS_BY_REC.get(rec, rec),
            "reason": rec_reason,
        }
        # 按角色分流【展示】，但不隐藏：疑似结局/干预变量单独成节，仍然完整可见
        if role in ("outcome", "outcome_suspect", "intervention"):
            res.suspect_vars.append(row)
        else:
            res.design_vars.append(row)

    usable = [d["column"] for d in res.design_vars
              if d["recommendation"] == "primary_match"]
    safe_only = [d["column"] for d in res.design_vars
                 if d["recommendation"] == "safe_but_uninformative"]
    cautions = [d["column"] for d in res.design_vars
                if d["recommendation"] == "caution_condition_dependent"]
    risky = [d["column"] for d in res.design_vars if d["recommendation"] == "risky"]
    flagged = [d["column"] for d in res.suspect_vars]
    dropped = [d["column"] for d in res.design_vars
               if d["recommendation"] == "not_usable"]

    # 排序：推荐优先；同档内 TVD 高者在前（信息量更大但仍安全）
    _PRIO = {"primary_match": 0, "safe_but_uninformative": 1, "risky": 2,
             "avoid": 3, "not_usable": 4, "never": 5}
    res.design_vars.sort(key=lambda d: (_PRIO.get(d["recommendation"], 9),
                                        -d["tvd_donor"]))
    res.suspect_vars.sort(key=lambda d: -d["tvd_donor"])

    # ── ⑤ 深度分层容量 ──
    res.depth_bins_ok = _depth_bin_capacity(ad, is_ctrl, cfg)

    # ── ⑥ 三项审计可行性 ──
    a1_notes = []
    if not res.counts_ok:
        a1_notes.append("原始计数不可用 → 检出率语义不可靠")
    if reference is None:
        a1_notes.append("未提供 --reference → R 无法计算，多数基因将判 "
                        "INSUFFICIENT (no reference)")
    else:
        a1_notes.append(f"参照表覆盖 {res.reference_genes} 个基因")
    res.feasibility["audit_1_floor"] = {
        "level": "ok" if res.counts_ok else "degraded",
        "status": "可运行" if res.counts_ok else "可运行但有前提",
        "notes": a1_notes}

    a2_notes = []
    # 是否存在取材/来源类注释列（含常量的——常量意味着"全部同一来源"，即已平衡）
    design_named = [c for c in ad.obs.columns
                    if c not in (condition_col, donor_col)
                    and _name_role(c, cfg) == "design"]
    max_tvd = max((d["tvd_donor"] for d in res.design_vars), default=0.0)

    if usable:
        a2_level = "ok"
        a2_status = f"可运行 · 推荐变量：{usable}"
        a2_notes.append(f"推荐变量（供者级 TVD 落在 "
                        f"[{cfg.match_var_tvd_low}, {cfg.hetero_max}] 之间）：{usable}")
    elif design_named and max_tvd < cfg.match_var_tvd_low and not risky:
        # 有取材注释，但各组在它上面已平衡 → 没有需要控制的混杂
        a2_level = "ok"
        a2_status = "无需匹配（两组在设计变量上已平衡）"
        a2_notes.append(f"检测到取材/来源类注释 {design_named}，但其供者级 TVD "
                        f"均 < {cfg.match_var_tvd_low} → 两组在来源上已平衡，"
                        f"matched comparison 不会改变结论")
    elif risky:
        a2_level = "degraded"
        a2_status = "可运行，但所有候选都与分组接近分离"
        a2_notes.append(
            f"全部候选供者级 TVD > {cfg.hetero_max}：{risky}。"
            f"匹配它们在控制来源差异的同时会抹掉真实的病例/对照差异"
            f"（把生物学一起匹配走），判定很可能直接是 NOT SEPARABLE")
        a2_notes.append("建议：改用【分层 + 供者级 TVD】作为可比性依据，"
                        "而非 matched comparison（见 audit_2_recommendation）")
    elif safe_only or cautions:
        a2_level = "degraded"
        a2_status = "可运行，但匹配收益有限"
        if safe_only:
            a2_notes.append(f"低信息变量：{safe_only}（与分组近乎独立，匹配几乎不控制任何东西）")
        if cautions:
            a2_notes.append(f"需人工甄别的变量：{cautions}"
                            f"（只在单一组内变化——可能是结局变量，也可能是取材变量）")
    else:
        a2_level, a2_status = "no", "不可运行 → 将报 INDETERMINATE"
        a2_notes.append(
            (f"没有任何可用设计变量（候选 {cands or '空'}；剔除：{dropped or '无'}；"
             f"按角色排除：{flagged or '无'}）"
             if cands else "obs 中没有可作为设计变量的类别型列")
            + " → 审计二将只能报 INDETERMINATE")
        a2_notes.append("原因：数据缺少 institution/donation_type/source/region 类注释；"
                        "这不是工具故障，是数据生态问题")
    res.feasibility["audit_2_composition"] = {
        "level": a2_level, "status": a2_status, "notes": a2_notes}

    # ── audit_2 总体建议（含可直接用于方法学描述的措辞）──
    if usable:
        a2r = {"method": "matched_comparison", "feasible": True,
               "recommended_variables": usable, "statement_for_paper": ""}
    elif a2_level == "ok":
        a2r = {"method": "no_matching_needed", "feasible": True,
               "recommended_variables": [],
               "statement_for_paper":
                   "Control and case groups were balanced on all available "
                   "procurement/source annotations (donor-level total-variation "
                   "distance below the pre-specified margin), so no matched "
                   "comparison was required."}
    elif risky:
        _top = max(res.design_vars, key=lambda d: d["tvd_donor"])
        a2r = {"method": "tvd_stratification", "feasible": False,
               "recommended_variables": [],
               "statement_for_paper":
                   "Control and case groups differed systematically in "
                   "procurement/source variables (donor-level total-variation "
                   f"distance up to {_top['tvd_donor']:.2f} for "
                   f"'{_top['column']}'); because no candidate variable could be "
                   "matched without also removing the biological contrast, we "
                   "report the imbalance explicitly and base comparability on the "
                   "donor-level total-variation distance rather than on a matched "
                   "comparison."}
    elif safe_only or cautions:
        a2r = {"method": "matched_comparison_marginal", "feasible": True,
               "recommended_variables": safe_only or cautions,
               "statement_for_paper":
                   "No candidate variable was both informative and non-confounded; "
                   "the matched comparison is reported with this caveat."}
    else:
        a2r = {"method": "tvd_only", "feasible": False,
               "recommended_variables": [],
               "statement_for_paper":
                   "No design/procurement annotation was available in this dataset; "
                   "a matched comparison could not be attempted, and comparability "
                   "was assessed only through the donor-level total-variation "
                   "distance where such variables existed."}
    a2r["excluded_by_role"] = flagged
    a2r["excluded_as_unusable"] = dropped
    a2r["needs_manual_review"] = cautions

    # ── 当存在"接近分离"变量时，给出三套可选策略（**不排序**，由使用者权衡）──
    if risky:
        _pref = usable or safe_only
        a2r["options"] = [
            {"id": "conservative",
             "label": "只用推荐变量匹配（最保守）",
             "config_hint": {"match_vars": _pref or []},
             "trade_off": f"保留可定量细胞数；但无法控制 {risky} 的混杂，"
                          f"须在方法学描述中把该局限写明"},
            {"id": "primary_plus_sensitivity",
             "label": "主分析 + 分层/敏感性分析（平衡）",
             # 键必须是 AuditConfig 真实字段（strata_cols 即审计三的分层依据），
             # 由 tests/test_check.py::test_all_config_hints_are_real_config_fields 钉住
             "config_hint": {"match_vars": _pref or [], "strata_cols": risky},
             "trade_off": f"主分析用 {_pref or '无'}；对 {risky} 分层或作为敏感性分析。"
                          f"兼顾混杂控制与生物学完整性，代价是需报告两组结果"},
            {"id": "no_matching",
             "label": "放弃匹配，以供者级 TVD 为可比性依据",
             "config_hint": {"run_matched_comparison": False},
             "trade_off": "承认设计共线、不做匹配；审计二以 TVD 判定，"
                          "在方法学描述中显式声明该局限。**这不是'退而求其次'，"
                          "当混杂是结构性的时它是最诚实的做法**"},
        ]
        a2r["options_note"] = ("以上三套策略的取舍取决于你对该变量是否代表生物学差异的"
                               "判断——工具不替你排序，也不写入配置。")
    res.audit_2_recommendation = a2r

    # ── 参照表覆盖（分母应是目标基因面板，不是数据集全部基因）──
    tg = [str(g) for g in (target_genes or [])]
    refset = set((reference or {}).keys())
    covered = [g for g in tg if g in refset]
    res.reference_coverage = {
        "n_reference_genes": len(refset),
        "n_target_genes": len(tg),
        "target_covered": covered,
        "target_missing": [g for g in tg if g not in refset],
        "pct_target_covered": (round(100.0 * len(covered) / len(tg), 1) if tg else None),
        "note": "只有被参照表覆盖的基因才能得到基于 R 的判定；未覆盖者将报 "
                "INSUFFICIENT (no reference)。分母应是目标基因面板，而非数据集全部基因。",
    }
    if reference is None and tg:
        # C6（v0.7.2，据第三方独立评估）：reference 未传入时，覆盖率 0 是“**未传入**”所致，
        # 不得断言“参照表只覆盖 0/N”（会造成“参照表缺少这些基因”的误读）。
        res.warnings.append(
            f"未提供参照表（reference=None）→ 覆盖率 {len(covered)}/{len(tg)} 为 0 系"
            f"**未传入**所致，而非参照表缺少这些基因；审计一无法计算 R，"
            f"请用 --reference 传入参照表后重跑。")
    elif tg and len(covered) < len(tg):
        res.warnings.append(
            f"参照表只覆盖 {len(covered)}/{len(tg)} 个目标基因；其余 "
            f"{len(tg) - len(covered)} 个将报 INSUFFICIENT (no reference)")

    # ── 参照表溯源（门控参数登记）：R 的可比性依赖门控口径一致 ──
    if reference_manifest:
        bgs = reference_manifest.get("backgrounds") or {}
        unknown = [k for k, v in bgs.items()
                   if str(v.get("provenance", "")).upper().startswith("NOT_")]
        res.reference_provenance = {
            "manifest_found": True,
            "n_backgrounds": len(bgs),
            "pct_ref_rule": reference_manifest.get("pct_ref_rule", ""),
            "gate_unknown_backgrounds": unknown,
            "known_gaps": reference_manifest.get("known_gaps", []),
            "formula": "R = pct_case / pct_ref",
        }
        if unknown:
            res.warnings.append(
                f"参照表中 {len(unknown)} 个背景的门控参数未登记：{unknown} → "
                f"R 的分母含口径未知的分量，跨数据集比较须谨慎")
    else:
        res.reference_provenance = {"manifest_found": False}
        res.warnings.append(
            "未找到参照表的 manifest（门控参数登记）→ R 的跨数据集可比性无法核验")

    db = res.depth_bins_ok
    a3_level = "ok" if db["n_bins_ok"] >= 1 else "degraded"
    res.feasibility["audit_3_robust"] = {
        "level": a3_level,
        "status": "可运行" if db["n_bins_ok"] >= 1 else "可运行但分层不足",
        "notes": [f"{db['n_bins_ok']}/{db['n_bins']} 个深度箱满足"
                  f"两组各 >= {db['min_cells_per_stratum']} 细胞"
                  + ("" if db["n_bins_ok"] else
                     "；分层检出率将无法产出，审计三退化为仅深度指标对比")]}

    # ── ⑦ 建议 ──
    if usable:
        res.recommendations.append(
            f"把选定的设计变量写入 config 的 match_vars（例如 {usable[:2]}），"
            f"再运行 composition —— 这一步由你决定，工具不会代选")
    over = [d["column"] for d in res.design_vars if d["status"] == "可用·需注意"]
    if over:
        res.recommendations.append(
            f"以下变量的供者级 TVD 已超 hetero_max={cfg.hetero_max}：{over}。"
            f"纳入其中任一并命中，判定将直接指向 NOT SEPARABLE —— 这通常是"
            f"**正确**的信号（该变量与分组共线），但若你认为它不代表取材/来源差异，"
            f"请显式记录理由。")
    # 方法学警告：候选列里往往混有结局/治疗后变量，纳入匹配会造成 collider bias
    res.recommendations.append(
        "匹配变量应是【取材/来源相关】的设计变量（机构、取材方式、分选来源、区域）。"
        "不要把结局或治疗后变量放进 match_vars —— 对结局变量匹配会引入 collider bias，"
        "使比较失去意义。"
        "本命令只做候选发现，不判断变量的方法学角色。")
    if res.n_donors_ctrl < 5 or res.n_donors_case < 5:
        res.recommendations.append(
            "供者数偏少：考虑将判定主要建立在检出率与事件数上，而非置换 P 值")
    if not res.counts_ok:
        res.recommendations.append(
            "用 --count-layer 指定原始计数层，或先 ad.X = ad.layers['counts']")
    return res
