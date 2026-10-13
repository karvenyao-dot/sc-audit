"""审计二 · v0.3 可识别性层（identifiability）的契约测试。

覆盖（每条对应 README / docs/methods.md 的一句可验证声明）：
  1. 花名册完备性：obs 每个非 QC 列必须落桶，未归类数 = 0，且计数自洽
  2. 变量分桶：病例组内恒定 → design-defining；病例组内有变异 → adjustable
  3. 共同支持域：只保留落在病例组设计剖面内的对照供者
  4. 阶梯锚点：R0（不匹配）的衰减比恒为 1.0
  5. 支持域下限：低于 matched_min_cells_warn → SUPPORT-BELOW-FLOOR 且抑制衰减比
  6. **单侧判定（2026-10-06 审计后修订）**：阈值落在观测区间之外 → 直接给方向判定
     （COMPARABLE / NOT SEPARABLE）**且必附强制警示**；两档并存 → THRESHOLD-SPANNING；
     关闭 `identifiability_one_sided_verdict` 时两端均报 INDETERMINATE（旧保守行为）
  7. 显式阶梯：API 传 ladder 被原样采纳
  8. 健壮性（2026-10-06 审计 P0/P1/P2 回归）：缺设计列不崩、0 基因命中不崩、
     空阶梯不崩、重复基因去重、`min_support_cells=0` 被尊重、多病例剖面不池化
  9. 外部数据冒烟回归（可选，h5ad 不存在则跳过）：不崩、判定合法、产物结构完整
     （具体预注册数字的复现在开发仓库的内部脚本里，不随公开树发布）
"""
import os

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from conftest import make_adata
from sc_audit.config import AuditConfig
from sc_audit.identifiability import (BUCKET_ADJUSTABLE, BUCKET_DESIGN_DEFINING,
                                      VERDICT_COMPARABLE, VERDICT_INDETERMINATE,
                                      VERDICT_SUPPORT_BELOW_FLOOR,
                                      VERDICT_THRESHOLD_SPANNING,
                                      audit_identifiability, build_variable_roster,
                                      common_support, _nondesign_columns)

# 外部数据冒烟回归（可选）：设置 SC_AUDIT_EXTERNAL_CASE_H5AD 后才会执行；否则跳过。
EXTERNAL_CASE = os.environ.get("SC_AUDIT_EXTERNAL_CASE_H5AD", "")
PANEL = ["TNNT2", "PLN", "ACTB", "MYL3"]


def _cfg(**kw):
    c = AuditConfig()
    for k, v in kw.items():
        setattr(c, k, v)
    return c


def _audit(ad, cfg=None, genes=PANEL, ladder=None):
    return audit_identifiability(ad, genes, "condition", "donor", "control",
                                 cfg=cfg or AuditConfig(), ladder=ladder)


def _set_detection(ad, gene, pcts_by_donor):
    """把某基因在指定供者上的检出率设为给定值（测试用，直接改矩阵）。"""
    gi = list(ad.var_names).index(gene)
    X = ad.X.toarray()
    rng = np.random.default_rng(7)
    for d, p in pcts_by_donor.items():
        m = (ad.obs["donor"] == d).values
        X[m, gi] = (rng.random(int(m.sum())) < p).astype(X.dtype)
    ad.X = sp.csr_matrix(X)
    return ad


# ──────────────────── ① 花名册完备性 ────────────────────

def test_roster_is_complete_no_silent_omission():
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=15)
    R, unclassified = build_variable_roster(ad, "condition", "donor", "control")
    nondesign = _nondesign_columns(ad, "condition", "donor", AuditConfig())
    assert unclassified == []
    assert len(R) + len(nondesign) == ad.obs.shape[1], \
        "花名册不完备：存在既未归类、也未被声明为非设计的列"
    assert set(R["bucket"]).isdisjoint({""}), "存在空桶"


def test_unbucketed_column_is_surfaced_not_skipped(monkeypatch):
    """列无法归类时必须显式浮出（unclassified 非空），不得静默跳过。"""
    ad = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=10)
    ad.obs["weird"] = [object() for _ in range(ad.n_obs)]

    real = __import__("sc_audit.identifiability", fromlist=["x"])

    def boom(*a, **k):
        raise TypeError("cannot bucket")

    monkeypatch.setattr(real, "_group_distribution_distance", boom)
    R, unclassified = real.build_variable_roster(ad, "condition", "donor", "control")
    assert any("weird" in u for u in unclassified)


# ──────────────────── ② 变量分桶 ────────────────────

def test_confounded_design_variables_are_design_defining():
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=15)
    R, _ = build_variable_roster(ad, "condition", "donor", "control")
    dd = dict(zip(R["column"], R["bucket"]))
    for c in ("institution", "donation_type", "source", "region"):
        assert dd[c] == BUCKET_DESIGN_DEFINING, f"{c} 应判 design-defining，实为 {dd[c]}"
    assert R.loc[R.column == "batch", "bucket"].iloc[0] == "excluded:donor-unique"


def test_varying_within_cases_is_adjustable():
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=15)
    vals = {}
    for i in range(5):
        vals[f"case_d{i}"] = "young" if i < 3 else "old"
    for i in range(5):
        vals[f"ctrl_d{i}"] = "young" if i % 2 == 0 else "old"
    ad.obs["age_range"] = ad.obs["donor"].map(vals)
    R, _ = build_variable_roster(ad, "condition", "donor", "control")
    assert R.loc[R.column == "age_range", "bucket"].iloc[0] == BUCKET_ADJUSTABLE
    assert R.loc[R.column == "age_range", "ent_case"].iloc[0] > 0


def test_expression_derived_score_is_excluded_from_bucketing():
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=12)
    ad.obs["score_Repolarization"] = np.linspace(0, 1, ad.n_obs)
    R, _ = build_variable_roster(ad, "condition", "donor", "control")
    assert R.loc[R.column == "score_Repolarization", "bucket"].iloc[0] == \
        "excluded:derived-from-expression"


# ──────────────────── ③ 共同支持域 ────────────────────

def test_common_support_keeps_only_same_profile_controls():
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=15)
    sup = common_support(ad, ["institution"], "condition", "donor", "control")
    assert sup["support_donors"] == ["ctrl_d0"], "只应保留与病例同机构的对照供者"
    assert sup["support_cells"] == 15
    assert sup["n_control_donors"] == 5


def test_support_is_vacuous_without_design_defining_variable():
    ad = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=10)
    sup = common_support(ad, [], "condition", "donor", "control")
    assert sup["vacuous"] is True and len(sup["support_donors"]) == 3


def test_multiple_case_profiles_are_not_pooled():
    """A1-F6：不得把多个病例剖面的对照并成一个支持域（防跨层比较）。

    注：在 `audit_identifiability` 主流程里，design-defining 变量在病例组内恒定
    （熵 = 0）⇒ 病例剖面按构造唯一，故多剖面分支是**防御性**的；
    真正的并集行为发生在 `common_support` 被传入"病例组内有变异的设计变量"时。
    """
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=10)
    vals = {f"case_d{i}": ("INST_A" if i < 2 else "INST_B") for i in range(4)}
    vals.update({f"ctrl_d{i}": ("INST_A" if i < 2 else "INST_B") for i in range(4)})
    ad.obs["institution"] = ad.obs["donor"].map(vals)
    sup = common_support(ad, ["institution"], "condition", "donor", "control")
    assert len(sup["case_profiles"]) == 2
    assert len(sup["support_by_profile"]) == 2
    assert len(sup["support_donors"]) == 2, "支持域必须只来自单一（最大）剖面"
    # 主流程里病例剖面按构造唯一（design-defining ⇒ 病例组内恒定）
    res = _audit(ad, _cfg(matched_min_cells_warn=0))
    assert len(res.case_profiles) == 1
    assert res.details["multi_case_profile"] is False


# ──────────────────── ④⑦ 阶梯 ────────────────────

def test_ladder_anchor_r0_has_ratio_one():
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=15,
                    gene_pcts={"TNNT2": (0.9, 0.2)})
    res = _audit(ad)
    r0 = res.ladder.loc[res.ladder.rule.str.startswith("R0")]
    assert len(r0) and np.allclose(r0["ratio"].dropna(), 1.0), "锚点规则衰减比必须为 1"


def test_explicit_ladder_is_honoured():
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=15,
                    gene_pcts={"TNNT2": (0.9, 0.2)})
    res = _audit(ad, ladder=[[], ["institution"], ["source"]])
    assert res.ladder_rules == ["R0 不匹配（锚点）", "R1 institution", "R2 source"]
    assert any("pre-specified" in n for n in res.notes)


def test_to_text_contains_all_four_sections():
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=15,
                    gene_pcts={"TNNT2": (0.9, 0.2)})
    txt = _audit(ad).to_text()
    for tag in ("① variable roster", "② common support", "③ prespecified matching ladder",
                "④ verdict"):
        assert tag in txt, f"文本报告缺少 {tag}"


# ──────────────────── ⑤⑥ 判定 ────────────────────

def test_below_cell_floor_suppresses_attenuation_and_reports_support_below_floor():
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=10,
                    gene_pcts={"TNNT2": (0.9, 0.2)})
    res = _audit(ad, _cfg(matched_min_cells_warn=50))
    assert res.support_cells == 10
    assert res.verdict.startswith(VERDICT_SUPPORT_BELOW_FLOOR)
    assert res.attenuated_suppressed is True
    assert any("NOT reported as evidence" in n for n in res.notes)


def test_one_sided_verdict_gives_direction_plus_mandatory_caveat():
    """单侧情形（阈值落在观测区间之外）：给方向判定 + **必附**警示。

    这是 2026-10-06 审计 A1-F1/F2 的修订点：旧实现在默认配置下
    COMPARABLE / NOT SEPARABLE 永不可达（一致性检查退化为空操作）。
    """
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=20,
                    gene_pcts={"TNNT2": (0.9, 0.4)})
    res = _audit(ad, _cfg(matched_min_cells_warn=0, attenuation_min=0.01))
    assert res.details["threshold_spanning"] is False
    assert res.details["n_low"] == 0 and res.details["n_high"] > 0
    assert res.verdict.startswith(VERDICT_COMPARABLE), \
        f"单侧（全档≥阈值）应报 COMPARABLE，实为 {res.verdict}"
    assert any("ONE-SIDED CAVEAT" in n for n in res.notes), "缺少强制警示"
    assert any("consistency check carries NO information" in n for n in res.notes)
    assert any("does NOT establish comparability" in n for n in res.notes)


def test_one_sided_verdict_can_be_disabled_for_the_conservative_reading():
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=20,
                    gene_pcts={"TNNT2": (0.9, 0.4)})
    res = _audit(ad, _cfg(matched_min_cells_warn=0, attenuation_min=0.01,
                          identifiability_one_sided_verdict=False))
    assert res.verdict.startswith(VERDICT_INDETERMINATE)
    assert any("one-sided outcomes are reported as INDETERMINATE" in n
               for n in res.notes)


def test_spanning_threshold_reports_threshold_spanning():
    """两档并存 → THRESHOLD-SPANNING（衰减比不是跨规则的单一估计量）。"""
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=30)
    _set_detection(ad, "TNNT2", {**{f"ctrl_d{i}": (0.9 if i < 3 else 0.1)
                                     for i in range(6)},
                                  **{f"case_d{i}": 0.3 for i in range(6)}})
    vals = {}
    for i in range(6):
        vals[f"ctrl_d{i}"] = "young" if i % 2 == 0 else "old"
    for i in range(6):
        vals[f"case_d{i}"] = "young" if i % 2 == 0 else "old"
    ad.obs["age_range"] = ad.obs["donor"].map(vals)
    res = _audit(ad, _cfg(matched_min_cells_warn=0, attenuation_min=0.5))
    assert res.details["n_low"] > 0 and res.details["n_high"] > 0, \
        f"本用例应产生两档并存的阶梯：{res.ladder_summary}"
    assert res.verdict.startswith(VERDICT_THRESHOLD_SPANNING), res.verdict
    assert any("per-rule medians" in n for n in res.notes)


def test_amplification_is_reported_when_ratio_exceeds_one():
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=25)
    _set_detection(ad, "TNNT2", {**{f"ctrl_d{i}": 0.05 for i in range(6)},
                                 **{f"case_d{i}": 0.9 for i in range(6)}})
    res = _audit(ad, _cfg(matched_min_cells_warn=0, attenuation_min=0.5))
    assert any("effect-size ratio > 1" in n or "scale note" in n for n in res.notes)


# ──────────────────── ⑧ 健壮性回归（审计 P0/P1/P2） ────────────────────

def test_missing_design_columns_do_not_crash():
    """P0-1：内置阶梯 R1 引用 obs 不存在的列 → 曾经 KeyError。"""
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=10)
    for c in ["institution", "donation_type", "source", "region", "batch"]:
        del ad.obs[c]
    res = _audit(ad, genes=["TNNT2", "ACTB"])
    assert res.support_cells > 0
    assert any("are absent from obs" in n for n in res.notes)
    # 缺列的规则不得退化成"不匹配"（那会把 ratio 写成 1.0）
    r1 = res.ladder_summary.loc[res.ladder_summary.rule.str.startswith("R1")].iloc[0]
    assert int(r1.n_cells) == 0, "缺列规则应判为不可评估（0 细胞）"
    assert r1.rule in res.excluded_from_consistency


def test_zero_target_genes_does_not_crash():
    """P0-2：目标基因全部不在 var_names → 曾经 AttributeError。"""
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=10)
    res = _audit(ad, _cfg(matched_min_cells_warn=0), genes=["NOPE1", "NOPE2"])
    assert res.verdict.startswith(VERDICT_INDETERMINATE)
    assert res.details["n_panel_genes_found"] == 0
    assert any("no target gene was found" in n for n in res.notes)
    assert list(res.ladder.columns), "空阶梯也必须带列名（否则下游 pivot 崩）"
    assert list(res.ladder_summary.columns)


def test_empty_gene_list_does_not_crash():
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=10)
    res = _audit(ad, _cfg(matched_min_cells_warn=0), genes=[])
    assert res.verdict.startswith(VERDICT_INDETERMINATE)


def test_empty_ladder_does_not_crash():
    """P0-2：ladder=[]（配置文件路径）→ 曾经 AttributeError。"""
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=10)
    res = _audit(ad, genes=["TNNT2", "ACTB"], ladder=[])
    assert res.ladder_rules == []
    assert list(res.ladder_summary.columns), "空阶梯也必须带列名"
    res2 = _audit(ad, _cfg(matched_min_cells_warn=0), genes=["TNNT2", "ACTB"], ladder=[])
    assert res2.verdict.startswith(VERDICT_INDETERMINATE)
    assert any("ladder is empty" in n for n in res2.notes)


def test_duplicate_genes_are_deduplicated():
    """P1-5：重复目标基因不得重复加权中位数。"""
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=15,
                    gene_pcts={"TNNT2": (0.9, 0.2)})
    res = _audit(ad, genes=["TNNT2", "TNNT2", "TNNT2", "ACTB"])
    assert sorted(res.ladder.gene.unique()) == ["ACTB", "TNNT2"]


def test_min_support_cells_zero_is_honoured():
    """P2-8：`0` 是合法值（不设下限），不得被 `or` 吞掉回落 50。"""
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=10,
                    gene_pcts={"TNNT2": (0.9, 0.2)})
    res = _audit(ad, _cfg(identifiability_min_support_cells=0))
    assert res.cell_floor == 0
    assert not res.verdict.startswith(VERDICT_SUPPORT_BELOW_FLOOR)


def test_chance_constancy_risk_is_reported_for_small_case_groups():
    ad = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=10,
                    gene_pcts={"TNNT2": (0.9, 0.2)})
    res = _audit(ad, _cfg(matched_min_cells_warn=0))
    note = next(n for n in res.notes if "chance-constancy risk" in n)
    assert "P=0.3" in note, f"风险值应在参照人群口径下计算（非退化 1.00）: {note}"
    assert "P=1.000" not in note, "退化值（病例组内占比恒为 1）不得出现"


def test_chance_constancy_risk_semantics():
    """B5（外部审查）：参照占比口径必须正确，且参照人群无变异时不可评估。"""
    from sc_audit.donor_utils import chance_constancy_risk
    assert abs(chance_constancy_risk(0.70, 12) - 0.0138) < 5e-4, "p=0.70,n=12 应为 ~0.0138"
    assert abs(chance_constancy_risk(0.5, 3) - 0.25) < 1e-9
    assert np.isnan(chance_constancy_risk(1.0, 12)), "p=1（参照人群无变异）应不可评估"
    assert np.isnan(chance_constancy_risk(0.0, 12))


def test_from_dict_warns_on_unknown_keys():
    """P2-11：拼错的锁定阈值不得被静默忽略。"""
    with pytest.warns(UserWarning, match="未知字段"):
        AuditConfig.from_dict({"floor_pct": 5.0, "floop_pct": 1.0})


def test_verdict_colour_mapping_treats_new_labels_as_bad():
    """P1-3：v0.3 负面判定不得被渲染成绿色。"""
    from sc_audit.report import _verdict_class
    for v in (VERDICT_SUPPORT_BELOW_FLOOR, VERDICT_THRESHOLD_SPANNING,
              "NOT SEPARABLE — x", "UNQUANTIFIABLE"):
        assert _verdict_class(v) == "bad", v
    assert _verdict_class(VERDICT_INDETERMINATE + " — x") == "warn"
    assert _verdict_class("QUANTIFIABLE") == "ok"


def test_single_group_condition_returns_indeterminate():
    """P2-10：condition 只有一个水平 → 显式 INDETERMINATE（不得 NaN 传播）。"""
    ad = make_adata(n_ctrl_donors=4, n_case_donors=0, cells_per_donor=10)
    res = _audit(ad, genes=["TNNT2", "ACTB"])
    assert res.verdict.startswith(VERDICT_INDETERMINATE)
    assert "does not contain both groups" in res.verdict
    assert any("--control-label" in n or "control label" in n for n in res.notes)


def test_all_missing_column_is_bucketed_separately():
    """P2-10：全缺失列不得被判成 design-defining 并进入适用范围声明。"""
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=10)
    ad.obs["blank_col"] = None
    R, _ = build_variable_roster(ad, "condition", "donor", "control")
    assert R.loc[R.column == "blank_col", "bucket"].iloc[0] == "excluded:all-missing"
    res = _audit(ad, _cfg(matched_min_cells_warn=0), genes=["TNNT2", "ACTB"])
    assert "blank_col" not in res.design_defining


# ──────────────────── ⑨ 真实数据回归（可选） ────────────────────

@pytest.mark.skipif(not (EXTERNAL_CASE and os.path.exists(EXTERNAL_CASE)),
                    reason="未设置 SC_AUDIT_EXTERNAL_CASE_H5AD（外部数据冒烟回归跳过）")
def test_external_dataset_smoke_regression():
    """外部数据冒烟回归：真实 h5ad 上不崩、判定落在已定义标签集合内、产物结构完整。

    ⚠️ 具体预注册数字（支持域供者/细胞数、阶梯中位衰减比等）**不在公开树内断言**；
    内部复现脚本不随公开树发布（设 `SC_AUDIT_EXTERNAL_CASE_H5AD` 后可在开发仓库内运行）。
    """
    import scanpy as sc
    from sc_audit.identifiability import (VERDICT_COMPARABLE, VERDICT_INDETERMINATE,
                                          VERDICT_NOT_SEPARABLE,
                                          VERDICT_SUPPORT_BELOW_FLOOR,
                                          VERDICT_THRESHOLD_SPANNING)
    ad = sc.read_h5ad(EXTERNAL_CASE)
    res = _audit(ad, genes=["MYL3", "PLN", "TNNT2", "ACTB", "GAPDH"])
    assert res.verdict.split(" — ")[0] in {
        VERDICT_COMPARABLE, VERDICT_NOT_SEPARABLE, VERDICT_INDETERMINATE,
        VERDICT_SUPPORT_BELOW_FLOOR, VERDICT_THRESHOLD_SPANNING}
    assert res.n_obs_columns > 0 and res.support_cells >= 0
    assert len(res.ladder_summary) >= 2, "阶梯至少应有锚点 + 一条匹配规则"


def test_all_verdicts_are_reachable():
    """死分支防护（A1-F1 的教训）：每个判定标签都必须能被**某个**场景真实产出。

    v0.3.0 曾出现"COMPARABLE / NOT SEPARABLE 永不可达"的死分支，
    原因是两个谓词被写成了同一个。此测试确保标签集合与实现一致。
    """
    from sc_audit.identifiability import (VERDICT_NOT_SEPARABLE, VERDICT_SUPPORT_BELOW_FLOOR)
    seen = set()

    # ① SUPPORT-BELOW-FLOOR：支持域 < 下限
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=10,
                    gene_pcts={"TNNT2": (0.9, 0.2)})
    seen.add(_audit(ad, _cfg(matched_min_cells_warn=50)).verdict.split(" — ")[0])

    # ② COMPARABLE：单侧（全部规则 ≥ 阈值）
    ad2 = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=20,
                     gene_pcts={"TNNT2": (0.9, 0.4)})
    seen.add(_audit(ad2, _cfg(matched_min_cells_warn=0,
                              attenuation_min=0.01)).verdict.split(" — ")[0])
    # ③ NOT SEPARABLE：单侧（全部规则 < 阈值）—— 把阈值抬到观测区间之上
    seen.add(_audit(ad2, _cfg(matched_min_cells_warn=0,
                              attenuation_min=99.0)).verdict.split(" — ")[0])
    # ④ THRESHOLD-SPANNING：两档并存
    ad3 = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=30)
    _set_detection(ad3, "TNNT2", {**{f"ctrl_d{i}": (0.9 if i < 3 else 0.1) for i in range(6)},
                                  **{f"case_d{i}": 0.3 for i in range(6)}})
    ad3.obs["age_range"] = ad3.obs["donor"].map(
        {**{f"ctrl_d{i}": ("young" if i % 2 == 0 else "old") for i in range(6)},
         **{f"case_d{i}": ("young" if i % 2 == 0 else "old") for i in range(6)}})
    seen.add(_audit(ad3, _cfg(matched_min_cells_warn=0,
                              attenuation_min=0.5)).verdict.split(" — ")[0])
    # ⑤ INDETERMINATE：无目标基因
    seen.add(_audit(ad2, _cfg(matched_min_cells_warn=0),
                    genes=["NOPE"]).verdict.split(" — ")[0])

    expected = {VERDICT_SUPPORT_BELOW_FLOOR, VERDICT_COMPARABLE, VERDICT_NOT_SEPARABLE,
                VERDICT_THRESHOLD_SPANNING, VERDICT_INDETERMINATE}
    assert seen == expected, f"存在不可达的判定标签：{expected - seen}"


def test_no_design_variable_yields_indeterminate_not_comparable():
    """可适用性修正（2026-10-06 跨数据集实跑）：缺元数据 ≠ 可比。

    当 obs 里没有任何设计/可调整变量（公共数据常见：只有 condition/donor），
    阶梯的非锚点规则会退化为"不匹配"，衰减比恒 = 1.000 → 旧逻辑会输出
    COMPARABLE，把"没有元数据"误读成"两组可比"。现改为 INDETERMINATE。
    """
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=15,
                    gene_pcts={"TNNT2": (0.9, 0.2)})
    for c in ["institution", "donation_type", "source", "region", "batch"]:
        if c in ad.obs.columns:
            del ad.obs[c]
    res = _audit(ad, genes=["TNNT2", "ACTB"])
    assert res.verdict.startswith(VERDICT_INDETERMINATE), res.verdict
    assert not res.verdict.startswith(VERDICT_COMPARABLE)
    assert res.details["ladder_uses_real_variable"] is False
    assert res.attenuated_suppressed is True
    assert any("ANNOTATION limitation" in n for n in res.notes)
    # 对照：有设计变量时不走这条分支
    ad2 = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=15,
                     gene_pcts={"TNNT2": (0.9, 0.2)})
    res2 = _audit(ad2, _cfg(matched_min_cells_warn=0), genes=["TNNT2", "ACTB"])
    assert res2.details["ladder_uses_real_variable"] is True


# ──────────────────── ⑩ 零信息列的守卫（v0.3.3） ────────────────────

def test_constant_zero_information_column_does_not_flip_to_comparable():
    """在该数据集内**恒定**（零信息）的列不得把真空支持域翻转成 COMPARABLE（v0.3.3）。

    机制：所有细胞同值的列会被归入 design-defining 桶 → 旧的 `vacuous` 守卫（仅在
    design-defining 为空时成立）失效 → 支持域 = 整个对照组、每条规则的衰减比恒 = 1
    → 落进单侧分支被报成 `COMPARABLE`，与本工具自己的文案
    （"do not read this state as 'COMPARABLE'"）自相矛盾。
    现：**支持域 = 整个对照组且无任何可评估规则缩域** → `INDETERMINATE`（注释局限）。
    """
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=20,
                    gene_pcts={"TNNT2": (0.9, 0.2)})
    for c in ["institution", "donation_type", "source", "region", "batch"]:
        if c in ad.obs.columns:
            del ad.obs[c]
    ad.obs["institution"] = "A"          # 全队列恒定 → 零信息
    ad.obs["region"] = "X"               # 换名同样命中
    res = _audit(ad, genes=["TNNT2", "ACTB"])
    assert res.verdict.startswith(VERDICT_INDETERMINATE), res.verdict
    assert not res.verdict.startswith(VERDICT_COMPARABLE)
    assert res.details["support_cells"] == res.details["n_control_cells"]
    assert res.attenuated_suppressed is True
    assert any("ANNOTATION limitation" in n for n in res.notes)


def test_real_stratum_column_still_reaches_support_below_floor():
    """真分层列（只有 1 个对照供者与病例同机构）仍走 SUPPORT-BELOW-FLOOR。

    守卫放宽**不得**改变"病例组在设计变量上恒定、支持域被压到地板以下"这条路径。
    """
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=20,
                    gene_pcts={"TNNT2": (0.9, 0.2)})
    donors = ad.obs["donor"].astype(str)
    ad.obs["institution"] = ["A" if d == "ctrl_d0" else "B" for d in donors]
    ad.obs.loc[ad.obs["condition"] == "disease", "institution"] = "A"
    res = _audit(ad, genes=["TNNT2", "ACTB"])
    assert res.verdict.startswith(VERDICT_SUPPORT_BELOW_FLOOR), res.verdict
    assert res.details["support_cells"] < res.details["n_control_cells"]
