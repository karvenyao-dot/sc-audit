"""预检子命令（`sc-audit check`）的功能测试。

核心契约（必须被钉住）：
  1. check 只【报告】，绝不写入/修改配置 —— 这是它不破坏"事前锁定"的前提
  2. 供者级唯一的列必须被判为不可用于匹配
  3. 缺少可用设计变量时，审计二必须被明确标为不可行（而非静默通过）
  4. 阻断项（分组/供者列缺失、计数不可用）必须出现在报告里，而不是抛异常中断
  5. 退出码语义：0 全可行 / 1 有审计不可行 / 2 存在阻断项
"""
import json
import os

import pytest

from sc_audit.check import check_readiness
from sc_audit.cli import main as cli_main
from sc_audit.config import AuditConfig

from conftest import make_adata


def _cfg(**kw):
    d = dict(n_perm=100, min_cells_per_stratum=5)
    d.update(kw)
    return AuditConfig(**d)


# ════════════════ 纯函数层 ════════════════

def test_check_reports_usable_and_donor_unique_columns():
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20)
    # 供者级唯一列（真实数据中 batch/lane/Channel 的形态）
    ad.obs["channel"] = ["c_" + str(d) for d in ad.obs["donor"].astype(str)]
    res = check_readiness(ad, "condition", "donor", "control",
                          cfg=_cfg(), target_genes=["TNNT2"])
    by = {d["column"]: d for d in res.design_vars}
    assert by["channel"]["donor_unique"] is True
    assert by["channel"]["status"] == "剔除"
    assert by["institution"]["status"].startswith("可用")


def test_check_flags_confounded_design_variable():
    """TVD 超阈值的变量必须被标为'可用·需注意'并说明后果。"""
    ad = make_adata(design_confounded=True)          # institution 与分组共线
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg())
    inst = {d["column"]: d for d in res.design_vars}["institution"]
    assert inst["tvd_donor"] > 0.5
    assert inst["status"] == "可用·需注意"
    assert "NOT SEPARABLE" in inst["reason"]


def test_check_marks_audit2_infeasible_without_design_vars():
    """没有任何可用设计变量 → 审计二必须被标为不可行，并说明是数据生态问题。"""
    ad = make_adata()
    for c in ("institution", "donation_type", "source", "region"):
        if c in ad.obs.columns:
            ad.obs[c] = ["v_" + str(d) for d in ad.obs["donor"].astype(str)]  # 全剔
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg(),
                          candidate_vars=["institution", "donation_type",
                                          "source", "region"])
    f2 = res.feasibility["audit_2_composition"]
    assert f2["level"] == "no"
    assert any("INDETERMINATE" in n for n in f2["notes"])
    assert any("数据生态" in n for n in f2["notes"])


def test_check_reports_blocker_instead_of_raising():
    """obs 解析失败必须作为 blocker 出现在报告里（而不是抛异常）。"""
    ad = make_adata()
    ad.obs = ad.obs.drop(columns=["condition"])
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg())
    assert res.blockers, "缺少 condition 列应产生 blocker"
    assert res.obs_ok is False
    assert all(v["level"] == "no" for v in res.feasibility.values())


def test_check_never_writes_config():
    """核心契约：check 不修改 cfg（不采纳任何变量）。"""
    ad = make_adata()
    cfg = _cfg()
    before = dict(cfg.to_dict())
    check_readiness(ad, "condition", "donor", "control", cfg=cfg)
    assert cfg.to_dict() == before, "check 修改了配置——违反'只报告不采纳'契约"


def test_check_depth_bin_capacity_counts_both_groups():
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=40)
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg())
    db = res.depth_bins_ok
    assert db["n_bins"] == len(AuditConfig().depth_edges) - 1
    assert sum(b["n_ctrl"] for b in db["bins"]) == int(
        (ad.obs["condition"] == "control").sum())
    assert db["n_bins_ok"] >= 1


def test_check_warns_on_multilevel_condition():
    ad = make_adata()
    idx = ad.obs.index[ad.obs["donor"] == "case_d0"]
    ad.obs.loc[idx, "condition"] = "DCM"
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg())
    assert res.case_scope == "all non-control"
    assert any("合并为 case" in w for w in res.warnings)


def test_check_warns_on_few_donors():
    ad = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=20)
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg())
    assert any("置换检验功效不足" in w for w in res.warnings)


# ════════════════ CLI 层 ════════════════

def test_check_cli_writes_both_reports(tmp_path):
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20)
    h5 = tmp_path / "a.h5ad"
    ad.write_h5ad(h5)
    out = tmp_path / "o"
    rc = cli_main(["check", "--input", str(h5), "--outdir", str(out),
                   "--quiet", "--no-count-check"])
    assert rc == 0
    assert (out / "check_report.txt").exists()
    j = json.loads((out / "check_report.json").read_text(encoding="utf-8"))
    assert j["schema"] == "sc-audit/readiness/1"
    assert j["groups"]["obs_ok"] is True
    assert "only REPORTS" in j["note"]


def test_check_cli_exit_code_1_when_audit2_infeasible(tmp_path):
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20)
    for c in ("institution", "donation_type", "source", "region"):
        if c in ad.obs.columns:
            ad.obs[c] = ["v_" + str(d) for d in ad.obs["donor"].astype(str)]
    h5 = tmp_path / "b.h5ad"
    ad.write_h5ad(h5)
    out = tmp_path / "o"
    rc = cli_main(["check", "--input", str(h5), "--outdir", str(out), "--quiet",
                   "--no-count-check", "--design-var-candidates",
                   "institution,donation_type,source,region"])
    assert rc == 1, "审计二不可行时应返回 1"


def test_check_cli_exit_code_2_on_blocker(tmp_path):
    ad = make_adata()
    ad.obs = ad.obs.drop(columns=["donor"])
    h5 = tmp_path / "c.h5ad"
    ad.write_h5ad(h5)
    rc = cli_main(["check", "--input", str(h5), "--outdir", str(tmp_path / "o"),
                   "--quiet", "--no-count-check"])
    assert rc == 2, "缺少 donor 列应返回 2（阻断）"


def test_check_cli_records_gating_when_preset_given(tmp_path):
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20)
    h5 = tmp_path / "d.h5ad"
    ad.write_h5ad(h5)
    out = tmp_path / "o"
    loose = tmp_path / "loose.yaml"
    loose.write_text("min_counts: 1\nmin_genes: 3\nmin_pos_markers: 1\n"
                     "pos_markers: [TTN]\nneg_markers: [PTPRC]\n", encoding="utf-8")
    cli_main(["check", "--input", str(h5), "--preset", "cardiac",
              "--config", str(loose), "--outdir", str(out), "--quiet",
              "--no-count-check"])
    j = json.loads((out / "check_report.json").read_text(encoding="utf-8"))
    g = j["input"]["gating"]
    assert g is not None and g["preset"] == "cardiac"
    assert "n_pass" in g and "missing_pos_markers" in g


# ════════════════ 变量角色与建议 ════════════════

def test_check_classifies_roles_from_column_names():
    from sc_audit.check import _name_role
    assert _name_role("institution") == "design"
    assert _name_role("donation_type") == "design"
    assert _name_role("source") == "design"
    assert _name_role("biopsy_site") == "design"
    assert _name_role("outcome_event") == "outcome"
    assert _name_role("mortality_30d") == "outcome"
    assert _name_role("treatment_flag") == "intervention"
    assert _name_role("treatment_combination") == "intervention"
    assert _name_role("Channel") == "technical"
    assert _name_role("seq_batch") == "technical"
    assert _name_role("age_range") == "confounder"
    assert _name_role("gender") == "confounder"
    assert _name_role("mystery_col") == "unknown"


def test_condition_dependent_does_not_hide_design_variable():
    """回归测试（核心修复）。

    "只在单一组内变化"同时匹配【结局变量】与【分组不平衡的取材变量】，
    二者在数据上无法区分。第一版据此把 institution 判为"疑似结局"并移出候选，
    等于把最该被匹配的变量藏起来 —— 假阳性比假阴性更糟。
    现在：词根定角色；统计信号只产出 condition_dependent 提示。
    """
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20)
    # 对照跨 2 个机构，病例全在 1 个机构 → cv=2, kv=1，TVD 恰为 0.5
    ctrl_d = [f"ctrl_d{i}" for i in range(4)]
    case_d = [f"case_d{i}" for i in range(4)]
    for i, d in enumerate(ctrl_d):
        ad.obs.loc[ad.obs["donor"] == d, "institution"] = "INST_A" if i < 2 else "INST_B"
    for d in case_d:
        ad.obs.loc[ad.obs["donor"] == d, "institution"] = "INST_A"
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg(),
                          candidate_vars=["institution"])
    cols = [d["column"] for d in res.design_vars]
    assert "institution" in cols, "取材变量不得因'条件依赖'被移出候选"
    inst = {d["column"]: d for d in res.design_vars}["institution"]
    assert inst["condition_dependent"] is True
    assert inst["recommendation"] == "caution_condition_dependent"
    assert "人工确认" in inst["reasoning"]
    assert not res.suspect_vars, "design 角色不应进入疑似结局区"


def test_check_separates_outcome_vars_but_keeps_them_visible():
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20)
    ad.obs["outcome_event"] = ["no"] * ad.n_obs
    idx = ad.obs.index[ad.obs["donor"].isin(["case_d0", "case_d1"])]
    ad.obs.loc[idx, "outcome_event"] = "yes"          # 结局变量：仅病例组有 yes
    ad.obs["treatment_flag"] = ["no"] * ad.n_obs
    ad.obs.loc[ad.obs.index[ad.obs["donor"] == "case_d2"], "treatment_flag"] = "yes"
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg(),
                          candidate_vars=["outcome_event", "treatment_flag"])
    sus = {d["column"]: d for d in res.suspect_vars}
    assert "outcome_event" in sus and sus["outcome_event"]["recommendation"] == "never"
    assert "collider" in sus["outcome_event"]["reasoning"]
    assert "treatment_flag" in sus and sus["treatment_flag"]["recommendation"] == "avoid"
    # 未隐藏：仍出现在报告里（suspect_vars 也是输出内容）
    j = res.to_json()
    assert len(j["suspect_design_variables"]) == 2


def test_check_ranks_primary_match_before_risky():
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=20)
    # source 与分组轻度相关（0.1~0.5 带内），institution 高度共线
    for i, d in enumerate([f"ctrl_d{j}" for j in range(5)]):
        ad.obs.loc[ad.obs["donor"] == d, "source"] = "SRC_A" if i < 3 else "SRC_B"
    for d in [f"case_d{j}" for j in range(5)]:
        ad.obs.loc[ad.obs["donor"] == d, "source"] = "SRC_A" if d != "case_d4" else "SRC_B"
    for i, d in enumerate([f"ctrl_d{j}" for j in range(5)]):
        ad.obs.loc[ad.obs["donor"] == d, "institution"] = "INST_B"
    for d in [f"case_d{j}" for j in range(5)]:
        ad.obs.loc[ad.obs["donor"] == d, "institution"] = "INST_A"
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg())
    order = [d["column"] for d in res.design_vars]
    assert "source" in order and "institution" in order
    assert order.index("source") < order.index("institution"), \
        f"推荐变量应排在 risky 之前: {order}"
    rec = {d["column"]: d["recommendation"] for d in res.design_vars}
    assert rec["source"] == "primary_match"
    assert rec["institution"] == "risky"


def test_audit2_recommendation_matched_comparison():
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=20)
    for i, d in enumerate([f"ctrl_d{j}" for j in range(5)]):
        ad.obs.loc[ad.obs["donor"] == d, "source"] = "SRC_A" if i < 3 else "SRC_B"
    for d in [f"case_d{j}" for j in range(5)]:
        ad.obs.loc[ad.obs["donor"] == d, "source"] = "SRC_A" if d != "case_d4" else "SRC_B"
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg())
    a2 = res.audit_2_recommendation
    assert a2["method"] == "matched_comparison"
    assert "source" in a2["recommended_variables"]


def test_audit2_recommendation_tvd_stratification_when_all_confounded():
    ad = make_adata(design_confounded=True, n_ctrl_donors=5, n_case_donors=5,
                    cells_per_donor=20)
    for c in ("institution", "source", "donation_type"):
        ad.obs[c] = "X"
    for d in [f"ctrl_d{j}" for j in range(5)]:
        ad.obs.loc[ad.obs["donor"] == d, "institution"] = "INST_B"
    for d in [f"case_d{j}" for j in range(5)]:
        ad.obs.loc[ad.obs["donor"] == d, "institution"] = "INST_A"
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg(),
                          candidate_vars=["institution"])
    a2 = res.audit_2_recommendation
    assert a2["method"] == "tvd_stratification"
    assert a2["feasible"] is False
    assert "total-variation distance" in a2["statement_for_paper"]


def test_audit2_recommendation_no_matching_needed_when_balanced():
    ad = make_adata(design_confounded=False, n_ctrl_donors=5, n_case_donors=5,
                    cells_per_donor=20)
    for d in ad.obs["donor"].astype(str).unique():
        ad.obs.loc[ad.obs["donor"].astype(str) == d, "source"] = "SRC_A"  # 常量
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg(),
                          candidate_vars=["source"])
    assert res.audit_2_recommendation["method"] == "no_matching_needed"
    assert res.feasibility["audit_2_composition"]["level"] == "ok"


def test_reference_coverage_uses_target_panel_as_denominator():
    from sc_audit.check import check_readiness as cr
    ad = make_adata()
    res = cr(ad, "condition", "donor", "control", cfg=_cfg(),
             reference={"TNNT2": 90.0}, target_genes=["TNNT2", "PLN", "MYL3"])
    rc = res.reference_coverage
    assert rc["n_target_genes"] == 3
    assert rc["target_covered"] == ["TNNT2"]
    assert set(rc["target_missing"]) == {"PLN", "MYL3"}
    assert rc["pct_target_covered"] == 33.3
    assert any("参照表只覆盖 1/3" in w for w in res.warnings)


# ════════════════ 显式声明与运行期保护 ════════════════

def test_role_overrides_take_priority_over_heuristics():
    """显式声明优先于词根：启发式认不出的结局变量由使用者标注。"""
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20)
    ad.obs["complication_30d"] = ["no"] * ad.n_obs
    ad.obs.loc[ad.obs.index[ad.obs["donor"] == "case_d0"], "complication_30d"] = "yes"
    ad.obs["mystery_x"] = ["a", "b"] * (ad.n_obs // 2) + ["a"] * (ad.n_obs % 2)
    cfg = _cfg(role_overrides={"complication_30d": "outcome",
                               "mystery_x": "ignore"})
    res = check_readiness(ad, "condition", "donor", "control", cfg=cfg,
                          candidate_vars=["complication_30d", "mystery_x"])
    sus = {d["column"]: d for d in res.suspect_vars}
    assert "complication_30d" in sus, "显式声明的结局变量应进入疑似区"
    assert sus["complication_30d"]["recommendation"] == "never"
    assert any("显式声明" in r for r in sus["complication_30d"]["role_reasons"])
    # ignore 应被移出候选
    allcols = [d["column"] for d in res.design_vars] + [d["column"] for d in res.suspect_vars]
    assert "mystery_x" not in allcols, "role_overrides=ignore 的列不应进入候选"


def test_treatment_alias_normalised_to_intervention():
    from sc_audit.check import _name_role
    cfg = _cfg(role_overrides={"some_col": "treatment"})
    assert _name_role("some_col", cfg) == "intervention"


def test_extra_keywords_extend_heuristics():
    """领域词根可扩展，用于启发式认不出的列名。"""
    from sc_audit.check import _name_role
    ad = make_adata()
    cfg = _cfg(extra_outcome_keywords=["troponin_peak"],
               extra_treatment_keywords=["colchicine"])
    assert _name_role("troponin_peak_x", cfg) == "outcome"
    assert _name_role("colchicine_arm", cfg) == "intervention"
    # 未扩展时识别不出
    assert _name_role("troponin_peak_x", _cfg()) == "unknown"


def test_options_block_present_when_risky_var_exists():
    """存在接近分离的变量时，必须给出三套可选策略（且不排序）。"""
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=20)
    for d in [f"ctrl_d{j}" for j in range(5)]:
        ad.obs.loc[ad.obs["donor"] == d, "institution"] = "INST_B"
    for d in [f"case_d{j}" for j in range(5)]:
        ad.obs.loc[ad.obs["donor"] == d, "institution"] = "INST_A"
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg(),
                          candidate_vars=["institution"])
    a2 = res.audit_2_recommendation
    assert "options" in a2 and len(a2["options"]) == 3
    ids = [o["id"] for o in a2["options"]]
    assert ids == ["conservative", "primary_plus_sensitivity", "no_matching"]
    for o in a2["options"]:
        assert o["trade_off"], f"{o['id']} 缺少取舍说明"
    assert "不替你排序" in a2["options_note"]


def test_runtime_guard_blocks_outcome_variable_in_match_vars(tmp_path):
    """运行期保护：match_vars 含结局变量 → 中止（防 collider bias）。"""
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20)
    ad.obs["outcome_event"] = ["no"] * ad.n_obs
    h5 = tmp_path / "g.h5ad"
    ad.write_h5ad(h5)
    cfgf = tmp_path / "c.yaml"
    cfgf.write_text("match_vars: [outcome_event, source]\nmin_counts: 1\nmin_genes: 3\n",
                    encoding="utf-8")
    import pytest as _pytest
    with _pytest.raises(SystemExit) as e:
        cli_main(["composition", "--input", str(h5), "--config", str(cfgf),
                  "--target-genes", "TNNT2", "--outdir", str(tmp_path / "o"),
                  "--quiet", "--no-count-check"])
    assert e.value.code == 2, "应中止"


def test_runtime_guard_can_be_explicitly_overridden(tmp_path):
    """显式放行后可运行（使用者承担判断责任）。"""
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20)
    ad.obs["outcome_event"] = ["no"] * ad.n_obs
    h5 = tmp_path / "h.h5ad"
    ad.write_h5ad(h5)
    cfgf = tmp_path / "c2.yaml"
    cfgf.write_text("match_vars: [source]\nmin_counts: 1\nmin_genes: 3\n",
                    encoding="utf-8")
    rc = cli_main(["composition", "--input", str(h5), "--config", str(cfgf),
                   "--target-genes", "TNNT2", "--allow-outcome-matching",
                   "--outdir", str(tmp_path / "o"), "--quiet", "--no-count-check"])
    assert rc == 0


def test_all_config_hints_are_real_config_fields():
    """持久防线（教训）。

    策略块曾建议 `run_matched_comparison` / `stratify_by` 两个键，
    而它们当时在代码中【根本不存在】——用户照做会得到被静默忽略的配置，
    即典型的"声称 ≠ 实现"。
    本测试把【输出中承诺的配置键】与【AuditConfig 真实字段】绑定：
    今后任何新增建议都必须先有实现。
    """
    from dataclasses import fields as _fields
    from sc_audit.config import AuditConfig
    real = {f.name for f in _fields(AuditConfig)}

    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=20)
    for d in [f"ctrl_d{j}" for j in range(5)]:
        ad.obs.loc[ad.obs["donor"] == d, "institution"] = "INST_B"
    for d in [f"case_d{j}" for j in range(5)]:
        ad.obs.loc[ad.obs["donor"] == d, "institution"] = "INST_A"
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg(),
                          candidate_vars=["institution"])
    opts = res.audit_2_recommendation.get("options") or []
    assert opts, "应存在策略块"
    for o in opts:
        for k in (o.get("config_hint") or {}):
            assert k in real, (
                f"策略 {o['id']} 建议了不存在的配置键 `{k}` —— "
                f"要么实现它，要么从 config_hint 移除（声称≠实现）")


def test_run_matched_comparison_false_disables_matching():
    """策略 3 必须真的能关闭匹配，而不只是文案。"""
    from sc_audit.composition import audit_control_composition
    ad = make_adata(design_confounded=True, n_ctrl_donors=5, n_case_donors=5,
                    cells_per_donor=20)
    res = audit_control_composition(ad, ["TNNT2"], "condition", "donor", "control",
                                    cfg=_cfg(run_matched_comparison=False,
                                             min_cells_per_stratum=1))
    assert res.match_attempted is False
    assert any("DISABLED" in n for n in res.notes), "应记录关闭原因供引用"
    # 对照组：默认 True 时应尝试匹配
    res2 = audit_control_composition(ad, ["TNNT2"], "condition", "donor", "control",
                                     cfg=_cfg(min_cells_per_stratum=1))
    assert res2.match_attempted is True


# ════════════════ 参照表溯源 / 门控参数登记════════════════

def test_reference_provenance_reported_from_manifest():
    """参照表带 manifest 时，必须报告门控未登记的背景并给出警告。"""
    ad = make_adata()
    man = {"backgrounds": {"bg_ok": {"provenance": "REPRODUCIBLE_FROM_SCRIPT_04"},
                           "bg_bad": {"provenance": "NOT_REPRODUCIBLE"}},
           "pct_ref_rule": "max over backgrounds",
           "known_gaps": ["某背景门控未登记"]}
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg(),
                          reference={"TNNT2": 90.0}, target_genes=["TNNT2"],
                          reference_manifest=man)
    rp = res.reference_provenance
    assert rp["manifest_found"] is True
    assert rp["gate_unknown_backgrounds"] == ["bg_bad"]
    assert rp["n_backgrounds"] == 2
    assert any("门控参数未登记" in w for w in res.warnings)


def test_missing_manifest_warns_r_not_verifiable():
    """没有 manifest → 必须明确说 R 的可比性无法核验（不静默通过）。"""
    ad = make_adata()
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg(),
                          reference={"TNNT2": 90.0}, target_genes=["TNNT2"])
    assert res.reference_provenance["manifest_found"] is False
    assert any("可比性无法核验" in w for w in res.warnings)


def test_small_donor_count_warning_states_the_minimum_attainable_p():
    """小样本警示必须给出正确的精确枚举下限（v0.3.3）。

    精确枚举下最小可达**双侧** P = 2/C(n,k)（两组等大）/ 1/C(n,k)（不等大，且要求
    极值划分唯一）。旧文案只写 1/C(n,k) 并断言"永远达不到显著"——对等大组且阳性
    细胞全部分离的情形是错的（实测 3v4 与 4v4 的下限都 < 0.05）。
    """
    ad = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=10)
    res = check_readiness(ad, "condition", "donor", "control", cfg=_cfg(),
                          target_genes=["TNNT2"])
    w = " ".join(res.warnings)
    assert "2/C(n,k)" in w and "1/C(n,k)" in w
    assert "永远达不到显著" not in w


def test_check_applies_count_layer_override(tmp_path):
    """`check --count-layer` 必须真的生效（v0.3.3）。

    历史缺陷：`_run_check` 只覆盖了 `n_perm`，`count_layer` / `ratio_basis` /
    `hetero_col` / `match_vars` 都被静默忽略 → 预检比正式运行**更严**，
    会把"计数在某个非标准 layer 里"误报成阻断项（而 `floor --count-layer` 能跑通）。
    """
    import json as _json
    import numpy as _np
    import scipy.sparse as sp
    ad = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=20)
    # 整数 counts 放进**非标准名**的 layer，再把 X 换成非整数 → 必须靠 --count-layer 才能通过
    ad.layers["mycounts"] = ad.X.copy()
    ad.X = sp.csr_matrix(_np.log1p(_np.asarray(ad.X.todense())))
    h5 = tmp_path / "d.h5ad"
    ad.write_h5ad(h5)

    rc1 = cli_main(["check", "--input", str(h5), "--target-genes", "TNNT2",
                    "--outdir", str(tmp_path / "c1"), "--quiet"])
    j1 = _json.loads((tmp_path / "c1" / "check_report.json").read_text())
    assert rc1 == 2 and j1.get("blockers"), "未指定 layer 时应报阻断（X 非计数）"

    rc2 = cli_main(["check", "--input", str(h5), "--target-genes", "TNNT2",
                    "--count-layer", "mycounts", "--outdir", str(tmp_path / "c2"),
                    "--quiet"])
    j2 = _json.loads((tmp_path / "c2" / "check_report.json").read_text())
    assert not j2.get("blockers"), (
        f"--count-layer 未生效（仍报阻断）：{j2.get('blockers')}")
    assert rc2 in (0, 1)
