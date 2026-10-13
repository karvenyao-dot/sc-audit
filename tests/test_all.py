"""sc-audit 测试套件。

修订（2026-10-01，据外部逻辑审计）：
  · 新增 B1/B2/B3/B4/B5 用例（原完全未覆盖）
  · 新增 S1（平衡多中心设计）、S3（k=0 不得为 ok）、S6（分子口径）用例
  · 新增 FDR、精确枚举、无参照路径、归一化 X 用例
  · 修正原 test_floor_verdict_below_5pct 的抽样脆弱性（改用更大 n）
"""
import numpy as np
import pandas as pd
import pytest

from sc_audit.config import AuditConfig, load_config
from sc_audit.gating import gate_cells, CARDIAC_PRESETS, preset_from_config
from sc_audit.io import (ensure_obs_fields, resolve_counts, compute_qc_metrics,
                         resolve_gene_names)
from sc_audit.floor import audit_detection_floor, perm_test, _bh_fdr
from sc_audit.composition import audit_control_composition, _group_distribution_distance
from sc_audit.robust import audit_robust_signal
from sc_audit.reference import compute_detection_ratio, build_reference_from_datasets

from conftest import make_adata, relabel


def _cfg(**kw):
    # 注意：**不覆写** min_cells_per_stratum / strat_ci_frac_ok，
    # 以便 S3c 的默认值（30 / 0.6）真正被集成测试执行到
    d = dict(n_perm=200, event_sparse=20, event_limited=50,
             exact_enum_max_units=40)
    d.update(kw)
    return AuditConfig(**d)


# ════════════════ 配置 ════════════════

def test_default_config_matches_paper_thresholds():
    c = AuditConfig()
    assert c.floor_pct == 5.0 and c.ratio_min == 0.3
    assert c.ratio_ok == 0.6 and c.pct_ok == 20.0
    assert c.event_sparse == 20 and c.event_limited == 50
    assert c.min_pos_markers == 2
    # 新增字段（原为硬编码）
    assert c.hetero_max == 0.5 and c.attenuation_min == 0.5
    assert c.depth_delta_max_pct == 20.0 and c.strat_min_delta_pp == 5.0
    assert c.fdr_method == "bh"


def test_load_preset_cardiac():
    c = load_config("cardiac")
    assert c.floor_pct == 5.0 and c.event_sparse == 20
    assert "TTN" in c.pos_markers and "PTPRC" in c.neg_markers
    assert getattr(c, "_source_path", "").endswith("cardiac.yaml")


def test_load_config_missing_raises():
    with pytest.raises(FileNotFoundError):
        load_config("no_such_preset_xyz")


def test_preset_from_config_is_single_source_of_truth():
    """门控阈值必须来自 AuditConfig，而非预设里另存一份。"""
    cfg = AuditConfig(min_counts=777, min_genes=321, min_pos_markers=3)
    pr = preset_from_config(cfg, "cardiac")
    assert pr["min_counts"] == 777 and pr["min_genes"] == 321
    assert pr["min_pos_markers"] == 3


# ════════════════ B1 · control_label 不得被覆盖 ════════════════

def test_explicit_control_label_is_respected():
    ad = make_adata()
    info = ensure_obs_fields(ad, "condition", "donor", control_label="disease")
    assert info["control_label"] == "disease"
    assert info["case_label"] == "control"
    assert info["n_control_cells"] == int((ad.obs["condition"] == "disease").sum())


def test_explicit_control_label_with_uncommon_names_does_not_raise():
    ad = make_adata()
    ad.obs["condition"] = ad.obs["condition"].map({"control": "NF", "disease": "DCM"})
    info = ensure_obs_fields(ad, "condition", "donor", control_label="NF")
    assert info["control_label"] == "NF" and info["case_label"] == "DCM"


def test_explicit_control_label_not_in_data_raises():
    ad = make_adata()
    with pytest.raises(ValueError, match="不在"):
        ensure_obs_fields(ad, "condition", "donor", control_label="healthy_old")


def test_multilevel_condition_warns_and_marks_case_as_all_noncontrol():
    ad = make_adata()
    idx = ad.obs.index[ad.obs["donor"] == "case_d0"]
    ad.obs.loc[idx, "condition"] = "DCM"
    with pytest.warns(UserWarning, match="合并为 case"):
        info = ensure_obs_fields(ad, "condition", "donor", control_label="control")
    assert info["case_label"] == "<all non-control>"
    assert info["n_condition_levels"] == 3


def test_donor_spanning_groups_raises():
    ad = make_adata()
    idx = ad.obs.index[ad.obs["donor"] == "ctrl_d0"][:5]
    ad.obs.loc[idx, "condition"] = "disease"
    with pytest.raises(ValueError, match="多个"):
        ensure_obs_fields(ad, "condition", "donor", control_label="control")


# ════════════════ B2 · batch 列不得破坏匹配 ════════════════

def test_sample_unique_batch_is_dropped_from_matching():
    """batch 为样本级唯一时，必须被剔除出匹配变量，否则匹配恒失败。"""
    ad = make_adata(design_confounded=False)     # 干净设计
    # N2：batch 是**供者级**唯一（每个供者一个取值），而不是细胞级唯一
    assert ad.obs["batch"].nunique() == ad.obs["donor"].nunique()
    assert ad.obs["batch"].nunique() < ad.n_obs

    # ① 默认：batch 不在 match_vars 候选中
    res = audit_control_composition(ad, ["TNNT2"], "condition", "donor", "control",
                                    cfg=_cfg())
    assert "batch" not in res.match_vars_used
    assert res.matched_cells > 0

    # ② 显式把 batch 放入 match_vars → 自动剔除（B2 的真正防护）
    cfg = _cfg(match_vars=["institution", "batch"], auto_filter_stratum_cols=True)
    res2 = audit_control_composition(ad, ["TNNT2"], "condition", "donor", "control",
                                     cfg=cfg)
    assert "batch" in res2.match_vars_dropped, res2.match_vars_dropped
    assert "batch" not in res2.match_vars_used
    assert res2.matched_cells > 0, "batch 未被剔除时匹配会恒为空"

    # ③ 关闭自动剔除 → 匹配必然失败（反证该防护是必要的）
    cfg3 = _cfg(match_vars=["institution", "batch"], auto_filter_stratum_cols=False)
    res3 = audit_control_composition(ad, ["TNNT2"], "condition", "donor", "control",
                                     cfg=cfg3)
    assert res3.matched_cells == 0, "样本级唯一列存在时匹配应当失败"


# ════════════════ N1 · 基因名回退必须接入审计 ════════════════

def test_gene_fallback_reaches_audit_end_to_end(tmp_path):
    """N1：var_names 为 ENSG 时，回退后审计必须真正检出基因（而非仅返回映射）。"""
    from sc_audit.cli import main as cli_main
    ad = make_adata(gene_pcts={"TNNT2": (0.9, 0.4)})
    symbols = list(ad.var_names)
    ad.var["gene_symbols"] = symbols
    ad.var_names = [f"ENSG{i:08d}" for i in range(ad.n_vars)]
    h5 = tmp_path / "ensg.h5ad"
    ad.write_h5ad(h5)

    out = tmp_path / "out"
    rc = cli_main(["floor", "--input", str(h5), "--target-genes", "TNNT2,ACTB",
                   "--condition-col", "condition", "--donor-col", "donor",
                   "--n-perm", "100", "--outdir", str(out), "--quiet",
                   "--no-count-check"])
    assert rc == 0
    df = pd.read_csv(out / "floor_audit.csv")
    assert set(df["gene"]) == {"TNNT2", "ACTB"}
    assert not (df["verdict"] == "NOT_IN_DATASET").any(), \
        f"基因名回退未接入审计: {df[['gene','verdict']].to_dict('records')}"


# ════════════════ N3 · yaml 与默认值必须一致 ════════════════

def test_yaml_preset_matches_defaults_for_locked_thresholds():
    """N3：cardiac.yaml 不得静默回退任何锁定阈值。

    v0.3.3：从“抽指定字段”改为**全字段**逐值比对（此前 52/66，靠注释声称
    “逐字段一致”是错宣称）。
    """
    y = load_config("cardiac").to_dict()
    d = AuditConfig().to_dict()
    assert set(y) == set(d), (f"cardiac.yaml 与默认值字段集不一致："
                              f"缺 {sorted(set(d) - set(y))}，多 {sorted(set(y) - set(d))}")
    for k in sorted(d):
        assert y[k] == d[k], f"{k}: yaml={y[k]!r} 默认={d[k]!r}"


# ════════════════ CLI 冒烟测试 ════════════════

def test_cli_smoke_run_generates_all_outputs(tmp_path):
    """CLI 层布线测试（B1 / N1 / N5 都是 CLI 层 bug，必须有冒烟覆盖）。"""
    from sc_audit.cli import main as cli_main
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=25)
    h5 = tmp_path / "smoke.h5ad"
    ad.write_h5ad(h5)
    out = tmp_path / "out"
    rc = cli_main(["run", "--input", str(h5),
                   "--target-genes", "TNNT2,ACTB,GAPDH",
                   "--condition-col", "condition", "--donor-col", "donor",
                   "--n-perm", "100", "--outdir", str(out), "--quiet",
                   "--no-count-check"])
    assert rc == 0
    for f in ("report.html", "floor_audit.csv", "composition.txt",
              "robustness.txt", "run_metadata.json"):
        assert (out / f).exists(), f"缺少产物 {f}"
    import json
    meta = json.loads((out / "run_metadata.json").read_text())
    assert meta["control_label"] == "control"
    assert "case_label" in meta and "case_scope" in meta
    assert "counts_source" in meta


def test_cli_cell_subset_failure_aborts(tmp_path):
    """N5：--cell-subset 表达式失败时必须中止，不得在全量数据上出结论。"""
    from sc_audit.cli import main as cli_main
    ad = make_adata()
    h5 = tmp_path / "s.h5ad"
    ad.write_h5ad(h5)
    with pytest.raises(SystemExit) as e:
        cli_main(["floor", "--input", str(h5), "--target-genes", "TNNT2",
                  "--cell-subset", "no_such_column == 1",
                  "--outdir", str(tmp_path / "o"), "--quiet", "--no-count-check"])
    assert e.value.code == 2


def test_cli_explicit_control_label_respected(tmp_path):
    """B1 的 CLI 级回归：显式 --control-label 必须生效（原 bug 正在此层）。"""
    import json
    from sc_audit.cli import main as cli_main
    ad = make_adata()
    h5 = tmp_path / "b1.h5ad"
    ad.write_h5ad(h5)
    out = tmp_path / "o"
    cli_main(["run", "--input", str(h5), "--target-genes", "TNNT2",
              "--control-label", "disease",
              "--n-perm", "100", "--outdir", str(out), "--quiet",
              "--no-count-check"])
    meta = json.loads((out / "run_metadata.json").read_text())
    assert meta["control_label"] == "disease"
    assert meta["case_label"] == "control"


# ════════════════ B5 · 原始计数校验 ════════════════

def test_noninteger_x_is_rejected_without_counts_layer():
    ad = make_adata(scale_x=True)
    with pytest.raises(ValueError, match="原始计数"):
        resolve_counts(ad, AuditConfig())


def test_noninteger_x_falls_back_to_counts_layer():
    import scipy.sparse as sp
    ad = make_adata(scale_x=True)
    raw = (ad.X > 0).astype(float) * 3          # 构造整数 counts 层
    ad.layers["counts"] = sp.csr_matrix(raw)
    # X 仍是 log-normalized：整体不满足整数 → 但 counts 层满足
    ad2, src = resolve_counts(ad, AuditConfig(count_layer="counts"))
    assert "counts" in src
    assert np.allclose(ad2.X.toarray() % 1, 0, atol=1e-6)


# ════════════════ 门控 ════════════════

def test_gating_passes_expected_cells():
    ad = make_adata()
    mask = gate_cells(ad, CARDIAC_PRESETS["cardiac"]["pos_markers"],
                      CARDIAC_PRESETS["cardiac"]["neg_markers"],
                      min_pos_markers=2, min_counts=0, min_genes=0)
    assert mask.sum() > 0.8 * ad.n_obs
    assert ad.uns["sc_audit_gating"]["n_pass"] == int(mask.sum())


def test_gating_rejects_marker_negative_cells():
    import scipy.sparse as sp
    ad = make_adata()
    X = ad.X.tolil()
    for r in range(10):
        X[r, ad.var_names.get_loc("PTPRC")] = 5
    ad.X = sp.csr_matrix(X)
    mask = gate_cells(ad, CARDIAC_PRESETS["cardiac"]["pos_markers"],
                      CARDIAC_PRESETS["cardiac"]["neg_markers"],
                      min_pos_markers=2, min_counts=0, min_genes=0)
    assert not mask[:10].any()


def test_gating_counts_missing_markers():
    ad = make_adata()
    gate_cells(ad, ["TTN", "NOT_A_GENE"], ["PTPRC"],
               min_pos_markers=1, min_counts=0, min_genes=0)
    assert "NOT_A_GENE" in ad.uns["sc_audit_gating"]["missing_pos_markers"]


# ════════════════ 审计一 · 判定与边界 ════════════════

def test_floor_verdict_below_5pct_is_unquantifiable():
    """检出率 <5% → UNQUANTIFIABLE（用大 n 降低抽样脆弱性）。"""
    ad = make_adata(n_ctrl_donors=20, n_case_donors=20, cells_per_donor=50,
                    gene_pcts={"BELOW": (0.03, 0.03), "ABOVE": (0.30, 0.30)})
    res = audit_detection_floor(ad, ["BELOW", "ABOVE"], "condition", "donor",
                                "control", cfg=_cfg(), reference={"ABOVE": 30.0})
    v = dict(zip(res.table["gene"], res.table["verdict"]))
    assert "UNQUANTIFIABLE" in v["BELOW"]
    assert "UNQUANTIFIABLE" not in v["ABOVE"]


def test_floor_no_reference_is_insufficient_not_quantifiable():
    """S6 修复：无参照时不得直接判 QUANTIFIABLE。"""
    ad = make_adata(n_ctrl_donors=10, n_case_donors=10, cells_per_donor=20,
                    gene_pcts={"G": (0.6, 0.6)})
    res = audit_detection_floor(ad, ["G"], "condition", "donor", "control", cfg=_cfg())
    row = res.table.iloc[0]
    assert np.isnan(row["ratio_R"])
    assert row["verdict"] == "INSUFFICIENT (no reference)"


def test_floor_verdict_uses_same_basis_as_ratio():
    """S6 + N4：verdict 与 R 共用同一个 ratio_basis 字段。"""
    # 病例组 0% / 对照组 40%：case 口径下 basis_pct = 0 → UNQUANTIFIABLE
    ad = make_adata(n_ctrl_donors=8, n_case_donors=8, cells_per_donor=25,
                    gene_pcts={"G": (0.40, 0.0)})
    res = audit_detection_floor(ad, ["G"], "condition", "donor", "control",
                                cfg=_cfg(ratio_basis="case"),
                                reference={"G": 40.0})
    row = res.table.iloc[0]
    assert row["case_pct"] == 0.0 and row["pct_study"] > 5.0
    assert "UNQUANTIFIABLE" in row["verdict"], row["verdict"]
    # 换 overall 口径 → 同一基因不再是 UNQUANTIFIABLE（证明口径确实生效）
    res2 = audit_detection_floor(ad, ["G"], "condition", "donor", "control",
                                 cfg=_cfg(ratio_basis="overall"),
                                 reference={"G": 40.0})
    assert res2.table.iloc[0]["basis_pct"] != row["basis_pct"]


def test_floor_event_strata_boundaries():
    cfg = _cfg(event_sparse=20, event_limited=50)
    ad = make_adata(n_ctrl_donors=10, n_case_donors=10, cells_per_donor=20,
                    gene_pcts={"G_LOW": (0.03, 0.03), "G_MID": (0.08, 0.08),
                               "G_HIGH": (0.60, 0.60)})
    res = audit_detection_floor(ad, ["G_LOW", "G_MID", "G_HIGH"], "condition",
                                "donor", "control", cfg=cfg)
    t = res.table.set_index("gene")
    assert t.loc["G_LOW", "stratum"] == "event_sparse"
    assert t.loc["G_MID", "stratum"] == "event_limited"
    assert t.loc["G_HIGH", "stratum"] == "event_sufficient"


def test_floor_event_sparse_gets_no_significance_test():
    ad = make_adata(gene_pcts={"RARE": (0.02, 0.02)})
    res = audit_detection_floor(ad, ["RARE"], "condition", "donor", "control", cfg=_cfg())
    assert np.isnan(res.table.iloc[0]["perm_P"])


def test_floor_binom_only_when_case_zero():
    """S4：二项检验仅在 case_pct == 0 时报告。"""
    ad = make_adata(n_ctrl_donors=10, n_case_donors=10, cells_per_donor=20,
                    gene_pcts={"ZERO": (0.30, 0.0), "NONZERO": (0.30, 0.30)})
    res = audit_detection_floor(ad, ["ZERO", "NONZERO"], "condition", "donor",
                                "control", cfg=_cfg())
    t = res.table.set_index("gene")
    assert not np.isnan(t.loc["ZERO", "binom_P_case0"])
    assert np.isnan(t.loc["NONZERO", "binom_P_case0"])


def test_floor_s7_definitions_are_consistent():
    """S7：unquantifiable property 与 summary 同定义。"""
    ad = make_adata(gene_pcts={"A": (0.01, 0.01), "B": (0.5, 0.5)})
    res = audit_detection_floor(ad, ["A", "B"], "condition", "donor", "control", cfg=_cfg())
    assert len(res.unquantifiable) == res.summary["n_unquantifiable"]


def test_floor_fdr_column_present_and_valid():
    """4.9：BH-FDR 必须输出 q 值，且 q >= p。"""
    ad = make_adata(n_ctrl_donors=10, n_case_donors=10, cells_per_donor=20,
                    gene_pcts={"G1": (0.5, 0.2), "G2": (0.5, 0.5), "G3": (0.6, 0.3)})
    res = audit_detection_floor(ad, ["G1", "G2", "G3"], "condition", "donor",
                                "control", cfg=_cfg(fdr_method="bh"))
    assert "q_value_BH" in res.table.columns
    sub = res.table.dropna(subset=["perm_P", "q_value_BH"])
    assert len(sub) > 0
    assert (sub["q_value_BH"] >= sub["perm_P"] - 1e-9).all()


def test_floor_exact_enumeration_is_deterministic():
    """4.3：供者数 <= 阈值时用精确枚举 → 同一数据不因种子而变。"""
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=20,
                    gene_pcts={"G": (0.5, 0.2)})
    r1 = perm_test(*_vec(ad, "G"), n_perm=200, seed=1, exact_max_units=40)
    r2 = perm_test(*_vec(ad, "G"), n_perm=200, seed=999, exact_max_units=40)
    assert r1["mode"] == "exact"
    assert r1["p"] == r2["p"]


def _vec(ad, gene):
    from sc_audit.io import gene_detection_vector
    v = gene_detection_vector(ad, gene).astype(float)
    cond = ad.obs["condition"].astype(str).values
    donor = ad.obs["donor"].astype(str).values
    return v, donor, cond == "control"


def test_bh_fdr_matches_manual():
    p = [0.01, 0.02, 0.03, 0.04, 0.05]
    q = _bh_fdr(p)
    assert abs(q[0] - 0.05) < 1e-9        # 0.01*5/1
    assert abs(q[-1] - 0.05) < 1e-9       # monotone capped


def test_floor_jackknife_covers_both_sides():
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=25,
                    gene_pcts={"G": (0.6, 0.15)})
    res = audit_detection_floor(ad, ["G"], "condition", "donor", "control", cfg=_cfg())
    jk = res.jackknife["G"]
    assert jk["sides"] == "both"
    assert jk["n_ctrl_donors"] == 6 and jk["n_case_donors"] == 6
    assert len(jk["table"]) == 12


def test_floor_missing_gene_reported_not_crashed():
    ad = make_adata()
    res = audit_detection_floor(ad, ["NOT_THERE"], "condition", "donor", "control",
                                cfg=_cfg())
    assert res.table.iloc[0]["verdict"] == "NOT_IN_DATASET"


# ════════════════ S1 · F_hetero → 分布距离 ════════════════

def test_balanced_multicenter_design_is_comparable():
    """三中心完全平衡：旧 F_hetero≈0.67 误判 NOT SEPARABLE；TVD=0 → COMPARABLE。"""
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, design_confounded=False)
    insts = ["A", "B", "C"] * 2
    ad = relabel(ad, "institution", {f"ctrl_d{i}": insts[i] for i in range(6)})
    ad = relabel(ad, "institution", {f"case_d{i}": insts[i] for i in range(6)})
    res = audit_control_composition(ad, ["TNNT2"], "condition", "donor", "control",
                                    cfg=_cfg())
    assert res.f_hetero == 0.0, res.f_hetero
    assert "COMPARABLE" in res.verdict and "NOT SEPARABLE" not in res.verdict


def test_fully_collinear_design_has_tvd_one():
    ad = make_adata(design_confounded=False)
    ad = relabel(ad, "institution", {f"ctrl_d{i}": "INST_B" for i in range(5)})
    res = audit_control_composition(ad, ["TNNT2"], "condition", "donor", "control",
                                    cfg=_cfg())
    assert res.f_hetero == 1.0
    assert "NOT SEPARABLE" in res.verdict


def test_distribution_distance_math():
    v = np.array(["A", "A", "B", "B"] + ["A", "A", "B", "B"])
    is_ctrl = np.array([True] * 4 + [False] * 4)
    tvd, cv, mode, tbl = _group_distribution_distance(v, is_ctrl)
    assert tvd == 0.0 and cv == 0.0
    v2 = np.array(["A"] * 4 + ["B"] * 4)
    tvd2, _, _, _ = _group_distribution_distance(v2, is_ctrl)
    assert tvd2 == 1.0


# ════════════════ 审计二 · 其余 ════════════════

def test_composition_detects_confounding(adata_confounded):
    res = audit_control_composition(adata_confounded, ["TNNT2"], "condition",
                                    "donor", "control", cfg=_cfg())
    assert "NOT SEPARABLE" in res.verdict
    assert res.f_hetero is not None and res.f_hetero > 0.5
    assert len(res.crosstabs) > 0 and len(res.crosstabs["institution"]) == 2


def test_composition_crosstabs_not_empty(adata_confounded):
    res = audit_control_composition(adata_confounded, ["TNNT2"], "condition",
                                    "donor", "control", cfg=_cfg())
    for k, ct in res.crosstabs.items():
        assert ct.size > 0, f"交叉表 {k} 为空"


def test_composition_clean_design_matches(adata_clean):
    res = audit_control_composition(adata_clean, ["TNNT2"], "condition", "donor",
                                    "control", cfg=_cfg())
    assert res.matched_cells > 0
    assert res.f_hetero == 0.0


def test_composition_attenuation_excludes_small_effects():
    """S2：full effect < 阈值 的基因不进衰减比汇总（验证排除行为，而非仅列存在）。"""
    # TNNT2 设 (0.9, 0.4) → full effect ≈ 50pp，确保它【保留】在汇总里；
    # 否则默认 (0.5, 0.5) 也可能 <5pp，导致 usable 为空而测试照过（P3）
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, design_confounded=True,
                    gene_pcts={"TNNT2": (0.9, 0.4)})
    # FLAT 两组检出率相同 → full effect ≈ 0 → 必被排除
    import scipy.sparse as sp
    X = ad.X.tolil()
    for g in ("MYH7",):
        X[:, ad.var_names.get_loc(g)] = 0.4        # 恒定 → 效应量 0
    ad.X = sp.csr_matrix(X)
    res = audit_control_composition(ad, ["TNNT2", "MYH7"], "condition", "donor",
                                    "control",
                                    cfg=_cfg(attenuation_min_full_effect_pp=5.0))
    assert res.attenuation is not None and len(res.attenuation) == 2
    assert "direction_preserved" in res.attenuation.columns
    excl = [n for n in res.notes if "excluded from the attenuation summary" in n]
    assert excl, "未报告被排除的基因"
    usable = res.attenuation.loc[
        res.attenuation["full_effect_pp"] >= 5.0, "attenuation_ratio"].dropna()
    assert len(usable) < len(res.attenuation), "至少应有一个基因被排除"
    assert len(usable) >= 1, "应至少有一个基因【保留】在衰减比汇总中"


# ════════════════ S3 · 审计三判定 ════════════════

def test_robust_zero_persisting_is_not_ok():
    """S3a：没有任何基因满足"分离持续"时，不得给绿色 ok 判定。"""
    ad = make_adata(n_ctrl_donors=8, n_case_donors=8, cells_per_donor=25,
                    gene_pcts={"FLAT": (0.5, 0.5)})     # 两组无差异
    res = audit_robust_signal(ad, ["FLAT"], "condition", "donor", "control",
                              housekeeping=["ACTB"], cfg=_cfg())
    assert res.verdict_class in ("warn", "bad"), res.verdict
    assert not res.verdict.startswith("DEPTH-COMPARABLE —")


def test_robust_sign_consistency_required():
    """S3b：方向翻转不得算作"分离持续"。"""
    import scipy.sparse as sp
    # 每供者 40 细胞、两组各 8 供者 → 320 细胞/组；深度分层按 n_genes 人工分两档
    ad = make_adata(n_ctrl_donors=8, n_case_donors=8, cells_per_donor=40,
                    gene_pcts={"G": (0.5, 0.5)})
    ng = ad.obs["n_genes"].values
    is_ctrl = (ad.obs["condition"].values == "control")
    lo, hi = ng < 650, ng >= 650
    # 强制每档两组各 >= 30 细胞，否则分层会被跳过（不再是"空通过"）
    idx_lo = np.where(lo)[0]; idx_hi = np.where(hi)[0]
    assert len(idx_lo) >= 60 and len(idx_hi) >= 60, (len(idx_lo), len(idx_hi))

    X = ad.X.tolil()
    gi = ad.var_names.get_loc("G")
    X[:, gi] = 0
    X[np.where(is_ctrl & lo)[0], gi] = 1      # 低深度：ctrl 检出
    X[np.where(~is_ctrl & hi)[0], gi] = 1     # 高深度：case 检出（方向翻转）
    ad.X = sp.csr_matrix(X)

    res = audit_robust_signal(ad, ["G"], "condition", "donor", "control",
                              housekeeping=["ACTB"],
                              cfg=_cfg(min_cells_per_stratum=30))
    assert "G" in res.stratified, "分层未生成，测试无效"
    d = res.stratified["G"]["delta_pp"]
    assert len(d) >= 2, f"分层数不足: {len(d)}"
    assert (d > 0).any() and (d < 0).any(), f"未构造出方向翻转: {d.tolist()}"
    assert "G" not in res.genes_persisting, "方向翻转的基因不得判为'分离持续'"


def test_robust_mito_is_separate_verdict_not_depth():
    ad = make_adata()
    obs = ad.obs
    is_ctrl = obs["condition"].values == "control"
    obs.loc[is_ctrl, "percent_mito"] = 5.0
    obs.loc[~is_ctrl, "percent_mito"] = 20.0        # +15 pp
    res = audit_robust_signal(ad, ["TNNT2"], "condition", "donor", "control",
                              housekeeping=["ACTB"], cfg=_cfg())
    assert "CELL-STATE DIFFERS" in res.cell_state_verdict
    assert "DEPTH-COMPARABLE" in res.depth_verdict


def test_robust_empty_group_records_error():
    ad = make_adata()
    ad.obs.loc[ad.obs["condition"] == "disease", "condition"] = "control"
    res = audit_robust_signal(ad, ["TNNT2"], "condition", "donor", "control",
                              housekeeping=["ACTB"], cfg=_cfg())
    assert res.errors or "INSUFFICIENT" in res.depth_verdict


def test_robust_housekeeping_table():
    ad = make_adata()
    res = audit_robust_signal(ad, ["TNNT2"], "condition", "donor", "control",
                              housekeeping=["ACTB", "GAPDH"],
                              reference={"ACTB": 50.0, "GAPDH": 50.0}, cfg=_cfg())
    assert res.housekeeping is not None and len(res.housekeeping) == 2


# ════════════════ 检出比 ════════════════

def test_detection_ratio_math():
    df = compute_detection_ratio({"A": 10.0, "B": 60.0}, {"A": 50.0, "B": 60.0})
    d = dict(zip(df["gene"], df["ratio_R"]))
    assert abs(d["A"] - 0.2) < 1e-9 and abs(d["B"] - 1.0) < 1e-9


def test_reference_build_max_vs_median():
    ds = {"s1": {"G": 10.0}, "s2": {"G": 30.0}, "s3": {"G": 20.0}}
    assert build_reference_from_datasets(ds, "max")["G"] == 30.0
    assert build_reference_from_datasets(ds, "median")["G"] == 20.0


def test_gene_name_fallback():
    ad = make_adata()
    ad.var["gene_symbols"] = list(ad.var_names)
    ad.var_names = [f"ENSG{i:08d}" for i in range(ad.n_vars)]
    info = resolve_gene_names(ad, ["TNNT2", "ACTB"])
    assert not info["missing"]
    assert info.get("fallback_column") == "gene_symbols"


# ════════════════ 端到端 ════════════════

def test_end_to_end_pipeline(tmp_path):
    ad = make_adata(gene_pcts={"MYL3": (0.02, 0.0), "TNNT2": (0.9, 0.5),
                               "ACTB": (0.7, 0.7)})
    cfg = _cfg()
    f = audit_detection_floor(ad, ["MYL3", "TNNT2", "ACTB"], "condition", "donor",
                              "control", cfg=cfg, reference={"MYL3": 60.0})
    c = audit_control_composition(ad, ["TNNT2"], "condition", "donor", "control", cfg=cfg)
    r = audit_robust_signal(ad, ["TNNT2"], "condition", "donor", "control",
                            housekeeping=["ACTB"], cfg=cfg)
    assert len(f.table) == 3 and c.verdict and r.verdict
    f.table.to_csv(tmp_path / "floor.csv", index=False)
    assert (tmp_path / "floor.csv").exists()


def test_report_renders_selfcontained(tmp_path):
    from sc_audit.report import render_report
    ad = make_adata(gene_pcts={"MYL3": (0.02, 0.0), "TNNT2": (0.9, 0.5)})
    cfg = _cfg()
    f = audit_detection_floor(ad, ["MYL3", "TNNT2"], "condition", "donor",
                              "control", cfg=cfg, reference={"MYL3": 60.0})
    path = render_report(floor=f, outdir=str(tmp_path), input_name="synthetic.h5ad")
    html = open(path, encoding="utf-8").read()
    assert "UNQUANTIFIABLE" in html and "data:image/png;base64" in html
    import re
    assert not [u for u in re.findall(r'src="([^"]+)"', html) if not u.startswith("data:")]


# ════════════════ P1 回归测试 ════════════════

def test_p1a_mito_survives_gene_fallback_with_counts_switch(tmp_path):
    """P1-a 回归：基因名回退必须早于 compute_qc_metrics。

    否则 compute_qc_metrics 用 var_names 的前缀（MT-）找线粒体基因时，
    var_names 还是 ENSG…  → percent_mito 算成 NaN，并在 force=True 时
    覆盖掉上游正确的值（修复叠加引入的回归）。
    本测试走【完整 run】且【不传 --no-count-check】，以覆盖 force=True 路径。
    """
    import json
    import scipy.sparse as sp
    from sc_audit.cli import main as cli_main
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=25,
                    gene_pcts={"TNNT2": (0.9, 0.4)})
    # X 设为**非整数**（模拟 log-normalized）+ 整数 counts 层 → 触发换源 →
    # qc_recomputed=True 路径。注：make_adata 现在默认产生整数计数（使 B5 的
    # 计数校验在集成测试中真正执行），故须显式把 X 变成非整数，否则
    # resolve_counts 会直接采用 X、不换源。
    _raw = np.asarray(ad.X.todense())
    ad.layers["counts"] = sp.csr_matrix(_raw.astype(int))
    ad.X = sp.csr_matrix(np.log1p(_raw))
    syms = list(ad.var_names)
    ad.var["gene_symbols"] = syms
    ad.var_names = [f"ENSG{i:08d}" for i in range(ad.n_vars)]
    h5 = tmp_path / "ensg_switch.h5ad"
    ad.write_h5ad(h5)
    out = tmp_path / "o"
    cli_main(["run", "--input", str(h5), "--target-genes", "TNNT2,ACTB",
              "--n-perm", "100", "--outdir", str(out), "--quiet"])
    rob = (out / "robustness.txt").read_text()
    assert "percent_mito_median" in rob, "robustness.txt 未输出 mito 行"
    seg = rob.split("percent_mito_median")[1].split("\n")[0]
    assert "nan" not in seg.lower(), f"percent_mito 被覆盖为 NaN（P1-a 回归）: {seg}"
    meta = json.loads((out / "run_metadata.json").read_text())
    assert meta["qc_recomputed"] is True, "换源后应重算 QC"
    assert meta["counts_source"] not in ("X", "X (assumed)")
    assert meta["n_vars"] > 0


def test_p1b_force_only_when_source_actually_switched(tmp_path):
    """P1-b：X 本就是计数时不得强制重算 QC（必须走 CLI 路径验证，P2-1）。

    原测试另起炉灶直接调 io 函数，验证的是默认参数而非 _prepare 里的 force 表达式。
    """
    import json
    import scipy.sparse as sp
    from sc_audit.cli import main as cli_main
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20)
    ad.X = sp.csr_matrix((ad.X > 0).astype(int))     # 纯整数 → 无需换源
    ad.obs["n_genes"] = 12345                        # 哨兵：若被重算会变小
    h5 = tmp_path / "int.h5ad"
    ad.write_h5ad(h5)
    out = tmp_path / "o"
    rc = cli_main(["run", "--input", str(h5), "--target-genes", "TNNT2,ACTB",
                   "--n-perm", "50", "--outdir", str(out), "--quiet"])
    assert rc == 0
    meta = json.loads((out / "run_metadata.json").read_text())
    assert meta["counts_source"] in ("X", "X (assumed)"), \
        f"整数 X 不应换源，实际 {meta['counts_source']}"
    assert meta["qc_recomputed"] is False, \
        f"未换源却重算了 QC（P1-b 回归）: {meta['qc_recomputed']}"
    # P3：docstring 说哨兵是"未重算"的旁证 —— 必须真的断言它
    rob = (out / "robustness.txt").read_text()
    assert "12345" in rob, \
        "哨兵 n_genes=12345 未出现在深度表中，说明 QC 被重算过（哨兵被覆盖）"


def test_p1c_housekeeping_uses_same_ratio_basis():
    """P1-c：管家基因 R 必须与目标基因用同一分子口径（cfg.ratio_basis）。"""
    ad = make_adata(n_ctrl_donors=8, n_case_donors=8, cells_per_donor=25,
                    gene_pcts={"ACTB": (0.3, 0.9)})     # 病例组更高
    c = audit_robust_signal(ad, ["TNNT2"], "condition", "donor", "control",
                            housekeeping=["ACTB"], reference={"ACTB": 50.0},
                            cfg=_cfg(ratio_basis="case"))
    v = ad.obs["condition"].values
    v_case = (v != "control").sum(); v_ctrl = (v == "control").sum()
    hk = c.housekeeping
    assert "basis" in hk.columns and hk["basis"].iloc[0] == "case"
    # case 口径下的检出率应明显高于 control 口径（构造上 0.9 vs 0.3）
    from sc_audit.io import gene_detection_vector
    pct_case = gene_detection_vector(ad, "ACTB")[v != "control"].mean() * 100
    assert abs(float(hk["pct_in_dataset"].iloc[0]) - pct_case) < 1.0, \
        f"管家基因用了错误口径: {hk.to_dict('records')}"


def test_p1d_upper_bound_column_names_match_semantics():
    """P1-d：两列 CP 上界的命名必须与语义一致。"""
    ad = make_adata(n_ctrl_donors=8, n_case_donors=8, cells_per_donor=25,
                    gene_pcts={"G": (0.30, 0.0)})
    res = audit_detection_floor(ad, ["G"], "condition", "donor", "control", cfg=_cfg())
    cols = set(res.table.columns)
    assert "case_pct_upper95" in cols, "病例组上界列名应为 case_pct_upper95"
    assert "overall_pct_upper95" in cols, "全体上界列名应为 overall_pct_upper95"
    assert "case0_upper95" not in cols and "pct_case_upper95" not in cols, \
        "旧的错位列名仍存在"
    row = res.table.iloc[0]
    # 病例组检出率 0 → 其单侧上界应远小于 overall 上界（后者含 30% 的对照组）
    assert row["case_pct_upper95"] < row["overall_pct_upper95"]


def test_housekeeping_excluded_from_persistence_count():
    """P3：管家基因是阳性对照，不得计入 genes_persisting。"""
    ad = make_adata(n_ctrl_donors=8, n_case_donors=8, cells_per_donor=40,
                    gene_pcts={"ACTB": (0.3, 0.9), "TNNT2": (0.9, 0.4)})
    res = audit_robust_signal(ad, ["TNNT2", "ACTB"], "condition", "donor", "control",
                              housekeeping=["ACTB"], cfg=_cfg(min_cells_per_stratum=30))
    assert "ACTB" not in res.genes_persisting, "管家基因不应计入分离持续"
    assert any("housekeeping" in n.lower() for n in res.notes)


# ════════════════ 多水平 + preset/subset 崩溃回归 ════════════════

def _loose_cfg_yaml(tmp_path):
    """一份阈值极宽松的 config，让合成数据也能通过门控。"""
    p = tmp_path / "loose.yaml"
    p.write_text(
        "min_counts: 1\nmin_genes: 3\nmin_pos_markers: 1\n"
        "pos_markers: [TTN]\nneg_markers: [PTPRC]\n", encoding="utf-8")
    return p


def test_multilevel_condition_with_preset_does_not_crash(tmp_path):
    """P1-2 回归：三水平 condition + --preset 不得因 case_label 哨兵而崩溃。

    原 bug：ensure_obs_fields 在多水平时把 case_label 设为哨兵
    "<all non-control>"；CLI 在门控后把它原样传回，
    io 又拿去和真实取值比较 → ValueError。
    本测试用【宽松 config】保证门控后仍有细胞，从而真正走到二次校验那一行。
    """
    import json
    from sc_audit.cli import main as cli_main
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=25)
    idx = ad.obs.index[ad.obs["donor"].isin(["case_d0", "case_d1", "case_d2"])]
    ad.obs.loc[idx, "condition"] = "DCM"
    h5 = tmp_path / "three.h5ad"
    ad.write_h5ad(h5)
    out = tmp_path / "o"
    rc = cli_main(["run", "--input", str(h5), "--target-genes", "TNNT2,ACTB",
                   "--preset", "cardiac", "--config", str(_loose_cfg_yaml(tmp_path)),
                   "--n-perm", "100", "--outdir", str(out), "--quiet",
                   "--no-count-check"])
    assert rc == 0, "三水平 + --preset（门控后仍有细胞）不应崩溃"
    assert (out / "report.html").exists()
    meta = json.loads((out / "run_metadata.json").read_text())
    assert meta["case_scope"] == "all non-control"
    assert meta["n_condition_levels"] >= 3


def test_gating_to_empty_aborts_with_actionable_message(tmp_path, capsys):
    """门控把细胞全过滤掉时，必须给出可诊断的中止（而非"标签不在 [] 中"）。"""
    from sc_audit.cli import main as cli_main
    ad = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=10)
    h5 = tmp_path / "tiny.h5ad"
    ad.write_h5ad(h5)
    with pytest.raises(SystemExit) as e:
        cli_main(["floor", "--input", str(h5), "--target-genes", "TNNT2",
                  "--preset", "cardiac",     # 默认阈值 500 counts —— 合成数据必然全灭
                  "--outdir", str(tmp_path / "o"), "--quiet", "--no-count-check"])
    assert e.value.code == 2
    err = capsys.readouterr().err
    assert "门控后没有细胞剩余" in err, f"错误信息不可诊断: {err}"
    assert "标签" not in err and "control_label" not in err


def test_multilevel_condition_with_cell_subset_does_not_crash(tmp_path):
    """P1-2 回归：三水平 + --cell-subset（同样会二次调用 ensure_obs_fields）。"""
    from sc_audit.cli import main as cli_main
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=25)
    idx = ad.obs.index[ad.obs["donor"] == "case_d0"]
    ad.obs.loc[idx, "condition"] = "DCM"
    h5 = tmp_path / "three2.h5ad"
    ad.write_h5ad(h5)
    out = tmp_path / "o2"
    rc = cli_main(["floor", "--input", str(h5), "--target-genes", "TNNT2",
                   "--cell-subset", "n_genes > 0", "--n-perm", "50",
                   "--outdir", str(out), "--quiet", "--no-count-check"])
    assert rc == 0, "三水平 + --cell-subset 不应崩溃"
    assert (out / "floor_audit.csv").exists()


def test_p2a_no_match_vars_branch_reports_correct_cause():
    """P2-a：所有匹配变量被剔除时，note 必须说明"未尝试匹配"而非"无人匹配"。"""
    ad = make_adata(design_confounded=False)
    # 把 match_vars 全设成供者级唯一的列 → 全部被剔除
    for c in ("institution", "donation_type", "source", "region"):
        # P3：用【供者级】唯一（真实场景），覆盖 donor_unique 分支
        ad.obs[c] = ["v_" + str(d) for d in ad.obs["donor"].astype(str)]
    cfg = _cfg(match_vars=["institution", "donation_type", "source", "region"])
    res = audit_control_composition(ad, ["TNNT2"], "condition", "donor", "control",
                                    cfg=cfg)
    assert res.match_attempted is False
    assert res.match_not_possible_on, "应记录无法匹配的变量"
    assert any("NO matching was ATTEMPTED" in n for n in res.notes), \
        f"note 未说明真实原因: {res.notes}"
    assert not any("NO control donor matched" in n for n in res.notes), \
        "不得给出'无人匹配'这种错误原因（实际根本没尝试）"


def test_p2_3_api_path_enforces_donor_nesting():
    """P2-3：Python API 直接调用也必须校验供者 ⊂ 组。"""
    ad = make_adata()
    idx = ad.obs.index[ad.obs["donor"] == "ctrl_d0"][:5]
    ad.obs.loc[idx, "condition"] = "disease"     # 制造跨组供者
    with pytest.raises(ValueError, match="多个"):
        audit_detection_floor(ad, ["TNNT2"], "condition", "donor", "control", cfg=_cfg())
    with pytest.raises(ValueError, match="多个"):
        audit_control_composition(ad, ["TNNT2"], "condition", "donor", "control", cfg=_cfg())
    with pytest.raises(ValueError, match="多个"):
        audit_robust_signal(ad, ["TNNT2"], "condition", "donor", "control", cfg=_cfg())


def test_fdr_exclusion_follows_housekeeping_cli(tmp_path):
    """P2-2：--housekeeping 必须同步到 FDR 的排除集 m。"""
    import json
    from sc_audit.cli import main as cli_main
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=25)
    h5 = tmp_path / "hk.h5ad"
    ad.write_h5ad(h5)
    out = tmp_path / "o"
    cli_main(["run", "--input", str(h5), "--target-genes", "TNNT2,ACTB,B2M",
              "--housekeeping", "ACTB,B2M", "--n-perm", "100",
              "--outdir", str(out), "--quiet", "--no-count-check"])
    meta = json.loads((out / "run_metadata.json").read_text())
    assert set(meta["config"]["fdr_exclude_from_m"]) == {"ACTB", "B2M"}, \
        "housekeeping 未同步到 FDR 排除集"


def test_p2a_partial_drop_is_flagged():
    """P2-a：**部分剔除**也必须记录，不能只在"全剔除"时提示。

    场景：institution 供者级唯一 → 被剔除；region 单水平（heart）
    → 剩下的 region 上所有对照都平凡匹配，衰减比≈1 但 verdict 仍 COMPARABLE。
    """
    ad = make_adata(design_confounded=False)
    ad.obs["institution"] = ["v_" + str(d) for d in ad.obs["donor"].astype(str)]
    ad.obs["donation_type"] = "biopsy"
    ad.obs["source"] = "Cells"
    ad.obs["region"] = "heart"          # 单水平
    cfg = _cfg(match_vars=["institution", "region"])
    res = audit_control_composition(ad, ["TNNT2"], "condition", "donor", "control",
                                    cfg=cfg)
    assert res.match_attempted is True, "仍然尝试了匹配（用 region）"
    assert res.match_not_possible_on == ["institution"], \
        f"部分剔除未被记录: {res.match_not_possible_on}"
    notes = " ".join(res.notes)
    assert "WEAKER than a same-design comparison" in notes, "缺少降级说明"
    assert "TRIVIALLY satisfied" in notes, "缺少'平凡匹配'说明"


def test_l5_jackknife_bh_crit_column():
    """L5：jackknife 增加 BH 临界值口径的计数列，且用主分析临界值（不各自做 BH）。"""
    ad = make_adata(n_ctrl_donors=8, n_case_donors=8, cells_per_donor=40,
                    gene_pcts={"G": (0.9, 0.3), "H": (0.85, 0.4)})
    res = audit_detection_floor(ad, ["G", "H", "ACTB"], "condition", "donor",
                                "control", cfg=_cfg(fdr_exclude_from_m=["ACTB"]))
    assert "jackknife_n_below_bh_crit" in res.table.columns
    assert "jackknife_n_sig" in res.table.columns, "未校正列应保留"
    # 被排除出 m 的基因无 BH 临界值 → NaN
    actb = res.table[res.table["gene"] == "ACTB"].iloc[0]
    assert actb["jackknife_n_below_bh_crit"] != actb["jackknife_n_below_bh_crit"], \
        "被 FDR 排除的基因该列应为 NaN"
    # 未校正口径 >= BH 口径（临界值更严）
    for _, r in res.table.iterrows():
        a, b = r["jackknife_n_sig"], r["jackknife_n_below_bh_crit"]
        if a == a and b == b:
            assert b <= a, f"{r['gene']}: BH 口径计数({b})不应超过未校正({a})"


def test_event_sparse_gene_gets_no_p_and_no_jackknife_note_matches():
    """事件稀疏基因：不计算全样本 P，也**不进入** jackknife；提示文案须与行为一致（v0.3.3）。

    旧文案写 "jackknife is still reported"，但实现里该分支 `continue` 跳过整段
    jackknife —— 出厂示例产物（`floor_jackknife.csv`）里确实没有事件稀疏基因。
    """
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=20,
                    gene_pcts={"MYL3": (0.02, 0.0), "TNNT2": (0.9, 0.5)})
    f = audit_detection_floor(ad, ["MYL3", "TNNT2", "ACTB"], "condition", "donor",
                              "control", cfg=_cfg())
    row = f.table.set_index("gene").loc["MYL3"]
    assert row["stratum"] == "event_sparse"
    note = str(row["perm_mode"])
    assert "event-sparse" in note
    assert "excluded from the jackknife" in note
    assert "still reported" not in note
    assert str(f.summary["perm_p_note"]).find("jackknife only") == -1
    assert f.jackknife is None or "MYL3" not in f.jackknife
    assert f.jackknife is None or "TNNT2" in f.jackknife
