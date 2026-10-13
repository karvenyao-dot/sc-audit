"""v0.7.2 回归：外部独立评估（2026-10-12）采信项 C4/C5/C6/C7/C8/C9/C11。

C4  审计四：门控通过率的**绝对量**与**按臂差**判据（原只比相对 margin → 0.3% 通过率仍 PASS）
C5  审计二：无 ≥2 水平设计变量 → `NOT ASSESSED`（原输出褒义的 COMPARABLE）
C6  `check`：reference 未传入时不得断言“参照表只覆盖 0/N”
C7  审计二：匹配集 < matched_min_cells_warn → 衰减比不作判据（与审计二b 一致）
C8  审计四：rank 分母用注释谱系总数（原显示 rank 1/0）
C9  审计一：表内注明 verdict 的基准组与阈值（stratum 与 verdict 基准不同）
C11 审计一：R2 未强制时显式登记（与 R1 对称）
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import anndata as ad

from conftest import make_adata
from sc_audit.config import AuditConfig

POS = ["TTN", "MYH7", "MYBPC3", "ACTN2", "TNNI3", "MYL2"]
NEG = ["PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"]
EXTRA = ["MALAT1", "RYR2", "CACNA1C", "KCNQ1", "SCN5A"]
GENES = POS + NEG + EXTRA


def _identity_ad(pos_p=0.85, neg_p=0.02, n_ctrl_don=2, n_case_don=2, cells=50,
                 with_ngenes=True, single_lineage=True):
    """构造审计四输入：detection 概率按臂可控（用于制造低通过率 / 按臂偏差）。"""
    rng = np.random.default_rng(0)
    obs = []
    for arm, nd in (("control", n_ctrl_don), ("disease", n_case_don)):
        for d in range(nd):
            for _i in range(cells):
                row = {"condition": arm, "donor": f"{arm[0]}{d}",
                       "lineage_names": "Cardiomyocytes", "source": "Cells"}
                if with_ngenes:
                    row["n_genes"] = 900
                obs.append(row)
    n = len(obs)
    X = np.zeros((n, len(GENES)), dtype=np.int32)
    for i, ob in enumerate(obs):
        p_pos = pos_p if isinstance(pos_p, float) else pos_p[ob["condition"]]
        p_neg = neg_p if isinstance(neg_p, float) else neg_p[ob["condition"]]
        for j, g in enumerate(GENES):
            p = p_pos if g in POS else (p_neg if g in NEG else 0.7)
            X[i, j] = 1 if rng.random() < p else 0
    return ad.AnnData(X=sp.csr_matrix(X),
                      obs=pd.DataFrame(obs, index=[f"c{i}" for i in range(n)]),
                      var=pd.DataFrame(index=GENES))


def _cfg():
    cfg = AuditConfig()
    cfg.identity_col = "lineage_names"
    cfg.identity_target = "Cardiomyocytes"
    return cfg


# ──────────────────────── C4 · gate pass-rate absolute + arm-delta ────────────────────────

def test_c4_low_absolute_gate_pass_not_pass():
    from sc_audit.identity import audit_cell_identity
    A = _identity_ad(pos_p=0.0)          # 阳性标志全部未检出 → 门控通过率 0%
    r = audit_cell_identity(A, "condition", "donor", "control", cfg=_cfg())
    assert r.target_gate_pass == 0.0
    assert r.verdict.startswith("GATE-NOT-USABLE"), r.verdict


def test_c4_arm_biased_gate_is_condition_collinear():
    from sc_audit.identity import audit_cell_identity
    A = _identity_ad(pos_p={"control": 0.9, "disease": 0.0})   # 通过率按臂极端偏
    r = audit_cell_identity(A, "condition", "donor", "control", cfg=_cfg())
    assert r.verdict.startswith("COMPOSITION-CONFOUNDED"), r.verdict
    assert any("condition-collinear" in n and "differs" in n for n in r.notes), r.notes


# ──────────────────────── C5 · NOT ASSESSED tier ────────────────────────

def test_c5_no_multi_level_design_variable_is_not_assessed_not_comparable():
    from sc_audit.composition import audit_control_composition
    A = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=20,
                   design_confounded=False)      # institution/donation_type/source/region 均单水平
    A.obs["condition"] = A.obs["condition"].astype(str)
    r = audit_control_composition(A, ["ACTB"], "condition", "donor", "control",
                                  cfg=AuditConfig())
    assert r.verdict.startswith("NOT ASSESSED"), r.verdict
    assert "COMPARABLE" not in r.verdict


def test_c5_multi_level_design_variable_still_assessed():
    from sc_audit.composition import audit_control_composition
    A = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=20,
                   design_confounded=True)       # 对照 institution 有两水平
    A.obs["condition"] = A.obs["condition"].astype(str)
    r = audit_control_composition(A, ["ACTB"], "condition", "donor", "control",
                                  cfg=AuditConfig())
    assert not r.verdict.startswith("NOT ASSESSED"), r.verdict


# ──────────────────────── C7 · attenuation suppressed when matched set tiny ────────────────────────

def test_c7_small_matched_set_drops_attenuation_from_criteria():
    from sc_audit.composition import audit_control_composition
    from scipy import sparse
    A = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=20,
                   design_confounded=True)
    A.obs["condition"] = A.obs["condition"].astype(str)
    # 让 TTN 的检出由 INST_A/INST_B 决定（INST_A=1, INST_B=0）：
    # 匹配子集只剩 INST_A 的 1 个对照供者 → matched effect = 0 → 衰减比 < 阈值。
    gi = list(A.var_names).index("TTN")
    X = A.X.toarray()
    X[:, gi] = (A.obs["institution"].astype(str).values == "INST_A").astype(float)
    A.X = sparse.csr_matrix(X)
    cfg = AuditConfig()
    r = audit_control_composition(A, ["TTN"], "condition", "donor", "control", cfg=cfg)
    assert r.matched_cells and r.matched_cells < cfg.matched_min_cells_warn
    assert r.attenuation_suppressed is True
    assert "median attenuation ratio=" not in r.verdict
    assert any("NOT used as a verdict criterion" in n for n in r.notes), r.notes


# ──────────────────────── C6 · check reference=None wording ────────────────────────

def test_c6_check_reference_none_says_not_passed_in():
    from sc_audit.check import check_readiness
    A = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=20)
    A.obs["condition"] = A.obs["condition"].astype(str)
    r = check_readiness(A, "condition", "donor", "control",
                        target_genes=["ACTB", "GAPDH"], reference=None)
    assert any("未提供参照表" in w for w in r.warnings), r.warnings
    assert not any("参照表只覆盖 0/" in w for w in r.warnings), r.warnings


# ──────────────────────── C8 · rank denominator uses annotation count ────────────────────────

def test_c8_rank_denominator_uses_annotation_count():
    from sc_audit.identity import audit_cell_identity
    A = _identity_ad(with_ngenes=False)   # 无 n_genes → complexity_by_lineage 不建立
    r = audit_cell_identity(A, "condition", "donor", "control", cfg=_cfg())
    assert r.n_lineages == 1
    assert "rank 1/1" in r.to_text(), r.to_text()


# ──────────────────────── C9 · verdict basis columns ────────────────────────

def test_c9_verdict_basis_columns_present():
    from sc_audit.floor import audit_detection_floor
    A = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=20)
    A.obs["condition"] = A.obs["condition"].astype(str)
    cfg = AuditConfig()
    t = audit_detection_floor(A, ["ACTB"], "condition", "donor", "control",
                              cfg=cfg).table
    assert "verdict_basis_group" in t.columns
    assert "verdict_floor_pct" in t.columns
    assert t.iloc[0]["verdict_basis_group"] == cfg.ratio_basis
    assert float(t.iloc[0]["verdict_floor_pct"]) == float(cfg.floor_pct)


# ──────────────────────── C11 · R2 not-enforced registration ────────────────────────

def test_c11_r2_not_enforced_is_registered():
    from sc_audit.floor import audit_detection_floor
    A = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=20)
    A.obs["condition"] = A.obs["condition"].astype(str)
    with pytest.warns(UserWarning, match="R2 .* NOT ENFORCED"):
        res = audit_detection_floor(A, ["ACTB"], "condition", "donor", "control",
                                    cfg=AuditConfig())
    assert res.summary.get("r2_enforcement_note")


def test_c11_r2_enforced_when_modality_declared():
    from sc_audit.floor import audit_detection_floor
    A = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=20)
    A.obs["condition"] = A.obs["condition"].astype(str)
    res = audit_detection_floor(A, ["ACTB"], "condition", "donor", "control",
                                cfg=AuditConfig(), modality="scRNA",
                                reference_modalities=["scRNA"])
    assert res.summary.get("r2_enforcement_note") is None
