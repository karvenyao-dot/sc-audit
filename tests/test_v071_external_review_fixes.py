"""v0.7.1 回归：第三方独立评估 独立评估（2026-10-12）经交叉验证后采信的修复项。

对应发现编号：
  C1  condition 的"标签副本"列被当 design-defining → 共同支持域被打到 0（假拒绝）
  C2  放行路径（--allow-cross-modality / --allow-depth-mismatch）只打标记、不限制判定
  C12 `check` 对 backed='r' 输入崩溃（绕过 io.x_matrix）
  C3  审计四可评估性备注硬编码（本文件覆盖 C1/C2/C12；C3 见 test_identity.py）
"""
import warnings

import numpy as np
import pytest

from conftest import make_adata
from sc_audit.config import AuditConfig


# ──────────────────────── C1 · redundant-with-condition ────────────────────────

def _ad_with_condition_duplicate():
    """case→'disease' / control→'healthy' 的 condition 副本列（与 condition 一对一）。"""
    A = make_adata(n_ctrl_donors=4, n_case_donors=3, cells_per_donor=20,
                   design_confounded=False)
    A.obs["condition"] = A.obs["condition"].astype(str)
    A.obs["condition_original"] = np.where(
        A.obs["condition"].values == "control", "healthy", "disease")
    return A


def test_c1_condition_duplicate_is_excluded_not_design_defining():
    from sc_audit.identifiability import audit_identifiability, BUCKET_REDUNDANT_CONDITION
    A = _ad_with_condition_duplicate()
    r = audit_identifiability(A, ["ACTB", "GAPDH"], "condition", "donor", "control",
                              cfg=AuditConfig())
    assert "condition_original" in r.excluded.get("redundant-with-condition", []), \
        r.roster.to_string()
    assert "condition_original" not in r.design_defining
    # 核心：不再因副本列把支持域打到 0
    assert r.support_cells > 0, "副本列仍被判 design-defining → 支持域被清空"
    assert r.support_cells == r.details["n_control_cells"]
    assert not str(r.verdict).startswith("SUPPORT-BELOW-FLOOR"), r.verdict


def test_c1_guard_does_not_misfire_on_real_design_variable():
    """对照组内有变异的 institution 是真正的 design-defining，不得被误排除。"""
    from sc_audit.identifiability import audit_identifiability, BUCKET_REDUNDANT_CONDITION
    A = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=20,
                   design_confounded=True)   # control: INST_A/INST_B 两水平
    A.obs["condition"] = A.obs["condition"].astype(str)
    r = audit_identifiability(A, ["ACTB", "GAPDH"], "condition", "donor", "control",
                              cfg=AuditConfig())
    assert "institution" in r.design_defining, r.roster.to_string()
    assert "institution" not in r.excluded.get("redundant-with-condition", [])


def test_c1_without_duplicate_behaviour_unchanged():
    from sc_audit.identifiability import audit_identifiability
    A = make_adata(n_ctrl_donors=4, n_case_donors=3, cells_per_donor=20,
                   design_confounded=False)
    A.obs["condition"] = A.obs["condition"].astype(str)
    r = audit_identifiability(A, ["ACTB", "GAPDH"], "condition", "donor", "control",
                              cfg=AuditConfig())
    assert not r.excluded.get("redundant-with-condition")


# ──────────────────────── C2 · allow-path must restrict verdict ────────────────────────

def _tiny_ad(pct=0.9):
    return make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=20,
                      gene_pcts={"ACTB": (pct, pct)})


def test_c2_cross_modal_reference_downgrades_verdict():
    from sc_audit.floor import audit_detection_floor
    A = _tiny_ad()
    A.obs["condition"] = A.obs["condition"].astype(str)
    with pytest.warns(UserWarning, match="跨模态"):
        t = audit_detection_floor(A, ["ACTB"], "condition", "donor", "control",
                                  cfg=AuditConfig(), reference={"ACTB": 90.0},
                                  modality="scRNA", reference_modalities=["snRNA"],
                                  allow_cross_modality=True).table
    assert str(t.iloc[0]["verdict"]).startswith("INSUFFICIENT (cross-modal reference)"), \
        t.iloc[0]["verdict"]
    assert bool(t.iloc[0]["verdict_contract_ok"]) is False


def test_c2_same_modality_not_downgraded_and_no_extra_column():
    from sc_audit.floor import audit_detection_floor
    A = _tiny_ad()
    A.obs["condition"] = A.obs["condition"].astype(str)
    t = audit_detection_floor(A, ["ACTB"], "condition", "donor", "control",
                              cfg=AuditConfig(), reference={"ACTB": 90.0},
                              modality="scRNA",
                              reference_modalities=["scRNA"]).table
    assert str(t.iloc[0]["verdict"]) == "QUANTIFIABLE"
    assert "verdict_contract_ok" not in t.columns   # 未受影响 → 输出零改动


def test_c2_depth_mismatch_downgrades_verdict():
    from sc_audit.floor import audit_detection_floor
    A = _tiny_ad()
    A.obs["condition"] = A.obs["condition"].astype(str)
    with pytest.warns(UserWarning, match="DEPTH-MISMATCH"):
        t = audit_detection_floor(A, ["ACTB"], "condition", "donor", "control",
                                  cfg=AuditConfig(), reference={"ACTB": 90.0},
                                  reference_depth=5000.0, case_depth=1000.0,
                                  allow_depth_mismatch=True).table
    assert str(t.iloc[0]["verdict"]).startswith("INSUFFICIENT (depth mismatch)"), \
        t.iloc[0]["verdict"]
    assert bool(t.iloc[0]["verdict_contract_ok"]) is False


# ──────────────────────── C12 · backed input must not crash ────────────────────────

def test_c12_backed_input_depth_bin_capacity(tmp_path):
    import anndata as ad
    A = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=20)
    A.obs["condition"] = A.obs["condition"].astype(str)
    del A.obs["n_genes"]                      # 强制走 X 路径（触发原崩溃）
    p = tmp_path / "backed.h5ad"
    A.write(p)
    B = ad.read_h5ad(p, backed="r")
    assert getattr(B, "isbacked", False)
    from sc_audit.check import _depth_bin_capacity
    is_ctrl = (B.obs["condition"].astype(str).values == "control")
    out = _depth_bin_capacity(B, is_ctrl, AuditConfig())
    assert out["n_bins"] >= 1 and "bins" in out
