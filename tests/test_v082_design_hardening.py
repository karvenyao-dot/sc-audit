"""v0.8.2 回归：据第三方复核（2026-10-12）加固 design 契约。

A.6  数据集名不匹配 → 报错（不再静默取 proposal 第一条）
A.7#10 export 校验 control/case 标签存在
A.7#11 proposal 含数据中不存在的列 → 记录（幻觉可测）
A.8  复验与 roster 同源（redundant_with_condition）
A.9  常量列（单水平）的 technical 不自动接受
B    proposal 出处块（attest）+ apply 记录 attested
"""
import json

import numpy as np
import pandas as pd
import pytest

from conftest import make_adata
from sc_audit.design_parse import (export_inventory, apply_proposal, attest_proposal,
                                   redundant_with_condition, _tier)


def _export(label="ds1", control="control"):
    A = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20,
                   design_confounded=False)
    A.obs["condition"] = A.obs["condition"].astype(str)
    return export_inventory(A, "condition", "donor", control, label=label)


def _prop(ds, role_map, extra_cols=()):
    cols = [{"column": c, "role": r, "confidence": 0.8, "evidence": "x",
             "abstain": False} for c, r in role_map.items()]
    cols += [{"column": c, "role": "technical", "confidence": 0.5, "evidence": "x",
              "abstain": False} for c in extra_cols]
    return {ds: {"contrast_column": "condition", "unit_column": "donor", "columns": cols}}


# ── A.6 ──
def test_a6_dataset_name_mismatch_raises(tmp_path):
    inv = _export(label="dsA")
    with pytest.raises(ValueError, match="不含数据集"):
        apply_proposal(inv, _prop("dsB", {"condition": "condition"}), str(tmp_path))


# ── A.7#10 ──
def test_a7_10_export_rejects_missing_control_label():
    A = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=20)
    A.obs["condition"] = A.obs["condition"].astype(str)
    with pytest.raises(ValueError, match="对照组为空"):
        export_inventory(A, "condition", "donor", "ZZZ", label="ds1")


# ── A.7#11 + B ──
def test_a7_11_hallucinated_columns_recorded(tmp_path):
    inv = _export()
    prop = _prop("ds1", {"condition": "condition", "donor": "unit"},
                 extra_cols=("ghost_col", "another_ghost"))
    s = apply_proposal(inv, prop, str(tmp_path))
    assert s["n_hallucinated"] == 2
    assert "ghost_col" in s["hallucinated_columns"]
    sheet = (tmp_path / "design_confirm_sheet.md").read_text()
    assert "不存在的列" in sheet


def test_b_inventory_sha_and_attest_roundtrip(tmp_path):
    inv = _export()
    assert len(inv["inventory_sha256"]) == 64
    pinv = tmp_path / "inventory.json"; pinv.write_text(json.dumps(inv))
    pp = tmp_path / "proposal.json"
    pp.write_text(json.dumps(_prop("ds1", {"condition": "condition", "donor": "unit"})))
    # apply 前：无出处
    s0 = apply_proposal(inv, json.loads(pp.read_text()), str(tmp_path))
    assert s0["attested"] is False
    # attest → apply 后：有出处
    gb = attest_proposal(str(pp), str(pinv), model="test-model-1")
    assert gb["inventory_sha256"] == inv["inventory_sha256"]
    s1 = apply_proposal(inv, json.loads(pp.read_text()), str(tmp_path))
    assert s1["attested"] is True and s1["generated_by"]["model"] == "test-model-1"
    assert "AI 生成、未经复核" in (tmp_path / "design_confirm_sheet.md").read_text()


# ── A.8 ──
def test_a8_redundant_uses_shared_predicate():
    assert redundant_with_condition(0.0, 0.0, 2) is True
    assert redundant_with_condition(0.0, 0.5, 2) is False
    assert redundant_with_condition(0.0, 0.0, 1) is False


def test_a8_donor_unique_confounder_is_flagged(tmp_path):
    """供者级 age（每供者互异、列名属 confounder）不应被判 donor_unique。"""
    A = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20,
                   design_confounded=False)
    A.obs["condition"] = A.obs["condition"].astype(str)
    donors = sorted(set(A.obs["donor"].astype(str)))
    A.obs["age"] = A.obs["donor"].astype(str).map(
        {d: str(20 + i) for i, d in enumerate(donors)}).astype(str)
    inv = export_inventory(A, "condition", "donor", "control", label="ds1")
    s = apply_proposal(inv, _prop("ds1", {"condition": "condition", "donor": "unit",
                                          "age": "donor_unique"}), str(tmp_path))
    row = pd.read_csv(tmp_path / "design_reverify.csv")
    r = row[row["column"] == "age"].iloc[0]
    assert r["reverify"] == "disagree"


# ── A.9 ──
def test_a9_constant_technical_not_auto_accepted():
    assert _tier("technical", {"n_levels": 1}) == "must_review"     # 常量列
    assert _tier("technical", {"n_levels": 3}) == "auto_accept"     # 非常量
    assert _tier("design_defining", {"n_levels": 1}) == "must_review"
