"""v0.8.0 回归：`design` 子命令（解耦版 AI 设计解析）的确定性契约。

export：从数据出 inventory（离线）；apply：读回 proposal → 确定性复验 + 确认单 + 注册表。
本模块**不调 LLM、不联网**。
"""
import json
import os

import numpy as np
import pandas as pd
import pytest

from conftest import make_adata
from sc_audit.design_parse import (export_inventory, write_inventory, apply_proposal,
                                   TIER_MUST_REVIEW, TIER_AUTO)


def _ad():
    A = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20,
                   design_confounded=False)
    A.obs["condition"] = A.obs["condition"].astype(str)
    A.obs["condition_original"] = np.where(
        A.obs["condition"].values == "control", "healthy", "disease")
    return A


def _export():
    A = _ad()
    return export_inventory(A, "condition", "donor", "control", label="test_ds")


def _proposal(role_map):
    return {"test_ds": {
        "contrast_column": "condition", "unit_column": "donor",
        "columns": [{"column": c, "role": r, "confidence": 0.8, "evidence": "x",
                     "abstain": False} for c, r in role_map.items()]}}


def test_export_has_expected_fields():
    inv = _export()
    assert inv["dataset"] == "test_ds"
    assert inv["n_donors"] == 8
    cols = {c["col"]: c for c in inv["columns"]}
    assert "donor_unique" in cols["donor"]
    assert cols["condition_original"]["donor_level_levels"] == 2
    assert cols["condition_original"]["ent_case"] <= 1e-9
    assert cols["condition_original"]["ent_ctrl"] <= 1e-9
    assert "baseline_roles" in inv and inv["baseline_roles"].get("condition_original")


def test_apply_agree_and_tiers(tmp_path):
    inv = _export()
    prop = _proposal({"condition": "condition", "donor": "unit",
                      "condition_original": "redundant_with_condition",
                      "batch": "technical", "region": "design_defining"})
    s = apply_proposal(inv, prop, str(tmp_path))
    # condition_original 是 condition 副本 → 结构复验 agree；batch 供者唯一 → technical 自动档
    df = pd.read_csv(tmp_path / "design_reverify.csv")
    rv = dict(zip(df["column"], df["reverify"]))
    assert rv["condition_original"] == "agree"
    assert rv["donor"] == "agree"
    assert s["n_disagree"] == 0
    assert (tmp_path / "design_confirm_sheet.md").exists()
    cfg = json.load(open(tmp_path / "design_config.json"))
    assert cfg["role_overrides"]["region"] == "design"
    assert "condition_original" in cfg["identifiability_nondesign_cols"]
    assert "不会自动写入" in cfg["_note"]


def test_apply_flags_structural_contradiction(tmp_path):
    """提案说 batch 是 design_defining，但 batch 每供者唯一 → 结构复验必须 disagree。"""
    inv = _export()
    prop = _proposal({"condition": "condition", "donor": "unit", "batch": "design_defining"})
    s = apply_proposal(inv, prop, str(tmp_path))
    df = pd.read_csv(tmp_path / "design_reverify.csv")
    rv = dict(zip(df["column"], df["reverify"]))
    assert rv["batch"] == "disagree"
    assert s["n_disagree"] == 1


def test_apply_registry_is_written_and_merged(tmp_path):
    inv = _export()
    prop = _proposal({"condition": "condition", "donor": "unit"})
    reg = str(tmp_path / "design_registry.json")
    apply_proposal(inv, prop, str(tmp_path), registry_path=reg)
    reg2 = str(tmp_path / "design_registry.json")
    assert os.path.exists(reg2)
    d = json.load(open(reg2))
    assert d["test_ds"]["condition"]["status"] == "proposed"
    # 第二次 apply 不应覆盖既有条目
    apply_proposal(inv, _proposal({"condition": "condition"}), str(tmp_path),
                   registry_path=reg2)
    d2 = json.load(open(reg2))
    assert d2["test_ds"]["condition"]["role"] == "condition"
