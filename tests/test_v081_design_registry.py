"""v0.8.1 回归：列注册表复用 + 人工确认 + 配置装配（design 闭环）。"""
import json

import numpy as np
import pandas as pd

from conftest import make_adata
from sc_audit.config import AuditConfig
from sc_audit.design_parse import (export_inventory, apply_proposal, confirm_registry,
                                   apply_design_config_to_cfg)


def _inv():
    A = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20,
                   design_confounded=False)
    A.obs["condition"] = A.obs["condition"].astype(str)
    return export_inventory(A, "condition", "donor", "control", label="ds1")


def _prop(role_map):
    return {"ds1": {"contrast_column": "condition", "unit_column": "donor",
                    "columns": [{"column": c, "role": r, "confidence": 0.8,
                                 "evidence": "x", "abstain": False}
                                for c, r in role_map.items()]}}


def test_registry_confirmed_overrides_ai_proposal(tmp_path):
    inv = _inv()
    reg = str(tmp_path / "design_registry.json")
    # 第一次 apply：batch 被 AI 判为 adjustable（结构复验会 disagree）
    apply_proposal(inv, _prop({"condition": "condition", "donor": "unit",
                               "batch": "adjustable"}), str(tmp_path), registry_path=reg)
    # 人工确认 batch=technical
    n = confirm_registry(reg, {"ds1": {"batch": "technical"}})
    assert n == 1
    # 第二次 apply：即使 AI 改了提案，注册表优先
    s = apply_proposal(inv, _prop({"condition": "condition", "donor": "unit",
                                   "batch": "design_defining"}), str(tmp_path),
                       registry_path=reg)
    df = pd.read_csv(tmp_path / "design_reverify.csv")
    row = df[df["column"] == "batch"].iloc[0]
    assert row["role"] == "technical" and row["tier"] == "confirmed"
    assert row["evidence"] == "registry(confirmed)"
    assert s["n_confirmed"] == 1
    # 注册表不被改回 proposed
    d = json.load(open(reg))
    assert d["ds1"]["batch"]["status"] == "confirmed"


def test_design_config_to_cfg():
    dc = {"role_overrides": {"region": "design", "fatal": "outcome"},
          "identifiability_nondesign_cols": ["condition_original"],
          "identity_source_col": "Sample_Source"}
    cfg = apply_design_config_to_cfg(AuditConfig(), dc)
    assert cfg.role_overrides["region"] == "design"
    assert cfg.role_overrides["fatal"] == "outcome"
    assert cfg.identifiability_nondesign_cols == ["condition_original"]
    assert cfg.identity_source_col == "Sample_Source"
