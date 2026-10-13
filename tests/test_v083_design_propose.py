"""v0.8.3 回归：design 内置调用器 `propose`（backend none/local/http）。

默认 none（不生成）；local 用本地命令；http 用 OpenAI 兼容端点（延迟导入，默认路径不碰网络）。
"""
import json

import numpy as np
import pytest

from conftest import make_adata
from sc_audit.design_parse import export_inventory, propose, _extract_json


def _inv():
    A = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=20,
                   design_confounded=False)
    A.obs["condition"] = A.obs["condition"].astype(str)
    return export_inventory(A, "condition", "donor", "control", label="ds1")


def test_backend_none_refuses():
    with pytest.raises(ValueError, match="none"):
        propose(_inv(), backend="none")


def test_backend_local_roundtrip(tmp_path):
    fixture = {"ds1": {"contrast_column": "condition", "unit_column": "donor",
                       "columns": [{"column": "condition", "role": "condition",
                                    "confidence": 0.9, "evidence": "x", "abstain": False},
                                   {"column": "donor", "role": "unit",
                                    "confidence": 0.9, "evidence": "x", "abstain": False}]}}
    fx = tmp_path / "fixture.json"
    fx.write_text(json.dumps(fixture))
    out = tmp_path / "proposal.json"
    r = propose(_inv(), backend="local", local_cmd=f"cat {fx}", out_path=str(out))
    assert out.exists()
    assert r["generated_by"]["backend"] == "local"
    assert r["generated_by"]["inventory_sha256"]
    assert "AI 生成、未经复核" in r["generated_by"]["disclaimer"]


def test_backend_http_requires_endpoint_and_model():
    with pytest.raises(ValueError, match="endpoint"):
        propose(_inv(), backend="http")


def test_extract_json_tolerates_fences_and_prose():
    assert _extract_json('前言 ```json\n{"a": 1}\n``` 后语') == {"a": 1}
    assert _extract_json('```{"b": 2}```') == {"b": 2}
