"""设计声明 + 能力闸门（v0.4.0 · Tier 1）契约测试。

覆盖：
  · 声明词表 / 默认 = 支持集 S0
  · 闸门语义：在 S0 内 → SUPPORTED(exit 0)；合法但超集 → DESIGN-NOT-SUPPORTED(exit 3)；
    非法取值（写错）→ 错误(exit 2)
  · 结构提示只警告、不改判定
  · CLI 集成：默认路径逐字节零影响；拒绝路径不产出任何判定产物
"""
import json
import os

import numpy as np
import pandas as pd
import pytest

import sc_audit
from sc_audit.design import (DesignDeclaration, DesignGate, check_design_support,
                             structural_hints, EXIT_DESIGN_NOT_SUPPORTED,
                             COMPARISON_TYPES, UNITS)
from conftest import make_adata


# ──────────────────────── ① 纯函数：声明与闸门 ────────────────────────

def test_default_declaration_is_supported_set():
    d = DesignDeclaration()
    assert d.comparison_type == "two_arm"
    assert d.unit == "donor"
    assert d.time_col is None
    assert d.is_default() is True
    g = check_design_support(d)
    assert isinstance(g, DesignGate)
    assert g.supported is True
    assert g.status == "SUPPORTED"
    assert g.exit_code == 0
    assert g.unsupported_reasons == []


@pytest.mark.parametrize("kw", [
    {"comparison_type": "paired"},
    {"comparison_type": "time_course"},
    {"comparison_type": "multi_arm"},
    {"comparison_type": "multi_state"},
    {"unit": "donor_timepoint"},
    {"unit": "sample"},
    {"unit": "cell"},
    {"time_col": "timepoint"},
])
def test_unsupported_declarations_are_refused(kw):
    """合法但超出支持集 → DESIGN-NOT-SUPPORTED + exit 3（不是 2）。"""
    d = DesignDeclaration(**kw)
    assert d.invalid_reasons() == []          # 合法
    assert d.is_default() is False
    g = check_design_support(d)
    assert g.supported is False
    assert g.status == "DESIGN-NOT-SUPPORTED"
    assert g.exit_code == EXIT_DESIGN_NOT_SUPPORTED == 3
    assert g.unsupported_reasons
    # 拒绝文案必须点明"与 INDETERMINATE 不同"
    assert "INDETERMINATE" in g.refusal_text()
    assert "DESIGN-NOT-SUPPORTED" in g.refusal_text()


@pytest.mark.parametrize("kw", [
    {"comparison_type": "two-arm"},     # 连字符拼错
    {"comparison_type": "twostage"},
    {"unit": "donors"},
    {"unit": ""},
])
def test_invalid_declarations_are_errors_not_unsupported(kw):
    """取值不在封闭词表内 = 调用方写错 → 非法（CLI 侧 exit 2），≠ 不支持。"""
    d = DesignDeclaration(**kw)
    assert d.invalid_reasons()             # 非空 = 非法


def test_gate_is_deterministic():
    d = DesignDeclaration(comparison_type="time_course")
    a = check_design_support(d)
    b = check_design_support(d)
    assert a.to_dict() == b.to_dict()


# ──────────────────────── ② 结构提示（只警告，不改判定） ────────────────────────

def test_structural_hint_on_multilevel_condition(adata_clean):
    ad = adata_clean.copy()
    # 造一个三水平 condition
    ad.obs["condition"] = ad.obs["condition"].astype(str)
    first_case = ad.obs.loc[ad.obs["condition"] == "disease", "donor"].iloc[0]
    ad.obs.loc[ad.obs["donor"] == first_case, "condition"] = "disease_subtype"
    hints = structural_hints(ad, "condition", "donor")
    assert any("个水平" in h for h in hints), hints


def test_no_structural_hint_on_two_level_condition(adata_clean):
    hints = structural_hints(adata_clean, "condition", "donor")
    assert hints == []


def test_structural_hint_does_not_change_gate(adata_clean):
    """提示是旁路信息：不因水平数 > 2 而把 two_arm 判成不支持。"""
    g = check_design_support(DesignDeclaration())   # two_arm / donor
    assert g.supported is True


# ──────────────────────── ③ CLI 集成 ────────────────────────

def _write_case(tmp_path, **kw):
    ad = make_adata(**kw)
    h5 = tmp_path / "case.h5ad"
    ad.write_h5ad(h5)
    return str(h5)


def _run_cli(args):
    from sc_audit.cli import main as cli_main
    return cli_main(args)


def _base_args(h5, outdir):
    return ["run", "--input", h5, "--target-genes", "KCNJ2,TNNT2,ACTB,GAPDH",
            "--outdir", str(outdir), "--no-count-check"]


def test_cli_default_path_is_supported_and_writes_metadata(tmp_path):
    h5 = _write_case(tmp_path)
    out = tmp_path / "out"
    rc = _run_cli(_base_args(h5, out))
    assert rc == 0
    meta = json.load(open(out / "run_metadata.json", encoding="utf-8"))
    assert meta["design"]["supported"] is True
    assert meta["design"]["status"] == "SUPPORTED"
    assert meta["design"]["declared"]["comparison_type"] == "two_arm"
    assert tuple(int(x) for x in meta["sc_audit_version"].split(".")[:2]) >= (0, 4)
    assert (out / "floor_audit.csv").exists()


def test_cli_default_matches_explicit_two_arm(tmp_path):
    """显式声明 S0 == 不声明（零影响）：判定产物逐字节一致。"""
    h5 = _write_case(tmp_path)
    out_a = tmp_path / "a"
    out_b = tmp_path / "b"
    assert _run_cli(_base_args(h5, out_a)) == 0
    assert _run_cli(_base_args(h5, out_b) + ["--comparison-type", "two_arm",
                                             "--unit", "donor"]) == 0
    ta = (out_a / "floor_audit.csv").read_bytes()
    tb = (out_b / "floor_audit.csv").read_bytes()
    assert ta == tb


@pytest.mark.parametrize("extra,codes", [
    (["--comparison-type", "time_course"], 3),
    (["--comparison-type", "multi_arm"], 3),
    (["--unit", "sample"], 3),
    (["--time-col", "timepoint"], 3),
])
def test_cli_unsupported_design_refuses_with_exit_3(tmp_path, extra, codes):
    h5 = _write_case(tmp_path)
    out = tmp_path / "out"
    with pytest.raises(SystemExit) as ei:
        _run_cli(_base_args(h5, out) + extra)
    assert ei.value.code == codes
    # 拒绝时**不得**产出任何判定产物
    assert not (out / "floor_audit.csv").exists()
    assert not (out / "identifiability.txt").exists()
    assert not (out / "run_metadata.json").exists()


def test_cli_check_also_gated(tmp_path):
    h5 = _write_case(tmp_path)
    out = tmp_path / "out"
    with pytest.raises(SystemExit) as ei:
        _run_cli(["check", "--input", h5, "--outdir", str(out),
                  "--comparison-type", "multi_state"])
    assert ei.value.code == 3
    assert not (out / "check_report.txt").exists()


def test_cli_invalid_enum_via_config_is_error_2(tmp_path):
    """配置文件里写非法取值（绕过 argparse choices）→ exit 2（错误，非不支持）。"""
    h5 = _write_case(tmp_path)
    out = tmp_path / "out"
    cfg = tmp_path / "bad.yaml"
    cfg.write_text("design_comparison_type: two-arm\n", encoding="utf-8")
    with pytest.raises(SystemExit) as ei:
        _run_cli(_base_args(h5, out) + ["--config", str(cfg)])
    assert ei.value.code == 2


def test_cli_multilevel_condition_still_runs_with_hint(tmp_path, capsys):
    """结构提示只警告：三水平 condition 在 two_arm 下仍照常出判定（exit 0）。"""
    ad = make_adata()
    ad.obs["condition"] = ad.obs["condition"].astype(str)
    first_case = ad.obs.loc[ad.obs["condition"] == "disease", "donor"].iloc[0]
    ad.obs.loc[ad.obs["donor"] == first_case, "condition"] = "disease_subtype"
    h5 = tmp_path / "ml.h5ad"
    ad.write_h5ad(h5)
    out = tmp_path / "out"
    rc = _run_cli(_base_args(str(h5), out))
    assert rc == 0
    assert (out / "floor_audit.csv").exists()
    err = capsys.readouterr().err
    assert "design-hint" in err
