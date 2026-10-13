"""参照契约（v0.5.0）：R1 五元组合并守卫 + R2 跨 modality 硬约束。

对应 `docs/REFERENCE_SCHEMA_v1.md` 第 III.2 的机器强制规则：
  R1  合并参照列时，五元组 (modality, intron_counting, gate_spec_hash,
      upstream_qc_hash, cell_source) 必须完全一致，否则拒绝合并并列出差异字段。
  R2  跨 modality 默认禁止，需 --allow-cross-modality 显式放行；
      放行后受影响基因的 R 标 cross_modality=True。
"""
import json
import os

import pytest

from sc_audit.config import AuditConfig
from sc_audit.floor import audit_detection_floor
from sc_audit.reference import (CrossModalityError, ReferenceContractError,
                                REFERENCE_FIVE_TUPLE, validate_reference_contract,
                                enforce_reference_contract, reference_modality_set)
from conftest import make_adata


# ──────────────────────── R2 · 跨 modality ────────────────────────

def _tiny_ad():
    A = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=8)
    A.obs['condition'] = A.obs['condition'].astype(str)
    return A


def test_r2_cross_modality_forbidden_by_default():
    A = _tiny_ad()
    with pytest.raises(CrossModalityError):
        audit_detection_floor(A, ['TNNT2'], 'condition', 'donor', 'control',
                              cfg=AuditConfig(), reference={'TNNT2': 90.0},
                              modality='scRNA', reference_modalities=['snRNA'])


def test_r2_allowed_marks_cross_modality_column():
    A = _tiny_ad()
    with pytest.warns(UserWarning, match='跨模态'):
        t = audit_detection_floor(A, ['TNNT2'], 'condition', 'donor', 'control',
                                  cfg=AuditConfig(), reference={'TNNT2': 90.0},
                                  modality='scRNA', reference_modalities=['snRNA'],
                                  allow_cross_modality=True).table
    assert t.iloc[0]['ref_modality_flag'] == 'cross-modal'
    assert bool(t.iloc[0]['cross_modality']) is True


def test_r2_same_modality_passes():
    A = _tiny_ad()
    t = audit_detection_floor(A, ['TNNT2'], 'condition', 'donor', 'control',
                              cfg=AuditConfig(), reference={'TNNT2': 90.0},
                              modality='scRNA',
                              reference_modalities=['scRNA', 'snRNA']).table
    assert t.iloc[0]['ref_modality_flag'] == 'same-modality'
    assert 'cross_modality' not in t.columns      # 非跨模态不加该列（保持默认输出不变）


def test_r2_unknown_when_undeclared():
    A = _tiny_ad()
    t = audit_detection_floor(A, ['TNNT2'], 'condition', 'donor', 'control',
                              cfg=AuditConfig(), reference={'TNNT2': 90.0}).table
    assert 'unknown' in str(t.iloc[0]['ref_modality_flag'])


# ──────────────────────── R1 · 五元组合并 ────────────────────────

def _write_ref(tmp_path, modality=None):
    p = tmp_path / "ref.csv"
    p.write_text("gene,pct_ref\nTNNT2,90.0\nACTB,95.0\n", encoding="utf-8")
    return str(p)


def _base_five(**over):
    d = {"modality": "scRNA", "intron_counting": "unknown",
         "gate_spec_hash": "aaa", "upstream_qc_hash": "bbb", "cell_source": "author_label"}
    d.update(over)
    return d


def test_r1_identical_five_tuple_passes(tmp_path):
    ref = _write_ref(tmp_path)
    man = {"column_contract": {
        "pct_ref_A": _base_five(), "pct_ref_B": _base_five()}}
    open(ref + ".manifest.json", "w").write(json.dumps(man))
    c = validate_reference_contract(ref, target_modality="scRNA")
    assert c.r1_enforced is True
    assert c.r1_violations == []
    enforce_reference_contract(c)      # 不抛


def test_r1_mismatch_refuses_and_lists_fields(tmp_path):
    ref = _write_ref(tmp_path)
    man = {"column_contract": {
        "pct_ref_A": _base_five(modality="scRNA"),
        "pct_ref_B": _base_five(modality="snRNA", gate_spec_hash="zzz")}}
    open(ref + ".manifest.json", "w").write(json.dumps(man))
    c = validate_reference_contract(ref, target_modality="scRNA")
    assert c.r1_enforced is True and c.r1_violations
    with pytest.raises(ReferenceContractError) as ei:
        enforce_reference_contract(c)
    msg = str(ei.value)
    assert "pct_ref_B" in msg and "modality" in msg and "gate_spec_hash" in msg


def test_r1_not_enforced_without_manifest_is_registered(tmp_path):
    ref = _write_ref(tmp_path)
    c = validate_reference_contract(ref, target_modality="scRNA")
    assert c.r1_enforced is False
    assert any("R1 未强制" in n for n in c.notes)     # 不静默
    enforce_reference_contract(c)                    # 无违规 → 不抛


def test_five_tuple_fields_are_frozen():
    assert REFERENCE_FIVE_TUPLE == ("modality", "intron_counting", "gate_spec_hash",
                                    "upstream_qc_hash", "cell_source")


def test_reference_modality_set_reads_csv_and_manifest(tmp_path):
    p = tmp_path / "r.csv"
    p.write_text("gene,pct_ref,modality\nTNNT2,90.0,scRNA\n", encoding="utf-8")
    assert reference_modality_set(str(p)) == ["scRNA"]
    open(str(p) + ".manifest.json", "w").write(
        json.dumps({"column_contract": {"pct_ref": _base_five(modality="snRNA")}}))
    assert set(reference_modality_set(str(p))) == {"scRNA", "snRNA"}


# ──────────────────────── CLI 层 ────────────────────────

def _case_h5(tmp_path):
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=10)
    h5 = tmp_path / "case.h5ad"
    ad.write_h5ad(h5)
    return str(h5)


def _cli(args):
    from sc_audit.cli import main as cli_main
    return cli_main(args)


def test_cli_cross_modality_exit2_without_flag(tmp_path):
    h5 = _case_h5(tmp_path)
    ref = tmp_path / "ref.csv"
    ref.write_text("gene,pct_ref,modality\nTNNT2,90.0,snRNA\nACTB,95.0,snRNA\n",
                   encoding="utf-8")
    with pytest.raises(SystemExit) as ei:
        _cli(["floor", "--input", h5, "--target-genes", "TNNT2,ACTB",
              "--reference", str(ref), "--modality", "scRNA",
              "--outdir", str(tmp_path / "o1"), "--no-count-check"])
    assert ei.value.code == 2


def test_cli_cross_modality_ok_with_flag(tmp_path):
    h5 = _case_h5(tmp_path)
    ref = tmp_path / "ref.csv"
    ref.write_text("gene,pct_ref,modality\nTNNT2,90.0,snRNA\nACTB,95.0,snRNA\n",
                   encoding="utf-8")
    out = tmp_path / "o2"
    rc = _cli(["floor", "--input", h5, "--target-genes", "TNNT2,ACTB",
               "--reference", str(ref), "--modality", "scRNA",
               "--allow-cross-modality", "--outdir", str(out), "--no-count-check"])
    assert rc == 0
    t = __import__("pandas").read_csv(out / "floor_audit.csv")
    assert "cross_modality" in t.columns
    assert bool(t["cross_modality"].iloc[0]) is True


def test_cli_r1_violation_exit2(tmp_path):
    h5 = _case_h5(tmp_path)
    ref = tmp_path / "ref.csv"
    ref.write_text("gene,pct_ref,modality\nTNNT2,90.0,scRNA\nACTB,95.0,scRNA\n",
                   encoding="utf-8")
    open(str(ref) + ".manifest.json", "w").write(json.dumps({"column_contract": {
        "pct_ref_A": _base_five(), "pct_ref_B": _base_five(gate_spec_hash="zzz")}}))
    with pytest.raises(SystemExit) as ei:
        _cli(["floor", "--input", h5, "--target-genes", "TNNT2,ACTB",
              "--reference", str(ref), "--modality", "scRNA",
              "--outdir", str(tmp_path / "o3"), "--no-count-check"])
    assert ei.value.code == 2
