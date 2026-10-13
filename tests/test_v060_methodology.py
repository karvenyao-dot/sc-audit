"""v0.6.0 方法学修订（A/B/C）契约测试。"""
import numpy as np
import pandas as pd
import pytest

from sc_audit.config import AuditConfig
from sc_audit.floor import audit_detection_floor
from sc_audit.reference import DepthMismatchError
from sc_audit.identifiability import audit_identifiability, positivity_overlap
from sc_audit.robust import audit_robust_signal
from sc_audit.identity import audit_cell_identity
from sc_audit.latent_batch import audit_latent_batch
from conftest import make_adata


def _ad(nc=6, nk=6, cpd=20, **kw):
    return make_adata(n_ctrl_donors=nc, n_case_donors=nk, cells_per_donor=cpd, **kw)


# ──────────────────────── A1 参照深度可比性 ────────────────────────

def test_a1_depth_mismatch_raises_by_default():
    A = _ad()
    with pytest.raises(DepthMismatchError):
        audit_detection_floor(A, ['TNNT2'], 'condition', 'donor', 'control',
                              cfg=AuditConfig(), reference={'TNNT2': 90.0},
                              reference_depth=1000.0, case_depth=3000.0)


def test_a1_depth_mismatch_allowed_adds_columns():
    A = _ad()
    with pytest.warns(UserWarning, match='REFERENCE-DEPTH-MISMATCH'):
        t = audit_detection_floor(A, ['TNNT2'], 'condition', 'donor', 'control',
                                  cfg=AuditConfig(), reference={'TNNT2': 90.0},
                                  reference_depth=1000.0, case_depth=3000.0,
                                  allow_depth_mismatch=True).table
    for c in ("case_depth", "ref_depth", "depth_delta_pct"):
        assert c in t.columns
    assert abs(float(t.iloc[0]['depth_delta_pct']) - 200.0) < 0.1


def test_a1_depth_match_ok_no_raise():
    A = _ad()
    t = audit_detection_floor(A, ['TNNT2'], 'condition', 'donor', 'control',
                              cfg=AuditConfig(), reference={'TNNT2': 90.0},
                              reference_depth=1000.0, case_depth=1100.0).table
    assert float(t.iloc[0]['depth_delta_pct']) == 10.0


def test_a1_no_depth_info_no_columns():
    """未提供参照深度 → 不加列（保持输出不变）。"""
    A = _ad()
    t = audit_detection_floor(A, ['TNNT2'], 'condition', 'donor', 'control',
                              cfg=AuditConfig(), reference={'TNNT2': 90.0}).table
    assert 'depth_delta_pct' not in t.columns


def test_a1_depth_matched_sensitivity_columns():
    A = _ad()
    t = audit_detection_floor(A, ['TNNT2'], 'condition', 'donor', 'control',
                              cfg=AuditConfig(), reference={'TNNT2': 90.0},
                              reference_depth=1000.0, case_depth=3000.0,
                              allow_depth_mismatch=True,
                              depth_matched_sensitivity=True).table
    assert 'case_pct_depth_matched' in t.columns
    assert 'verdict_depth_matched' in t.columns
    # 下采样（深度更低）后检出率应不高于原始
    if t.iloc[0]['case_pct_depth_matched'] == t.iloc[0]['case_pct_depth_matched']:
        assert t.iloc[0]['case_pct_depth_matched'] <= t.iloc[0]['case_pct'] + 1e-9


# ──────────────────────── A2 positivity / 协变量重叠 ────────────────────────

def _separated_ad():
    """病例全在 INST_A；对照大多在 INST_B，仅 1 个在 INST_A（≥50 细胞）。"""
    A = make_adata(n_ctrl_donors=4, n_case_donors=3, cells_per_donor=60,
                   design_confounded=False)
    ctrl = [f"ctrl_d{i}" for i in range(4)]
    for i, d in enumerate(ctrl):
        A.obs.loc[A.obs['donor'] == d, 'institution'] = 'INST_A' if i == 0 else 'INST_B'
    A.obs.loc[A.obs['condition'] == 'disease', 'institution'] = 'INST_A'
    return A


def test_a2_positivity_overlap_function():
    A = _separated_ad()
    p = positivity_overlap(A, ['institution'], 'condition', 'donor', 'control')
    assert p['overlap_case'] == 1.0          # 病例剖面在对照组中出现（1 个供者）
    assert p['overlap_ctrl'] == pytest.approx(0.25)   # 仅 1/4 对照有对应
    assert p['overlap_min'] == pytest.approx(0.25)


def test_a2_low_overlap_emits_mandatory_caveat():
    """A2 作并列诊断：低重叠时给强制警示（不覆盖阶梯判定）。"""
    A = _separated_ad()
    res = audit_identifiability(A, ['TNNT2'], 'condition', 'donor', 'control',
                                cfg=AuditConfig())
    ov = res.details.get('positivity_overlap', {})
    assert ov.get('overlap_min') is not None and ov['overlap_min'] < 0.5
    assert any("POSITIVITY-OVERLAP CAVEAT" in n for n in res.notes)
    # 不得因 positivity 就改写成 NOT SEPARABLE（阶梯判定保留，避免遮蔽衰减信息）
    assert "POSITIVITY" not in res.verdict


def test_a2_balanced_overlap_not_triggered():
    A = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=40,
                   design_confounded=False)   # institution 全 INST_A → 完全重叠
    res = audit_identifiability(A, ['TNNT2'], 'condition', 'donor', 'control',
                                cfg=AuditConfig())
    ov = res.details.get('positivity_overlap', {})
    assert ov.get('overlap_min') in (None, 1.0)


# ──────────────────────── B1 潜在批次层 ────────────────────────

def test_b1_runs_and_reports():
    A = _ad(nc=6, nk=6, cpd=25)
    r = audit_latent_batch(A, 'condition', 'donor', 'control', cfg=AuditConfig())
    assert r.verdict in ("LATENT-BATCH-RISK", "NO-LATENT-BATCH-SIGNAL", "INDETERMINATE")
    assert isinstance(r.to_text(), str) and "latent-batch scan" in r.to_text()


def test_b1_too_few_donors_is_indeterminate():
    A = _ad(nc=2, nk=2, cpd=20)
    r = audit_latent_batch(A, 'condition', 'donor', 'control', cfg=AuditConfig())
    assert r.verdict == "INDETERMINATE"


# ──────────────────────── B2 / B3 稳健信号 ────────────────────────

def test_b2_housekeeping_cv_and_applicability():
    A = _ad(nc=6, nk=6, cpd=40)
    ref = {'ACTB': 40.0, 'GAPDH': 95.0}      # 管家基线差异大
    r = audit_robust_signal(A, ['TNNT2'], 'condition', 'donor', 'control',
                            housekeeping=['ACTB', 'GAPDH'], reference=ref,
                            cfg=AuditConfig())
    assert r.housekeeping_cv == r.housekeeping_cv   # 非 NaN（>=2 管家）


def test_b3_depth_metric_sensitivity_dict():
    A = _ad(nc=6, nk=6, cpd=40)
    r = audit_robust_signal(A, ['TNNT2'], 'condition', 'donor', 'control',
                            housekeeping=['ACTB', 'GAPDH'], reference={'TNNT2': 80.0},
                            cfg=AuditConfig())
    d = r.depth_metric_sensitivity
    assert d.get('primary_metric') == 'n_genes'
    assert d.get('alt_metric') == 'total_counts'
    assert isinstance(d.get('consistent'), bool)


# ──────────────────────── C1–C5 透明项 ────────────────────────

def test_c4_c5_floor_summary_notes():
    A = _ad()
    r = audit_detection_floor(A, ['TNNT2'], 'condition', 'donor', 'control',
                              cfg=AuditConfig(), reference={'TNNT2': 90.0})
    assert "CELL-LEVEL" in r.summary['binom_p_note']
    assert "BH-FDR on a small target panel" in r.summary['fdr_small_panel_note']


def _ad_with_lineage():
    A = _ad(nc=5, nk=5, cpd=40)
    lin = np.where(A.obs['donor'].str.startswith('ctrl'), 'CM', 'CM')
    A.obs['lineage_names'] = lin
    A.obs.loc[A.obs.index[: int(len(A.obs) * 0.3)], 'lineage_names'] = 'FB'
    A.obs['predicted_doublet'] = False
    A.obs.loc[A.obs.index[: int(len(A.obs) * 0.05)], 'predicted_doublet'] = True
    return A


def test_c1_c2_c3_identity_notes():
    A = _ad_with_lineage()
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=AuditConfig())
    joined = " ".join(r.notes)
    assert "C1（ambient RNA）" in joined
    assert "C2（doublets）" in joined
    assert "C3（annotation dependence）" in joined


# ──────────────────────── v0.6.1 必修项（H1/H2/H3）────────────────────────

def test_h1_depth_metric_primary_field_removed():
    """H1：死配置字段已删（定义未用 → 移除，B3 仅作敏感性）。"""
    assert 'depth_metric_primary' not in AuditConfig().__dataclass_fields__


def test_h2_c2_has_no_silent_pass():
    """H2：C2 段不得静默吞错（应给显式 note）。"""
    import inspect
    import sc_audit.identity as I
    src = inspect.getsource(I)
    assert "解析失败" in src
    assert "except Exception:  # noqa: BLE001\n                pass" not in src


def test_h3_non_integer_counts_is_indeterminate():
    """H3：X 非整数计数 → 扫描不可解释，返回 INDETERMINATE 并声明。"""
    A = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=20, scale_x=True)
    r = audit_latent_batch(A, 'condition', 'donor', 'control', cfg=AuditConfig())
    assert r.verdict == "INDETERMINATE"
    assert any("integer counts" in n for n in r.notes)


def test_h3_oversize_matrix_is_skipped(monkeypatch):
    """H3：供者×基因超上限 → 跳过并声明（不尝试稠密化）。"""
    import sc_audit.latent_batch as LB
    monkeypatch.setattr(LB, "MAX_DENSE_ELEMENTS", 1)
    A = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=10)
    r = LB.audit_latent_batch(A, 'condition', 'donor', 'control', cfg=AuditConfig())
    assert r.verdict == "INDETERMINATE"
    assert any("SKIPPED" in n for n in r.notes)


# ──────────────────────── v0.6.2 M 组（M1–M5）────────────────────────

def test_m1_positivity_includes_adjustable():
    """M1：positivity 剖面 = design_defining ∪ adjustable。"""
    A = _ad(nc=6, nk=6, cpd=20)
    A.obs['age_range'] = A.obs['donor'].map(
        {d: ('old' if i % 2 else 'young') for i, d in enumerate(sorted(set(A.obs['donor'])))})
    p = positivity_overlap(A, ['institution'], 'condition', 'donor', 'control',
                           adjustable=['age_range'])
    assert 'age_range' in p['profile_vars'] and 'institution' in p['profile_vars']


def test_m2_reference_depth_prefers_median(tmp_path):
    """M2：参照深度优先中位口径。"""
    from sc_audit.reference import load_reference_depth
    import json
    ref = tmp_path / "r.csv"
    ref.write_text("gene,pct_ref\nTNNT2,90\n", encoding="utf-8")
    open(str(ref) + ".manifest.json", "w").write(
        json.dumps({"mean_umi": 1000, "median_umi": 5000}))
    assert load_reference_depth(str(ref)) == 5000.0


def test_m3_non_integer_depth_sensitivity_is_declared():
    """M3：计数非整数时下采样灵敏度跳过 → 显式 note（不得静默 NaN）。"""
    A = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=20, scale_x=True)
    r = audit_detection_floor(A, ['TNNT2'], 'condition', 'donor', 'control',
                              cfg=AuditConfig(), reference={'TNNT2': 90.0},
                              reference_depth=1000.0, case_depth=3000.0,
                              allow_depth_mismatch=True,
                              depth_matched_sensitivity=True)
    assert "not integers" in str(r.summary.get('reference_depth_note'))


def test_m4_check_enforces_depth_mismatch(tmp_path):
    """M4：`check` 也强制深度可比性。"""
    from sc_audit.cli import main as cli_main
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=10)
    h5 = tmp_path / "c.h5ad"
    ad.write_h5ad(h5)
    ref = tmp_path / "ref.csv"
    ref.write_text("gene,pct_ref,median_umi\nTNNT2,90,1\n", encoding="utf-8")
    with pytest.raises(SystemExit) as ei:
        cli_main(["check", "--input", str(h5), "--reference", str(ref),
                  "--outdir", str(tmp_path / "o")])
    assert ei.value.code == 2


def test_m5_verdict_helper_matches_chain():
    """M5：判定链单一实现（_verdict_from）——抽查分支字符串。"""
    from sc_audit.floor import _verdict_from
    cfg = AuditConfig()
    assert _verdict_from(float('nan'), 1.0, cfg, True).startswith("INDETERMINATE")
    assert _verdict_from(1.0, 1.0, cfg, True) == "UNQUANTIFIABLE (detection floor)"
    assert _verdict_from(50.0, 0.1, cfg, True) == "DROPOUT-LIMITED"
    assert _verdict_from(50.0, float('nan'), cfg, False) == "INSUFFICIENT (no reference)"
    assert _verdict_from(50.0, 0.9, cfg, True) == "QUANTIFIABLE"


# ──────────────────────── v0.7.0 L 组（L1–L4）────────────────────────

def test_l1_housekeeping_high_baseline_in_config():
    cfg = AuditConfig()
    assert 'housekeeping_high_baseline' in AuditConfig().__dataclass_fields__
    assert cfg.housekeeping_high_baseline == 0.8


def test_l23_latent_multiseed_stability():
    """L2（离散度选 HVG）+ L3（多初始化）：ari_by_seed 应给出 >=2 个种子。"""
    A = _big_ad(nc=6, nk=6, cpd=25)
    r = audit_latent_batch(A, 'condition', 'donor', 'control', cfg=AuditConfig())
    assert len(r.ari_by_seed) >= 2
    assert r.verdict in ("LATENT-BATCH-RISK", "NO-LATENT-BATCH-SIGNAL")


def _big_ad(nc=6, nk=6, cpd=25, **kw):
    """在 make_adata 基础上补足基因数（latent 扫描需 >=50 基因）。"""
    import anndata as ad
    import scipy.sparse as sp
    A = make_adata(n_ctrl_donors=nc, n_case_donors=nk, cells_per_donor=cpd, **kw)
    X = A.X.toarray() if sp.issparse(A.X) else np.asarray(A.X)
    extra = np.random.default_rng(0).integers(0, 4, size=(A.n_obs, 80)).astype('float32')
    X = np.hstack([X, extra])
    names = list(map(str, A.var_names)) + [f"x{i}" for i in range(80)]
    return ad.AnnData(X=sp.csr_matrix(X), obs=A.obs.copy(),
                      var=pd.DataFrame(index=names))


def test_l4_report_renders_latent_section(tmp_path):
    from sc_audit.report import render_report
    A = _big_ad(nc=6, nk=6, cpd=25)
    r = audit_latent_batch(A, 'condition', 'donor', 'control', cfg=AuditConfig())
    p = render_report(latent_batch=r, outdir=str(tmp_path))
    html = open(p, encoding='utf-8').read()
    assert "Latent-batch scan" in html
