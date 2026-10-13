"""审计四（细胞身份与捕获完整性）的单元测试。

构造带**已知真值**的合成数据，逐项验证八项检查的判定，并钉住三个关键失败模式：
  · 缺注释列 → INDETERMINATE（不静默通过）
  · 亚群只在单臂出现 → 必须检出
  · 门控规则与条件共线 → 必须检出
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import anndata as ad

from sc_audit.config import AuditConfig
from sc_audit.identity import audit_cell_identity, _detect_many, _tvd

POS = ['TTN', 'MYH7', 'MYBPC3', 'ACTN2', 'TNNI3', 'MYL2']
NEG = ['PTPRC', 'CD3E', 'DCN', 'COL1A1', 'PECAM1', 'VWF', 'LYZ']
EXTRA = ['MALAT1', 'RYR2', 'CACNA1C', 'KCNQ1', 'SCN5A']
GENES = POS + NEG + EXTRA


def _build(n_ctrl_don=3, n_case_don=3, cells=60, seed=0,
           cm_pct=0.5, other_lineages=('Endothelial cells', 'Myeloid cells',
                                       'Fibroblasts', 'T cells'),
           cm_complexity=1500, other_complexity=900,
           source_imbalance=False, collinear_neg=False,
           malat1_cm=0.9, include_annotation=True):
    """构造合成数据集。

    cm_pct            心肌细胞占比（门控通过后的比例由 pos/neg 决定）
    source_imbalance  True → 某来源的亚群只出现在对照组（复现"来源只在单臂出现"的结构）
    collinear_neg     True → 阴性标志在病例组升高（复现门控共线）
    """
    rng = np.random.default_rng(seed)
    obs, rows = [], []
    for arm, nd in (('control', n_ctrl_don), ('disease', n_case_don)):
        for d in range(nd):
            don = f"{arm[0]}{d}"
            for i in range(cells):
                is_cm = rng.random() < cm_pct
                lin = 'Cardiomyocytes' if is_cm else rng.choice(other_lineages)
                src = 'Cells'
                if source_imbalance and arm == 'control' and is_cm and rng.random() < 0.6:
                    src = 'CD45+'
                ng = int(rng.normal(cm_complexity if is_cm else other_complexity, 120))
                ng = max(ng, 150)
                obs.append({'condition': arm, 'donor': don, 'lineage_names': lin,
                            'source': src, 'n_genes': ng, '_cm': is_cm, '_arm': arm})
                rows.append(None)
    n = len(obs)
    X = np.zeros((n, len(GENES)), dtype=np.int32)
    for i, ob in enumerate(obs):
        cm = ob['_cm']
        for j, g in enumerate(GENES):
            if g in POS:
                p = 0.85 if cm else 0.05
            elif g in EXTRA:
                p = malat1_cm if (g == 'MALAT1' and cm) else (0.7 if cm else 0.7)
            else:
                p = 0.02
                if collinear_neg and g in ('PECAM1', 'VWF') and ob['_arm'] == 'disease':
                    p = 0.45
            # 高频基因给足计数，低频基因给 1
            X[i, j] = rng.poisson(3) + 1 if rng.random() < p else 0
    A = ad.AnnData(X=sp.csr_matrix(X),
                   obs=pd.DataFrame(obs, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=GENES))
    if not include_annotation:
        A.obs = A.obs.drop(columns=['lineage_names'])
    return A


def _cfg(**kw):
    c = AuditConfig()
    c.identity_col = 'lineage_names'
    c.identity_source_col = 'source'
    c.identity_target = 'Cardiomyocytes'
    c.identity_complexity_col = 'n_genes'
    for k, v in kw.items():
        setattr(c, k, v)
    return c


# ══════════════════ 基础 ══════════════════

def test_missing_annotation_is_indeterminate_not_silent_pass():
    """缺注释列必须 INDETERMINATE —— 绝不能静默返回 PASS。"""
    A = _build(include_annotation=False)
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    assert r.verdict.startswith('INDETERMINATE')
    assert any('lineage' in n for n in r.notes)


def test_clean_data_passes():
    """无任何缺陷时判 PASS。"""
    A = _build(source_imbalance=False, collinear_neg=False)
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    assert r.verdict.startswith('PASS'), r.notes


# ══════════════════ 检查 1/2/3 ══════════════════

def test_gate_margin_detects_lineage_specificity():
    A = _build()
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    assert r.gate_by_lineage is not None
    assert r.target_gate_pass > 80          # CM 应几乎全通过
    assert r.gate_margin > 5                # 且远高于其他谱系


def test_complexity_rank_is_one_when_target_is_lowest():
    """目标谱系的复杂度最低时必须排第 1。"""
    A = _build(cm_complexity=400, other_complexity=1200)
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    assert r.complexity_rank == 1
    assert any('LOWEST median genes/cell' in n for n in r.notes)


def test_nuclear_marker_rank_flags_nuclear_poor_capture():
    A = _build(malat1_cm=0.05)
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    assert r.nuclear_marker_rank == 1
    assert any('nuclear' in n for n in r.notes)


# ══════════════════ 检查 5（★ 亚群不平衡）══════════════════

def test_subset_imbalance_is_detected():
    """复现"某来源只出现在一个臂"的结构 —— 必须检出。"""
    A = _build(source_imbalance=True)
    A.obs['source'] = np.where(A.obs['source'].values == 'CD45+', 'CD45+', 'Cells')
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    assert r.imbalanced_subsets, "亚群不平衡未被检出"
    assert any('CD45+' in s for s in r.imbalanced_subsets)
    assert r.verdict.startswith('COMPOSITION-CONFOUNDED')


def test_no_imbalance_when_source_is_balanced():
    A = _build(source_imbalance=False)
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    assert not r.imbalanced_subsets


# ══════════════════ 检查 6（★ 门控共线）══════════════════

def test_gate_collinearity_is_detected():
    A = _build(collinear_neg=True)
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    assert r.collinear_markers, "门控共线未被检出"
    assert set(r.collinear_markers) & {'PECAM1', 'VWF'}
    assert any('condition-collinear' in n for n in r.notes)


def test_gate_pass_by_arm_reported():
    A = _build()
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    assert r.gate_pass_by_arm is not None
    assert set(r.gate_pass_by_arm.index) == {'control', 'case'}


# ══════════════════ 检查 7（模态）══════════════════

def test_cross_modality_is_flagged():
    A = _build()
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg(),
                            reference_modalities=['snRNA'] * 7, modality='scRNA')
    assert r.modality_matched is False
    assert any('cross-modal' in n for n in r.notes)


def test_same_modality_not_flagged():
    A = _build()
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg(),
                            reference_modalities=['snRNA', 'scRNA'], modality='scRNA')
    assert r.modality_matched is True


# ══════════════════ 边界与工具 ══════════════════

def test_missing_source_column_degrades_loudly():
    """缺 source 列时检查 5 跳过，但必须留下 note（不得静默）。"""
    A = _build()
    A.obs = A.obs.drop(columns=['source'])
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    assert any('source' in n for n in r.notes)


def test_single_lineage_degrades_but_still_runs_core_checks():
    """单谱系应只降级跨谱系对照（检查 1/2/3），**不得**连带关闭
    亚群不平衡与门控共线性这两个核心项。"""
    A = _build(source_imbalance=True)
    A.obs['lineage_names'] = 'Cardiomyocytes'
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    # note 已改为显式可评估性清单：
    # 逐项声明哪些检查降级、哪些完整评估，防止读者误读为"完整审计通过"
    assert any('ASSESSABILITY' in n for n in r.notes), "未声明可评估性"
    assert any('ASSESSABILITY (this run)' in n for n in r.notes), "未给出动态可评估清单"
    assert any('NOT evaluable' in n for n in r.notes), "未声明哪些检查未评估"
    # C3（v0.7.1）：硬编码的"其余检查全部已评估"断言必须已删除（会与实际跳过的检查矛盾）
    assert not any('FULLY EVALUATED' in n for n in r.notes), \
        "硬编码 FULLY EVALUATED 文案应已删除"
    assert r.gate_margin is None, "单谱系不应给出 margin"
    assert r.imbalanced_subsets, "单谱系下亚群不平衡检查被错误关闭"
    assert r.gate_pass_by_arm is not None, "单谱系下门控通过率检查被错误关闭"


def test_missing_positive_markers_is_indeterminate():
    A = _build()
    keep = [g for g in A.var_names if g not in POS]
    A = A[:, keep].copy()
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    assert r.verdict.startswith('INDETERMINATE')


def test_tvd_helper():
    a = pd.Series([10, 0, 0], index=['x', 'y', 'z'])
    b = pd.Series([0, 10, 0], index=['x', 'y', 'z'])
    assert _tvd(a, b) == pytest.approx(1.0)
    assert _tvd(a, a) == pytest.approx(0.0)


def test_detect_many_returns_only_present_genes():
    A = _build()
    d = _detect_many(A, ['TTN', 'NOT_A_GENE_XYZ'])
    assert 'TTN' in d and 'NOT_A_GENE_XYZ' not in d


def test_to_text_contains_verdict():
    A = _build(source_imbalance=True)
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=_cfg())
    t = r.to_text()
    assert 'Audit 4' in t and r.verdict in t


def test_check6_target_lineage_absent_from_one_arm_is_reported_not_nan():
    """目标谱系在某一臂**完全缺席**时，检查 6 必须显式说明"absent from arm"（v0.3.3）。

    历史行为：对空集合求均值 → `RuntimeWarning: Mean of empty slice` + 静默 NaN，
    读者反而看不出"该谱系在某一臂不存在"这一结构性问题（正是审计四的核心用例）。
    """
    A = _build()
    # 把目标谱系的细胞**全部**放进对照组：病例组该谱系为 0 个细胞
    lin = A.obs['lineage_names'].astype(str).copy()
    tgt = 'Cardiomyocytes' if 'Cardiomyocytes' in set(lin) else sorted(set(lin))[0]
    case_idx = A.obs.index[A.obs['condition'] != 'control']
    A.obs.loc[case_idx, 'lineage_names'] = 'NOT_TARGET'
    cfg = _cfg()
    cfg.identity_target = tgt
    r = audit_cell_identity(A, 'condition', 'donor', 'control', cfg=cfg)
    joined = " ".join(r.notes or [])
    assert "完全缺席" in joined, f"未报'目标谱系在单臂缺席'：{r.notes}"
    assert any('absent' in n or '完全缺席' in n for n in (r.notes or []))
    # 且不得留下 NaN 的逐臂行
    if r.negative_marker_by_arm is not None and len(r.negative_marker_by_arm):
        assert not r.negative_marker_by_arm.isna().any().any(), \
            f"逐臂表里仍有 NaN：{r.negative_marker_by_arm}"


def test_identifiability_rejects_matched_donors_and_run_accepts_ladder(tmp_path):
    """`identifiability --matched-donors` 必须显式拒绝；`run --ladder` 必须可用并写入元数据（v0.3.3）。"""
    import json
    from sc_audit.cli import main as cli_main
    from conftest import make_adata
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=20)
    h5 = tmp_path / "d.h5ad"
    ad.write_h5ad(h5)
    with pytest.raises(SystemExit) as e:
        cli_main(["identifiability", "--input", str(h5), "--target-genes", "TNNT2",
                  "--matched-donors", "ctrl_d0", "--outdir", str(tmp_path / "o1"),
                  "--quiet"])
    assert e.value.code == 2
    # run 接受 --ladder 并写入 run_metadata（含 perm_mode 留痕）
    out = tmp_path / "o2"
    rc = cli_main(["run", "--input", str(h5), "--target-genes", "TNNT2,ACTB",
                   "--ladder", "none;institution", "--n-perm", "100",
                   "--outdir", str(out), "--quiet", "--no-count-check"])
    assert rc == 0
    meta = json.loads((out / "run_metadata.json").read_text())
    assert meta["identifiability_ladder"] == [[], ["institution"]]
    assert "perm_mode" in meta


def test_check5_note_attribution_and_missing_target():
    """多谱系 + 未指定目标谱系时，检查 5 的说明必须是「未指定目标」，不是「obs 中无 source」（v0.3.3）。

    旧实现把两种完全不同的原因混进同一个 else 分支，于是即使 source 列**存在**，
    notes 也会写「obs 中无 'source'」——归因错误会误导使用者去改数据，而不是指定目标。
    """
    from sc_audit.identity import audit_cell_identity as _aud
    A = _build()                       # _build 默认含 source 列
    cfg = _cfg()
    A.obs['lineage_names'] = ['CM', 'FIB', 'EC'] * (A.n_obs // 3) + \
        ['CM'] * (A.n_obs - 3 * (A.n_obs // 3))
    cfg.identity_target = ''           # 显式空：多谱系下目标未定
    r = _aud(A, 'condition', 'donor', 'control', cfg=cfg)
    notes = ' '.join(r.notes or [])
    assert "obs 中无 'source'" not in notes, (
        '归因错误：source 列存在，却把检查 5 跳过归因于缺列：' + notes)
    assert ('未指定目标谱系' in notes) or ('未指定' in notes), (
        '未显式说明「未指定目标谱系」：' + notes)
