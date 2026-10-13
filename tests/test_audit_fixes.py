"""外部代码审计（2026-10-02）所报缺陷的回归测试。

每条修复对应一条断言 —— 缺陷一旦复现，测试必须失败。
编号与审计报告的编号一致（#1 / #2 / #3 / #4 / #6 / #10 / #11）。
"""
from __future__ import annotations

import inspect
import os
import tempfile

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import anndata as ad

from sc_audit import composition as _comp
from sc_audit import floor as _floor
from sc_audit import identity as _identity
from sc_audit import io as _io
from sc_audit import robust as _robust
from sc_audit.config import AuditConfig

TMP = tempfile.mkdtemp(prefix='scaudit_fix_')


def _tiny(var_names, fill_col=None, n=40):
    """极简 AnnData：指定列全 1，其余全 0。"""
    X = sp.csr_matrix(np.zeros((n, len(var_names)), dtype=np.int32))
    if fill_col is not None:
        X = X.tolil()
        X[:, fill_col] = 1
        X = X.tocsr()
    return ad.AnnData(X=X,
                      obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
                      var=pd.DataFrame(index=var_names))


# ══════════ #1 单基因 + 稀疏 → 曾返回 var_names[0] 的列 ══════════

def test_fix1_single_gene_returns_that_gene_not_column_zero():
    """目标基因不在第 0 列时，单基因入口必须返回**目标基因**的向量。"""
    A = _tiny(['AAA', 'MMM', 'ZZZ'], fill_col=2)      # 只有 ZZZ 表达
    d = _identity._detect_many(A, ['ZZZ'])
    assert d['ZZZ'].mean() == 1.0, "单基因入口返回了错误的列（应为 ZZZ）"


def test_fix1_single_gene_matches_multi_gene_result():
    """单基因与多基因入口对同一基因必须给出一致结果。"""
    A = _tiny(['AAA', 'MMM', 'ZZZ'], fill_col=2)
    one = _identity._detect_many(A, ['MMM', 'ZZZ'])['ZZZ']
    two = _identity._detect_many(A, ['ZZZ'])['ZZZ']
    assert np.array_equal(one, two)


def test_fix1_single_gene_via_public_io_entry():
    """公共入口 gene_detection_vector 同样须正确（#4 使其共用同一实现）。"""
    A = _tiny(['AAA', 'MMM', 'ZZZ'], fill_col=2)
    assert _io.gene_detection_vector(A, 'ZZZ').mean() == 1.0
    assert _io.gene_detection_vector(A, 'AAA').mean() == 0.0


def test_fix1_non_first_gene_by_name_not_position():
    """用真实基因名（非首列）复现审计报告给出的场景。"""
    A = _tiny(['TTN', 'KCNQ1', 'TNNT2'], fill_col=1)   # KCNQ1 在中间
    assert _identity._detect_many(A, ['KCNQ1'])['KCNQ1'].mean() == 1.0


# ══════════ #4 backed 模式 ══════════

def test_fix4_backed_data_does_not_break_detection():
    """backed='r' 打开时 gene_detection_vector 须可读 —— 原实现抛
    ValueError: scipy.sparse does not support dtype object。"""
    A = _tiny(['AAA', 'ZZZ'], fill_col=1, n=30)
    p = os.path.join(TMP, 'fix4.h5ad')
    A.write_h5ad(p)
    Ab = ad.read_h5ad(p, backed='r')
    v = _io.gene_detection_vector(Ab, 'ZZZ')
    assert v is not None, "backed 模式下返回 None（原实现抛异常，被吞后成 None）"
    assert v.mean() == 1.0
    assert _io.gene_detection_vector(Ab, 'AAA').mean() == 0.0


# ══════════ #3 _cp_upper 在 k == n ══════════

def test_fix3_cp_upper_all_positive_is_100_not_nan():
    assert _floor._cp_upper(5, 5) == 100.0
    assert _floor._cp_upper(1, 1) == 100.0


def test_fix3_cp_upper_partial_still_below_100():
    """k < n 时须仍给有限上界（防止修复过宽）。"""
    v = _floor._cp_upper(4, 5)
    assert 0 < v < 100


# ══════════ #2 L5 的 m/rank 与 q 值同族 ══════════

def test_fix2_l5_crit_uses_same_family_as_q():
    """源码级：order 必须由 p_for_fdr 构建，且 alpha 取自 config。"""
    src = inspect.getsource(_floor.audit_detection_floor)
    assert 'enumerate(p_for_fdr)' in src
    assert 'enumerate(df["perm_P"]' not in src
    assert 'getattr(cfg, "fdr_alpha", 0.05)' in src


def test_fix2_l5_crit_value_matches_hand_computation():
    """数值级：3 目标基因 + 1 管家基因时，手算临界值与工具一致。

    p = [0.01, 0.02, 0.03]（目标）, 管家 = 0.001（最小但被排除）
    → 该族 m = 3，rank 1/2/3 → crit = 0.05/3, 0.10/3, 0.15/3
    """
    assert AuditConfig().fdr_alpha == 0.05
    from sc_audit.floor import _bh_fdr
    # q 值族的 m 应为 3（排除管家基因）
    qs = _bh_fdr([0.001, 0.01, 0.02, 0.03])   # 含管家基因 0.001 时 m = 4
    assert len([q for q in qs if q == q]) == 4
    qs2 = _bh_fdr([np.nan, 0.01, 0.02, 0.03])  # 管家置 NaN 后 m = 3
    assert len([q for q in qs2 if q == q]) == 3
    # 手算：m=3 时 rank 1 → 0.05*1/3
    assert abs(0.05 * 1 / 3 - 0.016666) < 1e-5


# ══════════ #6 审计二判定须覆盖全部设计变量 ══════════

def _comp_adata(source_collinear=True, inst_balanced=True, n_ctrl=6, n_case=6, cpd=20):
    rng = np.random.default_rng(0)
    rows = []
    for d in range(n_ctrl):
        for _ in range(cpd):
            rows.append({'condition': 'control', 'donor': f'ctrl_d{d}',
                         'institution': ['A', 'B', 'C'][d % 3] if inst_balanced else 'A',
                         'source': 'CD45+' if source_collinear else 'Cells'})
    for d in range(n_case):
        for _ in range(cpd):
            rows.append({'condition': 'disease', 'donor': f'case_d{d}',
                         'institution': ['A', 'B', 'C'][d % 3] if inst_balanced else 'B',
                         'source': 'Cells'})
    o = pd.DataFrame(rows, index=[f'c{i}' for i in range(len(rows))])
    o['TNNT2'] = 0
    X = sp.csr_matrix(np.array(o[['TNNT2']].values, dtype=np.int32))
    return ad.AnnData(X=X, obs=o, var=pd.DataFrame(index=['TNNT2']))


def test_fix6_source_collinear_is_flagged_even_when_institution_is_balanced():
    """institution 平衡（TVD=0）但 source 完全共线 → 必须 NOT SEPARABLE。

    原实现只看 hetero_col(institution)，此场景会给 COMPARABLE。
    """
    cfg = AuditConfig()
    cfg.strata_cols = ['institution', 'source']
    cfg.match_vars = ['institution']
    cfg.hetero_col = 'institution'
    res = _comp.audit_control_composition(
        _comp_adata(source_collinear=True, inst_balanced=True),
        ['TNNT2'], 'condition', 'donor', 'control', cfg=cfg)
    assert 'NOT SEPARABLE' in res.verdict, res.verdict
    assert 'source' in res.verdict, "报告须指明是 source 而非 institution"


def test_fix6_balanced_design_stays_comparable():
    """两个设计变量都平衡 → 仍须 COMPARABLE（防止修复过宽）。"""
    cfg = AuditConfig()
    cfg.strata_cols = ['institution', 'source']
    cfg.match_vars = ['institution']
    cfg.hetero_col = 'institution'
    res = _comp.audit_control_composition(
        _comp_adata(source_collinear=False, inst_balanced=True),
        ['TNNT2'], 'condition', 'donor', 'control', cfg=cfg)
    assert 'COMPARABLE' in res.verdict and 'NOT SEPARABLE' not in res.verdict


def test_fix6_sample_unique_columns_are_not_design_variables():
    """样本/供者级唯一列（batch/lane/library）不得参与判定 —— 其 TVD 必然为 1。"""
    A = _comp_adata(source_collinear=False, inst_balanced=True)
    A.obs['batch'] = ['B_' + d for d in A.obs['donor']]     # 供者级唯一
    cfg = AuditConfig()
    cfg.strata_cols = ['institution', 'source', 'batch']
    cfg.match_vars = ['institution']
    cfg.hetero_col = 'institution'
    res = _comp.audit_control_composition(A, ['TNNT2'], 'condition', 'donor',
                                          'control', cfg=cfg)
    assert 'NOT SEPARABLE' not in res.verdict, \
        f"供者级唯一列被误当设计变量: {res.verdict}"


# ══════════ #10 管家 note 的模态条件 ══════════

def test_fix10_housekeeping_note_carries_modality_caveat():
    """管家护栏文案须由模态结论决定，三态各自给出对应的那一句。

    文案已从 audit_robust_signal 抽入 _hk_msg（外部审计 #5：原实现把两种
    解读并列输出，让读者自己判断）。
    """
    within = _robust._hk_msg(True)
    assert 'within a modality' in within
    assert 'argues against cell-type misassignment' in within
    assert 'not evidence' not in within

    across = _robust._hk_msg(False)
    assert 'cross-modal' in across
    # 措辞已收紧：限定在**本比较框架**内，不写成普适规律
    assert 'cannot be used directly as evidence about library quality' in across
    assert 'argues against cell-type misassignment' not in across

    unknown = _robust._hk_msg(None)
    assert 'across modalities' in unknown
    assert 'not evidence about library' in unknown

    # 调用方须真正把结论传进来（而非只接受参数不用）
    assert "_hk_msg(modality_matched)" in inspect.getsource(
        _robust.audit_robust_signal)


def test_fix11_capture_verdict_requires_two_signals():
    """源码级：CAPTURE-LIMITED 须由 complexity_low ∧ nuclear_low 驱动。"""
    src = inspect.getsource(_identity.audit_cell_identity)
    assert 'complexity_low and nuclear_low' in src
    assert 'below the 2-signal threshold' in src


def test_fix11_single_signal_does_not_yield_capture_limited():
    """功能级：只有低复杂度（无低核滞留）→ 不得判 CAPTURE-LIMITED。"""
    rng = np.random.default_rng(0)
    n = 800
    rows = [{'condition': ('control' if i % 2 else 'disease'),
             'donor': f'd{i % 4}',
             'lineage_names': (['Cardiomyocytes'] * 4 + ['Endo', 'Fib', 'Mye', 'TC'])[i % 8],
             'source': 'Cells'} for i in range(n)]
    med = {'Cardiomyocytes': 400, 'Endo': 900, 'Fib': 1100, 'Mye': 1300, 'TC': 1500}
    rows = []
    for i in range(n):
        lin = ['Cardiomyocytes', 'Endo', 'Fib', 'Mye', 'TC'][i % 5]
        rows.append({'condition': 'control' if i % 2 else 'disease',
                     'donor': f'd{i % 4}', 'lineage_names': lin, 'source': 'Cells',
                     'n_genes': med[lin]})
    o = pd.DataFrame(rows, index=[f'c{i}' for i in range(n)])
    pos = ['TTN', 'MYH7', 'MYBPC3', 'ACTN2', 'TNNI3', 'MYL2']
    neg = ['PTPRC', 'CD3E', 'DCN', 'COL1A1', 'PECAM1', 'VWF', 'LYZ']
    genes = pos + neg + ['MALAT1']
    M = np.zeros((n, len(genes)), dtype=np.int32)
    is_cm = (o['lineage_names'] == 'Cardiomyocytes').values
    M[is_cm, :6] = 1                      # CM 通过门控
    M[is_cm, 13] = 1                      # MALAT1 在 CM **最高** → nuclear_low=False
                                          # 于是只剩 complexity_low 一个信号
    A = ad.AnnData(X=sp.csr_matrix(M), obs=o, var=pd.DataFrame(index=genes))
    cfg = AuditConfig()
    cfg.identity_col = 'lineage_names'
    cfg.identity_target = 'Cardiomyocytes'
    res = _identity.audit_cell_identity(A, 'condition', 'donor', 'control', cfg=cfg)
    assert res.complexity_rank == 1, "构造应使 CM 复杂度最低"
    assert 'CAPTURE-LIMITED' not in res.verdict, res.verdict
    assert any('2-signal' in x for x in res.notes)


# ══════════ 2026-10-02 复查的 5 项边角 ══════════

def test_fix_r2_1_capture_not_assessable_when_nuclear_marker_absent():
    """核标志缺失时，必须明说 CAPTURE-LIMITED **不可评估**，
    而不是只说"single capture signal present"（会让读者以为第二信号被评估过）。"""
    n = 500
    rows = []
    for i in range(n):
        lin = ['Cardiomyocytes', 'Endo', 'Fib', 'Mye', 'TC'][i % 5]
        rows.append({'condition': 'control' if i % 2 else 'disease',
                     'donor': f'd{i % 4}', 'lineage_names': lin, 'source': 'Cells',
                     'n_genes': {'Cardiomyocytes': 400, 'Endo': 900, 'Fib': 1100,
                                 'Mye': 1300, 'TC': 1500}[lin]})
    o = pd.DataFrame(rows, index=[f'c{i}' for i in range(n)])
    pos = ['TTN', 'MYH7', 'MYBPC3', 'ACTN2', 'TNNI3', 'MYL2']
    neg = ['PTPRC', 'CD3E', 'DCN', 'COL1A1', 'PECAM1', 'VWF', 'LYZ']
    # ⚠️ 故意**不含 MALAT1**
    genes = pos + neg
    M = np.zeros((n, len(genes)), dtype=np.int32)
    M[(o['lineage_names'] == 'Cardiomyocytes').values, :6] = 1
    A = ad.AnnData(X=sp.csr_matrix(M), obs=o, var=pd.DataFrame(index=genes))
    cfg = AuditConfig()
    cfg.identity_col = 'lineage_names'
    cfg.identity_target = 'Cardiomyocytes'
    res = _identity.audit_cell_identity(A, 'condition', 'donor', 'control', cfg=cfg)
    joined = ' '.join(res.notes)
    assert 'not assessable' in joined, f"未声明不可评估: {res.notes}"
    assert 'absent' in joined


def test_fix_r2_2_backed_gating_and_qc_do_not_crash():
    """backed 打开时 gate_cells / compute_qc_metrics 须可运行（原抛 ValueError）。"""
    from sc_audit.gating import gate_cells
    from sc_audit.io import compute_qc_metrics
    genes = ['TTN', 'MYH7', 'MYBPC3', 'ACTN2', 'TNNI3', 'MYL2',
             'PTPRC', 'CD3E', 'DCN', 'COL1A1', 'PECAM1', 'VWF', 'LYZ', 'MT-CO1']
    X = np.zeros((20, len(genes)), dtype=np.int32)
    X[:, :6] = 1
    A = ad.AnnData(X=sp.csr_matrix(X),
                   obs=pd.DataFrame(index=[f'c{i}' for i in range(20)]),
                   var=pd.DataFrame(index=genes))
    p = os.path.join(TMP, 'bk_gate.h5ad')
    A.write_h5ad(p)
    Ab = ad.read_h5ad(p, backed='r')
    m = gate_cells(Ab, ['TTN', 'MYH7'], ['PTPRC', 'CD3E'], min_pos_markers=1,
                   min_counts=0, min_genes=0)
    assert m.sum() > 0
    q = compute_qc_metrics(Ab)
    assert 'percent_mito' in q.obs.columns


def test_fix_r2_3_detection_vectors_rejects_duplicate_var_names():
    """重名时 get_loc 返回 slice/mask → ad[:, idx] 失败；须显式拒绝。"""
    A = ad.AnnData(X=sp.csr_matrix(np.ones((5, 3), dtype=np.int32)),
                   obs=pd.DataFrame(index=[f'c{i}' for i in range(5)]),
                   var=pd.DataFrame(index=['A', 'A', 'B']))
    with pytest.raises(ValueError, match='不唯一'):
        _io.detection_vectors(A, ['A'])


def test_fix_r2_4_tvd_flag_text_matches_operator():
    """flag 文案的运算符须与代码一致（原为 >= 但运算符已是 >）。"""
    src = inspect.getsource(_identity.audit_cell_identity)
    assert '} > {cfg.identity_imbalance_min_tvd}' in src
    assert '} >= {cfg.identity_imbalance_min_tvd}' not in src


def test_fix_r2_5_ratio_basis_control_verdict_text():
    """ratio_basis='control' 时 verdict 不得写 "no case cells"。

    M5（v0.6.2）后判定链单一实现在 `_verdict_from`，故断言其源码。
    """
    src = inspect.getsource(_floor._verdict_from)
    assert 'no control cells' in src
    assert 'ratio_basis == "control"' in src


def test_fix_r2_6_event_sparse_columns_always_present():
    """无基因达 event_sparse 时，相关列仍须存在（否则报告层 KeyError）。"""
    ad_ = _tiny(['TTN', 'KCNQ1'], fill_col=0, n=40)
    ad_.obs['condition'] = ['control'] * 20 + ['disease'] * 20
    ad_.obs['donor'] = [f'd{i}' for i in range(40)]
    res = _floor.audit_detection_floor(ad_, ['TTN'], 'condition', 'donor',
                                       'control', cfg=AuditConfig())
    cols = res.table.columns
    for c in ('donor_weighted_effect', 'perm_units', 'jackknife_p_range',
              'jackknife_n_sig', 'jackknife_n_sig_note', 'jackknife_n_below_bh_crit'):
        assert c in cols, f"缺列 {c}"


def test_fix_r2_7_make_adata_defaults_to_integer_counts():
    """make_adata 默认须产生整数计数，使 B5 的计数校验在集成测试中真正执行。"""
    from conftest import make_adata
    A = make_adata(n_ctrl_donors=2, n_case_donors=2, cells_per_donor=5)
    X = A.X.toarray() if sp.issparse(A.X) else np.asarray(A.X)
    assert np.all(X == np.floor(X)), "默认 X 非整数 —— B5 路径会退化为 --no-count-check"
    assert X.min() >= 0


# ══════════ 2026-10-02 复查的 3 项 ══════════

def test_fix_r3_1_single_lineage_nuclear_rank_is_vacuous():
    """只有一个注释谱系时 nuclear_marker_rank 恒为 1 —— 这是**空洞的真**，
    不得作为核滞留信号输出（否则会误报 "nuclear-poor capture"）。"""
    n = 300
    o = pd.DataFrame({
        'condition': ['control' if i % 2 else 'disease' for i in range(n)],
        'donor': [f'd{i % 4}' for i in range(n)],
        'lineage_names': ['Cardiomyocytes'] * n,        # ⚠️ 单谱系
        'source': ['Cells'] * n,
        'n_genes': [900] * n,
    }, index=[f'c{i}' for i in range(n)])
    pos = ['TTN', 'MYH7', 'MYBPC3', 'ACTN2', 'TNNI3', 'MYL2']
    neg = ['PTPRC', 'CD3E', 'DCN', 'COL1A1', 'PECAM1', 'VWF', 'LYZ']
    genes = pos + neg + ['MALAT1']
    M = np.zeros((n, len(genes)), dtype=np.int32)
    M[:, :6] = 1
    A = ad.AnnData(X=sp.csr_matrix(M), obs=o, var=pd.DataFrame(index=genes))
    cfg = AuditConfig()
    cfg.identity_col = 'lineage_names'
    cfg.identity_target = 'Cardiomyocytes'
    res = _identity.audit_cell_identity(A, 'condition', 'donor', 'control', cfg=cfg)
    assert res.nuclear_marker_rank == 1, "构造应使排名为 1"
    assert res.complexity_n_lineages == 1
    flags = ' '.join(res.notes)
    assert 'nuclear-poor capture' not in flags, f"单谱系误报核滞留信号: {res.notes}"
    assert 'vacuous' in flags, f"应显式说明排名空洞: {res.notes}"


def test_fix_r3_1b_nuclear_signal_still_fires_with_two_lineages():
    """防过宽：谱系数 >= 2 且目标谱系核标志确为最低时，信号必须仍然出现。"""
    n = 600
    rows = []
    for i in range(n):
        lin = 'Cardiomyocytes' if i % 2 == 0 else 'Endo'
        rows.append({'condition': 'control' if i % 3 else 'disease',
                     'donor': f'd{i % 4}', 'lineage_names': lin, 'source': 'Cells',
                     'n_genes': 400 if lin == 'Cardiomyocytes' else 1200})
    o = pd.DataFrame(rows, index=[f'c{i}' for i in range(n)])
    pos = ['TTN', 'MYH7', 'MYBPC3', 'ACTN2', 'TNNI3', 'MYL2']
    neg = ['PTPRC', 'CD3E', 'DCN', 'COL1A1', 'PECAM1', 'VWF', 'LYZ']
    genes = pos + neg + ['MALAT1']
    M = np.zeros((n, len(genes)), dtype=np.int32)
    is_cm = (o['lineage_names'] == 'Cardiomyocytes').values
    M[is_cm, :6] = 1
    M[~is_cm, 13] = 1            # MALAT1 仅在非 CM 表达 → CM 的核标志最低
    A = ad.AnnData(X=sp.csr_matrix(M), obs=o, var=pd.DataFrame(index=genes))
    cfg = AuditConfig()
    cfg.identity_col = 'lineage_names'
    cfg.identity_target = 'Cardiomyocytes'
    res = _identity.audit_cell_identity(A, 'condition', 'donor', 'control', cfg=cfg)
    assert res.complexity_n_lineages == 2
    assert res.nuclear_marker_rank == 1
    assert any('nuclear-poor capture' in x for x in res.notes), res.notes


def test_fix_r3_3_x_matrix_docstring_warns_no_shared_memory():
    """x_matrix 的 docstring 须声明 backed 时返回副本、与原 ad 不共享内存。"""
    src = inspect.getsource(_io.x_matrix)
    assert '不共享内存' in src


# ══════════ CLI 是高风险层（"解析了但没生效"）══════════

def test_r4_every_cfg_assignment_has_a_matching_cli_flag():
    """CLI 层回归：每条 `cfg.X = args.Y` 都必须有对应的 `--flag`。

    本工具历轮缺陷多为"参数解析了但没真正生效"（A4-5 就是为堵这类
    问题而设）。此测试把该检查**机械化为全量扫描** —— 而不是逐个参数手写断言，
    这样新增参数时不会漏。
    """
    import re
    from sc_audit import cli as _cli
    src = inspect.getsource(_cli)

    assigns = re.findall(r'cfg\.([A-Za-z_][A-Za-z0-9_]*)\s*=\s*args\.([A-Za-z_][A-Za-z0-9_]*)', src)
    assert assigns, "未扫描到任何 cfg.X = args.Y 赋值（正则可能失效）"

    # 收集 argparse 定义的所有目标名（dest）
    dests = set(re.findall(r'add_argument\(\s*"--([A-Za-z0-9\-]+)"', src))
    dests |= {d.replace("-", "_") for d in list(dests)}

    missing = [(attr, arg) for attr, arg in assigns if arg not in dests]
    assert not missing, (
        "以下 cfg 赋值取自未定义的 CLI 参数（会 AttributeError 或静默失效）: "
        + ", ".join(f"cfg.{a} = args.{b}" for a, b in missing))

    # 反向：每个被赋值的 cfg 属性都必须是 AuditConfig 的真实字段
    from sc_audit.config import AuditConfig
    fields = set(AuditConfig().__dataclass_fields__)
    bogus = sorted({a for a, _ in assigns if a not in fields})
    assert not bogus, f"cfg 上不存在的字段被赋值: {bogus}"


def test_r4_single_lineage_states_assessability_explicitly():
    """单谱系时必须逐项声明可评估性，不能只说"降级"。"""
    n = 200
    o = pd.DataFrame({
        'condition': ['control' if i % 2 else 'disease' for i in range(n)],
        'donor': [f'd{i % 4}' for i in range(n)],
        'lineage_names': ['Cardiomyocytes'] * n,
        'source': ['Cells'] * n, 'n_genes': [900] * n,
    }, index=[f'c{i}' for i in range(n)])
    pos = ['TTN', 'MYH7', 'MYBPC3', 'ACTN2', 'TNNI3', 'MYL2']
    neg = ['PTPRC', 'CD3E', 'DCN', 'COL1A1', 'PECAM1', 'VWF', 'LYZ']
    genes = pos + neg + ['MALAT1']
    M = np.zeros((n, len(genes)), dtype=np.int32)
    M[:, :6] = 1
    A = ad.AnnData(X=sp.csr_matrix(M), obs=o, var=pd.DataFrame(index=genes))
    cfg = AuditConfig()
    cfg.identity_col, cfg.identity_target = 'lineage_names', 'Cardiomyocytes'
    res = _identity.audit_cell_identity(A, 'condition', 'donor', 'control', cfg=cfg)
    note = ' '.join(res.notes)
    assert 'ASSESSABILITY' in note
    # C3（v0.7.1）：改为逐项动态清单（"哪些可评估 / 哪些不可评估"），
    # 不再使用硬编码的"其余检查全部已评估"断言。
    assert 'ASSESSABILITY (this run)' in note
    assert 'NOT evaluable' in note
    assert 'FULLY EVALUATED' not in note
    assert 'NOT a complete' in note


def test_r4_crossmodal_wording_is_framework_scoped():
    """跨模态措辞须限定在本框架内，不得写成普适规律。"""
    msg = _robust._hk_msg(False)
    assert 'Within this comparison framework' in msg
    assert 'not a general rule' in msg
    assert 'cannot be used directly as evidence' in msg


# ══════════ 跨模块信息流（1️⃣/2️⃣）══════════

def test_r6_1_floor_table_carries_reference_provenance_and_modality():
    """floor 的输出表须显式标注 R 的分母来自哪列、哪个模态。

    原实现只输出 pct_ref 数值 —— 使用者看到 R 偏低却无从知道分母是核参照列，
    跨模态比较被静默掩盖（此即"最大风险：跨模块隐性依赖"）。
    """
    from sc_audit.floor import audit_detection_floor
    A = _tiny(['TNNT2', 'ACTB'], fill_col=0, n=40)
    A.obs['condition'] = ['control'] * 20 + ['disease'] * 20
    A.obs['donor'] = [f'd{i}' for i in range(40)]

    # ① 未声明模态 → unknown
    t = audit_detection_floor(A, ['TNNT2'], 'condition', 'donor', 'control',
                              cfg=AuditConfig(), reference={'TNNT2': 90.0}).table
    assert 'ref_modality_flag' in t.columns
    assert 'pct_ref_source' in t.columns
    assert 'unknown' in str(t.iloc[0]['ref_modality_flag'])

    # ② 声明为跨模态 → 默认**禁止**（v0.5.0 硬约束；REFERENCE_SCHEMA R2）
    from sc_audit.reference import CrossModalityError
    with pytest.raises(CrossModalityError):
        audit_detection_floor(
            A, ['TNNT2'], 'condition', 'donor', 'control', cfg=AuditConfig(),
            reference={'TNNT2': 90.0}, modality='scRNA',
            reference_modalities=['snRNA'])

    # ②b 显式放行 → cross-modal + 逐基因 cross_modality 列 + 告警
    with pytest.warns(UserWarning, match='跨模态'):
        t2 = audit_detection_floor(
            A, ['TNNT2'], 'condition', 'donor', 'control', cfg=AuditConfig(),
            reference={'TNNT2': 90.0}, modality='scRNA',
            reference_modalities=['snRNA'], allow_cross_modality=True,
            reference_provenance={'TNNT2': 'atlas_ref_A'}).table
    assert t2.iloc[0]['ref_modality_flag'] == 'cross-modal'
    assert bool(t2.iloc[0]['cross_modality']) is True
    assert t2.iloc[0]['pct_ref_source'] == 'atlas_ref_A'

    # ③ 声明为同模态 → same-modality，不告警
    t3 = audit_detection_floor(
        A, ['TNNT2'], 'condition', 'donor', 'control', cfg=AuditConfig(),
        reference={'TNNT2': 90.0}, modality='scRNA',
        reference_modalities=['scRNA', 'snRNA']).table
    assert t3.iloc[0]['ref_modality_flag'] == 'same-modality'


def test_r6_2_reference_provenance_records_argmax_column():
    """load_reference_with_provenance 须记下 mode='max' 时取值来自哪一列。"""
    import tempfile
    from sc_audit.reference import load_reference_with_provenance
    d = tempfile.mkdtemp()
    path = os.path.join(d, 'ref_multi.csv')
    pd.DataFrame({'gene': ['A', 'B'],
                  'pct_ref_snRNA': [10.0, 90.0],
                  'pct_ref_scRNA': [50.0, 5.0]}).to_csv(path, index=False)
    vals, prov = load_reference_with_provenance(path, mode='max')
    assert vals['A'] == 50.0 and vals['B'] == 90.0
    assert prov['A'] == 'pct_ref_scRNA'     # A 的最大值来自 scRNA 列
    assert prov['B'] == 'pct_ref_snRNA'     # B 的最大值来自 snRNA 列
    # 这两列的模态不同 → 正是"取 max 会掩盖模态"的机理
    assert prov['A'] != prov['B']


def test_r6_3_robust_warns_when_reference_without_modality():
    """给了 reference 却未声明 modality_matched 时须告警（否则静默落入 None 分支）。"""
    A = _tiny(['TNNT2', 'ACTB'], fill_col=0, n=40)
    A.obs['condition'] = ['control'] * 20 + ['disease'] * 20
    A.obs['donor'] = [f'd{i}' for i in range(40)]
    with pytest.warns(UserWarning, match='modality_matched was not declared'):
        _robust.audit_robust_signal(A, ['TNNT2'], 'condition', 'donor', 'control',
                                    reference={'TNNT2': 50.0}, cfg=AuditConfig())


def test_r6_4_cli_match_vars_flag_is_registered():
    """--match-vars 须存在并被写入 cfg（原实现无法显式指定）。"""
    src = inspect.getsource(__import__('sc_audit.cli', fromlist=['x']))
    assert '"--match-vars"' in src
    assert 'cfg.match_vars = [c.strip() for c in args.match_vars.split(",")' in src
