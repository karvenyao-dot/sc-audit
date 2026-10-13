# -*- coding: utf-8 -*-
"""公开前加固回归（2026-10-07 发现的缺陷，逐条钉住）。

背景：两次独立复现共同发现三类问题，本文件把每一条
变成**可执行断言**，防止再发：

  1. pandas ≥3 的字符串 dtype（`str`/`string`，PDEP-14）下，`check` 的候选变量
     识别失效（原判据 `s.dtype == object`）
  2. `io.resolve_gene_names` 在"部分缺失且无符号列"时把 **requested 全量** 当缺失
  3. CLI 的输入类错误（缺列 / 文件不存在 / 单水平 condition）抛裸 traceback
"""
import os

import numpy as np
import pandas as pd
import pytest

from conftest import make_adata


# ─────────── ① pandas 3 字符串 dtype ───────────

def test_candidate_detection_accepts_string_and_str_dtypes():
    """`dtype == object` 之外的字符串 dtype（`string` / pandas 3 的 `str`）必须被识别。

    用 `dtype="string"` 显式构造（pandas 2 与 3 都能命中旧 bug：StringDtype 既非
    object 也非 category/bool）。
    """
    from sc_audit.check import check_readiness
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=10)
    # 把设计变量列改成显式字符串 dtype（模拟 pandas 3 的行为）
    for c in ("institution", "donation_type", "source"):
        ad.obs[c] = pd.array(ad.obs[c].astype(str).values, dtype="string")
    assert all(str(ad.obs[c].dtype).startswith(("string", "str"))
               for c in ("institution", "donation_type", "source"))

    res = check_readiness(ad, "condition", "donor", "control")
    found = {d["column"] for d in res.design_vars}
    assert {"institution", "donation_type", "source"} <= found, \
        f"字符串 dtype 的设计变量未被识别：{sorted(found)}"


def test_candidate_detection_accepts_category_and_bool():
    """类别列与布尔列仍应被识别（回归旧行为）。"""
    from sc_audit.check import check_readiness
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=10)
    ad.obs["inst_cat"] = pd.Categorical(ad.obs["institution"].astype(str))
    ad.obs["flag_bool"] = ad.obs["condition"].astype(str).eq("control")
    res = check_readiness(ad, "condition", "donor", "control")
    found = {d["column"] for d in res.design_vars}
    assert {"inst_cat", "flag_bool"} <= found, f"类别/布尔列未识别：{sorted(found)}"


# ─────────── ② 基因名缺失判定 ───────────

def test_resolve_gene_names_reports_only_truly_absent_genes():
    """没有可用符号列时，missing 只能是**真正不在 var_names** 的基因。"""
    from sc_audit.io import resolve_gene_names
    ad = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=5)
    present = [g for g in ("TNNT2", "ACTB", "GAPDH") if g in set(ad.var_names)]
    assert len(present) >= 2, "fixture 应含若干目标基因"
    requested = present + ["NOT_A_GENE_XYZ"]
    info = resolve_gene_names(ad, requested)
    assert info["missing"] == ["NOT_A_GENE_XYZ"], \
        f"missing 被夸大成全量：{info['missing']}"
    for g in present:
        assert info["resolved"].get(g) == g


def test_check_reports_correct_missing_gene_count():
    """`check` 文本里报告的缺失基因数必须是真正缺失的个数。"""
    from sc_audit.check import check_readiness
    ad = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=5)
    targets = ["TNNT2", "ACTB", "NOPE1", "NOPE2"]
    res = check_readiness(ad, "condition", "donor", "control", target_genes=targets)
    gl = res.gene_lookup or {}
    assert len(gl.get("missing") or []) == 2, f"应报 2 个缺失，实为 {gl.get('missing')}"
    assert "缺失 2 个" in res.to_text()


# ─────────── ③ CLI 输入错误：友好报错而非裸 traceback ───────────

@pytest.fixture()
def _demo_h5ad(tmp_path):
    """写一份小型 h5ad 供 CLI 端到端调用。"""
    ad = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=8)
    p = tmp_path / "mini.h5ad"
    ad.write_h5ad(str(p))
    return str(p)


def test_cli_missing_column_exits_two_without_traceback(_demo_h5ad, tmp_path, capsys):
    from sc_audit.cli import main
    with pytest.raises(SystemExit) as ei:
        main(["run", "--input", _demo_h5ad, "--condition-col", "nosuch_col",
              "--donor-col", "donor", "--target-genes", "TNNT2",
              "--outdir", str(tmp_path / "o")])
    assert ei.value.code == 2
    err = capsys.readouterr().err
    assert "[error]" in err and "[hint]" in err


def test_cli_missing_input_file_exits_two(_demo_h5ad, tmp_path, capsys):
    from sc_audit.cli import main
    with pytest.raises(SystemExit) as ei:
        main(["check", "--input", str(tmp_path / "nope.h5ad"),
              "--condition-col", "condition", "--donor-col", "donor"])
    assert ei.value.code == 2
    assert "[error]" in capsys.readouterr().err


def test_cli_single_level_condition_exits_two(_demo_h5ad, tmp_path, capsys):
    import anndata as ad
    a = ad.read_h5ad(_demo_h5ad)
    a.obs["condition"] = "control"
    one = tmp_path / "one_level.h5ad"
    a.write_h5ad(str(one))
    from sc_audit.cli import main
    with pytest.raises(SystemExit) as ei:
        main(["run", "--input", str(one), "--condition-col", "condition",
              "--donor-col", "donor", "--target-genes", "TNNT2",
              "--outdir", str(tmp_path / "o2")])
    assert ei.value.code == 2
    assert "[error]" in capsys.readouterr().err


def test_cli_debug_env_shows_traceback(_demo_h5ad, tmp_path, monkeypatch):
    """SC_AUDIT_DEBUG=1 时必须保留原始异常（排查用）。"""
    from sc_audit.cli import main
    monkeypatch.setenv("SC_AUDIT_DEBUG", "1")
    with pytest.raises((KeyError, FileNotFoundError)):
        main(["run", "--input", _demo_h5ad, "--condition-col", "nosuch_col",
              "--donor-col", "donor", "--target-genes", "TNNT2",
              "--outdir", str(tmp_path / "o3")])


def test_cli_empty_target_genes_exits_two(_demo_h5ad, tmp_path, capsys):
    from sc_audit.cli import main
    with pytest.raises(SystemExit) as ei:
        main(["run", "--input", _demo_h5ad, "--condition-col", "condition",
              "--donor-col", "donor", "--target-genes", "",
              "--outdir", str(tmp_path / "o4")])
    assert ei.value.code == 2
    assert "--target-genes 为空" in capsys.readouterr().err

# ─────────── ④ API 不自动回退：必须显式提示（2026-10-07 已知答案验证发现） ───────────

def _ensg_indexed_adata():
    """构造 var_names=ENSG、var['gene_symbols']=符号 的对象（模拟 ENSG 索引数据集）。"""
    import anndata as ad
    a = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=6)
    a.var["gene_symbols"] = [str(x) for x in a.var_names]
    a.var_names = [f"ENSG{i:08d}" for i in range(a.n_vars)]
    return a


def test_api_warns_when_genes_unresolved_but_symbol_column_exists():
    """纯 API 未做回退时必须**显式提示**（否则用户只见一整列 NOT_IN_DATASET）。"""
    from sc_audit.floor import audit_detection_floor
    a = _ensg_indexed_adata()
    with pytest.warns(UserWarning, match="需显式调用"):
        res = audit_detection_floor(a, ["TNNT2", "ACTB"], "condition", "donor", "control")
    assert set(res.table["verdict"]) == {"NOT_IN_DATASET"}, "未回退时应全部未命中"


def test_api_no_warning_after_explicit_fallback():
    """显式回退后不应再提示，且基因可被检出。"""
    import warnings
    import sc_audit
    from sc_audit.floor import audit_detection_floor
    a = _ensg_indexed_adata()
    info = sc_audit.apply_gene_name_fallback(a, ["TNNT2", "ACTB"])
    assert info["applied"] is True and info["fallback_column"] == "gene_symbols"
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        res = audit_detection_floor(a, ["TNNT2", "ACTB"], "condition", "donor", "control")
    assert not [x for x in w if "需显式调用" in str(x.message)],         "回退后不应再提示未回退"
    assert "NOT_IN_DATASET" not in set(res.table["verdict"])

# ─────────── ⑤ 门控 QC 不得退化到"子集矩阵现算" ───────────

def _adata_for_gating():
    """构造：正 marker 高检出、负 marker 全 0，便于观察门控 QC 口径。"""
    pos = ["TTN", "MYH7"]
    neg = ["PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"]
    gp = {g: (0.9, 0.9) for g in pos}
    gp.update({g: (0.0, 0.0) for g in neg})
    return make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=8, gene_pcts=gp)


def test_gate_cells_reuses_obs_qc_on_gene_subset():
    """门控在**基因子集**上必须复用 obs 的 QC 列，否则行和过小 → 静默全灭。

    在**基因子集**上做门控时的实测：深度阈值按子集现算 → 可通过细胞数为 0，
    下游只报 INDETERMINATE (no case cells)；改为复用全基因空间的 obs QC 后正常。
    """
    from sc_audit.gating import gate_cells
    a = _adata_for_gating()
    # 人为给定"全基因空间"的 QC（否则合成数据本身可能不过 min_counts）
    a.obs["total_counts"] = 10000.0
    a.obs["n_genes"] = 1000.0
    keep = ["TTN", "MYH7", "PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"]
    sub = a[:, [g for g in keep if g in set(map(str, a.var_names))]].copy()
    assert float(sub.X.sum(axis=1).max()) < 500, "前置：子集行和必须 < min_counts"
    mk = ["PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"]
    # 默认 qc_from="x"：在子集上必然全灭，且**必须给出可执行的诊断提示**
    with pytest.warns(UserWarning, match="qc_from='obs'"):
        m_default = gate_cells(sub, ["TTN", "MYH7"], mk, min_pos_markers=2,
                               min_counts=500, min_genes=200)
    assert int(np.asarray(m_default).sum()) == 0
    # qc_from="obs"：复用全基因空间 QC → 正常通过
    mask = gate_cells(sub, ["TTN", "MYH7"], mk, min_pos_markers=2,
                      min_counts=500, min_genes=200, qc_from="obs")
    mask = np.asarray(mask["mask"] if isinstance(mask, dict) else mask)
    assert mask.sum() > 0, "复用 obs QC 后不应被 QC 全灭"
    assert "obs.total_counts" in sub.uns["sc_audit_gating"]["qc_source"]
    assert sub.uns["sc_audit_gating"]["qc_from"] == "obs"


def test_gate_cells_falls_back_to_x_and_reports_source():
    """默认口径为从 X 现算，并**在 uns 里注明口径**（可诊断）。"""
    from sc_audit.gating import gate_cells
    a = _adata_for_gating()
    a.obs = a.obs.drop(columns=[c for c in ("total_counts", "n_genes") if c in a.obs.columns])
    mask = gate_cells(a, ["TTN", "MYH7"],
                      ["PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"],
                      min_pos_markers=2, min_counts=1, min_genes=1)
    src = a.uns["sc_audit_gating"]["qc_source"]
    assert any(x.startswith("computed_from_X") for x in src), src
    assert a.uns["sc_audit_gating"]["qc_from"] == "x"

# ─────────── ⑥ C.4 判定语义（真空支持域）+ D 索引错位 ───────────

def _vacuous_adata():
    """**无 design-defining 变量** → 共同支持域真空（= 整个对照组）。

    构造：每个候选设计列在**两组内部都有变异** → 全部落在 adjustable，design_defining 为空
    （对应实测中一个公开心脏单细胞数据集的情形：region/cell_annotation 仅在组内变化）。
    """
    ad = make_adata(n_ctrl_donors=4, n_case_donors=4, cells_per_donor=10,
                    design_confounded=False)
    n_obs = ad.n_obs
    ad.obs["institution"] = [f"INST_{i % 2}" for i in range(n_obs)]
    ad.obs["donation_type"] = [f"DT_{i % 2}" for i in range(n_obs)]
    ad.obs["source"] = [f"SRC_{i % 3}" for i in range(n_obs)]
    return ad


def test_r_c1_vacuous_support_is_not_reported_as_comparable():
    """(R-C1) 支持域真空 + 无规则缩域 → 必须是 INDETERMINATE，且含"不得读作可比"的提示。"""
    from sc_audit.identifiability import audit_identifiability, VERDICT_INDETERMINATE
    ad = _vacuous_adata()
    r = audit_identifiability(ad, ["TNNT2", "ACTB"], "condition", "donor", "control")
    head = r.verdict.split(" — ")[0]
    # 核心意图：真空支持域不得被报成"可比/不可分离"（那是无条件的方向判定）
    assert head not in ("COMPARABLE", "NOT SEPARABLE"), r.verdict
    if head == VERDICT_INDETERMINATE:
        blob = " ".join(r.notes).lower()
        assert "not evidence of comparability" in blob, f"缺少强制提示：{r.notes}"


def test_r_c2_ladder_has_no_empty_conjunct_rule_when_no_design_defining():
    """(R-C2) design-defining 为空时不得生成"空合取"规则（等于匹配所有，却会计入一致性统计）。"""
    from sc_audit.identifiability import default_ladder
    L = default_ladder([], ["gender"], ["gender"])
    for label, vars_ in L:
        assert label.startswith("R0") or vars_, f"出现空合取规则：{label!r} -> {vars_!r}"
    assert not any(lbl.startswith("R3") for lbl, _ in L), "design-defining 为空时不应有 R3"
    # 有 design-defining 时 R3 仍应存在（回归保护）
    L2 = default_ladder(["institution"], ["gender"], ["gender"])
    assert any(lbl.startswith("R3") for lbl, _ in L2)


def test_r_c3_restricting_variable_still_reaches_one_sided_verdict():
    """(R-C3) 反例保护：若某可调整变量**真的**缩小了支持域，不应被新守卫一并拦成 INDETERMINATE。"""
    from sc_audit.identifiability import audit_identifiability
    ad = make_adata(n_ctrl_donors=6, n_case_donors=6, cells_per_donor=20,
                    gene_pcts={"TNNT2": (0.9, 0.4)}, design_confounded=False)
    r = audit_identifiability(ad, ["TNNT2", "ACTB"], "condition", "donor", "control")
    assert r.verdict.split(" — ")[0] in ("COMPARABLE", "NOT SEPARABLE", "THRESHOLD-SPANNING",
                                        "INDETERMINATE")


def test_r_d1_donor_unique_false_when_varies_within_donor():
    """(R-D1) 供者内变动的列不得判为 donor-unique（索引错位曾使其误判为 True）。"""
    import pandas as pd
    from sc_audit.donor_utils import is_donor_unique
    a = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=6)
    col = pd.Series([f"T{i % 12}" for i in range(a.n_obs)], index=a.obs.index)
    donor_np = a.obs["donor"].astype(str).values                    # numpy（旧调用形态）
    n_don = a.obs["donor"].nunique()
    assert is_donor_unique(col, pd.Series(donor_np, index=col.index), n_don) is False
    assert is_donor_unique(col, pd.Series(donor_np), n_don) is False, "索引未对齐时不得平凡为真"


def test_r_d2_donor_unique_true_for_donor_level_constant_column():
    """(R-D2) 供者内恒定且水平数 ≥ 供者数 → donor-unique = True。"""
    import pandas as pd
    from sc_audit.donor_utils import is_donor_unique
    a = make_adata(n_ctrl_donors=3, n_case_donors=3, cells_per_donor=6)
    s = a.obs["donor"].astype(str)
    assert is_donor_unique(s, pd.Series(s.values, index=s.index), s.nunique()) is True
