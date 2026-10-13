"""审计一 · 事件数阈值推导的**可执行**校核（吸收第三份外部审查 B1/B2）。

预注册口径的操作定义：病例组检出率 p̂ = x/n_case 的 Clopper–Pearson **双侧** 95%
区间的**相对半宽** (upper−lower)/2/p̂ 必须低于事先声明的容差，取满足该条件的最小 x。

本文件把"锁定 20"变成可复算断言，并**显式钉住**一个易混的口径错配：
±43.4% 与 ±50% 两列的 x_min 不同，n 依赖性也随之不同。

"""
import numpy as np
import pytest

from sc_audit.floor import (audit_detection_floor, derive_event_threshold,
                            event_threshold_table, relative_halfwidth_cp)
from conftest import make_adata

# 预注册容差表（复算示例 n_case = 200；x_min 随 n 变化）
TOLERANCE_TABLE_N200 = {0.50: 17, 0.434: 21, 0.40: 24, 0.33: 33, 0.25: 50}


def test_tolerance_table_matches_pre_registered_tolerances():
    got = event_threshold_table(200)
    for tol, x in TOLERANCE_TABLE_N200.items():
        assert got[tol] == x, f"±{tol*100:.1f}% 应得 x_min={x}，实得 {got[tol]}"


def test_xmin_dependence_on_n_case_uses_one_tolerance():
    """同一 ±43.4% 口径下的 n 依赖性：n=200 → 21；n=50 → 16。"""
    assert derive_event_threshold(200, 0.434) == 21
    assert derive_event_threshold(400, 0.434) == 22
    assert derive_event_threshold(50, 0.434) == 16
    # 对照：±50% 列更宽松（x_min 更小）
    assert derive_event_threshold(400, 0.50) == 17
    assert derive_event_threshold(50, 0.50) == 13


def test_relative_halfwidth_values():
    """n_case = 200 时 k = 19/20/21 的双侧相对半宽可复算。"""
    assert abs(relative_halfwidth_cp(20, 200) - 0.440) < 5e-4
    assert abs(relative_halfwidth_cp(19, 200) - 0.454) < 5e-4
    assert abs(relative_halfwidth_cp(21, 200) - 0.428) < 5e-4
    assert np.isnan(relative_halfwidth_cp(0, 200)), "k=0 时相对半宽不可定义"


def test_floor_table_exposes_two_sided_relative_halfwidth():
    """审计一必须逐基因输出预注册容差口径的量，使阈值推导可被读者当场核验（见 docs/methods.md §2.2）。"""
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=25,
                    gene_pcts={"TNNT2": (0.9, 0.2)})
    res = audit_detection_floor(ad, ["TNNT2", "ACTB"], "condition", "donor", "control")
    assert "cp_rel_halfwidth_case" in res.table.columns
    row = res.table.loc[res.table.gene == "TNNT2"].iloc[0]
    n_case = 5 * 25
    k = int(round(row["case_pct"] / 100 * n_case))
    assert abs(row["cp_rel_halfwidth_case"] - relative_halfwidth_cp(k, n_case)) < 1e-3


def test_derivation_is_not_wired_into_the_locked_threshold():
    """纪律：推导函数只用于校核，绝不自动改写事前锁定的 event_sparse。"""
    from sc_audit.config import AuditConfig
    assert AuditConfig().event_sparse == 20, "event_sparse 仍是事前锁定值 20"


def test_not_in_dataset_rows_keep_a_stable_schema():
    """2026-10-06 跨数据集实跑发现：基因全部未命中时 table 会丢核心列 → 下游 KeyError。

    场景：数据集 var_names 只有 Ensembl ID 且未做基因名回退。
    修法：退化路径也补齐核心列（值可为 NaN，但列不会时有时无）。
    """
    ad = make_adata()
    res = audit_detection_floor(ad, ["NOT_THERE", "ALSO_NOT_THERE"],
                                "condition", "donor", "control")
    assert set(res.table["verdict"]) == {"NOT_IN_DATASET"}
    for col in ("n_pos", "n_pos_case", "case_pct", "ctrl_pct", "ratio_R", "pct_ref",
                "cp_rel_halfwidth_case", "perm_P", "stratum", "gene", "verdict"):
        assert col in res.table.columns, f"退化路径缺少列 {col}（模式不稳定）"
    assert res.table["n_pos"].isna().all()
