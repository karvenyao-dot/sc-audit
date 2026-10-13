"""v0.7.3 回归：供者级连续协变量不再被误标为"供者唯一标识"。

发现的来源：design-parse 第三步把 GEO 文本里 obs 缺失的 kidney `age` 注入后，
`is_donor_unique` 因"每供者 1 值且取值互异"把它判为 `excluded:donor-unique`
（当成 id 列）。实为**不可平衡的供者级协变量**——两者对匹配的处置相同，
但归因不同。修复：confounder 角色的列 → 新桶 `excluded:donor-level-all-distinct`。
"""
import numpy as np
import pytest

from conftest import make_adata
from sc_audit.config import AuditConfig
from sc_audit.identifiability import (audit_identifiability,
                                      BUCKET_DONOR_SATURATED, BUCKET_DONOR_UNIQUE,
                                      build_variable_roster)


def _ad():
    A = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=20)
    A.obs["condition"] = A.obs["condition"].astype(str)
    # 供者级 age：每个供者一个互异取值（10 供者 → 10 个不同 age）
    donors = list(dict.fromkeys(A.obs["donor"].astype(str).tolist()))
    age = {d: str(20 + i) for i, d in enumerate(sorted(donors))}
    A.obs["age"] = A.obs["donor"].astype(str).map(age).astype(str)
    return A


def test_donor_level_confounder_not_labelled_donor_unique():
    A = _ad()
    R, _ = build_variable_roster(A, "condition", "donor", "control", cfg=AuditConfig())
    b = dict(zip(R["column"], R["bucket"]))
    assert b["age"] == BUCKET_DONOR_SATURATED, b
    assert b.get("batch") == BUCKET_DONOR_UNIQUE, b          # 技术/批次列仍为 donor-unique


def test_donor_level_confounder_emits_note():
    A = _ad()
    r = audit_identifiability(A, ["ACTB", "GAPDH"], "condition", "donor", "control",
                              cfg=AuditConfig())
    assert any("donor-level covariate" in n and "age" in n for n in r.notes), r.notes
    assert "donor-level-all-distinct" in r.excluded
