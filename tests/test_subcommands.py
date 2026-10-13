"""子命令冒烟测试 —— 每个子命令至少一条。

为什么单独成文件（直接教训）：
  `compare` 子命令因 `reference_mode=` 关键字被删除而抛 TypeError，
  但 115 项测试全绿——因为 **compare / composition / robust 三个子命令
  一条 CLI 测试都没有**。这类"改签名/删参数"式回归只能靠子命令级冒烟测试兜住。
"""
import json
import os

import pytest

from sc_audit.cli import main as cli_main

# 每个子命令的最小参数集（占位值在 fixture 中替换）
SUBCOMMANDS = [
    ("run", []),
    ("floor", []),
    ("composition", []),
    ("identifiability", []),
    ("robust", ["--reference", "@REF@"]),
    ("compare", ["--reference", "@REF@"]),
    ("check", []),
]


@pytest.fixture
def toy_inputs(tmp_path):
    """合成 h5ad + 参照表 csv。"""
    from conftest import make_adata
    ad = make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=25,
                    gene_pcts={"TNNT2": (0.9, 0.4), "PLN": (0.8, 0.35)})
    h5 = tmp_path / "toy.h5ad"
    ad.write_h5ad(h5)

    ref = tmp_path / "ref.csv"
    syms = list(ad.var_names)
    ref.write_text("gene,pct_ref\n" + "\n".join(f"{g},70.0" for g in syms),
                   encoding="utf-8")
    return h5, ref


@pytest.mark.parametrize("cmd,extra", SUBCOMMANDS, ids=[c for c, _ in SUBCOMMANDS])
def test_subcommand_smoke(tmp_path, toy_inputs, cmd, extra):
    """每个子命令都必须能跑通并产出预期文件。"""
    h5, ref = toy_inputs
    out = tmp_path / f"out_{cmd}"
    argv = [cmd, "--input", str(h5), "--target-genes", "TNNT2,PLN,ACTB",
            "--condition-col", "condition", "--donor-col", "donor",
            "--n-perm", "50", "--outdir", str(out), "--quiet", "--no-count-check"]
    argv += [str(ref) if e == "@REF@" else e for e in extra]

    rc = cli_main(argv)
    assert rc == 0, f"{cmd} 子命令返回 {rc}"

    expected = {
        "run": ["report.html", "floor_audit.csv", "composition.txt",
                "identifiability.txt", "identifiability_roster.csv",
                "identifiability_ladder_summary.csv", "robustness.txt",
                "run_metadata.json"],
        "floor": ["floor_audit.csv"],
        "composition": ["composition.txt"],
        "identifiability": ["identifiability.txt", "identifiability_roster.csv",
                            "identifiability_ladder.csv",
                            "identifiability_ladder_summary.csv"],
        "robust": ["robustness.txt"],
        "compare": ["detection_ratio.csv"],
        "check": ["check_report.json", "check_report.txt"],
    }[cmd]
    for f in expected:
        assert (out / f).exists(), f"{cmd} 未产出 {f}"


def test_compare_output_is_wellformed(tmp_path, toy_inputs):
    """compare 的产物结构正确（直接回归点）。"""
    import pandas as pd
    h5, ref = toy_inputs
    out = tmp_path / "cmp"
    rc = cli_main(["compare", "--input", str(h5), "--target-genes", "TNNT2,PLN",
                   "--reference", str(ref), "--outdir", str(out),
                   "--quiet", "--no-count-check"])
    assert rc == 0
    df = pd.read_csv(out / "detection_ratio.csv")
    assert {"gene", "pct_study", "pct_ref", "ratio_R"} <= set(df.columns)
    assert set(df["gene"]) == {"TNNT2", "PLN"}


def test_reference_reference_mode_is_not_passed(tmp_path, toy_inputs):
    """注册表 DEADPARAM 的强化：cli 不得再向 compute_detection_ratio 传 reference_mode。

    注意：必须检查**实际调用**（`reference_mode=`），不能只搜词——
    源码注释里会提到该参数名（说明它已被删除）。
    """
    import inspect
    import sc_audit.cli as cli
    src = inspect.getsource(cli)
    assert "reference_mode=" not in src, \
        "cli 仍在向 compute_detection_ratio 传已删除的 reference_mode（compare 会 TypeError）"
    # 同时确认签名里确实没有该参数（双向校验）
    from sc_audit.reference import compute_detection_ratio
    assert "reference_mode" not in inspect.signature(compute_detection_ratio).parameters


def test_reference_alt_without_reference_skips_flip_stats(tmp_path, toy_inputs, capsys):
    """P3：只给 --reference-alt 而无 --reference 时，不得把"缺参照"算成"口径敏感"。"""
    import pandas as pd
    h5, ref = toy_inputs
    out = tmp_path / "alt_only"
    rc = cli_main(["floor", "--input", str(h5), "--target-genes", "TNNT2",
                   "--reference-alt", str(ref), "--outdir", str(out),
                   "--quiet", "--no-count-check"])
    assert rc == 0
    df = pd.read_csv(out / "floor_audit.csv")
    assert "verdict_alt" not in df.columns, "无主参照时不应计算 verdict_alt"
    assert "缺少 --reference" in capsys.readouterr().err


# ══════════════════ 审计四 CLI（v0.2）══════════════════

@pytest.fixture
def identity_inputs(tmp_path):
    """带细胞类型注释 + 来源列的合成数据，且**刻意制造亚群不平衡**。"""
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import anndata as ad

    POS = ["TTN", "MYH7", "MYBPC3", "ACTN2", "TNNI3", "MYL2"]
    NEG = ["PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"]
    EXTRA = ["MALAT1", "RYR2", "MYL3", "ACTB", "GAPDH"]
    genes = POS + NEG + EXTRA
    rng = np.random.default_rng(7)
    rows = []
    for arm, nd in (("control", 4), ("disease", 4)):
        for d in range(nd):
            for i in range(40):
                is_cm = rng.random() < 0.5
                # ★ 某分选来源的细胞只出现在对照组
                src = "CD45+" if (is_cm and arm == "control" and rng.random() < 0.6) else "Cells"
                rows.append({"condition": arm, "donor": f"{arm[0]}{d}",
                             "lineage_names": "Cardiomyocytes" if is_cm else "Endothelial cells",
                             "source": src, "n_genes": int(rng.normal(900, 100)),
                             "_cm": is_cm, "_arm": arm})
    n = len(rows)
    X = np.zeros((n, len(genes)), dtype=np.int32)
    for i, r in enumerate(rows):
        for j, g in enumerate(genes):
            if g in POS:
                p = 0.85 if r["_cm"] else 0.05
            else:
                p = 0.02
            X[i, j] = (rng.poisson(3) + 1) if rng.random() < p else 0
    A = ad.AnnData(X=sp.csr_matrix(X),
                   obs=pd.DataFrame(rows, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=genes))
    h5 = tmp_path / "ident.h5ad"
    A.write_h5ad(h5)
    ref = tmp_path / "ref.csv"
    ref.write_text("gene,pct_ref\n" + "\n".join(f"{g},70.0" for g in genes), encoding="utf-8")
    return h5, ref


def test_identity_subcommand_runs_and_flags_imbalance(identity_inputs, tmp_path):
    """identity 子命令必须跑通、产出文件，并在 stdout 报出亚群不平衡。"""
    h5, _ = identity_inputs
    out = tmp_path / "ident_out"
    rc = cli_main(["identity", "--input", str(h5),
                   "--target-genes", "MYL3,ACTB",
                   "--condition-col", "condition", "--donor-col", "donor",
                   "--control-label", "control",
                   "--celltype-col", "lineage_names", "--source-col", "source",
                   "--target-lineage", "Cardiomyocytes", "--modality", "scRNA",
                   "--reference-modalities", "snRNA",
                   "--outdir", str(out), "--no-count-check"])
    assert rc == 0
    txt = (out / "identity.txt").read_text(encoding="utf-8")
    assert "Audit 4" in txt
    assert "subset imbalance" in txt or "CD45+" in txt
    assert "cross-modal" in txt


def test_identity_cli_args_reach_the_config(identity_inputs, tmp_path):
    """CLI 参数必须真的写进 cfg（不能只解析不使用）。"""
    h5, _ = identity_inputs
    out = tmp_path / "ident_out2"
    cli_main(["identity", "--input", str(h5), "--target-genes", "ACTB",
              "--condition-col", "condition", "--donor-col", "donor",
              "--control-label", "control",
              "--celltype-col", "unknown_column_xyz",
              "--outdir", str(out), "--quiet", "--no-count-check"])
    txt = (out / "identity.txt").read_text(encoding="utf-8")
    # 传了一个不存在的注释列 → 必须 INDETERMINATE，证明 --celltype-col 生效
    assert "INDETERMINATE" in txt


def test_identifiability_ladder_parsing(tmp_path, toy_inputs):
    """P2-9：`--ladder` 空值必须报错（不得静默忽略）；合法值被采纳。"""
    h5, _ = toy_inputs
    out = tmp_path / "out_ladder"
    base = ["identifiability", "--input", str(h5), "--target-genes", "TNNT2,PLN",
            "--condition-col", "condition", "--donor-col", "donor",
            "--outdir", str(out), "--quiet", "--no-count-check"]
    with pytest.raises(SystemExit) as e1:
        cli_main(base + ["--ladder", ""])
    assert e1.value.code == 2, "空 --ladder 应中止（exit 2）"
    with pytest.raises(SystemExit) as e2:
        cli_main(base + ["--ladder", "none"])
    assert e2.value.code == 2, "单规则应中止（无法做一致性检查）"
    assert cli_main(base + ["--ladder", "none;source"]) == 0
    txt = (out / "identifiability.txt").read_text(encoding="utf-8")
    assert "R1 source" in txt
