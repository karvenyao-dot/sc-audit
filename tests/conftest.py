"""测试用合成数据与公共 fixture。

修订（2026-10-01）：新增 batch 列（样本级唯一，用于 B2 测试）、
scale_x 选项（用于 B5 计数校验测试）、以及多水平 condition 构造 helper。
"""
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

try:
    import anndata as ad
    _HAS_ANNDATA = True
except ImportError:  # pragma: no cover
    _HAS_ANNDATA = False


def make_adata(n_ctrl_donors=5, n_case_donors=5, cells_per_donor=20,
               seed=0, gene_pcts=None, design_confounded=True,
               add_batch=True, scale_x=False):
    """构造带已知检出率的合成 AnnData。

    gene_pcts: {gene: (ctrl_pct, case_pct)}，取值 0-1。
    add_batch: 是否添加样本级唯一的 batch 列（测试 B2）。
    scale_x:   是否把 X 变成 log-normalized（测试 B5 的计数校验）。
    """
    if not _HAS_ANNDATA:
        pytest.skip("anndata 未安装")

    rng = np.random.default_rng(seed)
    donors, conds, insts, don_types, srcs = [], [], [], [], []
    for i in range(n_ctrl_donors):
        donors += [f"ctrl_d{i}"] * cells_per_donor
        conds += ["control"] * cells_per_donor
        if design_confounded:
            insts += ["INST_A" if i == 0 else "INST_B"] * cells_per_donor
            don_types += ["biopsy" if i == 0 else "DCD"] * cells_per_donor
            srcs += ["Cells" if i == 0 else "CD45+"] * cells_per_donor
        else:
            insts += ["INST_A"] * cells_per_donor
            don_types += ["biopsy"] * cells_per_donor
            srcs += ["Cells"] * cells_per_donor
    for j in range(n_case_donors):
        donors += [f"case_d{j}"] * cells_per_donor
        conds += ["disease"] * cells_per_donor
        insts += ["INST_A"] * cells_per_donor
        don_types += ["biopsy"] * cells_per_donor
        srcs += ["Cells"] * cells_per_donor

    n = len(donors)
    is_ctrl = np.array([c == "control" for c in conds])

    genes = ["TTN", "MYH7", "MYBPC3", "ACTN2", "TNNI3", "MYL2",
             "PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ",
             "MYL3", "KCNJ2", "TNNT2", "ACTB", "GAPDH", "MT-CO1"]
    extra = [g for g in (gene_pcts or {}) if g not in genes]
    genes = genes + extra

    default = {g: (0.5, 0.5) for g in genes}
    for g in genes[:6]:
        default[g] = (0.9, 0.9)
    for g in genes[6:13]:
        default[g] = (0.0, 0.0)
    if gene_pcts:
        default.update(gene_pcts)

    M = np.zeros((n, len(genes)), dtype=np.float32)
    for gi, g in enumerate(genes):
        pc, pk = default[g]
        p = np.where(is_ctrl, pc, pk)
        M[:, gi] = (rng.random(n) < p).astype(np.float32)
    # 背景列取**整数**计数：默认路径的矩阵由此整体为整数，使 B5 的
    # "X 必须是原始计数"校验在集成测试中真正被执行（原为常量 0.4 → 非整数，
    # 所有集成测试都靠 --no-count-check 绕过）。
    M = np.hstack([M, rng.integers(0, 3, size=(n, 8)).astype(np.float32)])

    if scale_x:
        # 模拟 log-normalized：非整数、且可能出现负值（scale 后）
        M = np.log1p(M * 10) - 0.5

    var = pd.DataFrame(index=genes + [f"bg{i}" for i in range(8)])
    obs_data = {
        "condition": conds, "donor": donors,
        "institution": insts, "donation_type": don_types, "source": srcs,
        "region": ["heart"] * n,
        "n_genes": rng.integers(400, 900, n),
        "total_counts": rng.integers(2000, 4000, n),
        "percent_mito": rng.uniform(5, 20, n),
    }
    if add_batch:
        # 供者级唯一 —— 真实数据中 batch / lane / library 的典型形态（N2 的直接场景）
        obs_data["batch"] = [f"B_{d}" for d in donors]
    obs = pd.DataFrame(obs_data, index=[f"cell{i}" for i in range(n)])
    var.index.name = None
    return ad.AnnData(X=sp.csr_matrix(M), obs=obs, var=var)


@pytest.fixture
def adata_confounded():
    return make_adata(design_confounded=True)


@pytest.fixture
def adata_clean():
    return make_adata(design_confounded=False)


def relabel(ad, col, mapping_by_donor):
    for d, v in mapping_by_donor.items():
        ad.obs.loc[ad.obs["donor"] == d, col] = v
    return ad
