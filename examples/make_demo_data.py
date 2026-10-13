#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成一份**合成**演示数据（不依赖任何外部数据集），用于试跑 sc-audit。

为什么需要它：真实单细胞数据要么体积大、要么带使用限制；公开仓库的示例应当**自包含**。
本脚本生成一份带"已知真值"的病例-对照 h5ad：
  · 12 个供者（6 对照 / 6 病例），每供者 25 个细胞
  · 设计变量：对照里 1 个供者与病例同机构（制造"部分共线"的可识别性情形）
  · 基因面板含：完全可检出（ACTB/GAPDH）、部分检出（TNNT2/PLN）、检出地板（MYL3）

用法： python3 examples/make_demo_data.py [输出路径，默认 examples/demo/demo.h5ad]
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import anndata as ad
except ImportError:
    raise SystemExit("需要 anndata：pip install anndata")


def make_demo(n_ctrl=6, n_case=6, per_donor=25, seed=0):
    rng = np.random.default_rng(seed)
    donors, conds, insts, don_types, srcs = [], [], [], [], []
    for i in range(n_ctrl):
        donors += [f"ctrl_d{i}"] * per_donor
        conds += ["control"] * per_donor
        insts += ["INST_A" if i == 0 else "INST_B"] * per_donor      # 仅 1 个对照与病例同机构
        don_types += ["biopsy" if i == 0 else "DCD"] * per_donor
        srcs += ["Cells" if i == 0 else "CD45+"] * per_donor
    for j in range(n_case):
        donors += [f"case_d{j}"] * per_donor
        conds += ["disease"] * per_donor
        insts += ["INST_A"] * per_donor
        don_types += ["biopsy"] * per_donor
        srcs += ["Cells"] * per_donor

    n = len(donors)
    is_ctrl = np.array([c == "control" for c in conds])
    genes = ["TTN", "MYH7", "TNNT2", "PLN", "ACTB", "GAPDH",
             "MYL3", "KCNJ2", "KCNQ1", "KCNH2", "SCN5A", "PECAM1", "PTPRC", "DCN"]
    pct = {"TTN": (0.95, 0.95), "MYH7": (0.9, 0.9), "TNNT2": (0.95, 0.30),
           "PLN": (0.85, 0.35), "ACTB": (0.9, 0.9), "GAPDH": (0.9, 0.9),
           "MYL3": (0.02, 0.0),        # 检出地板（病例组恰好全 0）
           "KCNJ2": (0.05, 0.02), "KCNQ1": (0.08, 0.03), "KCNH2": (0.10, 0.04),
           "SCN5A": (0.10, 0.05), "PECAM1": (0.01, 0.01), "PTPRC": (0.01, 0.01),
           "DCN": (0.01, 0.01)}
    M = np.zeros((n, len(genes)), dtype=np.float32)
    for gi, g in enumerate(genes):
        pc, pk = pct[g]
        p = np.where(is_ctrl, pc, pk)
        M[:, gi] = (rng.random(n) < p).astype(np.float32)
    M = np.hstack([M, rng.integers(0, 3, size=(n, 8)).astype(np.float32)])

    obs = pd.DataFrame({
        "condition": conds, "donor": donors, "institution": insts,
        "donation_type": don_types, "source": srcs, "region": ["heart"] * n,
        "n_genes": rng.integers(400, 900, n), "total_counts": rng.integers(2000, 4000, n),
        "percent_mito": rng.uniform(5, 20, n)},
        index=[f"cell{i}" for i in range(n)])
    var = pd.DataFrame(index=genes + [f"bg{i}" for i in range(8)])
    return ad.AnnData(X=sp.csr_matrix(M), obs=obs, var=var)


def main():
    import argparse

    ap = argparse.ArgumentParser(
        description="生成合成演示数据（不依赖任何外部数据集）")
    ap.add_argument("out", nargs="?", default=None,
                    help="输出路径（默认 examples/demo/demo.h5ad）")
    ap.add_argument("--out", dest="out_opt", default=None,
                    help="输出路径（与位置参数等价；显式命名避免误用）")
    ns = ap.parse_args()
    if ns.out and ns.out_opt and ns.out != ns.out_opt:
        ap.error(f"位置参数与 --out 都给且不一致：{ns.out!r} vs {ns.out_opt!r}")
    out = ns.out_opt or ns.out or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "demo", "demo.h5ad")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    A = make_demo()
    A.write_h5ad(out)
    print(f"已生成合成演示数据：{out}")
    print(f"  {A.n_obs} 细胞 / {A.obs['donor'].nunique()} 供者 / "
          f"{A.obs['condition'].nunique()} 组；基因 {A.n_vars}")
    print("\n下一步：\n  sc-audit check --input " + out +
          " --condition-col condition --donor-col donor\n"
          "  sc-audit run   --input " + out +
          " --target-genes MYL3,PLN,TNNT2,ACTB,GAPDH "
          "--condition-col condition --donor-col donor --outdir demo_out/")


if __name__ == "__main__":
    main()
