"""sc-audit 命令行入口。

修订记录（2026-10-01）：
  · B1  ensure_obs_fields 现在同时传 control_label 与 case_label
  · 新增 --case-label / --count-layer / --hetero-col
  · run 子命令补齐 --matched-donors / --housekeeping（原仅挂在子命令上）
  · B5  运行前调用 resolve_counts()，确保 X 为原始计数
  · 4.8 run_metadata.json 记录参照表路径 + SHA256 + 配置来源
  · compare 与管家基因 R 统一使用 cfg.ratio_basis（原用 overall，与默认不一致）
  · gating 阈值单一真值源：默认取自 AuditConfig（不再各自维护两套）
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import warnings
from typing import List, Optional

import pandas as pd

from . import __version__
from .config import load_config, AuditConfig
from .io import (load_data, ensure_obs_fields, compute_qc_metrics,
                 resolve_counts, resolve_gene_names, apply_gene_name_fallback)
from .gating import gate_cells, preset_from_config
from .floor import audit_detection_floor
from .composition import audit_control_composition
from .identifiability import audit_identifiability
from .robust import audit_robust_signal
from .identity import audit_cell_identity
from .reference import (load_reference_table, compute_detection_ratio,
                         reference_modalities, reference_modality_set,
                         load_reference_with_provenance, load_reference_depth,
                         validate_reference_contract, enforce_reference_contract,
                         CrossModalityError, ReferenceContractError, DepthMismatchError)
from .latent_batch import audit_latent_batch
from .check import check_readiness
from .report import render_report
from .design import (DesignDeclaration, DesignGate, check_design_support,
                     structural_hints, COMPARISON_TYPES, UNITS,
                     EXIT_DESIGN_NOT_SUPPORTED)


DEFAULT_TARGETS = "MYL3,PLN,TNNT2,ACTB,GAPDH"


def _sha256(path: str) -> str:
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return "unavailable"


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="sc-audit",
        description="Detection-floor and control-composition audit for "
                    "case-control single-cell studies.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""examples:
  sc-audit run --input mydata.h5ad --target-genes MYL3,PLN,TNNT2 --outdir out/
  sc-audit floor --input mydata.h5ad --target-genes MYL3 --ratio-basis case
  sc-audit run --input 10x_dir/ --preset cardiac --condition-col group --donor-col sample
""")
    p.add_argument("--version", action="version", version=f"sc-audit {__version__}")

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--input", default=None,
                        help=".h5ad / 10x mtx 目录 / 10x .h5（design apply 不需要）")
    common.add_argument("--condition-col", default="condition", help="obs 中分组列")
    common.add_argument("--donor-col", default="donor", help="obs 中供者列")
    common.add_argument("--control-label", default=None, help="对照组标签（默认自动识别）")
    common.add_argument("--case-label", default=None, help="病例组标签（默认自动识别）")
    common.add_argument("--target-genes", default=DEFAULT_TARGETS, help="逗号分隔的目标基因")
    common.add_argument("--config", default=None, help="配置预设名（cardiac）或 yaml 路径")
    common.add_argument("--preset", default=None,
                        help="细胞门控预设（cardiac / cardiomyocyte_atrial / t_cell）")
    common.add_argument("--cell-subset", default=None,
                        help="额外子集过滤（obs 表达式，如 'lineage_names==\"CM\"'）")
    common.add_argument("--reference", default=None,
                        help="外部参考检出率 CSV（列: gene, pct_ref[_*]）")
    common.add_argument("--n-perm", type=int, default=None, help="置换次数（默认 10000）")
    common.add_argument("--ratio-basis", default=None,
                        choices=["case", "overall", "control"],
                        help="检出比 R 的分子口径（默认 case）")
    common.add_argument("--verdict-basis", default=None,
                        choices=["case", "overall", "control"],
                        help="[已废弃] 判定与 R 现共用 --ratio-basis；此参数保留仅为兼容")
    common.add_argument("--count-layer", default=None, help="存放原始计数的 layer 名")
    common.add_argument("--hetero-col", default=None,
                        help="组间分布差异的基准变量（默认 institution）")
    common.add_argument("--reference-alt", default=None,
                        help="第二张参照表（参照口径敏感性分析，输出 verdict_alt 列）")
    common.add_argument("--reference-alt-mode", default=None,
                        choices=["max", "median"],
                        help="alt 参照表的聚合方式（默认与主表一致，避免把"
                             "“换文件”与“换聚合方式”两个因素混在一起）")
    common.add_argument("--no-count-check", action="store_true",
                        help="跳过原始计数校验（不推荐）")
    common.add_argument("--match-vars", default=None,
                        help="显式指定匹配变量（逗号分隔）；缺省由 config 自动选择")
    common.add_argument("--allow-outcome-matching", action="store_true",
                        help="允许把结局/治疗变量写进 match_vars（默认禁止："
                             "会造成 collider bias）")
    common.add_argument("--allow-cross-modality", action="store_true",
                        help="显式放行跨 modality 参照比较（默认禁止：不同 modality 的"
                             "检出率不可比）。放行后受影响基因的 R 标 cross_modality=True")
    common.add_argument("--allow-depth-mismatch", action="store_true",
                        help="A1：显式放行 case 与参照的深度不匹配（默认禁止：检出率随深度陡变）")
    common.add_argument("--depth-matched-sensitivity", action="store_true",
                        help="A1：把 case 下采样到参照深度，额外输出 depth-matched 的 R/判定列")
    common.add_argument("--with-latent-batch", action="store_true",
                        help="B1：额外跑一层数据驱动的潜在批次扫描（报警层，不改现有判定）")
    # ── 设计声明（v0.4.0 Tier 1）：只记录、不推断；超出支持集 → exit 3 ──
    common.add_argument("--comparison-type", default=None,
                        choices=list(COMPARISON_TYPES),
                        help="设计声明：比较结构（默认 two_arm）。其它值当前未实现，"
                             "声明后会**拒绝输出判定标签**（exit 3）")
    common.add_argument("--unit", default=None, choices=list(UNITS),
                        help="设计声明：分析/推断单位（默认 donor）")
    common.add_argument("--contrast-a", default=None,
                        help="设计声明：比较的 A 组标签（仅记录/留痕）")
    common.add_argument("--contrast-b", default=None,
                        help="设计声明：比较的 B 组标签（仅记录/留痕）")
    common.add_argument("--time-col", default=None,
                        help="设计声明：纵向时间列（当前未实现 → exit 3）")
    common.add_argument("--ref-time", default=None,
                        help="设计声明：纵向参考时间点（仅记录）")
    common.add_argument("--index-time", default=None,
                        help="设计声明：索引时间（仅记录，用于解释“事前/事后”）")
    common.add_argument("--index-event", default=None,
                        help="设计声明：索引事件（仅记录）")
    common.add_argument("--outdir", default="sc_audit_out")
    common.add_argument("--quiet", action="store_true")

    sub = p.add_subparsers(dest="cmd", required=True)
    subs = {}
    for name, help_ in [("run", "运行四项审计并生成报告"),
                        ("floor", "仅审计一：检出地板"),
                        ("composition", "仅审计二：对照构成"),
                        ("identifiability", "审计二 v0.3：可识别性 + 匹配规则阶梯"),
                        ("robust", "仅审计三：稳健信号对照"),
                        ("compare", "跨数据集检出比 R（需 --reference）"),
                        ("identity", "仅审计四：细胞身份与捕获完整性"),
                        ("check", "预检：这份数据能支撑哪些审计（只报告，不采纳）"),
                        ("design", "设计解析（解耦版）：export 出 inventory / apply 读回 proposal + 确定性复验 + 确认单")]:
        sp = sub.add_parser(name, parents=[common], help=help_)
        subs[name] = sp
        # 匹配与管家基因对 run/子命令都可传入（原仅子命令）
        if name in ("run", "composition", "identifiability"):
            sp.add_argument("--matched-donors", default=None,
                            help="与病例组同机构/同取材的对照供者，逗号分隔（默认自动检测）")
        if name in ("run", "identifiability"):
            sp.add_argument("--ladder", default=None,
                            help="事前锁定的匹配规则阶梯：用 ';' 分隔规则、用 '|' 分隔变量；"
                                 "空规则写作 'none'。例：'none;institution|source;source'")
        if name in ("run", "robust"):
            sp.add_argument("--housekeeping", default=None,
                            help="管家基因（默认 ACTB,GAPDH）")
        if name in ("run", "identity"):
            sp.add_argument("--celltype-col", default=None,
                              help="细胞类型注释列（默认读配置 identity_col）")
            sp.add_argument("--source-col", default=None,
                              help="分选/来源列；审计四检查 5/6 的核心")
            sp.add_argument("--target-lineage", default=None,
                              help="目标细胞类型名称（默认读 identity_target）")
        if name in ("run", "floor", "identity"):
            sp.add_argument("--modality", default=None,
                              help="被审计数据集的模态（scRNA/snRNA）；审计一 R2 + 检查 7")
            sp.add_argument("--reference-modalities", default=None,
                              help="参照表已有模态，逗号分隔；审计一 R2 + 检查 7")
        if name == "check":
            sp.add_argument("--design-var-candidates", default=None,
                            help="要评估的候选设计变量（逗号分隔；"
                                 "默认自动扫描 obs 中的类别型列）")
        if name == "design":
            sp.add_argument("--action", choices=["export", "apply", "confirm", "attest", "propose"], required=True,
                            help="export 出 inventory / propose 内置调用器 / apply 读回 / confirm 注册表 / attest 盖出处")
            sp.add_argument("--label", default=None, help="数据集标识（默认取文件名）")
            sp.add_argument("--geotext", default=None,
                            help="（export）GEO 自由文本文件路径，仅作为指针写入 inventory")
            sp.add_argument("--inventory", default=None, help="（apply）inventory.json")
            sp.add_argument("--proposal", default=None, help="（apply）proposal.json（外部 AI 生成）")
            sp.add_argument("--registry", default=None,
                            help="（apply/confirm）列注册表路径（默认 outdir/design_registry.json）")
            sp.add_argument("--decisions", default=None,
                            help="（confirm）人工裁决 JSON：{dataset:{col:role}}")
            sp.add_argument("--model", default=None, help="（attest）生成 proposal 的模型标识")
            sp.add_argument("--prompt", default=None, help="（attest）提示词文件路径（算 sha256）")
            sp.add_argument("--out", default=None, help="（attest/propose）输出文件（默认就地/outdir）")
            sp.add_argument("--backend", choices=["none", "local", "http"], default="none",
                            help="（propose）调用器后端：none（默认，不生成）/ local（本地命令）/ http（OpenAI 兼容端点）")
            sp.add_argument("--endpoint", default=None, help="（propose，http）如 https://api.openai.com/v1")
            sp.add_argument("--api-key-env", default="OPENAI_API_KEY",
                            help="（propose，http）存放 API key 的环境变量名")
            sp.add_argument("--local-cmd", default=None, help="（propose，local）本地命令（stdin 收提示词，stdout 出 JSON）")
            sp.add_argument("--temperature", type=float, default=0.0, help="（propose）采样温度（默认 0）")
            sp.add_argument("--timeout", type=int, default=180, help="（propose/attest）超时秒")
    return p


def _apply_design_args(cfg, args) -> None:
    """把 CLI 的设计声明写入 cfg（缺省不动 → 保持支持集 S0 默认）。"""
    _map = (("comparison_type", "design_comparison_type"),
            ("unit", "design_unit"),
            ("contrast_a", "design_contrast_a"),
            ("contrast_b", "design_contrast_b"),
            ("time_col", "design_time_col"),
            ("ref_time", "design_ref_time"),
            ("index_time", "design_index_time"),
            ("index_event", "design_index_event"))
    for arg_name, field_name in _map:
        val = getattr(args, arg_name, None)
        if val is not None:
            setattr(cfg, field_name, val)


def _design_gate_or_exit(cfg, n_condition_levels, quiet) -> DesignGate:
    """评估能力闸门（v0.4.0 Tier 1）。

    · 声明**非法**（取值不在封闭词表内） → [error] + SystemExit(2)
    · 声明合法但**超出支持集**            → 拒绝文案 + SystemExit(3)
    · 在支持集内                        → 返回 DesignGate（调用方继续）
    """
    decl = DesignDeclaration.from_config(cfg)
    bad = decl.invalid_reasons()
    if bad:
        print("[error] 设计声明非法（取值不在封闭词表内）：\n         - "
              + "\n         - ".join(bad)
              + f"\n         合法取值：comparison_type ∈ {list(COMPARISON_TYPES)}；"
                f"unit ∈ {list(UNITS)}",
              file=sys.stderr)
        raise SystemExit(2)
    gate = check_design_support(decl, n_condition_levels=n_condition_levels)
    if not gate.supported:
        print(gate.refusal_text(), file=sys.stderr)
        raise SystemExit(EXIT_DESIGN_NOT_SUPPORTED)
    if not quiet and not decl.is_default():
        print(gate.to_text())
    return gate


def _merged_ref_modalities(args) -> Optional[List[str]]:
    """参照模态集合 = --reference-modalities ∪ 参照表声明的模态（CSV `modality` 列 / manifest）。

    统一两处来源，使 floor 的 R2 判定（同/跨模态）与 CLI 层契约检查一致。
    """
    mods: List[str] = []
    if getattr(args, "reference_modalities", None):
        mods += [m.strip() for m in args.reference_modalities.split(",") if m.strip()]
    if getattr(args, "reference", None):
        try:
            mods += reference_modality_set(args.reference)
        except Exception:  # noqa: BLE001
            pass
    mods = sorted(set(mods))
    return mods or None


def _outcome_like_match_vars(ad, cfg) -> List[str]:
    """返回 match_vars 中角色为结局/治疗（或由配置显式声明）的列。

    判据与 check 一致：显式声明(role_overrides) > 领域词根 > 内置词根。
    仅看**列名/声明**，不使用分布统计（避免把取材变量误判成结局）。
    """
    from .check import _name_role
    _zh = {"outcome": "结局", "outcome_suspect": "疑似结局",
           "intervention": "治疗/干预"}
    bad: List[str] = []
    for c in (cfg.match_vars or []):
        if c not in ad.obs.columns:
            continue
        r = _name_role(c, cfg)
        if r in _zh:
            bad.append(f"{c}({_zh[r]})")
    return bad


def _prepare(args):
    cfg = load_config(args.config)
    _apply_design_args(cfg, args)
    if getattr(args, "allow_cross_modality", False):
        cfg.allow_cross_modality = True
    if getattr(args, "allow_depth_mismatch", False):
        cfg.allow_depth_mismatch = True
    if args.n_perm is not None:
        cfg.n_perm = args.n_perm
    if getattr(args, "verdict_basis", None):
        import warnings as _w
        _w.warn("--verdict-basis 已废弃：判定与 R 现共用 --ratio-basis；"
                "该参数被忽略", DeprecationWarning, stacklevel=2)
    if getattr(args, "ratio_basis", None):
        cfg.ratio_basis = args.ratio_basis
    # N4 修复：verdict 与 R 共用 cfg.ratio_basis（原 verdict_basis 字段已删除）
    if getattr(args, "count_layer", None):
        cfg.count_layer = args.count_layer
    if getattr(args, "hetero_col", None):
        cfg.hetero_col = args.hetero_col
    if getattr(args, "match_vars", None):
        # 显式覆盖：默认由 _pick_match_vars 从 cfg.match_vars 自动挑选，
        # 但使用者应能直接声明，不必改 config、也不依赖自动选择。
        cfg.match_vars = [c.strip() for c in args.match_vars.split(",") if c.strip()]
    # 审计二 v0.3：显式给出的匹配阶梯 = 事前锁定，写入 cfg 以留痕（run_metadata.json）
    if getattr(args, "ladder", None) is not None:
        _raw = (args.ladder or "").strip()
        if not _raw:
            print("[error] --ladder 给了空值：请写成 'none;institution|source' 形式"
                  "（'none' 代表不匹配的锚点规则），或直接省略该参数使用内置阶梯。",
                  file=sys.stderr)
            raise SystemExit(2)
        rules: List[List[str]] = []
        for chunk in _raw.split(";"):
            chunk = chunk.strip()
            if not chunk or chunk.lower() in ("none", "null", "-"):
                rules.append([])
            else:
                rules.append([c.strip() for c in chunk.split("|") if c.strip()])
        if len(rules) < 2:
            print("[error] --ladder 需至少给两条规则（含 R0 锚点与一条匹配规则），"
                  "否则无法做一致性检查。", file=sys.stderr)
            raise SystemExit(2)
        cfg.identifiability_ladder = rules

    ad = load_data(args.input)

    # 可识别性层不使用 `--matched-donors`（它的匹配由事前锁定的阶梯决定）。
    # 历史缺陷：该参数曾被 argparse 注册给 identifiability 但**从未生效**（静默 no-op）。
    # 现显式拒绝，避免“解析了但不生效”。
    if args.cmd == "identifiability" and getattr(args, "matched_donors", None):
        print("[error] `identifiability` 不接受 --matched-donors：它的匹配集合由 "
              "`--ladder`（事前锁定的规则阶梯）决定，没有单独的“指定对照供者”入口。\n"
              "       请在 `composition`/`run` 上用 --matched-donors，"
              "在 `identifiability`/`run` 上用 --ladder。", file=sys.stderr)
        raise SystemExit(2)

    # ── B5：确保原始计数 ──
    counts_source = "X (assumed)"
    switched = False
    if not args.no_count_check:
        ad, counts_source = resolve_counts(ad, cfg)
        switched = counts_source not in ("X", "X (assumed)")

    # ── 基因名回退（P1-a：必须早于 compute_qc_metrics）──
    # compute_qc_metrics 靠 var_names 的 "MT-" 前缀找线粒体基因；
    # 若此时 var_names 仍是 ENSG…，percent_mito 会被算成 NaN（并在 force=True 时
    # 覆盖掉上游正确的值）。故替换必须上移到这里。
    gene_info = {"missing": [], "applied": False}
    if args.target_genes:
        req = [g.strip() for g in args.target_genes.split(",") if g.strip()]
        # P3：探针集合 = 目标基因 ∪ 管家基因 ∪ 门控 marker。
        # 若只用 --target-genes 做探针，用户直接传 ENSG ID 时回退不触发，
        # 而门控 marker（TTN 等符号）全部缺失 → 所有细胞被过滤掉。
        probe = list(req)
        _hk_probe = (getattr(args, "housekeeping", None)
                     or ",".join(getattr(cfg, "housekeeping_genes", []) or []))
        probe += [g.strip() for g in _hk_probe.split(",") if g.strip()]
        if args.preset:
            _pp = preset_from_config(cfg, args.preset)
            probe += list(_pp.get("pos_markers", [])) + list(_pp.get("neg_markers", []))
        probe = list(dict.fromkeys(probe))          # 去重保序
        gene_info = apply_gene_name_fallback(ad, probe)
        if gene_info.get("applied") and not args.quiet:
            print(f"[genes] var_names replaced with '{gene_info['fallback_column']}' "
                  f"(original index kept in var['_orig_var_names'])")
        if gene_info.get("missing"):
            print(f"[warn] 以下基因不在数据集中：{gene_info['missing']}",
                  file=sys.stderr)

    info = ensure_obs_fields(ad, args.condition_col, args.donor_col,
                             args.control_label, args.case_label)

    # ── QC（P1-b：仅在【真正换源】且 config 允许时强制重算）──
    # 原写法 force=(counts_source != "X (assumed)") 恒为 True（resolve_counts 正常返回 "X"），
    # 会对 HVG 子集数据用子集口径覆盖上游的全转录组 n_genes，系统性偏小。
    qc_recomputed = bool(switched and cfg.recompute_qc_after_source_switch)
    ad = compute_qc_metrics(ad, force=qc_recomputed)

    gate_info = None
    if args.preset:
        pr = preset_from_config(cfg, args.preset)      # 单一真值源（cfg）
        mask = gate_cells(ad, pr["pos_markers"], pr["neg_markers"],
                          min_pos_markers=pr["min_pos_markers"],
                          min_counts=pr["min_counts"], min_genes=pr["min_genes"])
        ad = ad[mask].copy()
        gate_info = ad.uns.get("sc_audit_gating")
        # 门控把细胞全过滤掉时给出可诊断的中止，而不是让下游报"标签不存在"
        if ad.n_obs == 0:
            print(f"[error] --preset {args.preset} 门控后没有细胞剩余 "
                  f"(通过 {gate_info.get('n_pass')}/{gate_info.get('n_total')})。"
                  f"请放宽阈值（config 中的 min_counts / min_genes / min_pos_markers）"
                  f"或改用与本数据集匹配的 marker 面板。", file=sys.stderr)
            raise SystemExit(2)
        info = ensure_obs_fields(ad, args.condition_col, args.donor_col,
                                 info["control_label"], info.get("case_label"))

    if args.cell_subset:
        try:
            ad = ad[ad.obs.eval(args.cell_subset)].copy()
        except Exception as exc:
            # N5 修复：用户显式要求子集却失败时，必须中止而不是在全量数据上出结论
            print(f"[error] --cell-subset 评估失败：{type(exc).__name__}: {exc}\n"
                  f"        可用列：{list(ad.obs.columns)}\n"
                  f"        已中止，未对全量数据做任何审计。", file=sys.stderr)
            raise SystemExit(2)
        if ad.n_obs == 0:
            print("[error] --cell-subset 过滤后无剩余细胞，已中止。", file=sys.stderr)
            raise SystemExit(2)
        info = ensure_obs_fields(ad, args.condition_col, args.donor_col,
                                 info["control_label"], info.get("case_label"))

    ref = load_reference_table(args.reference) if args.reference else None
    # ── 参照契约（v0.5.0；R1 五元组合并 + R2 跨 modality）──
    # 违规即 [error] + exit 2（main 捕获 ValueError）。
    if args.reference:
        _rc = validate_reference_contract(
            args.reference, target_modality=getattr(args, "modality", None),
            allow_cross_modality=getattr(cfg, "allow_cross_modality", False))
        enforce_reference_contract(_rc)
        cfg._reference_contract = _rc.to_dict()
    # P2-b：alt 表默认用**与主表相同的聚合方式**，使翻转数可归因于"换文件"而非"换口径"
    alt_mode = getattr(args, "reference_alt_mode", None) or "max"
    ref_alt = (load_reference_table(args.reference_alt, mode=alt_mode)
               if getattr(args, "reference_alt", None) else None)
    genes = [g.strip() for g in args.target_genes.split(",") if g.strip()]
    # P2-2：用户显式指定 --housekeeping 时，阳性对照集必须同步到 FDR 的排除集 m；
    # 否则（如 B2M）会被排除出“分离持续”计数，却仍被算进 FDR 的假设数 m。
    if getattr(args, "housekeeping", None):
        cfg.fdr_exclude_from_m = [g.strip()
                                  for g in args.housekeeping.split(",") if g.strip()]

    # ── 运行期保护：match_vars 命中结局/治疗变量 → 中止（防 collider bias）──
    if getattr(cfg, "forbid_outcome_matching", True) and \
            not getattr(args, "allow_outcome_matching", False):
        bad = _outcome_like_match_vars(ad, cfg)
        if bad:
            print(f"[error] match_vars 含疑似结局/治疗变量：{bad}\n"
                  f"        对这些变量做匹配 = 条件化在结局上，会引入 collider bias，"
                  f"使暴露-结局关系无法解释。\n"
                  f"        若已确认语义（例如该列在两组中含义相同），"
                  f"用 --allow-outcome-matching 显式放行；\n"
                  f"        否则请从 match_vars 移除，或先用 `sc-audit check` 查看建议。",
                  file=sys.stderr)
            raise SystemExit(2)
    # 参照溯源：记下每个基因的参照值取自哪一列（参照溯源）
    _ref_prov = {}
    if args.reference:
        try:
            _, _ref_prov = load_reference_with_provenance(args.reference,
                                                          mode="max")
        except Exception:
            _ref_prov = {}
    return (cfg, ad, info, ref, genes, gate_info, counts_source, gene_info,
            ref_alt, qc_recomputed, alt_mode, _ref_prov)


def _run_check(args) -> int:
    """预检：只报告，不采纳。

    刻意绕过 _prepare()：_prepare 在 obs 解析失败时会抛异常，
    而预检命令的价值恰恰在于**把失败本身作为诊断内容报出来**。
    退出码：0 = 三项均可行；1 = 有审计不可行；2 = 存在阻断项。
    """
    cfg = load_config(args.config)
    _apply_design_args(cfg, args)
    if getattr(args, "allow_cross_modality", False):
        cfg.allow_cross_modality = True
    if getattr(args, "allow_depth_mismatch", False):
        cfg.allow_depth_mismatch = True
    if args.n_perm is not None:
        cfg.n_perm = args.n_perm
    # v0.3.3：预检必须与正式运行用**同一套**参数覆盖。旧实现只覆盖 n_perm，
    # 于是 `check --count-layer X` 被静默忽略 → 预检比正式运行更严，会误报“阻断”。
    if getattr(args, "count_layer", None):
        cfg.count_layer = args.count_layer
    if getattr(args, "ratio_basis", None):
        cfg.ratio_basis = args.ratio_basis
    if getattr(args, "hetero_col", None):
        cfg.hetero_col = args.hetero_col
    if getattr(args, "match_vars", None):
        cfg.match_vars = [c.strip() for c in args.match_vars.split(",") if c.strip()]
    ad = load_data(args.input)

    # 能力闸门（v0.4.0 Tier 1）：预检同样受约束——声明超出支持集时拒绝出"可行性"判定。
    try:
        _nlev = (int(ad.obs[args.condition_col].astype(str).nunique())
                 if args.condition_col in ad.obs.columns else None)
    except Exception:  # noqa: BLE001
        _nlev = None
    _design_gate_or_exit(cfg, _nlev, args.quiet)

    gate_info = None
    if args.preset:
        pr = preset_from_config(cfg, args.preset)
        mask = gate_cells(ad, pr["pos_markers"], pr["neg_markers"],
                          min_pos_markers=pr["min_pos_markers"],
                          min_counts=pr["min_counts"], min_genes=pr["min_genes"])
        ad = ad[mask].copy()
        gi = ad.uns.get("sc_audit_gating")
        if gi is not None:
            gate_info = dict(gi)
            gate_info["preset"] = args.preset

    if getattr(args, "cell_subset", None):
        try:
            ad = ad[ad.obs.eval(args.cell_subset)].copy()
        except Exception as exc:                      # noqa: BLE001
            print(f"[error] --cell-subset 评估失败：{type(exc).__name__}: {exc}",
                  file=sys.stderr)
            return 2

    ref = load_reference_table(args.reference) if args.reference else None
    # ── 参照契约（v0.5.0；R1 五元组合并 + R2 跨 modality）——预检同样受约束 ──
    if args.reference:
        enforce_reference_contract(validate_reference_contract(
            args.reference, target_modality=getattr(args, "modality", None),
            allow_cross_modality=getattr(cfg, "allow_cross_modality", False)))
        # M4（v0.6.2）：预检也检查参照深度可比性（与正式 run 一致）
        try:
            _rd = load_reference_depth(args.reference)
        except Exception:  # noqa: BLE001
            _rd = None
        _cd = None
        try:
            if "total_counts" in ad.obs.columns:
                _cd = float(pd.to_numeric(ad.obs["total_counts"],
                                          errors="coerce").median())
        except Exception:  # noqa: BLE001
            _cd = None
        if _rd and _cd:
            _ddp = (_cd - _rd) / _rd * 100.0
            if abs(_ddp) > float(getattr(cfg, "reference_depth_delta_max_pct", 50.0)) \
                    and not getattr(cfg, "allow_depth_mismatch", False):
                raise DepthMismatchError(
                    f"[check] case 深度（中位 total_counts={_cd:.0f}）与参照深度"
                    f"（median/mean_umi={_rd:.0f}）相差 {_ddp:+.1f}%，超事前锁定阈值"
                    f" {getattr(cfg, 'reference_depth_delta_max_pct', 50.0)}%。"
                    "检出率是深度的陡函数，跨深度算 R 会被系统性污染；"
                    "请换深度可比的参照表，或用 --allow-depth-mismatch 显式放行。")
    # 参照表含多模态列时主动告警：跨模态算 R 会系统性偏低，但仍可运行。
    if args.reference:
        try:
            _mods = reference_modalities(args.reference)
        except Exception:
            _mods = []
        if len(_mods) > 1:
            warnings.warn(
                f"参考表含 {len(_mods)} 个模态 {_mods}；跨模态计算的检出比 R "
                "系统性偏低，绝对值不可解释为生物学量。建议改用同模态参照表，"
                "或用 --modality / --reference-modalities 声明后由审计四标注。",
                stacklevel=2)
    genes = [g.strip() for g in args.target_genes.split(",") if g.strip()]
    cand = None
    if getattr(args, "design_var_candidates", None):
        cand = [c.strip() for c in args.design_var_candidates.split(",") if c.strip()]

    # 参照表旁若有 manifest（门控参数登记）则一并读取
    _man = None
    if args.reference:
        _mp = args.reference + ".manifest.json"
        if not os.path.exists(_mp):
            _mp = os.path.splitext(args.reference)[0] + ".manifest.json"
        if os.path.exists(_mp):
            with open(_mp, encoding="utf-8") as _fh:
                _man = json.load(_fh)
    rr = check_readiness(ad, args.condition_col, args.donor_col,
                         args.control_label, cfg=cfg, reference=ref,
                         candidate_vars=cand, target_genes=genes,
                         check_counts=not args.no_count_check,
                         reference_manifest=_man)
    rr.gating = gate_info

    txt = rr.to_text()
    with open(os.path.join(args.outdir, "check_report.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(txt + "\n")
    with open(os.path.join(args.outdir, "check_report.json"), "w",
              encoding="utf-8") as fh:
        json.dump(rr.to_json(), fh, ensure_ascii=False, indent=2)
    if not args.quiet:
        print(txt)

    if rr.blockers:
        return 2
    if any(v.get("level") == "no" for v in rr.feasibility.values()):
        return 1
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    """CLI 入口：把**输入侧**错误统一成 `[error] …` + exit 2（不再裸 traceback）。

    2026-10-07 发现：参数类错误走 `[error] …` + exit 2，而数据/文件类错误
    （obs 缺列 → KeyError、输入不存在 → FileNotFoundError、单水平 condition →
    ValueError）直接抛栈 + exit 1，风格不一致。此处统一；需要堆栈时设
    `SC_AUDIT_DEBUG=1`（开发/排查用）。
    """
    args = _build_parser().parse_args(argv)
    try:
        return _main_dispatch(args)
    except (FileNotFoundError, NotADirectoryError, KeyError, ValueError) as e:
        if os.environ.get("SC_AUDIT_DEBUG"):
            raise
        print(f"[error] {type(e).__name__}: {e}", file=sys.stderr)
        print("[hint] 先用 `sc-audit check --input <文件> --condition-col <列> "
              "--donor-col <列>` 预检输入；要看完整堆栈请设 SC_AUDIT_DEBUG=1。",
              file=sys.stderr)
        raise SystemExit(2)


def _run_design(args) -> int:
    """`design` 子命令：解耦式设计解析（**不调 LLM、不联网**）。"""
    from .design_parse import (export_inventory, write_inventory, apply_proposal,
                               confirm_registry)
    if args.action == "confirm":
        return _run_design_confirm(args)
    if args.action == "attest":
        return _run_design_attest(args)
    if args.action == "propose":
        return _run_design_propose(args)
    if args.action == "export":
        if not args.input:
            print("[error] design export 需要 --input", file=sys.stderr)
            return 2
        ad = load_data(args.input)
        if args.control_label is None:
            print("[error] design export 需要 --control-label（对照组的取值）", file=sys.stderr)
            return 2
        inv = export_inventory(ad, args.condition_col, args.donor_col,
                               args.control_label, label=args.label,
                               geotext=args.geotext)
        p = write_inventory(inv, args.outdir)
        if not args.quiet:
            print(f"-> {p}  ({len(inv['columns'])} 列, {inv['n_donors']} 供者)")
            print("   下一步：把 inventory.json 交给 AI（任一大模型/agent）生成 proposal.json；"
                  "再 `sc-audit design apply --inventory inventory.json --proposal proposal.json`")
        return 0
    # apply
    if not args.inventory or not args.proposal:
        print("[error] design apply 需要 --inventory 与 --proposal", file=sys.stderr)
        return 2
    with open(args.inventory, encoding="utf-8") as fh:
        inv = json.load(fh)
    with open(args.proposal, encoding="utf-8") as fh:
        prop = json.load(fh)
    reg = args.registry or os.path.join(args.outdir, "design_registry.json")
    s = apply_proposal(inv, prop, args.outdir, registry_path=reg)
    if not args.quiet:
        print(f"-> {os.path.join(args.outdir, 'design_confirm_sheet.md')}")
        print(f"-> {os.path.join(args.outdir, 'design_config.json')}  （待人工确认，不自动写入）")
        print(f"   {s['n_columns']} 列：必须逐条确认 {s['n_must_review']} / 批量 {s['n_batch']} / "
              f"自动 {s['n_auto']}；已确认 {s['n_confirmed']}；弃权 {s['n_abstain']}；"
              f"**复验不一致 {s['n_disagree']}**")
    return 0


def _run_design_confirm(args) -> int:
    from .design_parse import confirm_registry
    if not args.decisions:
        print("[error] design confirm 需要 --decisions {dataset:{col:role}}", file=sys.stderr)
        return 2
    reg = args.registry or os.path.join(args.outdir, "design_registry.json")
    with open(args.decisions, encoding="utf-8") as fh:
        dec = json.load(fh)
    n = confirm_registry(reg, dec)
    if not args.quiet:
        print(f"-> {reg}：{n} 条人工裁决已标 confirmed（下次 apply 不再重复问）")
    return 0


def _run_design_attest(args) -> int:
    from .design_parse import attest_proposal
    if not args.proposal or not args.inventory or not args.model:
        print("[error] design attest 需要 --proposal --inventory --model", file=sys.stderr)
        return 2
    gb = attest_proposal(args.proposal, args.inventory, args.model,
                         prompt_path=args.prompt, out_path=args.out)
    if not args.quiet:
        print(f"-> 已盖出处块：model={gb['model']} prompt_sha256={str(gb['prompt_sha256'])[:12]} "
              f"inventory_sha256={str(gb['inventory_sha256'])[:12]}（{'就地' if not args.out else args.out}）")
    return 0


def _run_design_propose(args) -> int:
    from .design_parse import propose
    if not args.inventory:
        print("[error] design propose 需要 --inventory", file=sys.stderr)
        return 2
    if args.backend == "none":
        print("[error] --backend none：本命令不生成提案；请自行生成 proposal.json 后用 apply；"
              "或用 --backend local|http。", file=sys.stderr)
        return 2
    with open(args.inventory, encoding="utf-8") as fh:
        inv = json.load(fh)
    api_key = os.environ.get(args.api_key_env) if args.api_key_env else None
    out = args.out or os.path.join(args.outdir, "proposal.json")
    r = propose(inv, backend=args.backend, model=args.model, endpoint=args.endpoint,
                api_key=api_key, local_cmd=args.local_cmd, temperature=args.temperature,
                timeout=args.timeout, prompt_path=args.prompt, out_path=out)
    if not args.quiet:
        print(f"-> {r['out']}（backend={args.backend}, model={r['generated_by']['model']}, "
              f"inventory_sha256={str(r['generated_by']['inventory_sha256'])[:12]}）；"
              "下一步：design --action apply")
    return 0


def _main_dispatch(args) -> int:
    os.makedirs(args.outdir, exist_ok=True)
    if args.cmd == "check":
        return _run_check(args)
    if args.cmd == "design":
        return _run_design(args)
    if not args.input:
        print("[error] 缺少 --input", file=sys.stderr)
        return 2
    (cfg, ad, info, ref, genes, gate_info, counts_src, gene_info,
     ref_alt, qc_recomputed, alt_mode, _ref_prov) = _prepare(args)
    # 空面板守卫（公开就绪度复核 2026-10-07）：`--target-genes ""` 会在审计一层
    # 抛未处理的 KeyError。此处与 `--ladder ""` 的处理风格一致：明确报错 + exit 2。
    if not genes and args.cmd in ("run", "floor", "robust", "identifiability", "compare"):
        print("[error] --target-genes 为空：请给至少一个目标基因（逗号分隔），"
              "例如 --target-genes MYL3,PLN。只想做预检请用 `sc-audit check`。",
              file=sys.stderr)
        raise SystemExit(2)
    ctrl = info["control_label"]

    # ── 能力闸门（v0.4.0 Tier 1）──
    # 声明超出支持集 → 拒绝出任何判定标签（exit 3），不做审计。
    # 结构提示**只警告，不改判定**（程序不推断，只索取）。
    gate = _design_gate_or_exit(cfg, info.get("n_condition_levels"), args.quiet)
    _design_hints = structural_hints(ad, args.condition_col, args.donor_col,
                                     cfg.design_time_col)
    if _design_hints and not args.quiet:
        for _h in _design_hints:
            print(f"[design-hint] {_h}", file=sys.stderr)

    if gate_info and not args.quiet:
        print(f"[gating] preset={args.preset}  passed {gate_info['n_pass']}/{gate_info['n_total']} "
              f"({gate_info['pct_pass']}%)")
    if not args.quiet and counts_src not in ("X", "X (assumed)"):
        print(f"[counts] using raw counts from {counts_src}")

    floor = comp = rob = ratio_tbl = None
    latent_batch_res = None
    # A1（v0.6.0）：case 深度（中位 total_counts）与参照深度（mean_umi），供深度可比性检查
    _case_depth = None
    try:
        if "total_counts" in ad.obs.columns:
            _case_depth = float(pd.to_numeric(ad.obs["total_counts"],
                                              errors="coerce").median())
    except Exception:  # noqa: BLE001
        _case_depth = None
    _ref_depth = None
    if args.reference:
        try:
            _ref_depth = load_reference_depth(args.reference)
        except Exception:  # noqa: BLE001
            _ref_depth = None

    if args.cmd in ("run", "floor"):
        floor = audit_detection_floor(
            ad, genes, args.condition_col, args.donor_col, ctrl, cfg=cfg,
            reference=ref,
            # 参照溯源：把模态与参照溯源一路传到输出表，
            # 使 R 的分母"来自哪列、哪个模态"在报告里自明
            modality=getattr(args, "modality", None),
            reference_modalities=_merged_ref_modalities(args),
            reference_provenance=_ref_prov,
            allow_cross_modality=getattr(cfg, "allow_cross_modality", False),
            reference_depth=_ref_depth, case_depth=_case_depth,
            allow_depth_mismatch=getattr(cfg, "allow_depth_mismatch", False),
            depth_matched_sensitivity=bool(getattr(args, "depth_matched_sensitivity",
                                                   False)))
        # L8：第二张参照表 → 参照口径敏感性（输出 verdict_alt 列）
        if ref_alt and not ref:
            # P3：没有主参照时主判定全为 INSUFFICIENT (no reference)，
            # 此时算"翻转数"会把"缺少参照"误算成"口径敏感"，虚高且误导
            print("[warn] 已给出 --reference-alt 但缺少 --reference；"
                  "主判定无参照，跳过参照口径敏感性统计。", file=sys.stderr)
        elif ref_alt:
            alt = audit_detection_floor(ad, genes, args.condition_col, args.donor_col,
                                        ctrl, cfg=cfg, reference=ref_alt,
                                        run_jackknife=False)
            amap = dict(zip(alt.table["gene"], alt.table["verdict"]))
            rmap = dict(zip(alt.table["gene"], alt.table["ratio_R"]))
            floor.table["verdict_alt"] = floor.table["gene"].map(amap)
            floor.table["ratio_R_alt"] = floor.table["gene"].map(rmap)
            n_flip = int((floor.table["verdict_alt"].notna() &
                          (floor.table["verdict"] != floor.table["verdict_alt"])).sum())
            floor.summary["reference_alt_flips"] = n_flip
            if not args.quiet:
                print(f"[reference-alt] {n_flip} verdict(s) change under the alternative "
                      f"reference table")
        floor.table.to_csv(os.path.join(args.outdir, "floor_audit.csv"), index=False)
        if floor.jackknife:
            rows = []
            for g, j in floor.jackknife.items():
                if "table" in j and len(j["table"]):
                    t = j["table"].copy(); t.insert(0, "gene", g); rows.append(t)
            if rows:
                pd.concat(rows, ignore_index=True).to_csv(
                    os.path.join(args.outdir, "floor_jackknife.csv"), index=False)
        if not args.quiet:
            print(floor.to_text())

    if args.cmd in ("run", "composition"):
        matched = None
        if getattr(args, "matched_donors", None):
            matched = [d.strip() for d in args.matched_donors.split(",") if d.strip()]
        comp = audit_control_composition(ad, genes, args.condition_col, args.donor_col,
                                         ctrl, matched_donors=matched, cfg=cfg)
        with open(os.path.join(args.outdir, "composition.txt"), "w", encoding="utf-8") as fh:
            fh.write(comp.to_text() + "\n")
        if not args.quiet:
            print("\n" + comp.to_text())

    # 审计二 · v0.3 可识别性层（与 composition 同源，但判定与产物独立）
    ident_res = None
    if args.cmd in ("run", "identifiability"):
        ident_res = audit_identifiability(
            ad, genes, args.condition_col, args.donor_col, ctrl, cfg=cfg,
            ladder=getattr(cfg, "identifiability_ladder", None))
        with open(os.path.join(args.outdir, "identifiability.txt"), "w",
                  encoding="utf-8") as fh:
            fh.write(ident_res.to_text() + "\n")
        ident_res.roster.to_csv(os.path.join(args.outdir, "identifiability_roster.csv"),
                                index=False)
        if len(ident_res.ladder_summary):
            ident_res.ladder_summary.to_csv(
                os.path.join(args.outdir, "identifiability_ladder_summary.csv"),
                index=False)
        if len(ident_res.ladder):
            ident_res.ladder.to_csv(
                os.path.join(args.outdir, "identifiability_ladder.csv"), index=False)
        if not args.quiet:
            print("\n" + ident_res.to_text())
        if args.cmd == "identifiability":
            return 0

    if args.cmd in ("run", "robust"):
        # P3：默认取 cfg.housekeeping_genes（单一真值源），
        # 不再硬编码 "ACTB,GAPDH"（避免与 fdr_exclude_from_m 默认值两份独立维护）
        hk = [g.strip() for g in
              (getattr(args, "housekeeping", None)
               or ",".join(getattr(cfg, "housekeeping_genes", []) or ["ACTB", "GAPDH"])
               ).split(",") if g.strip()]
        # 模态匹配结论：identity 在 robust 之后才跑，故此处由**同样的参数**
        # 独立算出（identity 内部也是这么算的），使管家护栏的解读文案
        # 能直接给出对应的那一句，而不是两种解读并列让读者自己判断。
        _mod_match = None
        if getattr(args, "modality", None) and getattr(args, "reference_modalities", None):
            _mods = [m.strip() for m in args.reference_modalities.split(",") if m.strip()]
            _mod_match = args.modality in _mods
        rob = audit_robust_signal(ad, genes, args.condition_col, args.donor_col,
                                  ctrl, housekeeping=hk, reference=ref, cfg=cfg,
                                  modality_matched=_mod_match)
        with open(os.path.join(args.outdir, "robustness.txt"), "w", encoding="utf-8") as fh:
            fh.write(rob.to_text() + "\n")
        if not args.quiet:
            print("\n" + rob.to_text())

    if args.cmd in ("run", "identity"):
        if getattr(args, "celltype_col", None):
            cfg.identity_col = args.celltype_col
        if getattr(args, "source_col", None):
            cfg.identity_source_col = args.source_col
        if getattr(args, "target_lineage", None):
            cfg.identity_target = args.target_lineage
        ref_mods = ([m.strip() for m in args.reference_modalities.split(",")]
                    if getattr(args, "reference_modalities", None) else None)
        ident = audit_cell_identity(ad, args.condition_col, args.donor_col, ctrl,
                                    cfg=cfg, reference_modalities=ref_mods,
                                    modality=getattr(args, "modality", None))
        with open(os.path.join(args.outdir, "identity.txt"), "w", encoding="utf-8") as fh:
            fh.write(ident.to_text() + "\n")
        if not args.quiet:
            print(ident.to_text())
        if args.cmd == "identity":
            return 0

    if args.cmd == "compare":
        if not ref:
            print("[error] compare 需要 --reference", file=sys.stderr)
            return 2
        from .io import gene_detection_vector
        cond = ad.obs[args.condition_col].astype(str).values
        is_ctrl = cond == ctrl
        study = {}
        for g in genes:
            v = gene_detection_vector(ad, g)
            if v is None:
                continue
            study[g] = float({
                "case": v[~is_ctrl].mean() * 100,
                "overall": v.mean() * 100,
                "control": v[is_ctrl].mean() * 100,
            }[cfg.ratio_basis])
        # P3（P1）：reference_mode 形参已在 reference.py 中删除
        # （聚合方式在 load_reference_table 构建参照时就已决定）；
        # 此处若继续传该关键字会 TypeError。compare 子命令原先没有测试，故未被发现。
        ratio_tbl = compute_detection_ratio(study, ref)
        ratio_tbl.to_csv(os.path.join(args.outdir, "detection_ratio.csv"), index=False)
        if not args.quiet:
            print(ratio_tbl.to_string(index=False))

    if getattr(args, "with_latent_batch", False):
        # B1（v0.6.0）：数据驱动的潜在批次扫描（报警层；不改现有判定）
        latent_batch_res = audit_latent_batch(ad, args.condition_col, args.donor_col,
                                              ctrl, cfg=cfg)
        with open(os.path.join(args.outdir, "latent_batch.txt"), "w",
                  encoding="utf-8") as fh:
            fh.write(latent_batch_res.to_text() + "\n")
        if not args.quiet:
            print("\n" + latent_batch_res.to_text())

    if args.cmd == "run":
        report = render_report(floor=floor, composition=comp, robust=rob,
                               reference_table=ratio_tbl, outdir=args.outdir,
                               title="sc-audit report", input_name=args.input,
                               info=info, identifiability=ident_res,
                               latent_batch=latent_batch_res)
        meta = {
            "input": args.input,
            "input_sha256": _sha256(args.input) if os.path.isfile(args.input) else "n/a (directory)",
            "condition_col": args.condition_col, "donor_col": args.donor_col,
            "control_label": ctrl, "case_label": info.get("case_label"),
            "case_scope": info.get("case_scope"),
            "n_condition_levels": info.get("n_condition_levels"),
            "condition_levels": info.get("condition_levels"),
            "n_cells": int(ad.n_obs), "n_donors": info["n_donors"],
            "gating": gate_info, "counts_source": counts_src,
            "qc_recomputed": bool(qc_recomputed),
            "n_vars": int(ad.n_vars),
            "n_vars_note": ("QC metrics reflect this gene set only; if the input is "
                            "an HVG/subset file they will be systematically smaller "
                            "than whole-transcriptome values"),
            "gene_lookup": gene_info,
            "config": cfg.to_dict(),
            "identifiability_ladder": (cfg.identifiability_ladder
                                        if getattr(cfg, "identifiability_ladder", None)
                                        else "built-in R0-R5"),
            "identifiability_verdict": (ident_res.verdict if ident_res is not None else None),
            "perm_mode": (floor.summary.get("perm_mode")
                          if floor is not None and isinstance(floor.summary, dict) else None),
            "config_source": getattr(cfg, "_source_path", "built-in defaults"),
            "reference_file": args.reference,
            "reference_sha256": _sha256(args.reference) if args.reference and
                                os.path.isfile(args.reference) else None,
            "reference_n_genes": (len(ref) if ref else 0),
            "reference_alt_file": getattr(args, "reference_alt", None),
            "reference_alt_sha256": (_sha256(args.reference_alt)
                                     if getattr(args, "reference_alt", None)
                                     and os.path.isfile(args.reference_alt) else None),
            "reference_alt_mode": alt_mode,
            "reference_alt_n_genes": (len(ref_alt) if ref_alt else 0),
            # 参照契约（v0.5.0；R1/R2）：附加字段，不改既有键含义
            "reference_contract": getattr(cfg, "_reference_contract", None),
            "latent_batch_verdict": (latent_batch_res.verdict
                                     if latent_batch_res is not None else None),
            # 设计声明与能力闸门（v0.4.0 Tier 1）：附加字段，不改既有键含义
            "design": {
                "declared": gate.declared,
                "supported": gate.supported,
                "status": gate.status,
                "supported_set": gate.supported_set,
                "structural_hints": _design_hints,
            },
            "sc_audit_version": __version__,
        }
        with open(os.path.join(args.outdir, "run_metadata.json"), "w", encoding="utf-8") as fh:
            json.dump(meta, fh, indent=2, ensure_ascii=False)
        if not args.quiet:
            print(f"\n-> {report}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
