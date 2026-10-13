"""防再发：断言"文档/修复表所声称的"必须在代码中**可机械验证**。

为什么单独成文件：
  前几轮的教训是——修复只写在文档的**叙述文字**里，而验证方式是
  `symbol in doc_text`；由于叙述文字本身包含那些符号，搜索永远为真，
  于是"联✅"全是假绿灯（N1–L8 全部如此）。
  本文件改为在**字段表 / 源码 / 默认值**层面断言，无法被叙述文字蒙混。

对应关系：
  N1  cli 真正替换 var_names         N2  匹配变量剔除判据（供者级）
  N3  yaml 与默认值逐字段一致         N4  verdict_basis 字段已删除
  N5  cell-subset 失败中止            N6  compute_qc_metrics(force=)
  L1  供者级 TVD                      L2  strat_ci_frac_ok 入 config
  L3  CP 上界单侧                     L4  fdr_exclude_from_m
  L6  exact_enum_max_units            L8  --reference-alt
"""
import inspect
import os

import pytest

from sc_audit.config import AuditConfig, load_config


def _fields():
    return set(AuditConfig.__dataclass_fields__)


# ════════════ 字段级（N3 / N4 / L2 / L4 / L6）════════════

def test_n4_verdict_basis_field_removed():
    """N4：verdict_basis 字段必须已删除（注释中提及历史名称不算）。"""
    assert "verdict_basis" not in _fields()
    assert "ratio_basis" in _fields()


def test_l6_field_renamed_to_exact_enum_max_units():
    """L6：字段名与语义一致——真正的门是组合数。"""
    f = _fields()
    assert "exact_enum_max_units" in f
    assert "exact_max_combinations" in f
    assert "exact_enum_max_donors" not in f, "旧字段名仍在"


def test_l4_fdr_exclusion_field_exists():
    assert "fdr_exclude_from_m" in _fields()
    assert set(AuditConfig().fdr_exclude_from_m) == {"ACTB", "GAPDH"}


def test_l2_ci_fraction_is_configurable():
    f = _fields()
    assert "strat_ci_frac_ok" in f
    assert AuditConfig().strat_ci_frac_ok == 0.6


def test_n3_yaml_matches_defaults_field_by_field():
    """N3：cardiac.yaml 不得静默回退任何锁定阈值。"""
    y, d = load_config("cardiac").to_dict(), AuditConfig().to_dict()
    locked = ("min_cells_per_stratum", "floor_pct", "event_sparse", "event_limited",
              "ratio_min", "ratio_ok", "pct_ok", "ratio_basis", "hetero_max",
              "hetero_col", "match_var_tvd_low", "attenuation_min", "attenuation_min_full_effect_pp",
              "depth_delta_max_pct", "strat_min_delta_pp", "strat_ci_frac_ok",
              "strat_ci_level", "housekeeping_ratio_ok", "mito_delta_warn_pp",
              "fdr_method", "fdr_exclude_from_m", "matched_min_cells_warn",
              "forbid_outcome_matching", "match_var_tvd_low",
              "run_matched_comparison",
              "exact_enum_max_units", "exact_max_combinations",
              "auto_filter_stratum_cols", "match_vars", "strata_cols",
              "pos_markers", "neg_markers", "min_pos_markers", "min_counts",
              "min_genes", "depth_edges")
    mism = {k: (y.get(k), d.get(k)) for k in locked if y.get(k) != d.get(k)}
    assert not mism, f"yaml 与默认值不一致: {mism}"


# ════════════ 源码级（N1 / N5 / N6 / L1 / L3 / L8 / 报告）════════════

def test_n1_cli_actually_replaces_var_names():
    """N1 + P1-a：替换逻辑必须存在、被 CLI 调用、且在 compute_qc_metrics **之前**。"""
    import sc_audit.cli as cli
    import sc_audit.io as io
    # N1：替换逻辑现在在 io.apply_gene_name_fallback（CLI 与 Python API 共用同一条路径）
    src_io = inspect.getsource(io.apply_gene_name_fallback)
    assert "ad.var_names = " in src_io, "N1：未真正替换 var_names"
    assert "var_names_make_unique()" in src_io
    # P1-a：CLI 必须调用它，且调用点在 compute_qc_metrics 之前
    src_cli = inspect.getsource(cli)
    assert "apply_gene_name_fallback(ad" in src_cli, "N1：CLI 未调用回退函数"
    i_genes = src_cli.index("apply_gene_name_fallback(ad")
    i_qc = src_cli.index("compute_qc_metrics(ad")
    assert i_genes < i_qc, (
        "P1-a 回归：基因名回退必须在 compute_qc_metrics 之前，"
        "否则 MT- 前缀匹配不到 → percent_mito 被算成/覆盖为 NaN")


def test_n5_cli_aborts_on_cell_subset_failure():
    import sc_audit.cli as cli
    src = inspect.getsource(cli)
    assert "SystemExit(2)" in src, "N5：cell-subset 失败未中止"
    assert "已中止" in src


def test_n6_qc_metrics_has_force():
    from sc_audit.io import compute_qc_metrics
    sig = inspect.signature(compute_qc_metrics)
    assert "force" in sig.parameters


def test_l1_donor_level_tvd_exists():
    import sc_audit.composition as comp
    src = inspect.getsource(comp)
    assert "donor=donor" in src, "L1：TVD 未按供者加权"
    assert "tvd_donor" in src
    sig = inspect.signature(comp._group_distribution_distance)
    assert "donor" in sig.parameters


def test_l3_cp_upper_is_one_sided():
    import sc_audit.floor as fl
    src = inspect.getsource(fl._cp_upper)
    assert "beta.ppf(level" in src, "L3：CP 上界未统一为单侧"
    assert "(1 - level) / 2" not in src


def test_l8_reference_alt_option():
    import sc_audit.cli as cli
    src = inspect.getsource(cli)
    assert "reference-alt" in src
    assert "verdict_alt" in src


def test_report_whitelist_includes_q_value():
    import sc_audit.report as rep
    src = inspect.getsource(rep)
    for col in ("q_value_BH", "donor_weighted_effect", "basis_pct", "perm_mode"):
        assert col in src, f"报告白名单缺少 {col}"


# ════════════ 文档 ↔ 代码 同步（防"声称 ≠ 证据"）════════════

DOC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "代码审计文档_sc-audit.md")


@pytest.mark.skipif(not os.path.exists(DOC), reason="审计文档未生成")
def test_doc_code_sync():
    """文档内嵌代码必须与仓库逐行一致（由开发仓库的维护脚本核验；该脚本不随包发布）。"""
    import subprocess
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script = os.path.join(repo, "tools", "check_doc_code_sync.py")
    if not os.path.exists(script):
        pytest.skip("自检脚本不存在")
    r = subprocess.run(["python3", script, DOC], cwd=repo,
                       capture_output=True, text=True)
    assert r.returncode == 0, f"DOC != CODE:\n{r.stdout}\n{r.stderr}"
    assert "DOC==CODE" in r.stdout


@pytest.mark.skipif(not os.path.exists(DOC), reason="审计文档未生成")
def test_doc_has_no_conflicting_test_counts():
    """文档声称的测试项数必须与实际一致（24/44/49 并存的教训）。

    注意正则要容忍 markdown 强调标记（**测试规模：** 49 项），
    否则会像第一次那样"匹配不到 → 空集合 → 假通过"。
    """
    import re
    from subprocess import run
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    doc = open(DOC, encoding="utf-8").read()

    # 只检查文档头部元信息块（全文搜索会误匹配代码块中的示例文本）
    head = "\n".join(doc.split("\n")[:20])
    counts = {int(m.group(1))
              for m in re.finditer(r"测试规模[：:]\s*\**\s*(\d+)\s*项", head)}
    assert counts, "文档中未找到 '测试规模：N 项' 声明（正则可能失效）"

    r = run(["python3", "-m", "pytest", "tests/", "--collect-only", "-q"],
            cwd=repo, capture_output=True, text=True)
    m = re.search(r"(\d+)\s+tests?\s+collected", r.stdout)
    assert m, f"无法解析实际测试项数:\n{r.stdout[-500:]}"
    actual = int(m.group(1))

    assert len(counts) <= 1, f"文档中出现多个互相矛盾的测试项数: {sorted(counts)}"
    claimed = list(counts)[0]
    assert claimed == actual, (
        f"文档声称 {claimed} 项，实际 {actual} 项 —— 请重新运行 "
        f"tools/build_audit_doc.py 生成文档")

# ════════════════════════════════════════════════════════════════
#  修复表声称 → 可执行断言 的注册表
# ════════════════════════════════════════════════════════════════
#
# 为什么需要它（核心意见）：
#   文档-代码同步脚本只保证【文档内嵌代码 == 仓库】，
#   但**不检查修复表（§7/§8/§9）里的 ✅ 是否真的兑现**。
#   结果就是 §9.2 声称"新增 match_not_possible_on" 而代码里根本没有——
#   文档与代码"一致地错误"。
#
# 机制：
#   · 每条修复表项在此登记一条断言（@claim("<ID>", "<声称>")）。
#   · test_claim_registry_covers_all_rounds 强制"曾声称过的 ID 都在表里"。
#   · 新增修复项却忘了加断言 → 该测试失败。
#   · 声称改了但代码没改 → 对应的 claim 断言失败。
#
# 因此：**声称 ≠ 证据 的问题无法再无声通过**。

import re as _re
from typing import Callable, Dict, Tuple

import sc_audit
import sc_audit.cli as _cli
import sc_audit.composition as _comp
import sc_audit.config as _cfg
import sc_audit.floor as _floor
import sc_audit.gating as _gating
import sc_audit.io as _io
import sc_audit.reference as _ref
import sc_audit.report as _rep
import sc_audit.robust as _rob

CLAIMS: Dict[str, Tuple[str, Callable[[], None]]] = {}

CLAIM_PREFIXES = ("B", "S", "N", "L", "P")


def claim(cid: str, desc: str):
    """把一个修复 ID 绑定到一条可执行断言。"""
    def deco(fn):
        CLAIMS[cid] = (desc, fn)
        return fn
    return deco


def _src(mod) -> str:
    return inspect.getsource(mod)


def _fields() -> set:
    return set(AuditConfig.__dataclass_fields__)


# ──────────────────────── B / S 类声称 ────────────────────────

@claim("B1", "显式 --control-label 不被自动推断覆盖（字段 + CLI 传参）")
def _():
    from sc_audit.io import ensure_obs_fields
    src = inspect.getsource(ensure_obs_fields)
    assert "control_label is None" in src, "对照标签未先判 None"
    assert "control_label not in uniq" in src, "未校验用户指定的标签存在"
    # CLI 必须把 case_label 一起传下去（原 bug 是只传了 3 个参数）
    assert "args.control_label, args.case_label" in _src(_cli)


@claim("B2", "匹配变量与展示变量分离 + 样本/供者级唯一列自动剔除")
def _():
    f = _fields()
    assert "match_vars" in f and "strata_cols" in f
    src = _src(_comp)
    assert "_pick_match_vars" in src and "auto_filter_stratum_cols" in src


@claim("B3", "多水平 condition 给出警告并记录 case_scope")
def _():
    src = inspect.getsource(_io.ensure_obs_fields)
    assert "case_scope" in src and "warnings" in src


@claim("B4", "供者跨组必须报错（入口 + API 双层）")
def _():
    assert hasattr(_io, "assert_donors_nested")
    for mod in (_floor, _comp, _rob):
        assert "assert_donors_nested" in _src(mod), f"{mod.__name__} 缺入口校验"


@claim("B5", "X 非原始计数时回退到 counts layer / raw，否则报错")
def _():
    assert hasattr(_io, "resolve_counts")
    assert "require_integer_counts" in _fields()


@claim("B6", "cardiac.yaml 随包分发（不在包外 ../configs/）")
def _():
    cfg_dir = os.path.join(os.path.dirname(sc_audit.__file__), "configs")
    assert os.path.exists(os.path.join(cfg_dir, "cardiac.yaml"))
    assert load_config("cardiac").min_cells_per_stratum == AuditConfig().min_cells_per_stratum
    assert "configs/*.yaml" in open(
        os.path.join(os.path.dirname(sc_audit.__file__), "..", "pyproject.toml"),
        encoding="utf-8").read()


@claim("S1", "F_hetero 已改为分布距离（TVD + Cramér's V），并在平衡设计下不误报")
def _():
    assert hasattr(_comp, "_group_distribution_distance")
    src = inspect.getsource(_comp._group_distribution_distance)
    assert "tvd" in src.lower() and "chi2" in src.lower()


@claim("S6", "verdict 与 R 共用同一分子口径（唯一字段 ratio_basis）")
def _():
    f = _fields()
    assert "ratio_basis" in f
    assert "verdict_basis" not in f, "N4 应已删除 verdict_basis"
    src = inspect.getsource(_floor.audit_detection_floor)
    assert "cfg.ratio_basis" in src


@claim("S7", "unquantifiable property 与 summary 用同一定义")
def _():
    src = inspect.getsource(_floor.DetectionFloorResult)
    assert "unquantifiable" in src


@claim("FDR", "BH-FDR 已实现，且阳性对照不参与假设集合 m")
def _():
    f = _fields()
    assert "fdr_method" in f and "fdr_exclude_from_m" in f
    src = inspect.getsource(_floor.audit_detection_floor)
    assert "q_value_BH" in src and "fdr_exclude_from_m" in src


@claim("CI", "Wilson CI + Clopper–Pearson 单侧上界")
def _():
    assert hasattr(_floor, "_cp_upper") and hasattr(_floor, "_wilson_ci")
    assert "(1 - level) / 2" not in inspect.getsource(_floor._cp_upper), "CP 应为单侧"


@claim("EXACT", "供者数/组合数可控时精确枚举（消除种子依赖）")
def _():
    assert "exact_enum_max_units" in _fields()
    assert "exact_max_combinations" in _fields()


@claim("JK-BOTH", "leave-one-donor 为双侧，且透传组合数阈值")
def _():
    sig = inspect.signature(_floor.leave_one_donor)
    assert "sides" in sig.parameters
    assert "exact_max_combinations" in sig.parameters


@claim("REF-HASH", "metadata 记录参照表路径与 SHA256（主表 + alt 表）")
def _():
    src = _src(_cli)
    assert "reference_sha256" in src and "reference_alt_sha256" in src
    assert "reference_alt_mode" in src


@claim("DEADPARAM", "reference.compute_detection_ratio 无死参数")
def _():
    sig = inspect.signature(_ref.compute_detection_ratio)
    assert "reference_mode" not in sig.parameters


# ──────────────────────── N / L 类声称 ────────────────────────

@claim("N1", "基因名回退真正替换 var_names（CLI + API 共用），且早于 QC")
def _():
    assert hasattr(_io, "apply_gene_name_fallback")
    assert "ad.var_names = " in inspect.getsource(_io.apply_gene_name_fallback)
    src = _src(_cli)
    assert "apply_gene_name_fallback(ad" in src
    assert src.index("apply_gene_name_fallback(ad") < src.index("compute_qc_metrics(ad"), \
        "P1-a：替换必须在 QC 之前"


@claim("N2", "匹配列剔除判据按【供者级】唯一（非细胞级）")
def _():
    src = inspect.getsource(_comp._pick_match_vars)
    assert "n_donors" in src or "donor_unique" in src


@claim("N3", "cardiac.yaml 与 AuditConfig 默认值逐字段一致")
def _():
    y, d = load_config("cardiac").to_dict(), AuditConfig().to_dict()
    locked = ("min_cells_per_stratum", "floor_pct", "event_sparse", "hetero_max",
              "attenuation_min", "depth_delta_max_pct", "strat_min_delta_pp",
              "strat_ci_frac_ok", "ratio_basis", "fdr_method", "fdr_exclude_from_m")
    mism = {k: (y.get(k), d.get(k)) for k in locked if y.get(k) != d.get(k)}
    assert not mism, f"yaml 与默认值不一致: {mism}"


@claim("N4", "verdict_basis 字段已删除，R 与判定共用 ratio_basis")
def _():
    assert "verdict_basis" not in _fields()
    assert "ratio_basis" in _fields()


@claim("N5", "--cell-subset 表达式失败必须中止（不得继续在全量数据上出结论）")
def _():
    assert "SystemExit(2)" in _src(_cli)


@claim("N6", "换源后强制重算 QC（force 参数存在）")
def _():
    assert "force" in inspect.signature(_io.compute_qc_metrics).parameters


@claim("L1", "TVD 判定用供者级加权（细胞级作参考）")
def _():
    sig = inspect.signature(_comp._group_distribution_distance)
    assert "donor" in sig.parameters
    assert "tvd_donor" in _src(_comp)


@claim("L2", "CI 重叠比例阈值入 config；verdict 不再声称 'all strata'")
def _():
    assert "strat_ci_frac_ok" in _fields()
    src = inspect.getsource(_rob.audit_robust_signal)
    assert "cfg.strat_ci_frac_ok" in src
    assert "in all depth strata" not in src


@claim("L3", "Clopper–Pearson 上界统一为单侧")
def _():
    assert "beta.ppf(level" in inspect.getsource(_floor._cp_upper)


@claim("L4", "FDR 排除阳性对照（fdr_exclude_from_m）")
def _():
    assert set(AuditConfig().fdr_exclude_from_m) == {"ACTB", "GAPDH"}


@claim("L6", "字段改名 exact_enum_max_units（旧名已删）")
def _():
    f = _fields()
    assert "exact_enum_max_units" in f and "exact_enum_max_donors" not in f


@claim("L8", "--reference-alt + --reference-alt-mode")
def _():
    src = _src(_cli)
    assert "reference-alt" in src and "verdict_alt" in src
    assert "reference_alt_mode" in src


@claim("REPORT-WL", "报告白名单含 q 值与供者加权效应量等新列")
def _():
    src = _src(_rep)
    for c in ("q_value_BH", "donor_weighted_effect", "basis_pct", "perm_mode",
              "case_pct_upper95", "overall_pct_upper95"):
        assert c in src, f"报告白名单缺 {c}"


# ──────────────────── P1-a~d / P2 / P3 类声称 ────────────────────

@claim("P1a", "基因名回退早于 compute_qc_metrics（防 percent_mito=NaN）")
def _():
    src = _src(_cli)
    assert src.index("apply_gene_name_fallback(ad") < src.index("compute_qc_metrics(ad")


@claim("P1b", "force 仅在真正换源且 config 允许时启用；字段被引用")
def _():
    src = _src(_cli)
    assert "switched = counts_source not in" in src
    assert "switched and cfg.recompute_qc_after_source_switch" in src
    assert "recompute_qc_after_source_switch" in _fields()
    assert "qc_recomputed" in src and "n_vars" in src


@claim("P1c", "管家基因 R 与目标基因同口径（输出 basis 列）")
def _():
    src = inspect.getsource(_rob.audit_robust_signal)
    assert "cfg.ratio_basis" in src, "管家基因仍用固定口径"
    assert '"basis" in' in src or "'basis'" in src or '"basis":' in src


@claim("P1d", "两列 CP 上界命名与语义一致")
def _():
    src = inspect.getsource(_floor.audit_detection_floor)
    assert "case_pct_upper95" in src and "overall_pct_upper95" in src
    assert "case0_upper95" not in src and "pct_case_upper95" not in src


@claim("P2a", "match_not_possible_on 字段存在 + 三分支 note")
def _():
    f = set(_comp.CompositionResult.__dataclass_fields__)
    assert "match_not_possible_on" in f, "§9.2 声称的字段不存在（P1 复发点）"
    assert "match_attempted" in f
    src = inspect.getsource(_comp.audit_control_composition)
    assert "NO matching was ATTEMPTED" in src, "缺少'未尝试匹配'分支"
    assert "matching WAS attempted" in src, "缺少'尝试过但无人匹配'分支"


@claim("P2b", "alt 参照表的聚合方式可指定，且 metadata 记录")
def _():
    src = _src(_cli)
    assert "reference-alt-mode" in src
    assert "reference_alt_mode" in src and "reference_alt_n_genes" in src


@claim("P2c", "jackknife 透传 exact_max_combinations")
def _():
    assert "exact_max_combinations" in inspect.signature(_floor.leave_one_donor).parameters


@claim("P2d", "TVD 分辨率 note 用每组供者数的较小者")
def _():
    src = inspect.getsource(_comp.audit_control_composition)
    assert "_n_min" in src and "TVD resolution note" in src


@claim("P3-HK", "管家基因不计入 genes_persisting")
def _():
    src = inspect.getsource(_rob.audit_robust_signal)
    assert "hk_upper" in src and "excluded from the persistence count" in src


@claim("P2-2", "--housekeeping 同步到 FDR 排除集")
def _():
    src = _src(_cli)
    assert "cfg.fdr_exclude_from_m = " in src


@claim("P2-3", "两个入口函数已导出（Python API 可用）")
def _():
    for name in ("apply_gene_name_fallback", "resolve_counts", "assert_donors_nested"):
        assert name in sc_audit.__all__, f"{name} 未导出"
        assert hasattr(sc_audit, name)


@claim("P1-2", "case_label 哨兵 '<all non-control>' 不再导致崩溃")
def _():
    src = inspect.getsource(_io.ensure_obs_fields)
    assert '"<all non-control>"' in src
    # 哨兵必须走"重新推断"分支，而不是拿去和真实取值比较
    assert "case_label == \"<all non-control>\"" in src


# ──────────────────────── 注册表自身的一致性 ────────────────────────

ALL_CLAIM_IDS = ("B1", "B2", "B3", "B4", "B5", "B6",
                 "S1", "S6", "S7", "FDR", "CI", "EXACT", "JK-BOTH",
                 "REF-HASH", "DEADPARAM",
                 "N1", "N2", "N3", "N4", "N5", "N6",
                 "L1", "L2", "L3", "L4", "L6", "L8", "REPORT-WL",
                 "P1a", "P1b", "P1c", "P1d",
                 "P2a", "P2b", "P2c", "P2d", "P3-HK",
                 "P2-2", "P2-3", "P1-2")



ALL_CLAIM_IDS = ALL_CLAIM_IDS + ("COMPARE", "SUBCOMMANDS", "HK-DEFAULT", "PROBE",
                                 "ALT-ONLY", "L5", "P2a-PARTIAL",
                                 "CHECK", "CHECK-ADVISORY",
                                 "CHECK-ROLE", "CHECK-RECOMMEND", "CHECK-REFCOV",
                                 "CHECK-OVERRIDE", "CHECK-OPTIONS", "CHECK-HINTS-REAL", "CHECK-REF-PROV")


def test_claim_registry_covers_all_rounds():
    """曾声称修过的每一条都必须在注册表里有断言（防止 §9.2 那类漏做）。"""
    missing = [c for c in ALL_CLAIM_IDS if c not in CLAIMS]
    assert not missing, f"以下声称没有对应断言，请补 @claim: {missing}"


# 记录实际执行过的声称（供 test_every_registered_claim_was_executed 兜底）
_CLAIMS_EXECUTED: set = set()


def pytest_generate_tests(metafunc):
    """在【收集期】展开声称注册表。

    ⚠️ 此处**不能**用 @pytest.mark.parametrize("cid", sorted(CLAIMS))。
    那样会在模块执行到该行时固化 CLAIMS 的快照 —— 而本文件中
    CHECK-* 系列声称定义在该行【之后】，于是它们虽被注册进 CLAIMS
    （所以 test_claim_registry_covers_all_rounds 通过），
    却永远不会被 parametrize 展开执行 —— 典型的"绿色但没跑"。
    pytest_generate_tests 在整段模块导入完毕后调用，故不会遗漏。
    """
    if "cid" in metafunc.fixturenames:
        metafunc.parametrize("cid", sorted(CLAIMS))


def test_claim_is_actually_delivered(cid):
    """逐条执行注册表中的断言：声称必须被兑现。"""
    _CLAIMS_EXECUTED.add(cid)
    desc, fn = CLAIMS[cid]
    try:
        fn()
    except AssertionError as e:
        pytest.fail(f"[{cid}] 声称『{desc}』未兑现: {e}")
    except Exception as e:  # noqa: BLE001
        pytest.fail(f"[{cid}] 声称『{desc}』检查时异常: {type(e).__name__}: {e}")


# ──────────────────── compare / P2-a 部分剔除 / L5 / P3 ────────────────────

@claim("COMPARE", "compare 子命令可运行（不得向已删参数的函数传关键字）")
def _():
    src = inspect.getsource(_cli)
    assert "reference_mode=" not in src, \
        "cli 仍传 reference_mode= → compare 会 TypeError（P1）"
    assert "compute_detection_ratio(study, ref)" in src


@claim("SUBCOMMANDS", "五个子命令各有冒烟测试文件")
def _():
    d = os.path.dirname(os.path.abspath(__file__))
    assert os.path.exists(os.path.join(d, "test_subcommands.py"))
    src = open(os.path.join(d, "test_subcommands.py"), encoding="utf-8").read()
    for c in ("run", "floor", "composition", "robust", "compare"):
        assert f'("{c}",' in src, f"缺少 {c} 子命令的冒烟测试"


@claim("HK-DEFAULT", "管家基因默认值取 cfg.housekeeping_genes（不再硬编码）")
def _():
    src = inspect.getsource(_cli)
    assert 'or "ACTB,GAPDH"' not in src, "仍硬编码 ACTB,GAPDH"
    assert "cfg, \"housekeeping_genes\"" in src or 'cfg, "housekeeping_genes"' in src


@claim("PROBE", "基因名回退探针 = 目标 ∪ 管家 ∪ 门控 marker")
def _():
    """按【行为】断言而非字面字符串（教训）。

    旧断言写死了 `_pp.get("pos_markers")`，而代码后来加了默认值
    `_pp.get("pos_markers", [])` —— 行为正确但断言失配。
    由于该声称当时从未被执行（parametrize 快照过早），这个失配被隐藏了很久。
    现改为：抽取"构造探针"的代码块，断言三个来源都在其中。
    """
    src = inspect.getsource(_cli)
    lines = src.splitlines()
    i = next(k for k, l in enumerate(lines) if "probe = list(" in l)
    j = next(k for k, l in enumerate(lines)
             if "apply_gene_name_fallback(ad, probe)" in l)
    block = "\n".join(lines[i:j])
    assert "pos_markers" in block, "探针未包含门控阳性 marker"
    assert "neg_markers" in block, "探针未包含门控阴性 marker"
    assert "housekeeping" in block, "探针未包含管家基因"
    assert "dict.fromkeys" in block, "探针未去重"


@claim("ALT-ONLY", "--reference-alt 无 --reference 时跳过翻转统计")
def _():
    src = inspect.getsource(_cli)
    assert "if ref_alt and not ref:" in src


@claim("L5", "jackknife 增加 BH 临界值口径列（不各自做 BH）")
def _():
    src = inspect.getsource(_floor.audit_detection_floor)
    assert "jackknife_n_below_bh_crit" in src
    # 修正（外部审计 #2）：α 已提入 config；且 m/rank 必须与 q 值**同族** ——
    # 由 p_for_fdr（已排除 fdr_exclude_from_m 的管家基因）构建，
    # 而非未排除的 df["perm_P"]。
    assert "alpha * rank / m_tests" in src, "临界值公式应为 alpha·rank/m"
    assert 'getattr(cfg, "fdr_alpha", 0.05)' in src, "alpha 应取自 config.fdr_alpha"
    assert "enumerate(p_for_fdr)" in src, \
        "order 应由 p_for_fdr 构建（与 q 值同族，排除管家基因）"
    assert 'enumerate(df["perm_P"]' not in src, \
        "order 不得由未排除管家基因的 df[perm_P] 构建"
    assert "jackknife_n_sig_note" in src, "未校正列应有语义注释"


@claim("P2a-PARTIAL", "P2-a 部分剔除分支也记录 match_not_possible_on")
def _():
    src = inspect.getsource(_comp.audit_control_composition)
    assert "design_dropped_b" in src
    assert "WEAKER than a same-design" in src
    assert "TRIVIALLY satisfied" in src


# ──────────────────── 预检子命令（P0）────────────────────

@claim("CHECK", "check 子命令存在且被 CLI 注册")
def _():
    import sc_audit.check as _chk
    assert hasattr(_chk, "check_readiness") and hasattr(_chk, "ReadinessResult")
    from sc_audit.cli import _build_parser
    p = _build_parser()
    subs = [a for a in p._subparsers._group_actions[0].choices]
    assert "check" in subs, f"check 未注册为子命令: {subs}"
    from sc_audit.cli import _run_check
    assert "check_report.json" in inspect.getsource(_run_check)


@claim("CHECK-ADVISORY", "check 只报告、不采纳：绝不写入 cfg.match_vars")
def _():
    import sc_audit.check as _chk
    src = inspect.getsource(_chk.check_readiness)
    # 不得对 cfg.match_vars / match_vars 赋值
    for bad in ("cfg.match_vars =", "cfg.match_vars=", "match_vars = ["):
        assert bad not in src, f"check 修改了匹配变量（违反只报告契约）: {bad}"
    # 不得调用会写配置的函数
    for bad in ("save_config", "yaml.dump", "write_text("):
        assert bad not in src, f"check 写入了文件/配置: {bad}"
    # 报告里必须明确声明"只报告"
    assert "only REPORTS" in inspect.getsource(_chk.ReadinessResult.to_json)


@claim("CHECK-ROLE", "变量角色分类：词根定角色，统计信号不覆盖（防假阳性）")
def _():
    import sc_audit.check as _chk
    src = inspect.getsource(_chk.classify_variable)
    assert "_name_role(col" in src, "角色未由词名决定"
    # 统计信号不得改写 role
    for bad in ('role = "outcome_suspect"', "role = 'outcome_suspect'"):
        assert bad not in src, "统计信号仍在覆盖角色（会把取材变量误判为结局）"
    assert "cond_dep" in src
    # 疑似区必须仍被输出（不隐藏）
    assert "suspect_design_variables" in inspect.getsource(_chk.ReadinessResult.to_json)


@claim("CHECK-RECOMMEND", "给出排序建议与 audit_2 总体建议（含方法学措辞）")
def _():
    import sc_audit.check as _chk
    src = inspect.getsource(_chk.recommend_role)
    for rec in ("primary_match", "risky", "safe_but_uninformative",
                "caution_condition_dependent", "not_usable", "never", "avoid"):
        assert rec in src, f"缺少建议档位 {rec}"
    src2 = inspect.getsource(_chk.check_readiness)
    assert "statement_for_paper" in src2
    for m in ("matched_comparison", "no_matching_needed", "tvd_stratification",
              "tvd_only", "matched_comparison_marginal"):
        assert m in src2, f"缺少 audit_2 建议方法 {m}"


@claim("CHECK-REFCOV", "参照覆盖以目标面板为分母（不是数据集全部基因）")
def _():
    import sc_audit.check as _chk
    src = inspect.getsource(_chk.check_readiness)
    assert "pct_target_covered" in src
    assert "target_missing" in src


@claim("CHECK-OVERRIDE", "变量角色可显式声明（priority > 词根），并有运行期保护")
def _():
    f = set(AuditConfig.__dataclass_fields__)
    for k in ("role_overrides", "extra_outcome_keywords",
              "extra_treatment_keywords", "forbid_outcome_matching"):
        assert k in f, f"缺少字段 {k}"
    import sc_audit.check as _chk
    src = inspect.getsource(_chk._name_role)
    assert "role_overrides" in src, "显式声明未接入角色判定"
    # CLI 必须有运行期保护 + 显式放行开关
    src_cli = inspect.getsource(_cli)
    assert "_outcome_like_match_vars" in src_cli
    assert "allow-outcome-matching" in src_cli
    assert "collider bias" in src_cli


@claim("CHECK-OPTIONS", "接近分离变量存在时给出三套可选策略（不排序、不自动采纳）")
def _():
    import sc_audit.check as _chk
    src = inspect.getsource(_chk.check_readiness)
    for oid in ("conservative", "primary_plus_sensitivity", "no_matching"):
        assert oid in src, f"缺少策略 {oid}"
    assert "options_note" in src
    assert "不替你排序" in src


@claim("CHECK-HINTS-REAL", "策略块建议的配置键必须真实存在（防'声称≠实现'复发）")
def _():
    import sc_audit.check as _chk
    src = inspect.getsource(_chk.check_readiness)
    for k in ("run_matched_comparison", "strata_cols", "match_vars"):
        assert k in src, f"策略块缺少真实键 {k}"
    assert "stratify_by" not in src, "策略块仍在建议不存在的键 stratify_by"
    f = set(AuditConfig.__dataclass_fields__)
    for k in ("run_matched_comparison", "strata_cols", "match_vars"):
        assert k in f, f"配置缺少字段 {k}"
    # composition 必须真的兑现 run_matched_comparison
    import sc_audit.composition as _comp
    src_comp = inspect.getsource(_comp.audit_control_composition)
    assert "run_matched_comparison" in src_comp, "composition 未兑现该开关"


@claim("CHECK-REF-PROV", "参照表门控参数必须可追溯；未登记则告警（R 可比性前提）")
def _():
    import sc_audit.check as _chk
    src = inspect.getsource(_chk.check_readiness)
    assert "reference_manifest" in src
    assert "gate_unknown_backgrounds" in src
    assert "可比性无法核验" in src
    # CLI 必须真的去读 manifest
    src_cli = inspect.getsource(_cli)
    assert "manifest.json" in src_cli


def test_every_registered_claim_was_executed(request):
    """兜底：注册的声称必须【真的跑过】。

    与 test_claim_registry_covers_all_rounds 的分工：
      · 那个检查「声称有没有断言」——注册表完整性
      · 本测试检查「断言有没有被执行」——parametrize 展开完整性

    曾有一个缺陷恰好通过了前者、逃过了后者（8 条 CHECK-* 声称
    注册齐全却从未运行），所以两道都要有。
    """
    # 过滤运行（-k / -m）下只跑了子集，此时不应判定"未执行"
    if request.config.option.keyword or request.config.option.markexpr:
        pytest.skip("过滤运行下不检查展开完整性")
    missing = sorted(set(CLAIMS) - _CLAIMS_EXECUTED)
    assert not missing, (
        f"以下声称注册了但从未被执行：{missing} —— "
        f"检查 parametrize 是否在声称定义之前就固化了快照")


# ══════════════ 审计四（v0.2 新增）：细胞身份与捕获完整性 ══════════════

@claim("A4-1", "审计四：缺注释列 → INDETERMINATE，绝不静默返回 PASS")
def _():
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import anndata as ad
    from sc_audit.config import AuditConfig
    from sc_audit.identity import audit_cell_identity
    g = ["TTN", "MYH7", "MYBPC3", "ACTN2", "TNNI3", "MYL2", "PTPRC", "CD3E",
         "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"]
    o = pd.DataFrame({"condition": ["control", "disease"] * 10,
                      "donor": [f"d{i}" for i in range(20)]})
    A = ad.AnnData(X=sp.csr_matrix(np.ones((20, len(g)), dtype=int)), obs=o,
                   var=pd.DataFrame(index=g))
    cfg = AuditConfig()
    cfg.identity_target = ""            # 无注释列、无目标
    r = audit_cell_identity(A, "condition", "donor", "control", cfg=cfg)
    assert r.verdict.startswith("INDETERMINATE"), r.verdict
    assert not r.verdict.startswith("PASS")


@claim("A4-2", "审计四：亚群不平衡（来源 × 条件，在细胞类型内部）必须检出")
def _():
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import anndata as ad
    from sc_audit.config import AuditConfig
    from sc_audit.identity import audit_cell_identity
    POS = ["TTN", "MYH7", "MYBPC3", "ACTN2", "TNNI3", "MYL2"]
    NEG = ["PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"]
    g = POS + NEG
    rng = np.random.default_rng(3)
    rows = []
    for arm in ("control", "disease"):
        for d in range(3):
            for i in range(40):
                cm = rng.random() < 0.6
                src = "CD45+" if (cm and arm == "control" and rng.random() < 0.7) else "Cells"
                rows.append({"condition": arm, "donor": f"{arm[0]}{d}", "source": src,
                             "lineage_names": "Cardiomyocytes" if cm else "Endothelial cells",
                             "n_genes": 900, "_cm": cm})
    X = np.zeros((len(rows), len(g)), dtype=int)
    for i, r in enumerate(rows):
        for j, gg in enumerate(g):
            p = (0.85 if gg in POS else 0.02) if r["_cm"] else 0.05
            X[i, j] = (rng.poisson(3) + 1) if rng.random() < p else 0
    A = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(rows), var=pd.DataFrame(index=g))
    cfg = AuditConfig()
    cfg.identity_col, cfg.identity_source_col = "lineage_names", "source"
    cfg.identity_target, cfg.identity_complexity_col = "Cardiomyocytes", "n_genes"
    r = audit_cell_identity(A, "condition", "donor", "control", cfg=cfg)
    assert r.imbalanced_subsets, "亚群不平衡未检出（这是审计四的核心项）"
    assert "COMPOSITION-CONFOUNDED" in r.verdict


@claim("A4-3", "审计四：门控规则的条件共线性必须检出（阴性标志随臂变化）")
def _():
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import anndata as ad
    from sc_audit.config import AuditConfig
    from sc_audit.identity import audit_cell_identity
    POS = ["TTN", "MYH7", "MYBPC3", "ACTN2", "TNNI3", "MYL2"]
    NEG = ["PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"]
    g = POS + NEG
    rng = np.random.default_rng(5)
    rows = []
    for arm in ("control", "disease"):
        for d in range(3):
            for i in range(40):
                rows.append({"condition": arm, "donor": f"{arm[0]}{d}", "source": "Cells",
                             "lineage_names": "Cardiomyocytes", "n_genes": 900})
    X = np.zeros((len(rows), len(g)), dtype=int)
    for i, r in enumerate(rows):
        for j, gg in enumerate(g):
            if gg in POS:
                p = 0.85
            elif gg in ("PECAM1", "VWF"):
                p = 0.45 if r["condition"] == "disease" else 0.02
            else:
                p = 0.02
            X[i, j] = (rng.poisson(3) + 1) if rng.random() < p else 0
    A = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(rows), var=pd.DataFrame(index=g))
    cfg = AuditConfig()
    cfg.identity_col, cfg.identity_source_col = "lineage_names", "source"
    cfg.identity_target, cfg.identity_complexity_col = "Cardiomyocytes", "n_genes"
    r = audit_cell_identity(A, "condition", "donor", "control", cfg=cfg)
    assert set(r.collinear_markers) & {"PECAM1", "VWF"}, r.collinear_markers
    assert r.gate_pass_by_arm is not None


@claim("A4-4", "审计四：无同模态参照时 R 必须标 cross-modal")
def _():
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import anndata as ad
    from sc_audit.config import AuditConfig
    from sc_audit.identity import audit_cell_identity
    g = ["TTN", "MYH7", "MYBPC3", "ACTN2", "TNNI3", "MYL2",
         "PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"]
    rows = [{"condition": ("control" if i < 20 else "disease"), "donor": f"d{i%6}",
             "source": "Cells", "lineage_names": "Cardiomyocytes", "n_genes": 900}
            for i in range(40)]
    A = ad.AnnData(X=sp.csr_matrix(np.ones((40, len(g)), dtype=int)),
                   obs=pd.DataFrame(rows), var=pd.DataFrame(index=g))
    cfg = AuditConfig()
    cfg.identity_col, cfg.identity_source_col = "lineage_names", "source"
    cfg.identity_target, cfg.identity_complexity_col = "Cardiomyocytes", "n_genes"
    aa = audit_cell_identity(A, "condition", "donor", "control", cfg=cfg,
                             reference_modalities=["snRNA"], modality="scRNA")
    assert aa.modality_matched is False
    assert any("cross-modal" in n for n in aa.notes)
    bb = audit_cell_identity(A, "condition", "donor", "control", cfg=cfg,
                             reference_modalities=["scRNA"], modality="scRNA")
    assert bb.modality_matched is True


@claim("A4-5", "审计四 CLI：identity 子命令存在，且 --celltype-col 等参数真正写入 cfg")
def _():
    src = _src(_cli)
    assert '"identity"' in src, "identity 未注册为子命令"
    for a in ("--celltype-col", "--source-col", "--target-lineage",
              "--modality", "--reference-modalities"):
        assert a in src, f"CLI 缺参数 {a}"
    # 参数必须真的写进 cfg（不能只解析不使用）
    assert "cfg.identity_col = args.celltype_col" in src
    assert "cfg.identity_source_col = args.source_col" in src
    assert "cfg.identity_target = args.target_lineage" in src
    assert hasattr(sc_audit, "audit_cell_identity")
    assert hasattr(sc_audit, "IdentityResult")


# ════════════════ 审计二 v0.3 可识别性层（2026-10-06 并入）════════════

@claim("ID-EXPORT", "identifiability API 已导出（包级公共 API）")
def _():
    import sc_audit
    for n in ("audit_identifiability", "IdentifiabilityResult",
              "build_variable_roster", "common_support", "build_matching_ladder",
              "default_ladder"):
        assert hasattr(sc_audit, n), f"包未导出 {n}"
        assert n in sc_audit.__all__, f"__all__ 缺 {n}"
    assert tuple(int(x) for x in sc_audit.__version__.split(".")[:2]) >= (0, 3), \
        f"版本不得低于 0.3（可识别性层引入版本）：{sc_audit.__version__}"


@claim("ID-FIELDS", "v0.3 阈值字段全部存在且为事前锁定默认值")
def _():
    f = _fields()
    for n in ("identifiability_ladder", "identifiability_min_support_cells",
              "identifiability_one_sided_verdict",
              "identifiability_suppress_below_floor", "identifiability_nondesign_cols",
              "identifiability_derived_cols", "identifiability_derived_prefixes"):
        assert n in f, f"缺字段 {n}"
    assert "identifiability_require_both_strata" not in f, \
        "旧字段仍在（审计 A1-F1 后应已被 identifiability_one_sided_verdict 取代）"
    d = AuditConfig()
    assert d.identifiability_ladder is None, "默认必须是内置阶梯（None）"
    assert d.identifiability_min_support_cells is None, \
        "默认必须复用 matched_min_cells_warn（单一真值源）"
    assert d.identifiability_one_sided_verdict is True, "单侧判定默认开启"
    assert d.matched_min_cells_warn == 50, "支持域下限保持 50（2026-10-06 确定）"


@claim("ID-PROTECT", "单侧判定：给方向 + 强制警示；两档并存报 THRESHOLD-SPANNING")
def _():
    import sc_audit.identifiability as _id
    src = _src(_id.audit_identifiability)
    assert "identifiability_one_sided_verdict" in src, "未实现单侧开关"
    assert "ONE-SIDED CAVEAT" in src, "缺少强制警示"
    assert "consistency check carries NO information" in src, "警示未说明一致性检查无信息"
    for lbl in ("VERDICT_SUPPORT_BELOW_FLOOR", "VERDICT_THRESHOLD_SPANNING",
                "VERDICT_COMPARABLE", "VERDICT_NOT_SEPARABLE",
                "VERDICT_INDETERMINATE"):
        assert hasattr(_id, lbl), f"缺少标签常量 {lbl}"
    assert "VERDICT_THRESHOLD_SPANNING" in src, "跨阈值标签未被使用"
    # 两档并存 ⇒ 必判 THRESHOLD-SPANNING（不得再退化成"永不可达"的死分支）
    assert "elif spanning:" in src


@claim("ID-COMPLETE", "花名册完备性：未归类列必须浮出，不得静默跳过")
def _():
    import sc_audit.identifiability as _id
    src = _src(_id.build_variable_roster)
    assert "unclassified" in src
    assert "except Exception" in src, "列不可归类时必须捕获并记录（而不是静默跳过）"
    src2 = _src(_id.audit_identifiability)
    assert "silently skipped" in src2, "报告须明写未静默跳过"


@claim("ID-SUPPRESS", "低于支持域下限 → 抑制衰减比（不出量化证据）")
def _():
    import sc_audit.identifiability as _id
    src = _src(_id.audit_identifiability)
    assert "attenuated_suppressed" in src
    assert "NOT reported as evidence" in src


@claim("ID-CLI", "identifiability 子命令已注册，产物写入 outdir，且阶梯可事前锁定")
def _():
    from sc_audit.cli import _build_parser
    p = _build_parser()
    subs = [a for a in p._subparsers._group_actions[0].choices]
    assert "identifiability" in subs, f"未注册子命令: {subs}"
    src = _src(_cli)
    assert "identifiability.txt" in src
    assert "identifiability_roster.csv" in src
    assert "identifiability_ladder_summary.csv" in src
    assert "--ladder" in src, "缺 --ladder（事前锁定阶梯入口）"
    assert "cfg.identifiability_ladder" in src, "阶梯必须写入 cfg 才能进 run_metadata"


@claim("ID-LADDER", "内置阶梯 R0–R5 的顺序与语义（锚点 → 现行口径 → 全量 → 设计变量 → 单变量 → 推荐窗）")
def _():
    from sc_audit.identifiability import default_ladder
    L = default_ladder(["institution", "region"], ["gender"], ["gender", "source"])
    labels = [k for k, _ in L]
    assert labels[0].startswith("R0") and L[0][1] == [], "R0 必须是锚点（不匹配）"
    assert labels[1].startswith("R1") and L[1][1] == ["institution", "donation_type",
                                                      "source", "region"], \
        "R1 必须是现行正文口径"
    assert labels[2].startswith("R2") and set(L[2][1]) == {"institution", "region", "gender"}
    assert labels[3].startswith("R3") and L[3][1] == ["institution", "region"]
    assert any(k == "R4 单变量 gender" for k, _ in L)
    assert labels[-1].startswith("R5") and L[-1][1] == ["gender", "source"]


@claim("ID-NUMBERS", "外部数据冒烟回归（h5ad 缺失时跳过）：不崩、判定合法、阶梯可算；预注册数字复现在 tools/regress_external_case.py")
def _():
    import os
    import numpy as np
    p = os.environ.get("SC_AUDIT_EXTERNAL_CASE_H5AD", "")
    if not (p and os.path.exists(p)):
        pytest.skip("未设置 SC_AUDIT_EXTERNAL_CASE_H5AD")
    import scanpy as sc
    from sc_audit.identifiability import audit_identifiability
    PANEL_14 = ["MYL3", "PLN", "TNNT2", "ACTB", "GAPDH"]
    ad = sc.read_h5ad(p)
    r = audit_identifiability(ad, PANEL_14, "condition", "donor", "control",
                              cfg=AuditConfig())
    assert r.support_cells >= 0 and r.n_obs_columns > 0
    meds = r.details["medians"]
    assert len(meds) >= 2, f"阶梯应至少两条可算规则，实为 {len(meds)}"
    assert r.verdict.split(" — ")[0] in {
        "COMPARABLE", "NOT SEPARABLE", "INDETERMINATE",
        "SUPPORT-BELOW-FLOOR", "THRESHOLD-SPANNING"}


@claim("ID-ROBUST", "健壮性：缺设计列 / 0 基因命中 / 空阶梯都不得崩溃（P0 回归）")
def _():
    import sc_audit.identifiability as _id
    src = _src(_id._match_donors)
    assert "avail" in src and "dropped" in src, "缺少缺列防护"
    src2 = _src(_id.build_matching_ladder)
    assert "LONG_COLS" in src2 and "SUM_COLS" in src2, \
        "退化情形必须返回**带列名**的空表（否则 Lad.rule 抛 AttributeError）"
    assert "dict.fromkeys" in src2, "目标基因未去重（P1-5）"


@claim("ID-DEDUP", "供者级工具已抽为单一实现（donor_utils），不再多处复制")
def _():
    import sc_audit.donor_utils as _du
    for n in ("donor_mode", "donor_level_values", "donor_group", "entropy",
              "is_donor_unique", "chance_constancy_risk"):
        assert hasattr(_du, n), f"donor_utils 缺 {n}"
    comp = _src(_comp._pick_match_vars)
    assert "is_donor_unique(" in comp, "composition 未改用共享实现（P2-12）"
    assert "per_donor <= 1" not in comp, "composition 仍保留重复的判据副本"


@claim("ID-VERDICT-COLOUR", "v0.3 负面判定不得渲染为绿色（P1-3）")
def _():
    from sc_audit.report import _verdict_class
    assert _verdict_class("SUPPORT-BELOW-FLOOR — x") == "bad"
    assert _verdict_class("THRESHOLD-SPANNING — x") == "bad"
    assert _verdict_class("INDETERMINATE — x") == "warn"


@claim("EVT-DERIV", "n_pos<20 的推导在代码里可复算，且推导函数不改写锁定阈值")
def _():
    import sc_audit
    from sc_audit.floor import derive_event_threshold, event_threshold_table
    for n in ("derive_event_threshold", "event_threshold_table", "relative_halfwidth_cp"):
        assert hasattr(sc_audit, n), f"未导出 {n}"
    # 预注册容差表必须能被代码恢复（x_min 随 n 变化，此处以 n=200 复算）
    assert event_threshold_table(200) == {0.50: 17, 0.434: 21, 0.40: 24,
                                          0.33: 33, 0.25: 50}
    # n 依赖性：同一 ±43.4% 口径下 x_min 随 n 单调上升
    assert derive_event_threshold(400, 0.434) == 22
    assert derive_event_threshold(50, 0.434) == 16
    src = _src(sc_audit.audit_detection_floor)
    assert "derive_event_threshold" not in src, \
        "推导函数不得进入判定路径（阈值必须保持事前锁定）"


@claim("VERDICT-REACH", "每个判定标签都有可达场景（防死分支）")
def _():
    d = os.path.dirname(os.path.abspath(__file__))
    src = open(os.path.join(d, "test_identifiability.py"), encoding="utf-8").read()
    assert "test_all_verdicts_are_reachable" in src
    assert "expected - seen" in src, "缺少不可达标签的失败断言"


@claim("ID-NOMETA", "无可用的设计/可调整变量时不得报 COMPARABLE（缺元数据 ≠ 可比）")
def _():
    import sc_audit.identifiability as _id
    src = _src(_id.audit_identifiability)
    assert "_uses_real_var" in src, "缺少\"是否用到真实变量\"的判据"
    assert "ANNOTATION limitation" in src, "缺少\"这是注释局限\"的显式说明"
    # 判据必须基于花名册的可用变量，而不是"列是否存在"
    assert "_usable_vars" in src and "res.design_defining" in src and "res.adjustable" in src
    # C.4（2026-10-07）：守卫由单一条件扩为"无可用变量 **或** 支持域真空且未缩域"
    assert "not _uses_real_var" in src and "_vacuous_no_restriction" in src, \
        "分支缺失 → 会退回 COMPARABLE"
    assert "_ANNOTATION_LIMITATION_NOTE" in src, "缺少统一的\"注释局限\"提示"


@claim("FLOOR-SCHEMA", "审计一返回表在退化路径（基因全部未命中）也保持稳定模式")
def _():
    import sc_audit.floor as _fl
    src = _src(_fl.audit_detection_floor)
    assert "模式稳定" in src, "缺少模式稳定性的实现说明"
    assert "\"n_pos\", \"n_pos_case\"" in src, "核心列未纳入补齐名单"
    assert "if _col not in df.columns:" in src
