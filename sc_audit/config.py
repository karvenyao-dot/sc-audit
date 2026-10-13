"""配置与阈值(全部为事前锁定值;默认值不应因结果而调整)。

修订记录（2026-10-01）：
  · 新增 hetero_col / hetero_max / attenuation_min —— 原为硬编码（审计二）
  · 新增 depth_delta_max_pct / strat_min_delta_pp / strat_ci_level /
    housekeeping_ratio_ok —— 原为硬编码（审计三）
  · 新增 fdr_method —— 多重比较校正（FDR 水平可锁定）
  · 新增 case_scope —— 明确多水平 condition 时 case 组的组成
  · 新增 count_layer / require_integer_counts —— 原始计数校验（B5）
  · 新增 match_vars / auto_filter_stratum_cols —— 匹配变量与展示变量分离（B2）

修订记录（2026-10-06，审计二 v0.3 可识别性层）：
  · 新增 7 个 identifiability_* 字段（见下方"审计二 · v0.3"块）
  · 同日多 Agent 审计后：`identifiability_require_both_strata` →
    **`identifiability_one_sided_verdict`**（旧写法使 COMPARABLE/NOT SEPARABLE 不可达）
  · `from_dict` 对未知键告警（P2-11，仅告警不报错，容忍历史配置文件）
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field, asdict
from typing import List, Optional, Dict, Any, Tuple

try:
    import yaml
    _HAS_YAML = True
except ImportError:  # pragma: no cover
    _HAS_YAML = False


# 打包内自带的预设配置目录(yaml 随包分发,见 pyproject.toml)
_PRESET_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "configs")


# ════════ 非设计列（单一真值源）════════
# 细胞级 QC / 分组 / 供者 / 细胞身份注释列：**不是设计变量**，不进变量花名册。
# 由 check.py 与 identifiability.py 共同导入（此前两份副本且不一致，见审计 A3/P2-12）。
DEFAULT_NONDESIGN_QC_COLS: Tuple[str, ...] = (
    "n_genes", "total_counts", "percent_mito", "n_counts",
    "n_genes_by_counts", "total_counts_mt", "pct_counts_mt",
    "doublet_score", "predicted_doublet", "log1p_total_counts",
    "log1p_n_genes_by_counts",
)


@dataclass
class AuditConfig:
    """审计阈值。默认值即事前锁定值,不应因结果而调整。"""

    # ════════ 审计一:检出地板 ════════
    floor_pct: float = 5.0            # pct_case < 5%  → 不可定量
    ratio_min: float = 0.3            # R < 0.3         → dropout-limited
    ratio_ok: float = 0.6             # R >= 0.6 且 pct >= pct_ok → 可定量
    pct_ok: float = 20.0
    # 检出比 R 的分子口径:'case'(病例组检出率,默认)
    #                       'overall'(数据集整体)| 'control'(对照组)
    ratio_basis: str = "case"
    # N4 修复:原 verdict_basis 与 ratio_basis 同名异义(后者在 floor 中曾是死参数)。
    # 现**只保留 ratio_basis 一个字段**,R 与判定共用同一分子口径。

    # ════════ 事件数分层(条件性敏感性分析)════════
    # event_sparse 的推导（外部审查 B2 要求把关系写在代码里）：
    #   病例组检出率 p̂ = x/n_case 用 Clopper–Pearson **精确双侧 95% 区间**估计，
    #   要求其**相对半宽** (upper−lower)/2/p̂ 低于事先声明的容差。锁定值 20
    #   对应容差 ±43.4%（即"要求相对精度优于 ±43.4%"）。
    #   ⚠️ x_min 随 n_case 变化（本机复算：n=50 → 16；n=100 → 19；n=200 → 21；
    #     n=400 → 22），故可移植的是**容差表**而不是数字本身；调用
    #     `sc_audit.floor.event_threshold_table(n)` 可自行复算。
    #   ⚠️ 该推导用的是**双侧**相对半宽，与本工具列 `case_pct_upper95`（**单侧**上界）
    #   不是同一个量（单侧 95% 上界 ≡ 双侧 90% 上界）；勿把两者混用，见 docs/methods.md §2.5。
    event_sparse: int = 20            # n_pos < 20      → 仅描述性统计
    event_limited: int = 50           # 20 <= n_pos < 50 → 附精确 CI 与 jackknife 范围
    #                                   n_pos >= 50     → 标准分析
    # 事件稀疏敏感性的替代阈值(输出 10/20/50 三档分层,不改主判定)
    event_sparse_grid: List[int] = field(default_factory=lambda: [10, 20, 50])

    # ════════ 置换检验 ════════
    n_perm: int = 10000
    seeds: List[int] = field(default_factory=lambda: [42])
    # 供者数 <= 该值时使用精确枚举(代替随机置换,消除种子依赖)
    # N4/L6 修复:实际的门是组合数;此上限仅作病态保护(勿设为小于实际供者数的值)
    exact_enum_max_units: int = 40
    exact_max_combinations: int = 5_000_000

    # ════════ 多重比较 ════════
    fdr_method: str = "bh"
    # BH 的族错误率;原实现在 floor.py 中硬编码 0.05,现提入 config
    fdr_alpha: float = 0.05            # 'bh' | 'none'（BH 的族错误率）
    # L4 修复:阳性对照基因(管家基因)不参与 FDR 的假设集合 m
    fdr_exclude_from_m: List[str] = field(default_factory=lambda: ["ACTB", "GAPDH"])

    # ════════ 细胞门控 ════════
    min_counts: int = 500
    min_genes: int = 200
    pos_markers: List[str] = field(default_factory=lambda: [
        "TTN", "MYH7", "MYBPC3", "ACTN2", "TNNI3", "MYL2"])
    neg_markers: List[str] = field(default_factory=lambda: [
        "PTPRC", "CD3E", "DCN", "COL1A1", "PECAM1", "VWF", "LYZ"])
    min_pos_markers: int = 2

    # ════════ 深度分层 ════════
    depth_edges: List[int] = field(default_factory=lambda: [
        0, 200, 400, 600, 800, 1200, 2000, 100000])
    min_cells_per_stratum: int = 30   # 原 5 导致检出率粒度 20pp 却用 5pp 判"持续"

    # ════════ 审计二:对照构成 ════════
    # 仅用于**交叉表展示**的变量(不参与匹配,避免 batch/lane 等样本级唯一值破坏匹配)
    strata_cols: List[str] = field(default_factory=lambda: [
        "institution", "donation_type", "source", "region", "batch", "lane", "library"])
    # 仅用于**匹配**的设计变量(按重要性排序;B2 修复)
    match_vars: List[str] = field(default_factory=lambda: [
        "institution", "donation_type", "source", "region"])
    # 自动剔除样本级/供者级唯一的列(B2)
    auto_filter_stratum_cols: bool = True
    # 是否执行"匹配供者"分析。False → 只报供者级 TVD(check 的策略 3:
    # 当设计变量与分组结构性共线时,放弃匹配是最诚实的做法)
    run_matched_comparison: bool = True
    # ── 变量角色的显式声明(优先级高于启发式;check 用)──
    # {列名: "design"|"confounder"|"technical"|"outcome"|"treatment"|"ignore"}
    # 用途:当启发式无法识别你的结局/治疗变量时(列名不含常见词根),由此显式声明。
    role_overrides: Dict[str, str] = field(default_factory=dict)
    # 领域词根扩展(追加到内置词根表之后)
    extra_outcome_keywords: List[str] = field(default_factory=list)
    extra_treatment_keywords: List[str] = field(default_factory=list)
    # ── 运行期保护:把结局/治疗变量写进 match_vars 会引入 collider bias ──
    # True → composition/run 在检测到这种组合时【中止】(除非 --allow-outcome-matching)
    forbid_outcome_matching: bool = True
    # 组间分布差异(总变差距离 TVD)超过该值 → NOT SEPARABLE(原硬编码 0.5)
    hetero_col: str = "institution"
    # hetero_max 的语义（外部审查 B3）：这是**事前锁定**的判据，不是从数据推导出来的值。
    # 语义：供者级 TVD（0 = 两组设计分布相同；1 = 完全分离）超过 0.5 即视为"设计变量与分组
    # 近乎共线"，因此生物学对比与技术/取材来源不可分离。0.5 是"两组的众数层开始互相缺失"
    # 的量级；它只用于**报警**，不产生任何估计值。已知的稳健性证据限于：本数据集上
    # ±10% 扰动（0.45/0.55）不改变判定；这在其它数据集上不构成"0.5 最优"的证明。
    hetero_max: float = 0.5
    # ── 以下仅供 check 的【排序建议】使用,不参与任何判定 ──
    # 匹配变量选择:TVD 低于此值 → 与分组近乎独立(安全但几乎不控制任何东西)
    match_var_tvd_low: float = 0.1
    # 衰减比中位数低于该值 → NOT SEPARABLE(原硬编码 0.5)
    attenuation_min: float = 0.5
    # 计算衰减比时只采用 full effect ≥ 该值(pp)的基因,避免噪声基因主导(S2 修复)
    attenuation_min_full_effect_pp: float = 5.0
    # 匹配子集小于该细胞数时强制提示低功效(原硬编码 50)
    matched_min_cells_warn: int = 50

    # ── 审计二 · v0.3 可识别性层(identifiability;2026-10-06 并入)──
    # 事前锁定的匹配规则阶梯;None → 使用内置阶梯 R0-R5(见 identifiability.default_ladder)。
    # 显式给出时 = 变量名列表的列表(空列表代表"不匹配(锚点)"),一旦给出即视为运行前锁定。
    identifiability_ladder: Optional[List[List[str]]] = None
    # 共同支持域(可识别性)下限;None → 复用 matched_min_cells_warn(单一真值源)
    identifiability_min_support_cells: Optional[int] = None
    # 保护①（2026-10-06 首次并入；同日审计 A1-F2 后改为下列形式）：
    # True  → 单侧情形（阈值落在观测档位区间之外）**直接给方向判定**
    #         （NOT SEPARABLE / COMPARABLE），但强制附"一致性检查无信息"警示；
    #         两档并存（阈值在区间内）→ THRESHOLD-SPANNING（不给单一结论）。
    # False → 两端均报 INDETERMINATE（旧保守行为，供对照）。
    identifiability_one_sided_verdict: bool = True
    # 低于支持域下限时抑制衰减比数字(不出量化证据)
    identifiability_suppress_below_floor: bool = True
    # 非设计列(细胞级 QC / 分组 / 供者 / 身份注释):不进花名册
    identifiability_nondesign_cols: List[str] = field(default_factory=list)
    # 由表达量派生的评分列(对其条件化 = collider)与前缀匹配
    identifiability_derived_cols: List[str] = field(default_factory=list)
    identifiability_derived_prefixes: List[str] = field(default_factory=lambda: ["score_"])

    # ════════ 审计三:稳健信号 ════════
    # |Δ genes/counts| ≥ 该值 → PARTIAL(原硬编码 20)
    depth_delta_max_pct: float = 20.0
    # 分层内"分离持续"的最小幅度(pp,且要求方向一致)(原硬编码 5)
    strat_min_delta_pp: float = 5.0
    strat_ci_level: float = 0.95      # Wilson CI 置信水平
    # L2 修复:CI 不重叠的箱比例达该值才判"分离持续"(原硬编码 0.6)
    strat_ci_frac_ok: float = 0.6
    housekeeping_ratio_ok: float = 0.6  # 管家基因 R 达标阈值(原硬编码 0.6)
    # 线粒体比例差异 ≥ 该值(pp 绝对值)→ 独立报告 CELL-STATE DIFFERS
    mito_delta_warn_pp: float = 5.0

    # ════════ 原始计数校验(B5)════════
    count_layer: Optional[str] = None   # 指定 counts 所在的 layer 名(None=自动探测)
    require_integer_counts: bool = True  # X 非非负整数时优先切到 counts,再否则报错
    # N6 修复:换源(X → counts layer / raw)后强制重算 obs 中的 QC 列
    recompute_qc_after_source_switch: bool = True

    # ════════ 设计声明（v0.4.0 · Tier 1）════════
    # 这组字段**只记录、不推断**：由调用方显式声明输入属于哪种比较设计；
    # 程序不尝试从数据猜。默认值 = 支持集 S0（二元、横断面、单位=供者），
    # 因此**不声明 = 与 v0.3.x 行为一致**。声明超出 S0 → 能力闸门拒绝出判定
    # （DESIGN-NOT-SUPPORTED，exit 3；见 sc_audit/design.py 与 docs/methods.md）。
    design_comparison_type: str = "two_arm"   # two_arm|paired|time_course|multi_arm|multi_state
    design_unit: str = "donor"                # donor|donor_timepoint|sample|cell
    design_contrast_a: Optional[str] = None    # 比较的 A 组标签（缺省由 condition/control 推）
    design_contrast_b: Optional[str] = None
    design_time_col: Optional[str] = None      # 纵向设计的时间列（未实现，声明即拒绝）
    design_ref_time: Optional[str] = None
    design_index_time: Optional[str] = None    # 索引时间（仅记录，用于解释“事前/事后”）
    design_index_event: Optional[str] = None   # 索引事件（仅记录）
    design_strata: List[str] = field(default_factory=list)  # 分层/匹配变量的显式声明（可选）

    # ════════ 参照契约（v0.5.0；REFERENCE_SCHEMA R1/R2）════════
    # R2 跨 modality：默认禁止（不同 modality 的检出率不可比）；
    # True = 显式放行（CLI --allow-cross-modality），放行后受影响基因标 cross_modality=True。
    # 仍不得用于 floor_pct / ratio_* 的可定量判定。
    allow_cross_modality: bool = False

    # ════════ 报告 ════════
    housekeeping_genes: List[str] = field(default_factory=lambda: ["ACTB", "GAPDH"])

    # ════════ v0.6.0 方法学修订（A/B/C）════════
    # A1 参照侧深度可比性：case 与参照的深度（中位 total_counts vs mean_umi）
    #   相对差超过此百分比 → 默认拒绝（--allow-depth-mismatch 放行）。事前锁定。
    reference_depth_delta_max_pct: float = 50.0
    # 显式放行深度不匹配（CLI --allow-depth-mismatch）；放行后 R 不得用于可定量判定
    allow_depth_mismatch: bool = False
    # A2 positivity / 协变量重叠：供者级设计剖面在另一组的支持比例低于此值 → NOT SEPARABLE
    #   （与细胞计数下限**并列**，不再混为一体）。事前锁定。
    positivity_overlap_min: float = 0.5
    # B2 管家护栏：各管家 R 的离散度（CV = std/mean）超过此值 → 告警“管家彼此不一致”
    housekeeping_cv_warn: float = 0.5
    # L1（v0.7.0）：管家基线“高基线”阈值（所有管家 pct ≥ 此值 → 声明护栏适用域局限）
    housekeeping_high_baseline: float = 0.8
    # B1 潜在批次层：病例-对照与供者技术簇的关联强度（调整后兰德指数 ARI）超过此值 → 报警
    latent_batch_ari_warn: float = 0.5

    # ────────────────── 审计四:细胞身份与捕获完整性(v0.2 新增)──────────────────
    # 由病例–对照单细胞数据的细胞身份取证实践逼出;前三项审计均无法发现它检出的问题。
    identity_col: str = "lineage_names"        # 细胞类型注释列(缺失则审计四报 INDETERMINATE)
    identity_source_col: str = "source"        # 分选/来源列(缺失则检查 5/6 跳过)
    identity_target: str = ""                  # 目标细胞类型;空则不设预期(检查 1/2 降级为报告)
    identity_complexity_col: str = "n_genes"   # 每细胞基因数所在 obs 列
    identity_nuclear_marker: str = "MALAT1"    # 核滞留转录本,作捕获完整性标志
    identity_long_genes: List[str] = field(default_factory=lambda: [
        "RYR2", "TTN", "CACNA1C", "KCNQ1", "SCN5A"])   # 长/低拷贝基因,作捕获完整性标志
    identity_complexity_bins: List[int] = field(default_factory=lambda: [200, 400, 600, 800, 1200])
    identity_bin_min_cells: int = 15           # 复杂度匹配箱的最小细胞数(双侧各需)
    identity_gate_margin: float = 5.0          # 目标谱系门控通过率 / 其他谱系中位数 的达标倍数
    identity_imbalance_min_cells: int = 20     # 亚群计数低于此不报不平衡(避免噪声)
    identity_imbalance_min_tvd: float = 0.5    # 亚群在供者间的分布 vs 目标细胞总体分布 的供者级 TVD 阈值
    identity_collinear_min_delta_pp: float = 10.0   # 阴性标志两臂差 ≥ 此值 → 门控共线警告
    # C4（v0.7.2，据第三方独立评估）：门控通过率的**绝对量**与**按臂差**判据。
    # 旧实现只比较“目标/其他谱系中位数”相对 margin，单谱系时 margin=None → 门控项永不触发；
    # 一个通过率 0.3%、且两臂差 3.4× 的门控规则仍会判 PASS（漏报）。
    identity_gate_pass_floor_pct: float = 1.0      # 目标自家门控通过率绝对下限(%);低于此 = 门控近乎失效
    identity_gate_arm_delta_max: float = 2.0       # 门控通过率的臂间比上限(max/min);超过 → 门控与分组共线
    identity_procurement_markers: List[str] = field(default_factory=lambda: [
        "PECAM1", "VWF", "DCN", "COL1A1"])     # 取材部位/方式标志(内皮/间质 ambient)


    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "AuditConfig":
        # P2-11（2026-10-06 审计）：未知键（含**拼错**的锁定阈值）不得静默忽略：
        # 历史键仅告警不报错（容忍旧配置文件），以保证事前锁定值不会被无声改变。
        known = {f for f in cls.__dataclass_fields__}
        unknown = sorted(k for k in (d or {}) if k not in known
                         and not str(k).startswith("_"))
        if unknown:
            import warnings
            warnings.warn(
                f"AuditConfig.from_dict: 忽略未知字段 {unknown}；"
                f"若本意是设定锁定阈值，请检查拼写（拼错会被静默忽略）。",
                UserWarning, stacklevel=2)
        return cls(**{k: v for k, v in (d or {}).items() if k in known})


def load_config(name_or_path: Optional[str] = None) -> AuditConfig:
    """加载配置。

    name_or_path:
      · None        → 返回默认阈值
      · 'cardiac'   → 加载随包分发的 configs/cardiac.yaml
      · 文件路径     → 加载该 yaml
    """
    if not name_or_path:
        return AuditConfig()

    path = name_or_path

    # 1 直接当作路径
    if os.path.exists(path):
        pass
    else:
        # 2 包内预设(唯一位置;包外旧的 ../configs/ 已删除,P3 修复死路径)
        cand = os.path.join(_PRESET_DIR, f"{name_or_path}.yaml")
        if not os.path.exists(cand):
            raise FileNotFoundError(
                f"配置未找到: {name_or_path!r}(已尝试 {cand})")
        path = cand

    if not _HAS_YAML:
        raise ImportError("读取 yaml 配置需要 pyyaml:pip install pyyaml")

    with open(path, encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    cfg = AuditConfig.from_dict(data)
    cfg._source_path = path  # 记录来源,供 metadata 使用
    return cfg
