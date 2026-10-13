"""设计声明与能力闸门（v0.4.0 · Tier 1）。

背景
-----
v0.3.x 及以前，工具**默认假设**输入是"二元、横断面、分析单位 = 供者"的病例–对照设计。
对纵向 / 治疗臂 / 多状态设计，若照旧计算，会**静默给出无效判定**（例如把复测单位当独立
观测、把多臂的对照含义搞错）。本模块**不实现**这些设计，只做两件事：

  ① **设计声明**（`DesignDeclaration`）：由调用方**显式声明**比较结构、分析单位、对照对、
     时间轴等。程序**不推断**——因为"比较的类型"是实验设计知识，不是数据能唯一决定的。
  ② **能力闸门**（`DesignGate`）：声明落在支持集 S0 内 → 正常执行；否则**拒绝输出任何判定
     标签**，返回 `DESIGN-NOT-SUPPORTED`。注意它与 `INDETERMINATE` **语义不同**：
     `INDETERMINATE` 是"支持的设计下、审计做了但无法给结论"；`DESIGN-NOT-SUPPORTED` 是
     "设计超出能力、**根本没做审计**"。

支持集（Tier 0，唯一实现、冻结）::

    S0 = {comparison_type: "two_arm", unit: "donor", time: None, arms: 2}

设计原则：**把复杂度接口化并推迟**——缺信息就要、不懂就拒绝；程序不去替代人的实验设计
判断。所有声明字段默认值即 S0，因此**不声明 = 与 v0.3.3 行为一致**（见 CHANGELOG 0.4.0）。
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional

# ── 声明的**封闭词表**：词表外的取值 = 声明非法（错误，exit 2），而非"不支持"（exit 3）──
COMPARISON_TYPES = ("two_arm", "paired", "time_course", "multi_arm", "multi_state")
UNITS = ("donor", "donor_timepoint", "sample", "cell")

# ── 当前实现的支持集（Tier 0）：只有这一个组合会真正执行审计 ──
SUPPORTED_COMPARISON_TYPES = ("two_arm",)
SUPPORTED_UNITS = ("donor",)

# 支持集的**可读描述**（用于拒绝文案与人可读输出；与上面的常量保持一致）
SUPPORTED_SET_TEXT = ('comparison_type="two_arm", unit="donor", '
                      'no time axis (time_col is None), 2 arms')

# 退出码：沿用 CLI 约定 —— 0 正常 / 1 有告警 / 2 错误 / **3 设计不支持**
EXIT_DESIGN_NOT_SUPPORTED = 3


@dataclass
class DesignDeclaration:
    """调用方显式声明的设计描述符。全部默认值 = 支持集 S0。"""

    comparison_type: str = "two_arm"
    unit: str = "donor"
    contrast_a: Optional[str] = None      # 比较的 A 组标签（缺省由 condition/control 推）
    contrast_b: Optional[str] = None
    time_col: Optional[str] = None        # 纵向设计的时间列
    ref_time: Optional[str] = None        # 纵向参考时间点
    index_time: Optional[str] = None      # 索引时间（**仅记录**，用于解释"事前/事后"）
    index_event: Optional[str] = None     # 索引事件（**仅记录**）
    strata: List[str] = field(default_factory=list)   # 分层/匹配变量的显式声明（可选）

    @classmethod
    def from_config(cls, cfg) -> "DesignDeclaration":
        """从 AuditConfig 读取声明（缺字段时回退到 S0 默认，保证向后兼容）。"""
        return cls(
            comparison_type=getattr(cfg, "design_comparison_type", "two_arm"),
            unit=getattr(cfg, "design_unit", "donor"),
            contrast_a=getattr(cfg, "design_contrast_a", None),
            contrast_b=getattr(cfg, "design_contrast_b", None),
            time_col=getattr(cfg, "design_time_col", None),
            ref_time=getattr(cfg, "design_ref_time", None),
            index_time=getattr(cfg, "design_index_time", None),
            index_event=getattr(cfg, "design_index_event", None),
            strata=list(getattr(cfg, "design_strata", []) or []),
        )

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def is_default(self) -> bool:
        """是否为全默认声明（= 支持集 S0，未显式声明任何东西）。

        默认声明下 CLI **不额外打印**设计块，以使默认路径的 stdout
        与 v0.3.x 保持一致。
        """
        return self.to_dict() == DesignDeclaration().to_dict()

    def invalid_reasons(self) -> List[str]:
        """声明**非法**（取值不在封闭词表内）的原因。空列表 = 合法。

        非法 ≠ 不支持：非法是调用方写错（exit 2），不支持是设计尚未实现（exit 3）。
        """
        errs: List[str] = []
        if self.comparison_type not in COMPARISON_TYPES:
            errs.append(f"comparison_type={self.comparison_type!r} 不在 {list(COMPARISON_TYPES)}")
        if self.unit not in UNITS:
            errs.append(f"unit={self.unit!r} 不在 {list(UNITS)}")
        return errs

    def unsupported_reasons(self) -> List[str]:
        """声明**合法但超出当前支持集**的原因。空列表 = 在 S0 内。"""
        reasons: List[str] = []
        if self.comparison_type not in SUPPORTED_COMPARISON_TYPES:
            reasons.append(
                f"comparison_type={self.comparison_type!r} 未实现"
                f"（当前仅支持 {list(SUPPORTED_COMPARISON_TYPES)}）")
        if self.unit not in SUPPORTED_UNITS:
            reasons.append(
                f"unit={self.unit!r} 未实现（当前仅支持 {list(SUPPORTED_UNITS)}）")
        if self.time_col:
            reasons.append(f"time_col={self.time_col!r}：纵向设计未实现")
        return reasons


@dataclass
class DesignGate:
    """能力闸门结果。"""

    declared: Dict[str, Any]
    supported: bool
    status: str                          # "SUPPORTED" | "DESIGN-NOT-SUPPORTED"
    unsupported_reasons: List[str] = field(default_factory=list)
    supported_set: str = SUPPORTED_SET_TEXT
    exit_code: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def to_text(self) -> str:
        lines = ["-- design declaration (v0.4.0 Tier 1) --"]
        for k, v in self.declared.items():
            lines.append(f"   {k:<16}= {v}")
        lines.append(f"   {'status':<16}= {self.status}")
        lines.append(f"   {'supported_set':<16}= {self.supported_set}")
        if not self.supported:
            lines.append("   unsupported:")
            for r in self.unsupported_reasons:
                lines.append(f"     - {r}")
        return "\n".join(lines)

    def refusal_text(self) -> str:
        """拒绝时的完整人可读文案（供 CLI 打印到 stderr）。"""
        return (
            "[refused] DESIGN-NOT-SUPPORTED\n"
            + self.to_text() + "\n"
            "  本工具尚未实现该设计的审计；为避免给出**无效判定**，"
            "未执行任何审计（不会输出 verdict 标签）。\n"
            "  注意：这与 INDETERMINATE 不同 —— 后者是支持的设计下"
            "审计已执行但无法给结论。\n"
            "  说明见 docs/methods.md「设计声明与支持集」；"
            "如需支持该设计，请把数据与设计信息交给维护者。\n"
        )


def check_design_support(decl: DesignDeclaration, *,
                         n_condition_levels: Optional[int] = None) -> DesignGate:
    """评估声明是否落在支持集内（确定性，无随机）。

    `n_condition_levels` 仅用于**结构提示**（见 `structural_hints`），不参与"是否支持"的判定
    ——`two_arm` 下"把其余非对照合并为 case"是既有且已声明的语义，不因水平数 > 2 而拒绝。
    """
    reasons = decl.unsupported_reasons()
    supported = not reasons
    gate = DesignGate(
        declared=decl.to_dict(),
        supported=supported,
        status="SUPPORTED" if supported else "DESIGN-NOT-SUPPORTED",
        unsupported_reasons=reasons,
        exit_code=0 if supported else EXIT_DESIGN_NOT_SUPPORTED,
    )
    return gate


def structural_hints(ad, condition_col: str, donor_col: str,
                     time_col: Optional[str] = None) -> List[str]:
    """**只读的结构提示**：数据形似某种更复杂的设计，但当前按 S0 处理。

    与闸门不同，提示**绝不改变任何判定标签**——它只把"你可能需要显式声明"这件事浮出来
    （本工具的一贯原则：程序不推断，只索取）。
    """
    hints: List[str] = []
    try:
        if condition_col in ad.obs.columns:
            nlev = int(ad.obs[condition_col].astype(str).nunique())
            if nlev > 2:
                hints.append(
                    f"condition 列有 {nlev} 个水平（形似多水平/多状态）：当前按 "
                    f"comparison_type='two_arm' 处理，即把对照以外的所有水平合并为 case。"
                    f"若这不是你的意图，请用 --comparison-type 声明。")
    except Exception:  # noqa: BLE001 —— 提示失败不得影响主流程
        pass
    try:
        if time_col and time_col in ad.obs.columns:
            hints.append(
                f"已声明 time_col={time_col!r}：纵向设计当前**不支持**（闸门会拒绝）。")
    except Exception:  # noqa: BLE001
        pass
    return hints


__all__ = [
    "COMPARISON_TYPES", "UNITS",
    "SUPPORTED_COMPARISON_TYPES", "SUPPORTED_UNITS", "SUPPORTED_SET_TEXT",
    "EXIT_DESIGN_NOT_SUPPORTED",
    "DesignDeclaration", "DesignGate",
    "check_design_support", "structural_hints",
]
