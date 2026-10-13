"""
sc-audit — 单细胞病例-对照研究的信号可信度审计工具
=====================================================

审计模块（见 docs/methods.md）:

  ① 检出地板审计 (detection-floor audit)
     判断目标基因在本数据集中是否具备【定量能力】。
     指标：检出率、检出比 R、绝对阳性细胞数（事件数）、二项检验、留一供者敏感性。

  ② 对照构成审计 (control-composition audit v0.3: identifiability + ladder)
     判断病例组与对照组是否【可比】。
     指标：条件 × 机构/取材/分选 交叉表、供者级总变差距离 TVD（细胞级作参考）、
           同队列匹配的效应量衰减。
     v0.3 新增：变量花名册完备性 → 病例组内熛分桶（design-defining / adjustable）
     → 共同支持域（可识别性） → 事前锁定的匹配规则阶梯 → 一致性检查。

  ③ 稳健信号对照 (robust-signal control)
     分离真信号与深度/技术假象，为审计提供【正面对照】。
     指标：深度匹配、深度分层后组间分离、管家基因内部对照。

用法:
    import sc_audit

    ad = sc_audit.load_data("mydata.h5ad")
    floor = sc_audit.audit_detection_floor(
        ad, genes=["MYL3", "PLN"], condition_col="condition",
        donor_col="donor", control_label="control",
    )
    print(floor.table)

命令行:
    sc-audit run --input mydata.h5ad --target-genes MYL3,PLN
    sc-audit floor --input mydata.h5ad --target-genes MYL3
"""

from .io import (load_data, ensure_obs_fields, compute_qc_metrics,
                 gene_detection_vector, resolve_counts, resolve_gene_names,
                 apply_gene_name_fallback, assert_donors_nested)
from .gating import gate_cells, CARDIAC_PRESETS
from .floor import (audit_detection_floor, DetectionFloorResult, derive_event_threshold,
                    event_threshold_table, relative_halfwidth_cp)
from .composition import audit_control_composition, CompositionResult
from .identifiability import (audit_identifiability, IdentifiabilityResult,
                              build_variable_roster, common_support,
                              positivity_overlap, build_matching_ladder, default_ladder)
from .robust import audit_robust_signal, RobustSignalResult
from .identity import audit_cell_identity, IdentityResult
from .reference import (compute_detection_ratio, reference_modalities,
                        reference_modality_set, load_reference_manifest,
                        load_reference_depth,
                        validate_reference_contract, enforce_reference_contract,
                        ReferenceContract, CrossModalityError,
                        ReferenceContractError, DepthMismatchError,
                        REFERENCE_FIVE_TUPLE)
from .latent_batch import audit_latent_batch, LatentBatchResult
from .check import check_readiness, ReadinessResult
from .config import AuditConfig, load_config
from .design import (DesignDeclaration, DesignGate, check_design_support,
                     structural_hints, COMPARISON_TYPES, UNITS,
                     SUPPORTED_COMPARISON_TYPES, SUPPORTED_UNITS,
                     EXIT_DESIGN_NOT_SUPPORTED)

__version__ = "0.8.4"

__all__ = [
    "load_data",
    "ensure_obs_fields",
    "compute_qc_metrics",
    "gene_detection_vector",
    # P2-3：Python API 用户所需的入口校验/预处理函数（与 CLI 共用同一实现）
    "resolve_counts",
    "resolve_gene_names",
    "apply_gene_name_fallback",
    "assert_donors_nested",
    "gate_cells",
    "CARDIAC_PRESETS",
    "audit_detection_floor",
    "DetectionFloorResult",
    "derive_event_threshold",
    "event_threshold_table",
    "relative_halfwidth_cp",
    "audit_control_composition",
    "CompositionResult",
    "audit_identifiability",
    "IdentifiabilityResult",
    "build_variable_roster",
    "common_support",
    "build_matching_ladder",
    "default_ladder",
    "audit_robust_signal",
    "RobustSignalResult",
    "audit_cell_identity",
    "IdentityResult",
    "compute_detection_ratio",
    # v0.5.0：参照契约（R1 五元组合并 / R2 跨 modality）
    "reference_modalities",
    "reference_modality_set",
    "load_reference_manifest",
    "validate_reference_contract",
    "enforce_reference_contract",
    "ReferenceContract",
    "CrossModalityError",
    "ReferenceContractError",
    "DepthMismatchError",
    "REFERENCE_FIVE_TUPLE",
    "load_reference_depth",
    "positivity_overlap",
    "audit_latent_batch",
    "LatentBatchResult",
    "check_readiness",
    "ReadinessResult",
    "AuditConfig",
    "load_config",
    # v0.4.0：设计声明与能力闸门（Tier 1）
    "DesignDeclaration",
    "DesignGate",
    "check_design_support",
    "structural_hints",
    "COMPARISON_TYPES",
    "UNITS",
    "SUPPORTED_COMPARISON_TYPES",
    "SUPPORTED_UNITS",
    "EXIT_DESIGN_NOT_SUPPORTED",
]
