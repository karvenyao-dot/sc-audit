# sc-audit 文档导航

**用途：** 说明仓库里每一份文档是什么、以哪份为准。引用任何规则或数字前先看这里。

## 1. 使用与参考（现行）

| 文件 | 内容 |
|---|---|
| `README.md` | 项目说明、四项审计口径、子命令概览 |
| `GETTING_STARTED.md` | 安装 → 三步跑通 → 读结果 → 命令速查 → API → FAQ |
| `CHANGELOG.md` | 版本变更记录 |
| `docs/methods.md` | 方法学：操作定义、阈值、判定真值表（§7 为审计二 v0.3 层） |
| `docs/v0.3_threshold_contract.json` | 审计二 v0.3 层的阈值契约（机器可读） |
| `docs/audit-2-diagnostic-guide.md` | 审计二诊断指南与方法学措辞建议 |
| `examples/` | 合成演示数据生成器 + 示例输出（`bash examples/run_demo.sh` 可复跑） |

## 2. 代码内自证

| 位置 | 内容 |
|---|---|
| `tests/` | 契约测试（含"声称 ↔ 代码"字段级断言）|
| `sc_audit/config.py` | 全部事前锁定阈值（默认值不应因结果调整）|
| `.github/workflows/ci.yml` | CI：多 Python 版本跑 `pytest` |

> **关于测试跳过项**：公开树内共 **409** 项测试（collected），实测 **405 passed + 4 skipped**；跳过的 4 项为 2 项需外部数据集（环境变量门控，如
> `SC_AUDIT_EXTERNAL_CASE_H5AD`）与 2 项断言依赖**开发仓库的审计文档**（不随包发布）而跳过。其余全部执行。

## 3. 阈值速查（现行，事前锁定）

| 阈值 | 值 | 出处 |
|---|---|---|
| `floor_pct` | 5% | `config.py` |
| `event_sparse` / `event_limited` | 20 / 50 | 推导见 `docs/methods.md` §2.2 与 `floor.derive_event_threshold()` |
| `hetero_max` | 0.5（供者级 TVD） | `config.py`（事前锁定，非从数据推导） |
| `attenuation_min` / `_full_effect_pp` | 0.5 / 5 pp | `config.py` |
| `matched_min_cells_warn` | 50 | `config.py`（v0.3 共同支持域下限一并复用） |
| `identity_collinear_min_delta_pp` | 10 pp | `config.py` |

## 4. 不随本树发布的内部文档

以下名称在 `README.md` / `CHANGELOG.md` / `docs/methods.md` 的历史记录或说明中被提及，但**不随本发布树分发**（属开发仓库的内部审计与提案材料）：`docs/REFERENCE_SCHEMA_v1.md`（参照表契约设计稿）、`文献_逐层面映射_*.md`、`方法学修订提案_*.md`、`代码逻辑与漏洞评审_*.md`、`代码审计文档_*.md`、`tools/check_doc_code_sync.py`、`tools/regress_external_case.py`、`tools/build_audit_doc.py`。被提到不等于本树包含该文件。
