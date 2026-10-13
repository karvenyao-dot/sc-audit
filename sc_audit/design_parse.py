"""`design` 子命令的确定性核心（解耦版 AI 设计解析，v0.8.0）。

定位（与 `tools`/docs 一致）：
  · 本模块**不调用任何语言模型、不联网**。它只定义"契约"：
      `export` —— 从数据产出 **inventory.json**（列清单 + 结构统计 + 确定性基线角色 + 可选 GEO 文本指针）；
      `apply`  —— 读回 **proposal.json**（由外部任一大模型/agent 生成），
                  做**确定性复验**、生成**人工确认单**与**配置声明**，并合并**列注册表**。
  · AI 在工具**之外**：谁生成 proposal 都行（本项目的 agent、以后换的模型、甚至手工）。
  · `apply` 只依赖 inventory.json（不需要原始 h5ad）→ **可完全离线**。
  · **绝不自动写配置**：`design_config.json` 是"待人工确认的声明"，不是已生效的配置。

v0.8.3 追加 **`propose`（可选内置调用器）**：把 inventory→proposal 变成一条命令（`--backend none|local|http`）。
⚠️ 这是**唯一**会触碰网络/子进程的路径，且**默认关闭**（`none`）；network/subprocess 均**延迟导入**，
故模块顶层与默认路径（export/apply/confirm/attest）**仍然零网络依赖**。
"""
from __future__ import annotations

import datetime as _dt
import json
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .donor_utils import (donor_group, donor_level_values, entropy as _entropy,
                          is_donor_unique)
from .identifiability import (BUCKET_ADJUSTABLE, BUCKET_DESIGN_DEFINING,
                              BUCKET_DONOR_SATURATED, BUCKET_DONOR_UNIQUE,
                              BUCKET_REDUNDANT_CONDITION, build_variable_roster,
                              redundant_with_condition)
from .check import _name_role
from .config import AuditConfig

# 提案角色词表 → sc-audit 配置角色
PROPOSAL_ROLES = ("condition", "unit", "redundant_with_condition", "design_defining",
                  "adjustable", "outcome", "intervention", "technical",
                  "donor_unique", "unknown")
_ROLE_TO_OVERRIDE = {"design_defining": "design", "adjustable": "confounder",
                     "outcome": "outcome", "intervention": "intervention",
                     "technical": "technical"}
# 风险分层（2026-10-12 用户确认"分不同层级"）
TIER_MUST_REVIEW = ("condition", "unit", "redundant_with_condition",
                    "design_defining", "adjustable", "unknown")
TIER_BATCH = ("outcome", "intervention")
TIER_AUTO = ("technical", "donor_unique")


def _tier(role: str, col: Optional[Dict] = None) -> str:
    if role in TIER_MUST_REVIEW:
        return "must_review"
    if role in TIER_BATCH:
        return "batch_confirm"
    # A.9（v0.8.2）：常量列（单水平）的 `technical` 与"scope 常量"结构不可区分
    # → 不得自动接受（必须人工）。这是对"自动接受无结构护栏"的修补。
    if role == "technical" and col is not None and int(col.get("n_levels", 2)) <= 1:
        return "must_review"
    return "auto_accept"


def _inventory_sha(inv: Dict) -> str:
    """inventory 的内容哈希（供 proposal 绑定输入快照；A/B 契约）。"""
    import hashlib
    core = {k: inv.get(k) for k in ("dataset", "condition_col", "donor_col",
                                    "control_label", "n_cells", "n_genes",
                                    "n_donors", "columns")}
    b = json.dumps(core, ensure_ascii=False, sort_keys=True).encode("utf-8")
    return hashlib.sha256(b).hexdigest()


# ══════════════════ export ══════════════════

def export_inventory(ad, condition_col: str, donor_col: str, control_label: str,
                     label: Optional[str] = None, geotext: Optional[str] = None,
                     cfg: Optional[AuditConfig] = None) -> Dict:
    """把数据集的 obs 结构导出为 inventory（**离线、确定性**）。"""
    cfg = cfg or AuditConfig()
    o = ad.obs
    cond = o[condition_col].astype(str)
    is_ctrl = (cond == str(control_label)).values
    donor = o[donor_col].astype(str)
    n_don = int(donor.nunique())
    # A.7#10（v0.8.2）：control/case 标签必须在数据里存在；否则各列的“两臂恒定/熵”全为
    # 常量，与真实结构不可区分（静默成功会造成对“常量”的误读）——直接报错。
    if int(is_ctrl.sum()) == 0:
        raise ValueError(
            f"--control-label='{control_label}' 在数据里没有任何匹配细胞（对照组为空）；"
            f"该数据集 condition 取值为 {sorted(cond.unique())[:8]}")
    if int((~is_ctrl).sum()) == 0:
        raise ValueError(
            f"除 control='{control_label}' 外没有病例细胞（病例组为空）；"
            f"condition 取值为 {sorted(cond.unique())[:8]}")

    cols = []
    for c in o.columns:
        s = o[c]
        sc = s.astype(str)
        per_donor_const = bool((sc.groupby(donor).nunique() <= 1).all())
        dl = donor_level_values(o, c, donor)
        dgrp = donor_group(is_ctrl, donor).reindex(dl.index)
        cols.append(dict(
            col=c, dtype=str(s.dtype), n_nan=int(s.isna().sum()),
            n_levels=int(sc.nunique()),
            per_donor_const=per_donor_const,
            donor_level_levels=int(dl.nunique()),
            donor_unique=bool(is_donor_unique(s, donor, n_don)),
            ent_case=round(float(_entropy(dl[dgrp == "case"])), 4),
            ent_ctrl=round(float(_entropy(dl[dgrp == "ctrl"])), 4),
            n_levels_case=int(dl[dgrp == "case"].nunique()),
            n_levels_ctrl=int(dl[dgrp == "ctrl"].nunique()),
            sample_values=[str(v) for v in pd.Series(sc.unique())[:8]],
        ))

    # 确定性基线角色（离线可用；AI 提案可选地覆盖/补位）
    baseline = {}
    if len(o.columns):
        R, _un = build_variable_roster(ad, condition_col, donor_col, control_label, cfg=cfg)
        baseline = {r["column"]: r["bucket"] for _, r in R.iterrows()}

    inv = dict(tool="sc-audit", schema="design_inventory/1",
               dataset=label or os.path.basename(str(getattr(ad, "filename", "dataset"))),
               condition_col=condition_col, donor_col=donor_col,
               control_label=str(control_label),
               n_cells=int(o.shape[0]), n_genes=int(ad.n_vars), n_donors=n_don,
               geotext=geotext, columns=cols, baseline_roles=baseline)
    inv["inventory_sha256"] = _inventory_sha(inv)
    return inv


def write_inventory(inv: Dict, outdir: str) -> str:
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, "inventory.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(inv, f, ensure_ascii=False, indent=1)
    return p


# ══════════════════ apply（确定性复验）══════════════════

def _reverify(col: Dict, role: str, inv: Dict) -> Tuple[str, str]:
    """把提案角色与**结构性判据**对照。返回 (verdict, note)。

    verdict ∈ agree / disagree / unverifiable。结构性判据**不交给 AI**。
    """
    c = col["col"]
    e0 = 1e-9
    if role == "condition":
        ok = (c == inv["condition_col"])
        return (("agree" if ok else "disagree"),
                "该列即声明的 condition 列" if ok else
                f"声明的 condition 列是 '{inv['condition_col']}'")
    if role == "unit":
        ok = (c == inv["donor_col"])
        return (("agree" if ok else "disagree"),
                "该列即声明的 donor 列" if ok else
                f"声明的 donor 列是 '{inv['donor_col']}'")
    if role == "redundant_with_condition":
        # A.8（v0.8.2）：与 roster 共用**同一判据** redundant_with_condition()
        ok = redundant_with_condition(col["ent_case"], col["ent_ctrl"], col["n_levels"])
        return (("agree" if ok else "disagree"),
                "两臂各自恒定且 ≥2 水平（condition 的重编码）" if ok else
                "结构上不满足 condition 副本判据（非两臂恒定/水平<2）")
    if role == "donor_unique":
        # A.8（v0.8.2）：与 roster 一致——confounder 角色的列即使供者内恒定且水平=供者数，
        # 也不是标识符（而属 v0.7.3 的“供者级全异协变量”桶）→ 不能判 donor_unique。
        if bool(col["donor_unique"]) and _name_role(c) == "confounder":
            return ("disagree",
                    "供者级协变量但每供者取值互异（不可匹配，也**非**标识符；应为策略：人工）")
        ok = bool(col["donor_unique"])
        return (("agree" if ok else "disagree"),
                "供者内恒定且水平数≥供者数" if ok else "非供者唯一")
    if role in ("design_defining", "adjustable"):
        if role == "design_defining":
            ok = col["ent_case"] <= e0
            return (("agree" if ok else "disagree"),
                    "病例组内恒定" if ok else "病例组内**有**变化（应为 adjustable）")
        ok = col["ent_case"] > e0
        return (("agree" if ok else "disagree"),
                "病例组内有变化" if ok else "病例组内恒定（应为 design_defining）")
    # outcome / intervention / technical / unknown：结构无法验证 → 纯语义
    return ("unverifiable", "结构性判据不适用（语义判断）")


def apply_proposal(inv: Dict, proposal: Dict, outdir: str,
                   registry_path: Optional[str] = None) -> Dict:
    """读回 proposal → 确定性复验 → 确认单 + 配置声明 + 注册表。"""
    os.makedirs(outdir, exist_ok=True)
    ds = inv.get("dataset")
    # A.6（v0.8.2）：数据集名不匹配时**直接报错**，绝不静默取 proposal 的第一条
    #（否则会把**别人的**设计提案按到本数据集上，且无任何异常迹象）。
    if ds not in proposal:
        raise ValueError(
            f"proposal 不含数据集 '{ds}'（现有键：{sorted(proposal.keys())}）；"
            "为避免静默换用他人提案，此处直接报错。请确认 --inventory 的 dataset "
            "与 --proposal 的键一致。")
    prop = proposal.get(ds) or {}
    pcols = {c["column"]: c for c in prop.get("columns", [])}
    inv_cols = {c["col"] for c in inv["columns"]}
    # A.7#11（v0.8.2）：记录 proposal 里**数据中不存在**的列（幻觉可测）。
    hallucinated = sorted(set(pcols) - inv_cols)
    generated_by = (prop.get("generated_by") if isinstance(prop, dict) else None) \
        or (proposal.get("generated_by") if isinstance(proposal, dict) else None)

    # 注册表（人工裁决缓存）：**已 confirmed 的列以注册表为准**，重跑不再重复问 AI。
    reg = {}
    if registry_path and os.path.exists(registry_path):
        try:
            reg = json.load(open(registry_path))
        except Exception:
            reg = {}
    confirmed = {c: e.get("role") for c, e in (reg.get(ds) or {}).items()
                 if isinstance(e, dict) and e.get("status") == "confirmed"}

    rows = []
    for c in inv["columns"]:
        name = c["col"]
        if name in confirmed:                       # 注册表复用：人工确认优先
            role = str(confirmed[name])
            rv, note = _reverify(c, role, inv)
            rows.append(dict(column=name, role=role, tier="confirmed", confidence=1.0,
                             reverify=rv, reverify_note=note, abstain=False,
                             evidence="registry(confirmed)",
                             baseline=inv.get("baseline_roles", {}).get(name, "")))
            continue
        p = pcols.get(name)
        role = str(p.get("role", "unknown")) if p else "unknown"
        conf = float(p.get("confidence", 0.0)) if p else 0.0
        ev = str(p.get("evidence", "")) if p else ""
        abstain = bool(p.get("abstain", False)) if p else True
        if abstain:
            role = "unknown"
        rv, note = _reverify(c, role, inv)
        rows.append(dict(column=name, role=role, tier=_tier(role, c), confidence=conf,
                         reverify=rv, reverify_note=note, abstain=abstain,
                         evidence=ev, baseline=inv.get("baseline_roles", {}).get(name, "")))

    df = pd.DataFrame(rows)

    # 配置声明（待人工确认；**不自动写**）
    role_overrides = {r["column"]: _ROLE_TO_OVERRIDE[r["role"]]
                      for _, r in df.iterrows() if r["role"] in _ROLE_TO_OVERRIDE}
    nondesign = [r["column"] for _, r in df.iterrows()
                 if r["role"] in ("redundant_with_condition", "donor_unique")]
    # identity_source_col：若某 adjustable/design 列名像来源列，提示为候选
    src_candidates = [r["column"] for _, r in df.iterrows()
                      if r["role"] in ("adjustable", "design_defining")
                      and any(k in r["column"].lower()
                              for k in ("source", "site", "institution", "center",
                                        "centre", "region", "procure", "sample_source"))]
    config = dict(condition_col=inv["condition_col"], donor_col=inv["donor_col"],
                  role_overrides=role_overrides,
                  identifiability_nondesign_cols=nondesign,
                  identity_source_col=(src_candidates[0] if src_candidates else None),
                  _note="本文件是**待人工确认**的设计声明，sc-audit 不会自动写入配置。")

    # 确认单（opt-out：人只动不同意的行）
    lines = [f"# 设计确认单 · {ds}", "", f"- 数据集：`{ds}`　细胞 {inv['n_cells']} / 基因 {inv['n_genes']} / 供者 {inv['n_donors']}",
             f"- condition=`{inv['condition_col']}`　donor=`{inv['donor_col']}`　control=`{inv['control_label']}`",
             "- 使用：**只修改你不同意的行**（默认全采纳）；'复验'列 = 确定性结构判据是否同意该角色。"]
    # B（v0.8.2）：proposal 出处（可溯源/钉版）
    if generated_by:
        _gb = generated_by
        _hash_ok = (inv.get("inventory_sha256") and _gb.get("inventory_sha256") == inv.get("inventory_sha256"))
        lines.append(f"- 提案出处：model=`{_gb.get('model','?')}`　prompt_sha256=`{str(_gb.get('prompt_sha256'))[:12]}`　"
                     f"inventory_sha256 匹配={_hash_ok}　ts={_gb.get('ts','?')}")
        lines.append(f"- ⚠️ **AI 生成、未经复核**：{_gb.get('disclaimer','')}")
    else:
        lines.append("- ⚠️ **本 proposal 无 generated_by 出处块**（未 attest）：无法证明模型/提示词/输入快照；"
                     "建议用 `design --action attest` 补上。")
    if hallucinated:
        lines.append(f"- ⚠️ **proposal 含数据中不存在的列 {len(hallucinated)} 个**（已忽略）：{hallucinated}")
    lines.append("")
    for tname, tkey in (("⓪ 已确认（注册表复用，不再重复问）", "confirmed"),
                        ("① 必须逐条确认", "must_review"),
                        ("② 批量确认", "batch_confirm"),
                        ("③ 自动接受（建议抽检）", "auto_accept")):
        sub = df[df["tier"] == tkey]
        if not len(sub):
            continue
        lines.append(f"## {tname}（{len(sub)} 列）")
        lines.append("")
        lines.append("| 列 | 提案角色 | 确信度 | 确定复验 | 依据/备注 |")
        lines.append("|---|---|---|---|---|")
        for _, r in sub.iterrows():
            lines.append(f"| `{r['column']}` | {r['role']} | {r['confidence']:.2f} | "
                         f"{r['reverify']} | {r['reverify_note']} |")
        lines.append("")
    disagree = df[df["reverify"] == "disagree"]
    lines.append(f"## ⚠️ 复验不一致（{len(disagree)} 列）")
    lines.append("")
    if len(disagree):
        lines.append("以下列的提案角色与**确定性结构判据矛盾**，必须人工裁决：")
        for _, r in disagree.iterrows():
            lines.append(f"- `{r['column']}`：提案 `{r['role']}` ↔ {r['reverify_note']}")
    else:
        lines.append("（无）")
    lines.append("")

    with open(os.path.join(outdir, "design_confirm_sheet.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    with open(os.path.join(outdir, "design_config.json"), "w", encoding="utf-8") as f:
        json.dump(config, f, ensure_ascii=False, indent=1)
    df.to_csv(os.path.join(outdir, "design_reverify.csv"), index=False)

    # 注册表（per-dataset 缓存人工裁决；重跑不再重复问）
    entry = reg.setdefault(ds, {})
    now = _dt.datetime.now().isoformat(timespec="seconds")
    for _, r in df.iterrows():
        if r["column"] not in entry:
            entry[r["column"]] = dict(role=r["role"], status="proposed",
                                      reverify=r["reverify"], ts=now)
    if registry_path:
        os.makedirs(os.path.dirname(registry_path) or ".", exist_ok=True)
        with open(registry_path, "w", encoding="utf-8") as f:
            json.dump(reg, f, ensure_ascii=False, indent=1)

    return dict(n_columns=len(df),
                n_must_review=int((df["tier"] == "must_review").sum()),
                n_batch=int((df["tier"] == "batch_confirm").sum()),
                n_auto=int((df["tier"] == "auto_accept").sum()),
                n_confirmed=int((df["tier"] == "confirmed").sum()),
                n_disagree=int(len(disagree)),
                n_abstain=int(df["abstain"].sum()),
                abstain_rate=round(float(df["abstain"].mean()), 4) if len(df) else 0.0,
                n_hallucinated=len(hallucinated),
                hallucinated_columns=hallucinated,
                attested=bool(generated_by), generated_by=generated_by,
                config=config)


def attest_proposal(proposal_path: str, inventory_path: str, model: str,
                    prompt_path: Optional[str] = None, out_path: Optional[str] = None,
                    notes: Optional[str] = None) -> Dict:
    """给 proposal 盖上**出处块**（可溯源/钉版；**不调 LLM、不联网**）。

    草案 §4.4「可溯源 + 钉版」的落地：记录 model / prompt 哈希 / inventory 哈希 /
    输入快照 / 时间 / 工具版本 + “AI 生成、未经复核”声明。重写 JSON（默认就地）。"""
    import hashlib
    prop = json.load(open(proposal_path, encoding="utf-8"))
    inv = json.load(open(inventory_path, encoding="utf-8"))
    ph = None
    if prompt_path and os.path.exists(prompt_path):
        ph = hashlib.sha256(open(prompt_path, "rb").read()).hexdigest()
    try:
        import sc_audit as _pkg
        _ver = getattr(_pkg, "__version__", "?")
    except Exception:
        _ver = "?"
    gb = dict(model=str(model), prompt_sha256=ph,
              inventory_sha256=inv.get("inventory_sha256"),
              dataset=inv.get("dataset"),
              ts=_dt.datetime.now().isoformat(timespec="seconds"),
              tool_version=_ver, attested=True,
              disclaimer="AI 生成、未经复核；须经人工确认后方可采信")
    if notes:
        gb["notes"] = str(notes)
    ds = inv.get("dataset")
    if ds in prop and isinstance(prop[ds], dict):
        prop[ds]["generated_by"] = gb        # 出处块跟数据集条目走
    else:
        prop["generated_by"] = gb
    out = out_path or proposal_path
    with open(out, "w", encoding="utf-8") as f:
        json.dump(prop, f, ensure_ascii=False, indent=1)
    return gb


def confirm_registry(registry_path: str, decisions: Dict[str, Dict[str, str]]) -> int:
    """人工确认入口：把 decisions={dataset:{col:role}} 写回注册表并标 confirmed。

    确认后的列在下一次 `apply` 中**不再重复问 AI**（注册表优先）。"""
    reg = {}
    if os.path.exists(registry_path):
        try:
            reg = json.load(open(registry_path))
        except Exception:
            reg = {}
    now = _dt.datetime.now().isoformat(timespec="seconds")
    n = 0
    for ds, cols in decisions.items():
        for col, role in cols.items():
            e = reg.setdefault(ds, {}).setdefault(col, {})
            e["role"] = str(role)
            e["status"] = "confirmed"
            e["ts"] = now
            n += 1
    os.makedirs(os.path.dirname(registry_path) or ".", exist_ok=True)
    with open(registry_path, "w", encoding="utf-8") as f:
        json.dump(reg, f, ensure_ascii=False, indent=1)
    return n


def apply_design_config_to_cfg(cfg, design_config: Dict) -> object:
    """把（**人工确认后**的）`design_config.json` 写入 AuditConfig。

    注意：本函数**不自动调用**——必须在人工确认后由调用方显式执行（闭环第二步）。
    """
    if design_config.get("role_overrides"):
        cfg.role_overrides = dict(design_config["role_overrides"])
    if design_config.get("identifiability_nondesign_cols"):
        cfg.identifiability_nondesign_cols = list(
            design_config["identifiability_nondesign_cols"])
    if design_config.get("identity_source_col"):
        cfg.identity_source_col = str(design_config["identity_source_col"])
    return cfg


# ══════════════════ propose（可选内置调用器，v0.8.3）══════════════════
# 说明：这是**唯一**会让工具触碰网络/子进程的路径，且**默认关闭**（backend="none"）。
# 延迟导入 requests/urllib/subprocess → 模块顶层与默认路径**仍然零网络依赖**。
# "谁调模型"由使用者显式选择；产物一律带 generated_by 出处块（可溯源）。

PROPOSE_SYSTEM = (
    "你是单细胞病例-对照研究里的『设计变量解析器』。只输出一个 JSON 对象，不要任何解释文字。"
    "对给定 inventory 的**每一个 obs 列**给出一个角色，词表恰好为："
    "condition（病例/对照对比列本身）、unit（分析单位=供者）、"
    "redundant_with_condition（与 condition 的重编码/更细切分，不得当设计变量）、"
    "design_defining（病例组内恒定的设计/范围变量）、adjustable（病例组内有变化的混杂）、"
    "outcome（结局/collider）、intervention（治疗/干预）、"
    "donor_unique（与供者 1:1 的样本/文库 id）、technical（QC/处理/批次/注释等非设计噪音）、"
    "unknown（无可引用证据→弃权）。规则：每个 role 必须给出 evidence（逐字引用 inventory 的取值/统计"
    "或 GEO 原文片段）；引不出就 abstain=true 且 role=\"unknown\"；不许编造 inventory 里没有的列。"
)


def build_propose_prompt(inv: Dict, prompt_text: Optional[str] = None,
                         geotext_text: Optional[str] = None) -> str:
    """构造喂给模型的提示词（若给了 `prompt_text` 则用其覆盖）。"""
    if prompt_text:
        body = prompt_text
    else:
        body = ("请为下列数据集产出设计提案 JSON。\n\n"
                "输出 schema：\n"
                '{"<dataset>": {"contrast_column": "<列名或 null>", "unit_column": "<列名>", '
                '"arm_labels": ["..."], "columns": [{"column": "...", "role": "...", '
                '"confidence": 0.0, "evidence": "...", "abstain": false}]}}\n\n')
    body += "\n【inventory】\n" + json.dumps(inv, ensure_ascii=False)
    if geotext_text:
        body += "\n\n【GEO 自由文本（增强模式）】\n" + geotext_text[:20000]
    return body


def _extract_json(text: str) -> Dict:
    """从模型输出里抽出第一个 JSON 对象（容忍 markdown 围栏/前后废话）。"""
    t = str(text).strip()
    if t.startswith("```"):
        t = t.strip("`")
        t = t.split("\n", 1)[-1] if "\n" in t else t
    i, j = t.find("{"), t.rfind("}")
    if i < 0 or j <= i:
        raise ValueError("模型输出里未找到 JSON 对象")
    return json.loads(t[i:j + 1])


def _generated_by(inv: Dict, model: str, prompt_text: str) -> Dict:
    import hashlib
    try:
        import sc_audit as _pkg
        _ver = getattr(_pkg, "__version__", "?")
    except Exception:
        _ver = "?"
    return dict(model=str(model),
                prompt_sha256=hashlib.sha256(prompt_text.encode("utf-8")).hexdigest(),
                inventory_sha256=inv.get("inventory_sha256"),
                dataset=inv.get("dataset"),
                ts=_dt.datetime.now().isoformat(timespec="seconds"),
                tool_version=_ver, backend=None, attested=True,
                disclaimer="AI 生成、未经复核；须经人工确认后方可采信")


def propose_via_http(inv: Dict, model: str, endpoint: str, api_key: Optional[str],
                     temperature: float = 0.0, timeout: int = 180,
                     prompt_text: Optional[str] = None,
                     geotext_text: Optional[str] = None) -> Tuple[Dict, str]:
    """调 OpenAI 兼容的 `/chat/completions`（**延迟导入 urllib**）。"""
    import urllib.request  # 延迟导入：默认路径不触发
    prompt = build_propose_prompt(inv, prompt_text, geotext_text)
    body = json.dumps({"model": model, "temperature": float(temperature),
                       "messages": [{"role": "system", "content": PROPOSE_SYSTEM},
                                    {"role": "user", "content": prompt}],
                       "response_format": {"type": "json_object"}}).encode("utf-8")
    headers = {"Content-Type": "application/json"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"
    req = urllib.request.Request(endpoint.rstrip("/") + "/chat/completions",
                                 data=body, headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = json.loads(r.read().decode("utf-8"))
    text = data["choices"][0]["message"]["content"]
    return _extract_json(text), prompt


def propose_via_local(inv: Dict, cmd: str, timeout: int = 900,
                      prompt_text: Optional[str] = None,
                      geotext_text: Optional[str] = None) -> Tuple[Dict, str]:
    """把提示词喂给**本地命令**（如 llama.cpp/ollama 的 CLI），取 stdout 的 JSON。"""
    import subprocess  # 延迟导入
    prompt = build_propose_prompt(inv, prompt_text, geotext_text)
    r = subprocess.run(cmd, shell=True, input=prompt.encode("utf-8"),
                       capture_output=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"本地命令失败 rc={r.returncode}: {r.stderr.decode()[:500]}")
    return _extract_json(r.stdout.decode("utf-8")), prompt


def propose(inv: Dict, backend: str = "none", model: Optional[str] = None,
            endpoint: Optional[str] = None, api_key: Optional[str] = None,
            local_cmd: Optional[str] = None, temperature: float = 0.0,
            timeout: int = 180, prompt_path: Optional[str] = None,
            out_path: Optional[str] = None) -> Dict:
    """内置调用器：inventory → proposal.json（带出处块）。backend ∈ none/local/http。

    `none`（默认）→ 报错（请自行生成 proposal 后用 `apply`）。
    """
    backend = str(backend or "none").lower()
    if backend == "none":
        raise ValueError("backend='none'：本命令不生成提案。请自行生成 proposal.json 后用 "
                         "`design --action apply`；或用 --backend local|http。")
    prompt_text = None
    if prompt_path and os.path.exists(prompt_path):
        prompt_text = open(prompt_path, encoding="utf-8").read()
    geotext_text = None
    if inv.get("geotext") and os.path.exists(str(inv["geotext"])):
        geotext_text = open(str(inv["geotext"]), encoding="utf-8").read()
    if backend == "http":
        if not endpoint or not model:
            raise ValueError("--backend http 需要 --endpoint 与 --model")
        prop, used_prompt = propose_via_http(inv, model, endpoint, api_key, temperature,
                                             timeout, prompt_text, geotext_text)
    elif backend == "local":
        if not local_cmd:
            raise ValueError("--backend local 需要 --local-cmd")
        prop, used_prompt = propose_via_local(inv, local_cmd, timeout, prompt_text,
                                              geotext_text)
    else:
        raise ValueError(f"未知 backend: {backend}")
    ds = inv.get("dataset")
    if ds not in prop:
        raise ValueError(f"模型输出不含数据集 '{ds}'（键：{sorted(prop.keys())}）")
    prop[ds]["generated_by"] = _generated_by(inv, model or backend, used_prompt)
    prop[ds]["generated_by"]["backend"] = backend
    out = out_path or os.path.join(os.path.dirname(inv.get("_path", "")), "proposal.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(prop, f, ensure_ascii=False, indent=1)
    return dict(out=out, dataset=ds, generated_by=prop[ds]["generated_by"])
