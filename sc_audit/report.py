"""审计报告生成（自包含 HTML，图以内嵌 base64 呈现）。"""
from __future__ import annotations

import base64
import io
import os
from datetime import datetime
from typing import Dict, Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


_CSS = """
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
     margin:0;padding:28px 34px;color:#1a1a1a;line-height:1.55;background:#fafafa}
h1{font-size:23px;margin:0 0 4px}
h2{font-size:17px;margin:28px 0 10px;padding-bottom:5px;border-bottom:2px solid #e3e3e3}
h3{font-size:14px;margin:18px 0 6px;color:#333}
.meta{color:#666;font-size:13px;margin-bottom:18px}
table{border-collapse:collapse;font-size:12.5px;background:#fff;margin:8px 0}
th,td{border:1px solid #dcdcdc;padding:5px 9px;text-align:left}
th{background:#f2f2f2;font-weight:600}
tr:nth-child(even){background:#fbfbfb}
.verdict{display:inline-block;padding:3px 10px;border-radius:11px;font-size:12px;font-weight:600}
.bad{background:#fdecea;color:#a32116;border:1px solid #f5c6c0}
.warn{background:#fff5e6;color:#8a5300;border:1px solid #ffe0b3}
.ok{background:#eaf6ec;color:#1b5e20;border:1px solid #c3e6c8}
.note{background:#f7f7f7;border-left:3px solid #999;padding:8px 12px;margin:10px 0;font-size:12.5px}
.figure{margin:14px 0;max-width:100%}
code{background:#f0f0f0;padding:1px 5px;border-radius:3px;font-size:12px}
"""


def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _df_html(df: pd.DataFrame, cls: str = "") -> str:
    if df is None or not len(df):
        return "<p><em>（无数据）</em></p>"
    return df.to_html(index=False, border=0, classes=cls, escape=False)


def _verdict_class(v: str) -> str:
    """判定 → 配色。bad=红 / warn=橙 / ok=绿。

    2026-10-06 审计 A2-F3：v0.3 的负面判定（SUPPORT-BELOW-FLOOR / THRESHOLD-SPANNING）
    原先未收录，落到默认 `ok`（**绿色**）—— 把"设计不可识别"渲染得比 "NOT SEPARABLE" 更健康。
    同时保留旧标签 NON-OVERLAPPING / RULE-DEPENDENT 以兼容历史报告。
    """
    v = (v or "").upper()
    if ("UNQUANTIFIABLE" in v or "NOT SEPARABLE" in v
            or "SUPPORT-BELOW-FLOOR" in v or "NON-OVERLAPPING" in v
            or "THRESHOLD-SPANNING" in v or "RULE-DEPENDENT" in v):
        return "bad"
    if "DROPOUT" in v or "INSUFFICIENT" in v or "PARTIAL" in v or "INDETERMINATE" in v:
        return "warn"
    return "ok"


def _floor_figure(res) -> Optional[str]:
    df = res.table
    if df is None or "pct_study" not in df.columns:
        return None
    d = df.dropna(subset=["pct_study"]).sort_values("pct_study")
    if not len(d):
        return None
    fig, ax = plt.subplots(figsize=(9, max(3.0, 0.34 * len(d) + 1.4)))
    colors = ["#C44E52" if s == "event_sparse" else
              "#DD8452" if s == "event_limited" else "#4C72B0" for s in d["stratum"]]
    ax.barh(np.arange(len(d)), d["pct_study"], color=colors, edgecolor="white")
    if "pct_ref" in d.columns:
        ax.plot(d["pct_ref"].values, np.arange(len(d)), "D", ms=5, color="#2E7D32",
                label="best reference")
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels([f"{g}  (n={int(n)})" for g, n in zip(d["gene"], d["n_pos"])], fontsize=9)
    ax.axvline(5, color="red", ls="--", lw=1)
    ax.set_xlabel("Detection rate in target cells (%)")
    ax.set_title("Audit 1 · detection rates and absolute event counts\n"
                 "(red = event-sparse, n_pos < 20)", fontsize=11, loc="left")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    return _fig_to_b64(fig)


def _composition_figure(res) -> Optional[str]:
    if not res.crosstabs:
        return None
    panels = [(k, v) for k, v in res.crosstabs.items() if v.shape[0] >= 2][:3]
    if not panels:
        return None
    fig, axes = plt.subplots(1, len(panels), figsize=(4.1 * len(panels), 3.2), squeeze=False)
    pal = ["#1F4E79", "#D68910", "#196F3D", "#C0392B", "#8172B3"]
    for i, (name, ct) in enumerate(panels):
        ax = axes[0][i]
        left = np.zeros(ct.shape[0])
        for j, col in enumerate(ct.columns):
            vals = ct[col].values.astype(float)
            ax.barh(np.arange(ct.shape[0]), vals, left=left, height=0.6,
                    color=pal[j % len(pal)], edgecolor="white", label=str(col))
            left += vals
        ax.set_yticks(np.arange(ct.shape[0]))
        ax.set_yticklabels([str(x) for x in ct.index], fontsize=9)
        ax.set_xlabel("cells", fontsize=9)
        ax.set_title(f"condition × {name}", fontsize=10, loc="left")
        ax.legend(fontsize=7.5, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Audit 2 · condition vs design variables", fontsize=11, x=0.005, ha="left")
    return _fig_to_b64(fig)


def _robust_figure(res) -> Optional[str]:
    if not res.stratified:
        return None
    genes = list(res.stratified)[:4]
    fig, axes = plt.subplots(1, len(genes), figsize=(3.3 * len(genes), 3.4), squeeze=False, sharey=True)
    for i, g in enumerate(genes):
        ax = axes[0][i]
        t = res.stratified[g]
        x = np.arange(len(t)); w = 0.38
        ax.bar(x - w/2, t["ctrl_pct"], w, color="#4C72B0", label="control")
        ax.bar(x + w/2, t["case_pct"], w, color="#C44E52", label="case")
        ax.set_xticks(x); ax.set_xticklabels(t["stratum"], rotation=45, ha="right", fontsize=8)
        ax.set_title(g, fontsize=10)
        if i == 0:
            ax.set_ylabel("detection rate (%)", fontsize=9)
            ax.legend(fontsize=8, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Audit 3 · detection stratified by cell depth", fontsize=11, x=0.005, ha="left")
    return _fig_to_b64(fig)


def render_report(floor=None, composition=None, robust=None, reference_table=None,
                  outdir: str = ".", title: str = "sc-audit report",
                  input_name: str = "", extra_notes=None, info=None,
                  identifiability=None, latent_batch=None) -> str:
    """生成自包含 HTML 报告，返回文件路径。"""
    os.makedirs(outdir, exist_ok=True)
    H = []
    H.append(f"<html><head><meta charset='utf-8'><title>{title}</title>")
    H.append(f"<style>{_CSS}</style></head><body>")
    H.append(f"<h1>{title}</h1>")
    meta = [f"generated {datetime.now():%Y-%m-%d %H:%M}"]
    if input_name:
        meta.insert(0, f"input: <code>{os.path.basename(input_name)}</code>")
    if floor is not None:
        meta.append(f"{floor.n_cells} cells · {floor.n_donors} donors · "
                    f"control {floor.n_control} / case {floor.n_case}")
        s = floor.summary
        meta.append(f"R basis: <code>{s.get('ratio_basis')}</code> · "
                    f"permutation: {s.get('perm_mode')} · "
                    f"FDR: {s.get('fdr_method')} (m={s.get('fdr_n_tests', 0)})")
        if not s.get("has_reference"):
            meta.append("<b>no reference table supplied</b> — R-dependent verdicts "
                        "fall back to INSUFFICIENT (no reference)")
    H.append(f"<div class='meta'>" + " ｜ ".join(meta) + "</div>")

    if info and info.get("n_condition_levels", 1) > 1:
        H.append(f"<div class='note'><b>Multiple non-control levels detected</b> — "
                 f"all audits pool every non-control level into the case group "
                 f"(case_label=<code>{info.get('case_label')}</code>, "
                 f"scope={info.get('case_scope')}). If this is not intended, "
                 f"subset the data first.</div>")

    # ---- 审计一 ----
    if floor is not None:
        H.append("<h2>① Audit 1 · Detection floor</h2>")
        img = _floor_figure(floor)
        if img:
            H.append(f"<img class='figure' src='{img}'>")
        cols = [c for c in ["gene", "n_pos", "n_pos_case", "pct_study", "ctrl_pct",
                            "case_pct", "basis_pct", "pct_ref", "ratio_R",
                            "effect_size_pp", "donor_weighted_effect",
                            "binom_P_case0", "case_pct_upper95", "overall_pct_upper95",
                            "wilson_lo", "wilson_hi",
                            "stratum", "perm_P", "q_value_BH", "perm_mode",
                            "jackknife_p_range", "jackknife_n_sig",
                            "jackknife_n_below_bh_crit", "verdict",
                            "verdict_alt", "ratio_R_alt"] if c in floor.table.columns]
        H.append(_df_html(floor.table[cols]))
        if "q_value_BH" in floor.table.columns:
            n_q = int((floor.table["q_value_BH"] < 0.05).sum()) if len(floor.table) else 0
            H.append(f"<div class='note'>Multiplicity: BH-FDR applied to "
                     f"m = {floor.summary.get('fdr_n_tests', 0)} testable gene(s) "
                     f"(excluded from m: "
                     f"{', '.join(floor.summary.get('fdr_excluded_from_m', [])) or 'none'}); "
                     f"<b>{n_q}</b> gene(s) have q &lt; 0.05. Raw permutation P and "
                     f"FDR-adjusted q are reported side by side.</div>")
        if "verdict_alt" in floor.table.columns:
            H.append(f"<div class='note'>Reference-basis sensitivity: "
                     f"<b>{floor.summary.get('reference_alt_flips', 0)}</b> verdict(s) "
                     f"change under the alternative reference table "
                     f"(columns <code>verdict_alt</code> / <code>ratio_R_alt</code>).</div>")
        n_sp = floor.summary.get("n_event_sparse", 0)
        n_uq = floor.summary.get("n_unquantifiable", 0)
        H.append(f"<div class='note'><b>{n_uq}</b> gene(s) fall below the "
                 f"{floor.summary.get('floor_pct_threshold')}% detection floor and "
                 f"<b>{n_sp}</b> are event-sparse (n_pos &lt; "
                 f"{floor.summary.get('event_sparse_threshold')}): group-level comparisons "
                 f"for these genes should not be reported as findings.</div>")
        # L4（v0.7.0）：渲染 v0.6.x 新增的透明项（深度可比性 / 口径说明）
        if floor.summary.get("reference_depth_note"):
            H.append(f"<div class='note'>⚠ {floor.summary.get('reference_depth_note')}</div>")
        if floor.summary.get("binom_p_note"):
            H.append(f"<div class='note'>{floor.summary.get('binom_p_note')}</div>")
        if floor.summary.get("fdr_small_panel_note"):
            H.append(f"<div class='note'>{floor.summary.get('fdr_small_panel_note')}</div>")

    # ---- 审计二 ----
    if composition is not None:
        H.append("<h2>② Audit 2 · Control composition</h2>")
        H.append(f"<p><span class='verdict {_verdict_class(composition.verdict)}'>"
                 f"{composition.verdict}</span></p>")
        img = _composition_figure(composition)
        if img:
            H.append(f"<img class='figure' src='{img}'>")
        for k, ct in composition.crosstabs.items():
            H.append(f"<h3>condition × {k}</h3>")
            disp = ct.copy()
            disp.index.name = "condition"
            H.append(_df_html(disp.reset_index()))
        if composition.attenuation is not None:
            H.append("<h3>Effect-size attenuation (matched comparison)</h3>")
            H.append(_df_html(composition.attenuation))
        for n in composition.notes:
            H.append(f"<div class='note'>{n}</div>")

    # ---- 审计二 · v0.3 可识别性层 ----
    if identifiability is not None:
        H.append("<h2>② b Audit 2 v0.3 · Identifiability &amp; matching ladder</h2>")
        H.append(f"<p><span class='verdict {_verdict_class(identifiability.verdict)}'>"
                 f"{identifiability.verdict}</span></p>")
        # P1-4（2026-10-06 审计）：v0.2 单规则口径与 v0.3 层同时产出时，
        # 两个判定可能不一致（真实数据：NOT SEPARABLE vs SUPPORT-BELOW-FLOOR）——
        # 必须显式交代二者关系，否则读者无法判断以谁为准。
        if composition is not None:
            H.append(f"<div class='note'><b>How this relates to Audit 2 above.</b> "
                     f"The section above is the single-rule convention "
                     f"(<code>{composition.verdict.split(' — ')[0]}</code>), which matches on "
                     f"the configured variables once and reports the resulting attenuation. "
                     f"This section is the <b>v0.3 identifiability layer</b> "
                     f"(<code>{identifiability.verdict.split(' — ')[0]}</code>), which first asks "
                     f"whether a matched estimate exists at all and then how the attenuation "
                     f"behaves across a prespecified ladder of rules. When the two differ, the "
                     f"identifiability layer is the <b>recommended</b> reading; the single-rule "
                     f"verdict is retained for continuity with earlier versions and must not be "
                     f"quoted as the confounding share.</div>")
        if len(identifiability.roster):
            H.append("<h3>① Variable roster (every non-QC obs column must be bucketed)</h3>")
            disp = identifiability.roster.sort_values(["bucket", "tvd_donor"],
                                                      ascending=[True, False])
            H.append(_df_html(disp.reset_index(drop=True)))
        H.append(f"<div class='note'>completeness: {identifiability.n_obs_columns} "
                 f"columns = non-design {identifiability.n_nondesign_columns} + bucketed "
                 f"{len(identifiability.roster)} + unclassified "
                 f"{len(identifiability.unclassified)}<br>"
                 f"<b>design-defining</b> (constant within cases → scope statement, not "
                 f"matchable): {', '.join(identifiability.design_defining) or 'none'}<br>"
                 f"<b>adjustable</b> (varies within cases → true confounders): "
                 f"{', '.join(identifiability.adjustable) or 'none'}</div>")
        if identifiability.unclassified:
            H.append(f"<div class='note'>⚠ unclassified: "
                     f"{identifiability.unclassified}</div>")
        H.append("<h3>② Common support (identifiability before matching)</h3>")
        H.append(f"<div class='note'>case design profile(s): "
                 f"{identifiability.case_profiles}<br>control donors inside the profile: "
                 f"{identifiability.support_donors} → "
                 f"<b>{identifiability.support_cells}</b> cells "
                 f"(of {identifiability.n_control_donors} control donors) · "
                 f"pre-specified floor = {identifiability.cell_floor} cells</div>")
        _po = (identifiability.details or {}).get("positivity_overlap")
        if _po:
            H.append(f"<div class='note'>A2 donor-level positivity overlap "
                     f"(design-defining ∪ adjustable): case {_po.get('overlap_case')} · "
                     f"ctrl {_po.get('overlap_ctrl')} · min {_po.get('overlap_min')}"
                     f"{(' · vars: ' + ', '.join(_po.get('profile_vars'))) if _po.get('profile_vars') else ''}</div>")
        for c, t in identifiability.balance.items():
            H.append(f"<h3>support-domain balance · {c}</h3>")
            H.append(_df_html(t.reset_index().rename(columns={"index": c})))
        H.append("<h3>③ Prespecified matching ladder</h3>")
        if len(identifiability.ladder_summary):
            H.append(_df_html(identifiability.ladder_summary))
        if len(identifiability.ladder):
            piv = identifiability.ladder.pivot_table(index="rule", columns="gene",
                                                     values="ratio", aggfunc="first")
            H.append("<h4>attenuation ratio (matched / full effect), rules × genes "
                     "(R0 = anchor, ratio ≡ 1 by definition)</h4>")
            H.append(_df_html(piv.reset_index()))
        if identifiability.attenuated_suppressed:
            H.append("<div class='note'>⚠ attenuation ratios are <b>suppressed</b>: every "
                     "matched subset is decided by too few donors to constitute "
                     "evidence of the confounding share.</div>")
        for n in identifiability.notes:
            H.append(f"<div class='note'>{n}</div>")

    # ---- 审计三 ----
    if robust is not None:
        H.append("<h2>③ Audit 3 · Robust-signal control</h2>")
        cls = getattr(robust, "verdict_class", None) or _verdict_class(robust.verdict)
        H.append(f"<p><span class='verdict {cls}'>"
                 f"{robust.verdict}</span></p>")
        H.append("<h3>Depth comparability</h3>")
        H.append(_df_html(robust.depth_table.reset_index().rename(columns={"index": "metric"})))
        if getattr(robust, "errors", None):
            for e in robust.errors:
                H.append(f"<div class='note'>⚠ {e}</div>")
        if getattr(robust, "cell_state_verdict", ""):
            H.append(f"<div class='note'>{robust.cell_state_verdict}</div>")
        img = _robust_figure(robust)
        if img:
            H.append(f"<img class='figure' src='{img}'>")
        if robust.housekeeping is not None:
            H.append("<h3>Housekeeping internal control</h3>")
            H.append(_df_html(robust.housekeeping))
        for n in robust.notes:
            H.append(f"<div class='note'>{n}</div>")

    # ---- B1 潜在批次层（v0.6.0；L4：v0.7.0 渲染）----
    if latent_batch is not None:
        H.append("<h2>⑤ Latent-batch scan (exploratory, alarm-only)</h2>")
        cls = _verdict_class(getattr(latent_batch, "verdict", ""))
        H.append(f"<p><span class='verdict {cls}'>{latent_batch.verdict}</span></p>")
        _ari = getattr(latent_batch, "ari", float('nan'))
        H.append(f"<div class='note'>donors {latent_batch.n_donors} · genes "
                 f"{latent_batch.n_genes_used} · HVG {latent_batch.n_hvg} · "
                 f"ARI(cluster, condition) = {_ari:.3f} · seed {latent_batch.seed}</div>")
        for n in getattr(latent_batch, "notes", []):
            H.append(f"<div class='note'>{n}</div>")

    # ---- 跨数据集检出比 ----
    if reference_table is not None and len(reference_table):
        H.append("<h2>④ Cross-dataset detection ratio</h2>")
        H.append(_df_html(reference_table))
        H.append("<div class='note'>R = pct_study / pct_ref (best reference context). "
                 "Absolute detection rates are not compared across datasets.</div>")

    if extra_notes:
        H.append("<h2>Notes</h2><div class='note'>" + "<br>".join(map(str, extra_notes)) + "</div>")

    H.append("<h2>Interpretation guide</h2><div class='note'>"
             "· <b>UNQUANTIFIABLE</b>: detection rate below the floor — do not report "
             "group-level comparisons for this gene.<br>"
             "· <b>event-sparse</b> is a <i>separate axis</i>: such genes get no permutation "
             "P (and thus no q) because the absolute number of positive cells is too small; "
             "a gene can be event-sparse while still receiving a verdict from its detection "
             "rate. Check the <code>stratum</code> column independently of "
             "<code>verdict</code>.<br>"
             "· <b>INSUFFICIENT (no reference)</b>: no reference table was supplied, so R "
             "cannot be computed and the gene is not called quantifiable by default.<br>"
             "· <b>NOT SEPARABLE</b>: condition collinear with design variables — biological "
             "signal cannot be distinguished from institution/procurement/batch effects.<br>"
             "· <b>DEPTH-COMPARABLE</b>: separation survives depth stratification, so it is "
             "not a depth artefact (this is a positive control, not an endorsement of the "
             "biological claim).<br>"
             "· <b>SUPPORT-BELOW-FLOOR</b> (Audit 2 v0.3): the common support is below the "
             "prespecified cell floor and the case group is constant on the design-defining "
             "variable(s) — no adjusted between-group estimate exists, and the attenuation "
             "ratio is not reported as evidence.<br>"
             "· <b>THRESHOLD-SPANNING</b> (Audit 2 v0.3): the across-ladder median attenuation "
             "ratio moves across the prespecified threshold, so it is not a single estimand "
             "across rules.<br>"
             "· <b>COMPARABLE / NOT SEPARABLE</b> (Audit 2 v0.3, one-sided): all evaluable rules "
             "fall on the same side of the threshold; the verdict carries a mandatory caveat "
             "because the consistency check then has no information.<br>"
             "· <b>INDETERMINATE</b> (Audit 2 v0.3): nothing could be evaluated (no target gene, "
             "no evaluable rule, or the case group spans multiple design profiles).<br>"
             "· <b>LATENT-BATCH-RISK</b> (B1): the case/control split is collinear with donor "
             "expression clusters — compatible with an unregistered batch. ALARM ONLY; do not "
             "match on those clusters (they may encode biology).<br>"
             "· <b>POSITIVITY-OVERLAP CAVEAT</b> (A2): donor-level covariate overlap is "
             "insufficient; the ladder estimate applies only to the overlap population.</div>")
    H.append("</body></html>")

    out = os.path.join(outdir, "report.html")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(H))
    return out
