"""
Differential Expression Metabolites (DEM) – multi-comparison volcano plots.
Supports t-test, VIP (from PLSDA), and Spearman (ordered 3-group) methods.
"""

from __future__ import annotations

import json
import uuid

import numpy as np
import pandas as pd
from flask import Blueprint, flash, redirect, render_template, request, url_for

from ..auth import login_required
from ..services.session_store import load_json, load_metab, save_json, save_plotnine, session_path

dem_bp = Blueprint("dem", __name__, url_prefix="/dem")


def _cols_for_group(groups_dict, g):
    """Return list of (group, sample) MultiIndex tuples for a group name."""
    return [(g, s) for s in groups_dict.get(g, [])]


def _calc_ttest(metab, g1, g2, xcut, ycut):
    """t-test + fold-change volcano. Returns volcano_df, img_name."""
    from scipy import stats

    cols1 = _cols_for_group(load_json("groups"), g1)
    cols2 = _cols_for_group(load_json("groups"), g2)
    df_num = metab[cols1 + cols2].apply(pd.to_numeric, errors="coerce").dropna(how="all")
    if df_num.empty:
        raise ValueError("No valid numeric data for selected groups")

    mean1 = df_num[cols1].mean(axis=1)
    mean2 = df_num[cols2].mean(axis=1)
    log2fc = np.log2((mean1 + 1e-12) / (mean2 + 1e-12))
    log2fc = log2fc.replace([np.inf, -np.inf], np.nan).fillna(0)

    pvals = []
    for idx in df_num.index:
        a = df_num.loc[idx, cols1].dropna()
        b = df_num.loc[idx, cols2].dropna()
        if len(a) > 1 and len(b) > 1 and a.std() > 0 and b.std() > 0:
            _, p = stats.ttest_ind(a, b)
        else:
            p = 1.0
        pvals.append(p)
    pvals = pd.Series(pvals, index=df_num.index)
    pvals = pvals.replace(0, 1e-300)
    # ttest_ind 在数据方差为 0 时返回 NaN，需填充为 1.0（无显著性）
    pvals = pvals.fillna(1.0)
    neglog10p = -np.log10(pvals)
    neglog10p = neglog10p.replace([np.inf, -np.inf], 0)

    volcano_df = pd.DataFrame({
        "index": df_num.index,
        "log2FC": log2fc.values,
        "negLog10P": neglog10p.values,
        "Metabolite name": metab[('_', 'Metabolite name')] if ('_', 'Metabolite name') in metab.columns else metab.index,
    })

    def classify(r):
        if r["log2FC"] >= xcut and r["negLog10P"] >= ycut:
            return "up"
        elif r["log2FC"] <= -xcut and r["negLog10P"] >= ycut:
            return "dn"
        return "no"
    volcano_df["sig"] = volcano_df.apply(classify, axis=1)

    from ...plot import plot_volcano as _plot_volcano
    p = _plot_volcano(
        volcano_df, x="log2FC", y="negLog10P", fill="sig",
        xcut=xcut, ycut=ycut,
    )
    return volcano_df, p


def _calc_vip(metab, g1, g2, xcut, ycut):
    """VIP + fold-change volcano. Returns volcano_df, plotnine object."""
    vip_records = load_json("vip_df")
    if vip_records is None:
        raise ValueError("VIP data not found. Please run PLSDA first.")

    vip_map = {row["index"]: row["VIP"] for row in vip_records}

    cols1 = _cols_for_group(load_json("groups"), g1)
    cols2 = _cols_for_group(load_json("groups"), g2)
    df_num = metab[cols1 + cols2].apply(pd.to_numeric, errors="coerce").dropna(how="all")
    if df_num.empty:
        raise ValueError("No valid numeric data for selected groups")

    mean1 = df_num[cols1].mean(axis=1)
    mean2 = df_num[cols2].mean(axis=1)
    log2fc = np.log2((mean1 + 1e-12) / (mean2 + 1e-12))
    log2fc = log2fc.replace([np.inf, -np.inf], np.nan).fillna(0)

    vip_vals = pd.Series([vip_map.get(i, 0) for i in df_num.index], index=df_num.index)

    volcano_df = pd.DataFrame({
        "index": df_num.index,
        "log2FC": log2fc.values,
        "VIP": vip_vals.values,
        "Metabolite name": metab[('_', 'Metabolite name')] if ('_', 'Metabolite name') in metab.columns else metab.index,
    })

    def classify(r):
        if r["log2FC"] >= xcut and r["VIP"] >= ycut:
            return "up"
        elif r["log2FC"] <= -xcut and r["VIP"] >= ycut:
            return "dn"
        return "no"
    volcano_df["sig"] = volcano_df.apply(classify, axis=1)

    from ...plot import plot_volcano as _plot_volcano
    p = _plot_volcano(
        volcano_df, x="log2FC", y="VIP", fill="sig",
        xcut=xcut, ycut=ycut, ylab="VIP",
    )
    return volcano_df, p


def _calc_spearman(metab, ordered_groups, corr_thd, p_thd):
    """Spearman correlation (ordered 3 groups). Returns volcano_df, plotnine object."""
    df_spearman, plot = metab.spearman(
        groups=ordered_groups,
        corr_thd=corr_thd,
        p_thd=p_thd,
        save_fig_to=None,
    )
    # df_spearman columns: kid, corr, pval, -log_pval, monot
    # Rename to match volcano format
    volcano_df = pd.DataFrame({
        "index": df_spearman.index,
        "log2FC": df_spearman["corr"].values,
        "negLog10P": df_spearman["-log_pval"].values,
        "Metabolite name": metab[('_', 'Metabolite name')] if ('_', 'Metabolite name') in metab.columns else metab.index,
    })
    volcano_df["sig"] = df_spearman["monot"]

    # Re-plot with consistent styling
    from ...plot import plot_volcano as _plot_volcano
    p = _plot_volcano(
        volcano_df, x="log2FC", y="negLog10P", fill="sig",
        xcut=corr_thd, ycut=-np.log10(p_thd),
        xlab=r"$r_s$", ylab=r"-$\log_{10}(p)$",
    )
    return volcano_df, p


@dem_bp.route("/volcano", methods=["GET", "POST"])
@login_required
def volcano():
    metab = load_metab()
    groups = load_json("groups")
    if metab is None:
        flash("Please upload data first", "warning")
        return redirect(url_for("data.upload"))
    if not groups or len(groups) < 2:
        flash("Need at least two groups for comparison", "warning")
        return redirect(url_for("data.preview"))

    group_names = list(groups.keys())
    comparisons = load_json("comparisons") or []

    if request.method == "POST":
        action = request.form.get("action", "add")

        if action == "delete":
            cid = request.form.get("cid")
            comparisons = [c for c in comparisons if c.get("id") != cid]
            save_json("comparisons", comparisons)
            flash("Comparison removed", "success")
            return redirect(request.url)

        # ── Add new comparison ──
        name = request.form.get("name", "Untitled").strip()
        method = request.form.get("method", "ttest")
        xcut = float(request.form.get("xcut", 1.0))
        ycut = float(request.form.get("ycut", 2.0))

        try:
            if method == "spearman":
                g1 = request.form.get("sg1", "")
                g2 = request.form.get("sg2", "")
                g3 = request.form.get("sg3", "")
                if not all([g1, g2, g3]):
                    raise ValueError("Please select all three ordered groups for Spearman")
                if len({g1, g2, g3}) != 3:
                    raise ValueError("The three groups must be distinct")
                volcano_df, p = _calc_spearman(metab, [g1, g2, g3], xcut, 10 ** (-ycut))
                cmp = {
                    "id": uuid.uuid4().hex[:8],
                    "name": name,
                    "method": "spearman",
                    "spearman_groups": [g1, g2, g3],
                    "xcut": xcut,
                    "ycut": ycut,
                }
            elif method == "vip":
                g1 = request.form.get("g1", group_names[0])
                g2 = request.form.get("g2", group_names[1])
                volcano_df, p = _calc_vip(metab, g1, g2, xcut, ycut)
                cmp = {
                    "id": uuid.uuid4().hex[:8],
                    "name": name,
                    "method": "vip",
                    "numerator": g1,
                    "denominator": g2,
                    "xcut": xcut,
                    "ycut": ycut,
                }
            else:  # ttest
                g1 = request.form.get("g1", group_names[0])
                g2 = request.form.get("g2", group_names[1])
                volcano_df, p = _calc_ttest(metab, g1, g2, xcut, ycut)
                cmp = {
                    "id": uuid.uuid4().hex[:8],
                    "name": name,
                    "method": "ttest",
                    "numerator": g1,
                    "denominator": g2,
                    "xcut": xcut,
                    "ycut": ycut,
                }

            up = volcano_df[volcano_df["sig"] == "up"]
            dn = volcano_df[volcano_df["sig"] == "dn"]

            # Save plot image
            img_url = save_plotnine(f"volcano_{cmp['id']}", p)

            # 确保 DataFrame 字典化后无 numpy scalar / NaN 等 JSON 不兼容类型
            clean_df = volcano_df.copy()
            for col in clean_df.select_dtypes(include=[np.number]).columns:
                clean_df[col] = clean_df[col].astype(float)
            clean_df = clean_df.replace({np.nan: None, np.inf: None, -np.inf: None})

            cmp.update({
                "img_url": img_url,
                "volcano_df": clean_df.to_dict(orient="records"),
                "up_count": len(up),
                "dn_count": len(dn),
            })
            comparisons.append(cmp)
            save_json("comparisons", comparisons)
            flash(f"Comparison '{name}' added", "success")
        except Exception as exc:
            flash(f"Analysis error: {exc}", "error")

        return redirect(request.url)

    return render_template(
        "dem_volcano.html",
        group_names=group_names,
        comparisons=comparisons,
    )
