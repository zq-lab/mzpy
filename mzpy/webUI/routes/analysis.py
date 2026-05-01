"""
PCA / PLSDA / VIP filtering blueprint.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from flask import Blueprint, flash, redirect, render_template, request, url_for

from ..auth import login_required
from ..services.plot_adapter import plot_pca, plot_plsda
from ..services.session_store import load_json, load_metab, save_json, save_metab

analysis_bp = Blueprint("analysis", __name__, url_prefix="/analysis")


def _get_groups_and_metab():
    metab = load_metab()
    groups = load_json("groups")
    if metab is None:
        return None, None, redirect(url_for("data.upload"))
    if not groups:
        return None, None, redirect(url_for("data.preview"))
    return metab, groups, None


@analysis_bp.route("/pca", methods=["GET"])
@login_required
def pca():
    metab, groups, redir = _get_groups_and_metab()
    if redir:
        return redir

    # metab 列是 MultiIndex (group, sample)，需要构造对应的元组列表
    selected_cols = [(g, s) for g, samples in groups.items() for s in samples]
    grp_list = [g for g, samples in groups.items() for _ in samples]

    try:
        df_num = metab[selected_cols].apply(pd.to_numeric, errors="coerce")
        df_num = df_num.dropna(how="all").T  # 转置为 sample × metabolite
    except Exception as exc:
        flash(f"Data subset error: {exc}", "error")
        return redirect(url_for("data.preview"))

    img_url = plot_pca(df_num, groups=grp_list, add_ellipse=True)
    return render_template("analysis_pca.html", img_url=img_url, groups=groups)


@analysis_bp.route("/plsda", methods=["GET", "POST"])
@login_required
def plsda():
    metab, groups, redir = _get_groups_and_metab()
    if redir:
        return redir

    selected_cols = [(g, s) for g, samples in groups.items() for s in samples]
    grp_list = [g for g, samples in groups.items() for _ in samples]

    try:
        df_num = metab[selected_cols].apply(pd.to_numeric, errors="coerce")
        df_num = df_num.dropna(how="all").T  # 转置为 sample × metabolite
    except Exception as exc:
        flash(f"Data subset error: {exc}", "error")
        return redirect(url_for("data.preview"))

    if request.method == "POST":
        vip_threshold = float(request.form.get("vip_threshold", 1.0))
        vip_df = load_json("vip_df")
        if vip_df is None:
            flash("Please run PLSDA first", "warning")
            return redirect(request.url)

        vip_df = pd.DataFrame(vip_df)
        keep = vip_df[vip_df["VIP"] >= vip_threshold]["index"].tolist()
        metab_filt = metab.loc[keep].copy()
        save_metab(metab_filt)
        save_json("filter_info", {
            "type": "VIP",
            "threshold": vip_threshold,
            "before": int(metab.shape[0]),
            "after": int(metab_filt.shape[0]),
        })
        flash(f"Filtered by VIP ≥ {vip_threshold}: {metab.shape[0]} → {metab_filt.shape[0]}", "success")
        return redirect(url_for("dem.volcano"))

    try:
        from sklearn.cross_decomposition import PLSRegression

        y = pd.Categorical(grp_list).codes
        pls = PLSRegression(n_components=2)
        X = df_num.values
        pls.fit(X, y)
        scores = pls.x_scores_

        img_url = plot_plsda(scores, grp_list, data_transfer=None)

        W0 = pls.x_weights_ / np.sqrt(np.sum(pls.x_weights_ ** 2, axis=0))
        q = pls.y_loadings_
        SS = np.diag(q.T @ q).reshape(1, -1)
        VIP = np.sqrt(df_num.shape[1] * (W0 ** 2) @ SS.T / np.sum(SS))
        vip_df = pd.DataFrame({
            "index": df_num.columns,
            "VIP": VIP.ravel(),
            "Metabolite name": metab[('_', 'Metabolite name')] if ('_', 'Metabolite name') in metab.columns else metab.index,
        })
        save_json("vip_df", vip_df.to_dict(orient="records"))

        return render_template(
            "analysis_plsda.html",
            img_url=img_url,
            vip_table=vip_df.sort_values("VIP", ascending=False).head(50).to_html(
                classes="dataframe", border=0, index=False
            ),
        )
    except Exception as exc:
        flash(f"PLSDA error: {exc}", "error")
        return redirect(url_for("data.preview"))
