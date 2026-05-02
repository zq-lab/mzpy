"""
Data upload and preview blueprint.
"""

from __future__ import annotations

import pandas as pd
from flask import Blueprint, flash, jsonify, redirect, render_template, request, url_for

from ..auth import login_required
from ..services.session_store import (
    load_json,
    load_metab,
    save_json,
    save_metab,
    save_upload,
)
from ...metab import Metab, read_csv, read_msd_ali

data_bp = Blueprint("data", __name__, url_prefix="/data")


@data_bp.route("/upload", methods=["GET", "POST"])
@login_required
def upload():
    if request.method == "POST":
        files = request.files.getlist("alignment_file")
        files = [f for f in files if f and f.filename != ""]
        if not files:
            flash("No data file selected", "error")
            return redirect(request.url)

        is_msdial = request.form.get("is_msdial") == "on"

        try:
            if is_msdial:
                # ── MS-Dial mode ──
                if len(files) == 1:
                    fpath = save_upload(files[0], files[0].filename)
                    metab = read_msd_ali(str(fpath), drop_null_ms=True)
                    filename = files[0].filename
                else:
                    parts = []
                    for f in files:
                        fpath = save_upload(f, f.filename)
                        parts.append(read_msd_ali(str(fpath), drop_null_ms=True))
                    metab = pd.concat(parts)
                    metab = Metab(metab)
                    metab = metab.wash()
                    metab[('_', 'kid')] = metab[('_', 'Metabolite name')].str.split("|", n=1).str[0]
                    if ('_', 'MS/MS spectrum') in metab.columns:
                        metab = metab.drop(columns=[('_', 'MS/MS spectrum')])
                    filename = ", ".join(f.filename for f in files)
            else:
                # ── Generic CSV/Excel mode ──
                group_file = request.files.get("group_file")
                if not group_file or group_file.filename == "":
                    flash("Group info table is required for non-MS-Dial data", "error")
                    return redirect(request.url)

                gpath = save_upload(group_file, group_file.filename)

                if len(files) == 1:
                    fpath = save_upload(files[0], files[0].filename)
                    metab = read_csv(str(fpath), str(gpath))
                    filename = files[0].filename
                else:
                    parts = []
                    for f in files:
                        fpath = save_upload(f, f.filename)
                        parts.append(read_csv(str(fpath), str(gpath)))
                    metab = pd.concat(parts)
                    metab = Metab(metab)
                    filename = ", ".join(f.filename for f in files)

            save_metab(metab)
            save_json("upload_info", {
                "filename": filename,
                "rows": int(metab.shape[0]),
                "cols": int(metab.shape[1]),
            })
            flash(f"Loaded {metab.shape[0]} metabolites × {metab.shape[1]} columns", "success")
            return redirect(url_for("data.preview"))
        except Exception as exc:
            flash(f"Parse error: {exc}", "error")
            return redirect(request.url)

    return render_template("data_upload.html")


@data_bp.route("/preview", methods=["GET", "POST"])
@login_required
def preview():
    metab = load_metab()
    if metab is None:
        flash("Please upload an alignment file first", "warning")
        return redirect(url_for("data.upload"))

    info = load_json("upload_info") or {}

    # 自动从 Metab 的 MultiIndex 列中提取分组信息
    groups = {}
    samples = []
    for c in metab.columns:
        if c[0] != '_':
            groups.setdefault(c[0], []).append(c[1])
            samples.append({"name": c[1], "group": c[0]})
    save_json("groups", groups)

    # ── AJAX: 更新分组 ──
    if request.method == "POST" and request.form.get("action") == "update_groups":
        mapping = {}
        for key, val in request.form.items():
            if key.startswith("sample_group_"):
                sample = key.replace("sample_group_", "")
                mapping[sample] = val
        new_metab = metab.reassign_groups(mapping)
        save_metab(new_metab)
        flash("Group assignments updated", "success")
        return jsonify({"ok": True})

    # ── POST: 应用过滤 ──
    if request.method == "POST" and request.form.get("action") == "filter":
        before = int(metab.shape[0])

        # 化学注释过滤
        filter_chem = request.form.get("filter_chem") == "on"
        filter_msms = False
        total_score = 0.0
        if filter_chem:
            filter_msms = request.form.get("filter_msms") == "on"
            if filter_msms:
                if ('_', 'MS/MS matched') in metab.columns:
                    metab = metab[metab[('_', 'MS/MS matched')] == True]
                else:
                    flash("Column 'MS/MS matched' not found, skipping MS/MS filter", "warning")

            try:
                total_score = float(request.form.get("total_score", 0.0))
            except ValueError:
                total_score = 0.0

            if total_score > 0:
                if ('_', 'Total score') in metab.columns:
                    metab = metab[metab[('_', 'Total score')] > total_score]
                else:
                    flash("Column 'Total score' not found, skipping score filter", "warning")

        # 最小表达量过滤
        try:
            min_expr = float(request.form.get("min_expression", 0.0))
        except ValueError:
            min_expr = 0.0
        if min_expr > 0:
            metab = metab.filter_min_expression(min_expr)

        save_metab(metab)
        info["rows"] = int(metab.shape[0])
        save_json("upload_info", info)
        save_json("filter_info", {
            "type": "combined",
            "filter_chem": filter_chem,
            "filter_msms": filter_msms,
            "total_score": total_score,
            "min_expression": min_expr,
            "before": before,
            "after": int(metab.shape[0]),
        })
        flash(f"Filtered: {before} → {metab.shape[0]} metabolites", "success")
        return redirect(request.url)

    # 当前所有组名（用于下拉框）
    all_groups = sorted(groups.keys())

    preview_df = metab.head(20).to_html(
        classes="dataframe", border=0, index=False,
        max_cols=50, max_rows=20,
    )

    return render_template(
        "data_preview.html",
        info=info,
        preview=preview_df,
        groups=groups,
        samples=samples,
        all_groups=all_groups,
    )
