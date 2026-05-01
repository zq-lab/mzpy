"""
Data upload and preview blueprint.
"""

from __future__ import annotations

from flask import Blueprint, flash, redirect, render_template, request, url_for

from ..auth import login_required
from ..services.session_store import (
    load_json,
    load_metab,
    save_json,
    save_metab,
    save_upload,
)
from ...metab import read_msd_ali

data_bp = Blueprint("data", __name__, url_prefix="/data")


@data_bp.route("/upload", methods=["GET", "POST"])
@login_required
def upload():
    if request.method == "POST":
        file = request.files.get("alignment_file")
        if not file or file.filename == "":
            flash("No file selected", "error")
            return redirect(request.url)

        fpath = save_upload(file, file.filename)

        try:
            metab = read_msd_ali(str(fpath), drop_null_ms=True)
            save_metab(metab)
            save_json("upload_info", {
                "filename": file.filename,
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
    for c in metab.columns:
        if c[0] != '_':
            groups.setdefault(c[0], []).append(c[1])
    save_json("groups", groups)

    # ── 化学注释过滤 ──
    if request.method == "POST":
        filter_chem = request.form.get("filter_chem") == "on"
        if filter_chem:
            before = int(metab.shape[0])

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

            save_metab(metab)
            info["rows"] = int(metab.shape[0])
            save_json("upload_info", info)
            save_json("filter_info", {
                "type": "chemical_annotation",
                "filter_msms": filter_msms,
                "total_score": total_score,
                "before": before,
                "after": int(metab.shape[0]),
            })
            flash(f"Filtered by annotation: {before} → {metab.shape[0]} metabolites", "success")
            return redirect(request.url)

    preview_df = metab.head(20).to_html(
        classes="dataframe", border=0, index=False,
        max_cols=50, max_rows=20,
    )

    return render_template(
        "data_preview.html",
        info=info,
        preview=preview_df,
        groups=groups,
    )
