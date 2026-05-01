"""
Venn diagram blueprint.
"""

from __future__ import annotations

from flask import Blueprint, flash, redirect, render_template, request, url_for

from ..auth import login_required
from ..services.plot_adapter import plot_venn
from ..services.session_store import load_json, save_json

venn_bp = Blueprint("venn", __name__, url_prefix="/venn")


@venn_bp.route("/diagram", methods=["GET", "POST"])
@login_required
def diagram():
    up_ids = load_json("up_ids") or []
    dn_ids = load_json("dn_ids") or []

    if not up_ids and not dn_ids:
        flash("Please run differential analysis first", "warning")
        return redirect(url_for("dem.volcano"))

    comparisons = load_json("dem_params") or {}
    g1 = comparisons.get("g1", "Group A")
    g2 = comparisons.get("g2", "Group B")

    data = {
        f"{g2} vs {g1} (up)": up_ids,
        f"{g2} vs {g1} (dn)": dn_ids,
    }

    if request.method == "POST":
        subset = request.form.get("subset", "up")
        if subset == "up":
            selected = up_ids
        elif subset == "dn":
            selected = dn_ids
        else:
            selected = list(set(up_ids) & set(dn_ids))

        from ..services.session_store import load_metab
        metab = load_metab()
        kegg_ids = []
        if metab is not None:
            for idx in selected:
                name = metab.loc[idx, "Metabolite name"] if "Metabolite name" in metab.columns else str(idx)
                import re
                m = re.search(r'(C\d{5})', str(name))
                if m:
                    kegg_ids.append(m.group(1))
        save_json("selected_kegg_ids", kegg_ids)
        save_json("selected_subset", subset)
        flash(f"Selected {len(kegg_ids)} KEGG compound IDs", "success")
        return redirect(url_for("enrichment.fella"))

    img_url = plot_venn(data, color_type="qualitative")
    return render_template(
        "venn_diagram.html",
        img_url=img_url,
        up_count=len(up_ids),
        dn_count=len(dn_ids),
    )
