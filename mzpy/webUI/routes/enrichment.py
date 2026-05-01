"""
Fella enrichment analysis blueprint.
"""

from __future__ import annotations

from flask import Blueprint, flash, redirect, render_template, request, url_for

from ..auth import login_required
from ..services.session_store import load_json, save_json

enrich_bp = Blueprint("enrichment", __name__, url_prefix="/enrichment")


@enrich_bp.route("/fella", methods=["GET", "POST"])
@login_required
def fella():
    kegg_ids = load_json("selected_kegg_ids") or []
    if not kegg_ids:
        flash("Please select a subset from the Venn diagram first", "warning")
        return redirect(url_for("venn.diagram"))

    if request.method == "POST":
        organism = request.form.get("organism", "hsa").strip()
        h5_path = request.form.get("h5_path", "").strip()
        method = request.form.get("method", "diffusion")

        if not h5_path:
            flash("Please provide the path to the FELLA HDF5 database file", "error")
            return redirect(request.url)

        try:
            from ...fella import Fella
        except ImportError as exc:
            flash(f"Fella import error: {exc}", "error")
            return redirect(request.url)

        try:
            fella_obj = Fella(organism=organism, h5_path=h5_path)
            fella_obj.enrich(compounds=kegg_ids, methods=method)

            res_df = fella_obj.generate_results_table(method=method)
            enz_df = fella_obj.generate_enzymes_table(method=method)
            fig = fella_obj.plot(method=method)

            from ..services.session_store import save_plot
            img_url = save_plot("fella_network", fig)

            save_json("enrichment_info", {
                "organism": organism,
                "method": method,
                "input_count": len(kegg_ids),
                "result_count": int(res_df.shape[0]) if res_df is not None else 0,
            })

            return render_template(
                "enrichment_fella.html",
                img_url=img_url,
                result_table=(res_df.head(100).to_html(classes="dataframe", border=0, index=False)
                              if res_df is not None else None),
                enzyme_table=(enz_df.head(100).to_html(classes="dataframe", border=0, index=False)
                              if enz_df is not None else None),
                organism=organism,
                method=method,
            )
        except Exception as exc:
            flash(f"Enrichment error: {exc}", "error")
            return redirect(request.url)

    return render_template("enrichment_fella.html", kegg_count=len(kegg_ids))
