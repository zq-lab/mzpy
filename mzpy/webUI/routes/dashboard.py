"""
Dashboard / home page blueprint.
"""

from flask import Blueprint, render_template

from ..auth import login_required

dash_bp = Blueprint("dashboard", __name__)


@dash_bp.route("/")
@login_required
def index():
    return render_template("dashboard.html")
