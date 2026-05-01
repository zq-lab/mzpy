"""
Authentication blueprint.
"""

from __future__ import annotations

import functools

from flask import (
    Blueprint,
    flash,
    redirect,
    render_template,
    request,
    session,
    url_for,
)

from . import config as cfg

auth_bp = Blueprint("auth", __name__, url_prefix="/auth")


def login_required(view):
    """Decorator: require login when a password is set."""
    @functools.wraps(view)
    def wrapped(*args, **kwargs):
        if cfg.passwd_is_set() and not session.get("logged_in"):
            return redirect(url_for("auth.login", next=request.path))
        return view(*args, **kwargs)
    return wrapped


@auth_bp.route("/login", methods=["GET", "POST"])
def login():
    if not cfg.passwd_is_set():
        return redirect(url_for("dashboard.index"))

    if request.method == "POST":
        password = request.form.get("password", "")
        if cfg.verify_passwd(password):
            session["logged_in"] = True
            nxt = request.args.get("next") or url_for("dashboard.index")
            return redirect(nxt)
        flash("Invalid password", "error")

    return render_template("login.html")


@auth_bp.route("/logout")
def logout():
    session.pop("logged_in", None)
    return redirect(url_for("auth.login"))
