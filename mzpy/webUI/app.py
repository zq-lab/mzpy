"""
Flask application factory for the mzpy Web UI.
"""

from __future__ import annotations

import os
from pathlib import Path

from flask import Flask, send_from_directory

from . import config as cfg
from .auth import auth_bp
from .routes.dashboard import dash_bp
from .routes.data import data_bp
from .routes.analysis import analysis_bp
from .routes.dem import dem_bp
from .routes.venn import venn_bp
from .routes.enrichment import enrich_bp


def create_app() -> Flask:
    app = Flask(
        __name__,
        template_folder=str(Path(__file__).parent / "templates"),
        static_folder=str(Path(__file__).parent / "static"),
    )
    app.config["SECRET_KEY"] = cfg.get_or_create_secret()
    app.config["MAX_CONTENT_LENGTH"] = cfg.MAX_CONTENT_LENGTH
    app.config["UPLOAD_FOLDER"] = str(cfg.UPLOAD_FOLDER)

    # Ensure plot output dir exists
    os.makedirs(cfg.PLOT_FOLDER, exist_ok=True)

    # Register blueprints
    app.register_blueprint(auth_bp)
    app.register_blueprint(dash_bp)
    app.register_blueprint(data_bp)
    app.register_blueprint(analysis_bp)
    app.register_blueprint(dem_bp)
    app.register_blueprint(venn_bp)
    app.register_blueprint(enrich_bp)

    @app.context_processor
    def inject_globals():
        return {
            "passwd_is_set": cfg.passwd_is_set(),
        }

    # Serve plot PNGs from the session plot folder (outside static/)
    @app.route("/static/plots/<path:filename>")
    def serve_plot(filename):
        return send_from_directory(cfg.PLOT_FOLDER, filename)

    return app
