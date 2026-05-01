"""
Server-side session storage for metabolomics analysis pipelines.
"""

from __future__ import annotations

import json
import os
import pickle
import shutil
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from ..config import PLOT_FOLDER, SESSION_ROOT, UPLOAD_FOLDER


def _session_dir(sid: str) -> Path:
    d = SESSION_ROOT / sid
    d.mkdir(parents=True, exist_ok=True)
    return d


def get_session_id() -> str:
    """Return the current Flask session ID (create if needed)."""
    from flask import session
    if "sid" not in session:
        session["sid"] = uuid.uuid4().hex
    return session["sid"]


def session_path(filename: str = "") -> Path:
    sid = get_session_id()
    return _session_dir(sid) / filename


def save_upload(file_storage, filename: str) -> Path:
    """Save an uploaded file to the session upload folder."""
    dest = session_path("upload")
    dest.mkdir(parents=True, exist_ok=True)
    fpath = dest / filename
    file_storage.save(str(fpath))
    return fpath


def list_uploads() -> List[Path]:
    up = session_path("upload")
    if not up.exists():
        return []
    return sorted(up.iterdir())


def save_metab(metab_obj) -> None:
    """Pickle a Metab object."""
    fpath = session_path("metab.pkl")
    with open(fpath, "wb") as f:
        pickle.dump(metab_obj, f)


def load_metab() -> Any:
    """Load the pickled Metab object if present."""
    fpath = session_path("metab.pkl")
    if not fpath.exists():
        return None
    with open(fpath, "rb") as f:
        return pickle.load(f)


class _NumpyJSONEncoder(json.JSONEncoder):
    """Handle numpy scalars / arrays and pandas NA when serialising to JSON."""

    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.integer, np.floating, np.bool_)):
            return obj.item()
        if pd.isna(obj):
            return None
        return super().default(obj)


def save_json(name: str, data: dict) -> None:
    fpath = session_path(f"{name}.json")
    with open(fpath, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2, cls=_NumpyJSONEncoder)


def load_json(name: str) -> Optional[dict]:
    fpath = session_path(f"{name}.json")
    if not fpath.exists():
        return None
    with open(fpath, "r", encoding="utf-8") as f:
        return json.load(f)


def save_csv(name: str, df: pd.DataFrame) -> Path:
    fpath = session_path(f"{name}.csv")
    df.to_csv(fpath, index=False)
    return fpath


def load_csv(name: str) -> Optional[pd.DataFrame]:
    fpath = session_path(f"{name}.csv")
    if not fpath.exists():
        return None
    return pd.read_csv(fpath)


def save_plot(name: str, fig) -> str:
    """Save a matplotlib figure and return a web-accessible path."""
    sid = get_session_id()
    dest = PLOT_FOLDER / sid
    dest.mkdir(parents=True, exist_ok=True)
    fpath = dest / f"{name}.png"
    fig.savefig(fpath, dpi=300, bbox_inches="tight")
    return f"/static/plots/{sid}/{name}.png"


def save_plotnine(name: str, p) -> str:
    """Save a plotnine figure and return a web-accessible path."""
    sid = get_session_id()
    dest = PLOT_FOLDER / sid
    dest.mkdir(parents=True, exist_ok=True)
    fpath = dest / f"{name}.png"
    p.save(str(fpath), dpi=300)
    return f"/static/plots/{sid}/{name}.png"


def cleanup_session() -> None:
    """Remove all files for the current session."""
    sid = get_session_id()
    d = _session_dir(sid)
    if d.exists():
        shutil.rmtree(d)
    pd = PLOT_FOLDER / sid
    if pd.exists():
        shutil.rmtree(pd)
