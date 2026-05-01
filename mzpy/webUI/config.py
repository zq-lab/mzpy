"""
Configuration and password management for the mzpy Web UI.
"""

from __future__ import annotations

import json
import os
import secrets
import tempfile
from pathlib import Path

from werkzeug.security import check_password_hash, generate_password_hash

# ── Paths ──────────────────────────────────────────────────────
CONFIG_DIR = Path.home() / ".config" / "mzpy"
PASSWD_FILE = CONFIG_DIR / "webui_passwd.json"
SECRET_FILE = CONFIG_DIR / "webui_secret.key"

# Server-side session storage (file-based, isolated per session)
SESSION_ROOT = Path(tempfile.gettempdir()) / "mzpy_webui_sessions"
SESSION_ROOT.mkdir(parents=True, exist_ok=True)

UPLOAD_FOLDER = SESSION_ROOT / "uploads"
UPLOAD_FOLDER.mkdir(parents=True, exist_ok=True)

PLOT_FOLDER = SESSION_ROOT / "plots"
PLOT_FOLDER.mkdir(parents=True, exist_ok=True)

MAX_CONTENT_LENGTH = 50 * 1024 * 1024  # 50 MB


# ── Password helpers ───────────────────────────────────────────

def save_passwd(password: str) -> None:
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    hashed = generate_password_hash(password)
    PASSWD_FILE.write_text(json.dumps({"hash": hashed}), encoding="utf-8")
    # Restrict permissions
    os.chmod(PASSWD_FILE, 0o600)


def remove_passwd() -> None:
    if PASSWD_FILE.exists():
        PASSWD_FILE.unlink()


def load_passwd_hash() -> str | None:
    if not PASSWD_FILE.exists():
        return None
    data = json.loads(PASSWD_FILE.read_text(encoding="utf-8"))
    return data.get("hash")


def verify_passwd(password: str) -> bool:
    h = load_passwd_hash()
    if h is None:
        return True
    return check_password_hash(h, password)


def passwd_is_set() -> bool:
    return PASSWD_FILE.exists()


# ── Flask secret key ───────────────────────────────────────────

def get_or_create_secret() -> str:
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    if SECRET_FILE.exists():
        return SECRET_FILE.read_text(encoding="utf-8").strip()
    key = secrets.token_hex(32)
    SECRET_FILE.write_text(key, encoding="utf-8")
    os.chmod(SECRET_FILE, 0o600)
    return key
