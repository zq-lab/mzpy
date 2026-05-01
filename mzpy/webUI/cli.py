"""
CLI entry point for the mzpy Web UI.

Usage:
    mzpy                  Start the local server (127.0.0.1:5000)
    mzpy run --host 0.0.0.0 --port 8080
    mzpy set passwd       Interactively set an access password
    mzpy remove passwd    Remove the access password
"""

from __future__ import annotations

import argparse
import getpass
import os
import sys
import threading
import webbrowser

from .config import PASSWD_FILE, save_passwd, remove_passwd


def _open_browser(url: str, delay: float = 1.2) -> None:
    """Open the default web browser after a short delay (allow server to start)."""
    def _open() -> None:
        webbrowser.open(url)
    threading.Timer(delay, _open).start()


def _cmd_run(args):
    from .app import create_app
    app = create_app()
    url = f"http://{args.host}:{args.port}"
    print(f"Starting mzpy Web UI at {url}")
    print("Press Ctrl+C to stop")
    _open_browser(url)
    app.run(host=args.host, port=args.port, debug=args.debug)


def _cmd_set_passwd(_args):
    if os.path.exists(PASSWD_FILE):
        print(f"A password already exists at {PASSWD_FILE}")
        overwrite = input("Overwrite? [y/N]: ").strip().lower()
        if overwrite not in ("y", "yes"):
            print("Cancelled.")
            return

    pw = getpass.getpass("New password: ")
    if not pw:
        print("Password cannot be empty.")
        sys.exit(1)
    pw2 = getpass.getpass("Confirm password: ")
    if pw != pw2:
        print("Passwords do not match.")
        sys.exit(1)

    save_passwd(pw)
    print(f"Password saved to {PASSWD_FILE}")


def _cmd_remove_passwd(_args):
    if not os.path.exists(PASSWD_FILE):
        print("No password is currently set.")
        return
    remove_passwd()
    print("Password removed. The server will no longer require authentication.")


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="mzpy",
        description="mzpy Web UI – metabolomics analysis in the browser",
    )
    subparsers = parser.add_subparsers(dest="command")

    # run
    run_parser = subparsers.add_parser("run", help="Start the Flask server")
    run_parser.add_argument("--host", default="127.0.0.1", help="Bind host")
    run_parser.add_argument("--port", type=int, default=5000, help="Bind port")
    run_parser.add_argument("--debug", action="store_true", help="Enable debug mode")

    # set passwd
    set_parser = subparsers.add_parser("set", help="Manage settings")
    set_sub = set_parser.add_subparsers(dest="set_command")
    set_sub.add_parser("passwd", help="Set an access password")

    # remove passwd
    rm_parser = subparsers.add_parser("remove", help="Remove settings")
    rm_sub = rm_parser.add_subparsers(dest="rm_command")
    rm_sub.add_parser("passwd", help="Remove the access password")

    args = parser.parse_args()

    if args.command is None:
        # Default: start server
        from .app import create_app
        app = create_app()
        url = "http://127.0.0.1:5000"
        print(f"Starting mzpy Web UI at {url}")
        print("Press Ctrl+C to stop")
        _open_browser(url)
        app.run(host="127.0.0.1", port=5000, debug=False)
        return

    if args.command == "run":
        _cmd_run(args)
    elif args.command == "set" and args.set_command == "passwd":
        _cmd_set_passwd(args)
    elif args.command == "remove" and args.rm_command == "passwd":
        _cmd_remove_passwd(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
