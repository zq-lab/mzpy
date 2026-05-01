
'''
国内源在安装时，依赖配置上容易出问题，最好不修改pip源的global设置
'''

#!/usr/bin/env python3
import platform
import shutil
import subprocess
import sys

if platform.system() == "Linux":
    has_pkg_config = shutil.which("pkg-config") is not None
    cairo_found = False
    if has_pkg_config:
        result = subprocess.run(
            ["pkg-config", "--exists", "cairo"],
            capture_output=True,
        )
        cairo_found = result.returncode == 0

    if not cairo_found:
        print(
            "\n"
            "============================================================\n"
            "  ERROR: System dependency 'cairo' is missing!\n"
            "  mzpy depends on 'pycairo', which requires the cairo\n"
            "  system library to compile.\n"
            "\n"
            "  Please install both packages first:\n"
            "\n"
            "      sudo apt install libcairo2-dev python3-cairo\n"
            "============================================================\n",
            file=sys.stderr,
        )
        sys.exit(1)

from setuptools import setup

setup()
