"""Runs the browser engine's Node tests (tests/browser) when Node and Pyodide are available."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent / "browser"
ENGINE = Path(__file__).resolve().parent.parent / "site" / "app" / "engine"

pytestmark = pytest.mark.skipif(
    shutil.which("node") is None
    or not (HERE / "node_modules" / "pyodide").is_dir()
    or not all((ENGINE / f"mafft-{name}.wasm").is_file()
               for name in ("disttbfast", "tbfast", "dndpre", "dvtditr")),
    reason="needs node, `npm install` in tests/browser and the built engine files",
)


def test_browser_engine_node_tests() -> None:
    proc = subprocess.run(
        ["node", "--test", "mafft.test.mjs", "engine.test.mjs"], cwd=HERE,
        env={**os.environ, "PYTHON": sys.executable}, capture_output=True, text=True, timeout=1200,
    )
    assert proc.returncode == 0, proc.stdout[-4000:] + proc.stderr[-2000:]
