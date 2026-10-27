from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from tests.synthetic import make_family, write_fasta

needs_mafft = pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft not installed")


@pytest.fixture
def family_base() -> str:
    return make_family()[0]


@pytest.fixture
def family(tmp_path: Path) -> Path:
    return write_fasta(tmp_path / "family.fasta", make_family()[1])
