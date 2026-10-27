from __future__ import annotations

from pathlib import Path

import pytest

from viraprime.check import load_assays
from viraprime.models import Assay, PipelineError

F = "ACGTACGTACGTACGTAC"
R = "TTGGCCAATTGGCCAATT"


def _write(tmp_path: Path, data: bytes | str) -> Path:
    p = tmp_path / "assays.tsv"
    if isinstance(data, str):
        data = data.encode()
    p.write_bytes(data)
    return p


# ---------- load_assays ----------

def test_load_assays_header_comments_blank(tmp_path: Path) -> None:
    p = _write(tmp_path, f"name\tforward\treverse\n# c\n\nA1\t{F.lower()}\t{R}\n")
    assert load_assays(p, 5) == [Assay("A1", F, R)]


def test_load_assays_u_to_t(tmp_path: Path) -> None:
    p = _write(tmp_path, f"A1\tACGUACGUACGUACGUAC\t{R}\n")
    assert load_assays(p, 5)[0].forward == F


def test_load_assays_excel_export(tmp_path: Path) -> None:
    data = b"\xef\xbb\xbfname\tforward\treverse\r\nA1\t" + F.encode() + b"\t" + R.encode() + b"\t\t\r\n"
    assert load_assays(_write(tmp_path, data), 5) == [Assay("A1", F, R)]


def test_load_assays_strips_field_spaces(tmp_path: Path) -> None:
    p = _write(tmp_path, f" A1 \t {F} \t{R}\n")
    assert load_assays(p, 5) == [Assay("A1", F, R)]


def test_load_assays_missing_file(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="not found"):
        load_assays(tmp_path / "nope.tsv", 5)


def test_load_assays_header_only(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="no primer pairs"):
        load_assays(_write(tmp_path, "name\tforward\treverse\n"), 5)


def test_load_assays_wrong_field_count(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="line 1: expected 3 tab-separated columns"):
        load_assays(_write(tmp_path, f"A1\t{F}\n"), 5)


def test_load_assays_invalid_name(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="invalid name"):
        load_assays(_write(tmp_path, f"a b\t{F}\t{R}\n"), 5)


def test_load_assays_duplicate_names(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="duplicate.*A1"):
        load_assays(_write(tmp_path, f"A1\t{F}\t{R}\nA1\t{F}\t{R}\n"), 5)


def test_load_assays_invalid_character(tmp_path: Path) -> None:
    text = f"name\tforward\treverse\nA1\tACGTXCGTACGTACGTAC\t{R}\n"
    with pytest.raises(PipelineError, match=r"forward.*position 5.*'X'"):
        load_assays(_write(tmp_path, text), 5)


def test_load_assays_internal_space(tmp_path: Path) -> None:
    text = f"name\tforward\treverse\nA1\tACGT ACGTACGTACGTAC\t{R}\n"
    with pytest.raises(PipelineError, match="position 5"):
        load_assays(_write(tmp_path, text), 5)


def test_load_assays_too_short_for_window(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="shorter than --three-prime-window"):
        load_assays(_write(tmp_path, f"A1\tACGT\t{R}\n"), 5)


def test_load_assays_typo_on_first_line_is_not_a_header(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match=r"line 1: forward.*position 5.*'X'"):
        load_assays(_write(tmp_path, f"A1\tACGTXCGTACGTACGTAC\t{R}\n"), 5)
