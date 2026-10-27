"""Checking primer pairs the user already has against every aligned sequence."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from viraprime.models import Assay, PipelineError
from viraprime.regions import NAME_RE
from viraprime.scoring import IUPAC


def _clean(field: str) -> str:
    return field.strip().upper().replace("U", "T")


def _looks_like_header(fwd: str) -> bool:
    """A column title like 'forward' is mostly not A/C/G/T; a primer, even with a typo, is."""
    return sum(ch in "ACGT" for ch in fwd) * 2 < len(fwd)


def _check_oligo(seq: str, what: str, where: str, window: int) -> None:
    for pos, ch in enumerate(seq, start=1):
        if ch not in IUPAC:
            raise PipelineError(f"{where}: {what} primer has an invalid character at position {pos}: '{ch}'")
    if len(seq) < window:
        raise PipelineError(
            f"{where}: {what} primer ({len(seq)} bases) is shorter than --three-prime-window ({window})"
        )


def load_assays(path: Path, window: int) -> list[Assay]:
    """Read a tab-separated name/forward/reverse file (optional header, # comments)."""
    path = Path(path)
    if not path.is_file():
        raise PipelineError(f"Primer file not found: {path}")
    assays: list[Assay] = []
    seen_data = False
    for lineno, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = [f.strip() for f in line.split("\t")]
        while fields and not fields[-1]:
            fields.pop()
        where = f"{path.name} line {lineno}"
        if len(fields) != 3:
            raise PipelineError(f"{where}: expected 3 tab-separated columns (name, forward, reverse)")
        name, fwd, rev = fields[0], _clean(fields[1]), _clean(fields[2])
        if not seen_data:
            seen_data = True
            if _looks_like_header(fwd):
                continue  # header line
        if not NAME_RE.fullmatch(name):
            raise PipelineError(f"{where}: invalid name '{name}' (use letters, digits, _ . -)")
        _check_oligo(fwd, "forward", where, window)
        _check_oligo(rev, "reverse", where, window)
        assays.append(Assay(name, fwd, rev))
    if not assays:
        raise PipelineError(f"Primer file contains no primer pairs: {path}")
    dups = sorted(n for n, c in Counter(a.name for a in assays).items() if c > 1)
    if dups:
        raise PipelineError(f"{path.name}: duplicate pair names: {', '.join(dups)}")
    return assays
