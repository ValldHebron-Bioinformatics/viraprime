"""Synthetic sequence family used by the tests and the README example."""

from __future__ import annotations

import random
from pathlib import Path

BASES = "ACGT"
LENGTH = 1200
N_SEQS = 10


def _substitute(rng: random.Random, base: str) -> str:
    return rng.choice([b for b in BASES if b != base])


def make_family(seed: int = 42) -> tuple[str, list[tuple[str, str]]]:
    """Return (base sequence, records) for a 10-member family.

    1-based layout of the base:
      * 61-540 identical in every sequence (conserved block)
      * 1-60 and 541-600: 3 SNPs per variant
      * 701-1100: every base substituted with probability 0.25 (hypervariable)
      * seq03 has a 3 bp deletion at 30-32, seq05 a 2 bp insertion after 570
    """
    rng = random.Random(seed)
    base = "".join(rng.choice(BASES) for _ in range(LENGTH))
    snp_positions = list(range(0, 60)) + list(range(540, 600))  # 0-based
    records: list[tuple[str, str]] = [("seq01", base)]
    for i in range(2, N_SEQS + 1):
        seq = list(base)
        for pos in rng.sample(snp_positions, 3):
            seq[pos] = _substitute(rng, seq[pos])
        for pos in range(700, 1100):
            if rng.random() < 0.25:
                seq[pos] = _substitute(rng, seq[pos])
        s = "".join(seq)
        if i == 3:
            s = s[:29] + s[32:]
        if i == 5:
            s = s[:570] + "TT" + s[570:]
        records.append((f"seq{i:02d}", s))
    return base, records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f">{rid}\n{seq}\n" for rid, seq in records))
    return path


def make_aligned_family(seed: int = 42) -> list[tuple[str, str]]:
    """make_family as a gapped alignment, so tests can use --aligned without MAFFT.

    seq03's 3-base deletion (base positions 30-32) becomes gaps, and every
    sequence except seq05 gets a 2-column gap where seq05 has its insertion.
    """
    rows = []
    for rid, seq in make_family(seed)[1]:
        if rid == "seq03":
            seq = seq[:29] + "---" + seq[29:]
        if rid != "seq05":
            seq = seq[:570] + "--" + seq[570:]
        rows.append((rid, seq))
    return rows
