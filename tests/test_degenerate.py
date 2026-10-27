from __future__ import annotations

from consensus_primers.degenerate import degenerate_pair, fold, tm_range
from consensus_primers.models import Alignment, ScoringParams

STRICT = ScoringParams(max_mismatches=0, window=5)
FWD = "ACGTTGCAAGCTTCGGATCC"
REV = "TTGACCGTAGGCATCGATGC"
FCOLS = list(range(20))
RCOLS = list(range(40, 60))


def row(fwd: str = FWD, rev: str = REV) -> str:
    return fwd + "A" * 20 + rev


def sub(seq: str, pos: int, base: str) -> str:
    return seq[:pos] + base + seq[pos + 1:]


def test_snp_becomes_iupac() -> None:
    aln = Alignment(["a", "b", "c", "d"], [row(), row(), row(sub(FWD, 3, "C")), row(sub(FWD, 3, "C"))])
    # FWD[3] is T; adding C gives Y
    f, r, cov = degenerate_pair(FWD, FCOLS, REV, RCOLS, aln, STRICT, max_fold=8)
    assert f[3] == "Y" and f[:3] + f[4:] == FWD[:3] + FWD[4:]
    assert r == REV
    assert cov.coverage_pct == 100.0


def test_reverse_site_degenerated() -> None:
    aln = Alignment(["a", "b"], [row(), row(rev=sub(REV, 10, "A"))])  # REV[10] is G -> R
    f, r, cov = degenerate_pair(FWD, FCOLS, REV, RCOLS, aln, STRICT, max_fold=8)
    assert r[10] == "R" and f == FWD and cov.coverage_pct == 100.0


def test_fold_cap_respected() -> None:
    variants = [row()] + [row(sub(FWD, p, "A" if FWD[p] != "A" else "C")) for p in (1, 4, 7, 10)]
    aln = Alignment([f"s{i}" for i in range(5)], variants)
    f, _, cov = degenerate_pair(FWD, FCOLS, REV, RCOLS, aln, STRICT, max_fold=4)
    assert fold(f) <= 4
    assert cov.coverage_pct < 100.0
    assert len(cov.covered) == 3


def test_three_prime_terminal_never_degenerate() -> None:
    aln = Alignment(["a", "b"], [row(), row(sub(FWD, 19, "A"))])
    f, r, cov = degenerate_pair(FWD, FCOLS, REV, RCOLS, aln, STRICT, max_fold=8)
    assert f == FWD and r == REV
    assert cov.covered == ["a"]


def test_reverse_three_prime_terminal_never_degenerate() -> None:
    aln = Alignment(["a", "b"], [row(), row(rev=sub(REV, 0, "G"))])
    _, r, _ = degenerate_pair(FWD, FCOLS, REV, RCOLS, aln, STRICT, max_fold=8)
    assert r == REV


def test_fold() -> None:
    assert fold("ARYN") == 16
    assert fold("ACGT") == 1


def test_tm_range_orders() -> None:
    lo, hi = tm_range("ACGTRACGTACGTACGTAC")
    assert lo < hi
    lo2, hi2 = tm_range("ACGTAACGTACGTACGTAC")
    assert lo2 == hi2
