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


# ---------- locate / check_assay ----------

import logging
import random
import time
from itertools import product

import primer3
from Bio.Seq import reverse_complement

from viraprime.check import TM_FOLD_LIMIT, check_assay, locate
from viraprime.conservation import build_consensus, score_columns
from viraprime.degenerate import tm_range
from viraprime.models import Alignment, Consensus, Location, ScoringParams
from viraprime.scoring import IUPAC, score_pair

PARAMS = ScoringParams(max_mismatches=2, window=5)


def _family(mutate: dict[int, dict[int, str]] | None = None) -> tuple[Alignment, Consensus]:
    """6 rows over a seeded 140 bp base, SNPs at 60 and 120, plus one mostly-gap column at 70."""
    rng = random.Random(3)
    base = [rng.choice("ACGT") for _ in range(140)]
    rows = []
    for k in range(6):
        seq = list(base)
        if k == 1:
            seq[60] = "A" if seq[60] != "A" else "C"
        if k == 2:
            seq[120] = "G" if seq[120] != "G" else "T"
        for pos, b in (mutate or {}).get(k, {}).items():
            seq[pos] = b
        seq.insert(70, "T" if k == 0 else "-")
        rows.append("".join(seq))
    aln = Alignment([f"s{k}" for k in range(6)], rows)
    return aln, build_consensus(score_columns(aln))


def _gc(seq: str) -> float:
    return 100 * sum(b in "GC" for b in seq) / len(seq)


def test_locate_exact() -> None:
    _, cons = _family()
    assert locate(cons.seq[30:50], cons.seq) == Location(30, 0)


def test_locate_best_with_mismatches() -> None:
    _, cons = _family()
    site = list(cons.seq[30:50])
    site[5] = "A" if site[5] != "A" else "C"
    assert locate("".join(site), cons.seq) == Location(30, 1)


def test_locate_degenerate() -> None:
    _, cons = _family()
    site = list(cons.seq[30:50])
    i = next(i for i, b in enumerate(site) if b in "AG")
    site[i] = "R"
    assert locate("".join(site), cons.seq) == Location(30, 0)


def test_locate_ambiguous() -> None:
    assert locate("ACGTAC", "TTACGTACGGACGTACTT") is None


def test_locate_too_long() -> None:
    assert locate("A" * 30, "ACGT") is None


def test_check_assay_matches_score_pair() -> None:
    aln, cons = _family()
    fwd, rev_site = cons.seq[10:30], cons.seq[80:100]
    res = check_assay(Assay("P", fwd, reverse_complement(rev_site)), aln, cons, PARAMS, 95)
    assert res.located is True
    assert res.coverage == score_pair(
        fwd, cons.cons_to_col[10:30], rev_site, cons.cons_to_col[80:100], aln, PARAMS
    )
    assert res.fwd == Location(10, 0) and res.rev == Location(80, 0)
    assert res.fwd_tm == pytest.approx(primer3.calc_tm(fwd)) and res.fwd_gc == pytest.approx(_gc(fwd))
    assert res.below_threshold is False


def test_check_assay_wrong_order() -> None:
    aln, cons = _family()
    res = check_assay(Assay("P", cons.seq[80:100], reverse_complement(cons.seq[10:30])), aln, cons, PARAMS, 95)
    assert res.located is False
    assert res.note == "reverse site is not downstream of forward site"
    assert res.coverage.coverage_pct == 0.0
    assert res.coverage.not_covered == {sid: "not located" for sid in aln.ids}
    assert res.below_threshold is True


def test_check_assay_overlap() -> None:
    aln, cons = _family()
    res = check_assay(Assay("P", cons.seq[10:30], reverse_complement(cons.seq[29:49])), aln, cons, PARAMS, 95)
    assert res.located is False
    assert res.note == "reverse site is not downstream of forward site"


def test_check_assay_degenerate_tm() -> None:
    aln, cons = _family()
    fwd = list(cons.seq[10:30])
    i = next(i for i, b in enumerate(fwd) if b in "AG")
    fwd[i] = "R"
    fwd = "".join(fwd)
    res = check_assay(Assay("P", fwd, reverse_complement(cons.seq[80:100])), aln, cons, PARAMS, 95)
    assert res.degenerate is True
    assert res.fwd_tm_range == tm_range(fwd)
    assert res.fwd_tm == pytest.approx(sum(res.fwd_tm_range) / 2)
    expansions = ["".join(v) for v in product(*(sorted(IUPAC[b]) for b in fwd))]
    assert res.fwd_gc == pytest.approx(sum(_gc(e) for e in expansions) / len(expansions))


def test_check_assay_high_fold_skips_tm() -> None:
    aln, cons = _family()
    fwd = "N" * 12 + cons.seq[22:30]
    t0 = time.monotonic()
    res = check_assay(Assay("P", fwd, reverse_complement(cons.seq[80:100])), aln, cons, PARAMS, 95)
    assert time.monotonic() - t0 < 5
    assert TM_FOLD_LIMIT == 256
    assert res.fwd_tm is None and res.fwd_tm_range is None
    assert "Tm not computed: fwd fold 16777216 > 256" in res.note


def test_check_assay_wrong_orientation_warns(caplog: pytest.LogCaptureFixture) -> None:
    aln, cons = _family()
    with caplog.at_level(logging.WARNING, logger="viraprime"):
        check_assay(Assay("P", cons.seq[10:30], cons.seq[80:100]), aln, cons, PARAMS, 95)
    assert any("wrong orientation" in r.getMessage() for r in caplog.records)


def test_check_assay_below_threshold() -> None:
    _, cons0 = _family()
    three_prime = cons0.seq[29]
    aln, cons = _family({3: {29: "A" if three_prime != "A" else "C"}})
    res = check_assay(Assay("P", cons.seq[10:30], reverse_complement(cons.seq[80:100])), aln, cons, PARAMS, 100)
    assert res.below_threshold is True
    assert "s3" in res.coverage.not_covered
