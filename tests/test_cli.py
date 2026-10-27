from __future__ import annotations

from pathlib import Path

import pytest

from consensus_primers.cli import main
from tests.conftest import needs_mafft
from tests.synthetic import make_family

AMP = ["--amplicon-min", "100", "--amplicon-max", "300"]


def test_amplicon_args_required() -> None:
    with pytest.raises(SystemExit):
        main(["-i", "x.fa", "-o", "out"])


def test_input_and_aligned_exclusive() -> None:
    with pytest.raises(SystemExit):
        main(["-i", "x.fa", "--aligned", "y.fa", "-o", "out", *AMP])


def test_bad_tm_order_returns_1(family: Path, tmp_path: Path) -> None:
    rc = main(["-i", str(family), "-o", str(tmp_path / "o"), *AMP, "--tm-min", "62", "--tm-opt", "60"])
    assert rc == 1


def test_window_longer_than_primer_returns_1(family: Path, tmp_path: Path) -> None:
    rc = main(["-i", str(family), "-o", str(tmp_path / "o"), *AMP, "--three-prime-window", "30"])
    assert rc == 1


def test_missing_input_returns_1(tmp_path: Path) -> None:
    assert main(["-i", str(tmp_path / "missing.fa"), "-o", str(tmp_path / "o"), *AMP]) == 1


def test_unknown_reference_returns_1(family: Path, tmp_path: Path) -> None:
    rc = main(["-i", str(family), "-o", str(tmp_path / "o"), *AMP, "--reference-id", "nope"])
    assert rc == 1


def test_bad_region_returns_1(family: Path, tmp_path: Path) -> None:
    rc = main(["-i", str(family), "-o", str(tmp_path / "o"), *AMP, "--region", "500-100"])
    assert rc == 1


@needs_mafft
def test_region_outside_reference_returns_1(family: Path, tmp_path: Path) -> None:
    rc = main(["-i", str(family), "-o", str(tmp_path / "o"), *AMP, "--region", "1000-5000"])
    assert rc == 1


@needs_mafft
def test_reference_id_with_description(tmp_path: Path) -> None:
    _, recs = make_family()
    fa = tmp_path / "desc.fa"
    fa.write_text("".join(f">{rid.replace('seq', 'seqX')} some description\n{s}\n" for rid, s in recs))
    rc = main(["-i", str(fa), "-o", str(tmp_path / "o"), *AMP, "--reference-id", "seqX02",
               "--region", "100-500"])
    assert rc == 0
