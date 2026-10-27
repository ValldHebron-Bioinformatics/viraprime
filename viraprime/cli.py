"""Command-line interface and pipeline orchestration."""

from __future__ import annotations

import argparse
import logging
import shlex
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import NoReturn

from Bio.Seq import reverse_complement

from viraprime.align import run_aligner
from viraprime.conservation import build_consensus, score_columns
from viraprime.coords import CoordMap
from viraprime.degenerate import degenerate_pair, tm_range
from viraprime.design import run_primer3
from viraprime.io import (
    alignment_from_records,
    load_fasta,
    write_alignment,
    write_report,
    write_tsv,
)
from viraprime.lengths import filter_by_length, ungapped_length
from viraprime.models import (
    Alignment,
    CandidatePair,
    ColumnStats,
    Consensus,
    DesignParams,
    LengthFilterResult,
    PipelineError,
    Record,
    Region,
    RegionResult,
    RunReport,
    ScoredPair,
    ScoringParams,
)
from viraprime.ranking import failure_summary, rank
from viraprime.regions import resolve_regions
from viraprime.scoring import score_pair

log = logging.getLogger("viraprime")

EXIT_OK, EXIT_ERROR, EXIT_BELOW_THRESHOLD = 0, 1, 2


class _ListHandler(logging.Handler):
    def __init__(self) -> None:
        super().__init__()
        self.lines: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.lines.append(self.format(record))


class _Parser(argparse.ArgumentParser):
    """ArgumentParser whose usage errors exit 1 (exit 2 means 'below threshold')."""

    def error(self, message: str) -> NoReturn:
        self.print_usage(sys.stderr)
        self.exit(EXIT_ERROR, f"{self.prog}: error: {message}\n")


def build_parser() -> argparse.ArgumentParser:
    p = _Parser(
        prog="viraprime",
        description="Design consensus PCR primers from a multi-FASTA of related sequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-i", "--input", type=Path, help="multi-FASTA of unaligned sequences")
    src.add_argument("--aligned", type=Path, help="pre-aligned FASTA (skips the alignment step)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--aligner", choices=["mafft", "muscle"], default="mafft")
    p.add_argument("--reference-id", help="reference for coordinates (default: first sequence)")

    reg = p.add_mutually_exclusive_group()
    reg.add_argument("--region", action="append", metavar="[NAME:]START-END",
                     help="target region, 1-based inclusive reference coordinates; repeatable "
                          "(default: whole reference)")
    reg.add_argument("--regions-file", type=Path, help="TSV with columns name, start, end")

    p.add_argument("--amplicon-min", type=int, required=True, help="minimum amplicon length (bp)")
    p.add_argument("--amplicon-max", type=int, required=True, help="maximum amplicon length (bp)")
    p.add_argument("--tm-min", type=float, default=57.0)
    p.add_argument("--tm-opt", type=float, default=60.0)
    p.add_argument("--tm-max", type=float, default=63.0)
    p.add_argument("--max-tm-diff", type=float, default=2.0, help="max Tm difference within a pair")
    p.add_argument("--primer-len-min", type=int, default=18)
    p.add_argument("--primer-len-max", type=int, default=25)
    p.add_argument("--gc-min", type=float, default=40.0)
    p.add_argument("--gc-max", type=float, default=60.0)

    p.add_argument("--length-tolerance", type=float, default=0.10,
                   help="allowed deviation from the median length, as a fraction")
    p.add_argument("--outlier-action", choices=["warn", "drop", "abort"], default="drop")

    p.add_argument("--max-mismatches", type=int, default=2, help="max mismatches per primer")
    p.add_argument("--three-prime-window", type=int, default=5,
                   help="3' bases that must match exactly")
    p.add_argument("--min-coverage", type=float, default=95.0, help="minimum coverage (%%)")
    p.add_argument("--allow-degenerate", action="store_true", help="allow IUPAC bases in primers")
    p.add_argument("--max-degeneracy", type=int, default=8, help="max fold per degenerate primer")
    p.add_argument("--candidates", type=int, default=200, help="primer3 pairs to re-score")
    p.add_argument("--num-pairs", type=int, default=5, help="pairs to report per region")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging and tracebacks")
    return p


def check_params(a: argparse.Namespace) -> None:
    problems: list[str] = []
    if not a.tm_min <= a.tm_opt <= a.tm_max:
        problems.append("--tm-min <= --tm-opt <= --tm-max is required")
    if not 1 <= a.primer_len_min <= a.primer_len_max:
        problems.append("1 <= --primer-len-min <= --primer-len-max is required")
    if not 0 <= a.gc_min <= a.gc_max <= 100:
        problems.append("0 <= --gc-min <= --gc-max <= 100 is required")
    if not 1 <= a.amplicon_min <= a.amplicon_max:
        problems.append("1 <= --amplicon-min <= --amplicon-max is required")
    if a.length_tolerance <= 0:
        problems.append("--length-tolerance must be > 0")
    if not 1 <= a.three_prime_window <= a.primer_len_min:
        problems.append("--three-prime-window must be between 1 and --primer-len-min")
    if a.max_mismatches < 0:
        problems.append("--max-mismatches must be >= 0")
    if not 0 <= a.min_coverage <= 100:
        problems.append("--min-coverage must be between 0 and 100")
    if a.max_tm_diff < 0:
        problems.append("--max-tm-diff must be >= 0")
    for name in ("max_degeneracy", "candidates", "num_pairs"):
        if getattr(a, name) < 1:
            problems.append(f"--{name.replace('_', '-')} must be >= 1")
    if problems:
        raise PipelineError("Invalid parameters: " + "; ".join(problems))


def _score_candidate(
    cand: CandidatePair,
    cons: Consensus,
    aln: Alignment,
    sparams: ScoringParams,
    allow_degenerate: bool,
    max_fold: int,
) -> ScoredPair:
    f, r = cand.forward, cand.reverse
    fwd_cols = cons.cons_to_col[f.start: f.start + f.length]
    rev_cols = cons.cons_to_col[r.start: r.start + r.length]
    fwd_site, rev_site = f.seq, reverse_complement(r.seq)
    cov = score_pair(fwd_site, fwd_cols, rev_site, rev_cols, aln, sparams)
    fwd_seq, rev_seq, degenerate = f.seq, r.seq, False
    if allow_degenerate and cov.not_covered:
        new_f, new_r, new_cov = degenerate_pair(
            fwd_site, fwd_cols, rev_site, rev_cols, aln, sparams, max_fold
        )
        if (new_f, new_r) != (fwd_site, rev_site):
            fwd_seq, rev_seq, cov, degenerate = new_f, reverse_complement(new_r), new_cov, True
    return ScoredPair(
        pair=cand, fwd_seq=fwd_seq, rev_seq=rev_seq, fwd_cols=fwd_cols, rev_cols=rev_cols,
        coverage=cov, degenerate=degenerate,
        fwd_tm_range=tm_range(fwd_seq) if degenerate else None,
        rev_tm_range=tm_range(rev_seq) if degenerate else None,
    )


def _design_region(
    region: Region,
    included: tuple[int, int],
    cons: Consensus,
    aln: Alignment,
    a: argparse.Namespace,
    dparams: DesignParams,
    sparams: ScoringParams,
) -> RegionResult:
    try:
        design = run_primer3(cons, included, dparams)
    except PipelineError as exc:
        raise PipelineError(f"Region '{region.name}': {exc}") from exc
    scored = [
        _score_candidate(c, cons, aln, sparams, a.allow_degenerate, a.max_degeneracy)
        for c in design.pairs
    ]
    ranked = rank(scored, a.min_coverage, a.num_pairs)
    n_passing = sum(1 for s in scored if s.coverage.coverage_pct >= a.min_coverage)
    failures = failure_summary(scored) if ranked.below_threshold else []
    if not design.pairs:
        log.warning("Region %s: primer3 found no candidate pairs", region.name)
    elif ranked.below_threshold:
        log.warning(
            "Region %s: no pair reached %g%% coverage (best %.1f%%)",
            region.name, a.min_coverage, ranked.pairs[0].coverage.coverage_pct,
        )
    else:
        log.info(
            "Region %s: %d of %d candidate pairs meet %g%% coverage",
            region.name, n_passing, len(design.pairs), a.min_coverage,
        )
    return RegionResult(region, included, len(design.pairs), n_passing, ranked, design.explain, failures)


@dataclass
class _Inputs:
    records: list[Record]
    lengths: list[int]
    lf: LengthFilterResult
    reference_id: str
    ref_len: int


@dataclass
class _Aligned:
    aln: Alignment
    stats: list[ColumnStats]
    cons: Consensus
    cmap: CoordMap


def _load_inputs(a: argparse.Namespace) -> _Inputs:
    """Load and validate the input, filter by length and pick the reference."""
    if a.input is not None:
        records = load_fasta(a.input)
    else:
        records = load_fasta(a.aligned, allow_gaps=True)
        alignment_from_records(records)  # validates equal row lengths
    if a.reference_id is not None and a.reference_id not in {r.id for r in records}:
        raise PipelineError(f"--reference-id '{a.reference_id}' is not in the input")
    lengths = [ungapped_length(r.seq) for r in records]
    lf = filter_by_length(records, a.length_tolerance, a.outlier_action, a.reference_id)
    # default reference: first sequence left after length filtering
    reference_id = a.reference_id or lf.kept[0].id
    ref_len = ungapped_length(next(r.seq for r in lf.kept if r.id == reference_id))
    return _Inputs(records, lengths, lf, reference_id, ref_len)


def _align_inputs(a: argparse.Namespace, inp: _Inputs) -> _Aligned:
    """Align (or take the given alignment), write it, and build consensus and coordinates."""
    kept = inp.lf.kept
    aln = run_aligner(kept, a.aligner) if a.input is not None else alignment_from_records(kept)
    a.outdir.mkdir(parents=True, exist_ok=True)
    write_alignment(aln, a.outdir / "alignment.fasta")
    stats = score_columns(aln)
    cons = build_consensus(stats)
    return _Aligned(aln, stats, cons, CoordMap(aln, inp.reference_id, cons))


def run(a: argparse.Namespace, command: str, log_lines: list[str]) -> int:
    check_params(a)
    inp = _load_inputs(a)
    regions = resolve_regions(a.region, a.regions_file, inp.ref_len)
    al = _align_inputs(a, inp)
    records, lengths, lf, reference_id = inp.records, inp.lengths, inp.lf, inp.reference_id
    aln, stats, cons, cmap = al.aln, al.stats, al.cons, al.cmap
    # validate every region (incl. its consensus length) before any primer design
    included: dict[str, tuple[int, int]] = {}
    for region in regions:
        cmap.validate_region(region, a.amplicon_min)
        included[region.name] = cmap.region_to_included(region)
        if included[region.name][1] < a.amplicon_min:
            raise PipelineError(
                f"Region '{region.name}' spans only {included[region.name][1]} consensus "
                f"positions (the reference has insertions there), shorter than "
                f"--amplicon-min ({a.amplicon_min})"
            )

    dparams = DesignParams(
        a.tm_min, a.tm_opt, a.tm_max, a.max_tm_diff, a.primer_len_min, a.primer_len_max,
        a.gc_min, a.gc_max, a.amplicon_min, a.amplicon_max, a.candidates,
    )
    sparams = ScoringParams(a.max_mismatches, a.three_prime_window)
    results = [
        _design_region(r, included[r.name], cons, aln, a, dparams, sparams) for r in regions
    ]

    write_tsv(results, cmap, a.outdir / "primers.tsv")
    failed = [res.region.name for res in results if res.rank.below_threshold]
    if failed:
        log.warning("Regions without a pair meeting --min-coverage: %s", ", ".join(failed))
    params = {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(a).items()}
    report = RunReport(
        command=command, params=params, n_input=len(records),
        median_len=float(statistics.median(lengths)), min_len=min(lengths), max_len=max(lengths),
        length=lf, reference_id=reference_id, aln_length=aln.length, cons_length=len(cons.seq),
        n_conserved=sum(1 for s in stats if s.majority_freq == 1.0),
        mean_conservation=statistics.fmean(s.majority_freq for s in stats),
        min_coverage=a.min_coverage, regions=results, log_lines=log_lines,
    )
    write_report(report, a.outdir / "report.txt")
    log.info("Wrote results to %s", a.outdir)
    return EXIT_BELOW_THRESHOLD if failed else EXIT_OK


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)

    fmt = logging.Formatter("%(levelname)s: %(message)s")
    stderr_handler = logging.StreamHandler(sys.stderr)
    memory = _ListHandler()
    for h in (stderr_handler, memory):
        h.setFormatter(fmt)
        log.addHandler(h)
    log.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    log.propagate = False
    try:
        return run(args, "viraprime " + shlex.join(argv), memory.lines)
    except PipelineError as exc:
        log.error("%s", exc)
        if args.verbose:
            log.exception("Traceback:")
        return EXIT_ERROR
    finally:
        log.removeHandler(stderr_handler)
        log.removeHandler(memory)


if __name__ == "__main__":
    sys.exit(main())
