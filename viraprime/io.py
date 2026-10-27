"""File input/output: FASTA loading/validation and output writers."""

from __future__ import annotations

import csv
import logging
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

from Bio import SeqIO

from viraprime.models import Alignment, PipelineError, Record, RegionResult, RunReport

if TYPE_CHECKING:
    from viraprime.coords import CoordMap

log = logging.getLogger(__name__)

IUPAC_NT = set("ACGTRYSWKMBDHVN")


def load_fasta(path: Path, allow_gaps: bool = False) -> list[Record]:
    """Load a multi-FASTA, upper-case it, convert U->T and validate it."""
    path = Path(path)
    if not path.is_file():
        raise PipelineError(f"Input file not found: {path}")
    text = path.read_text().replace("\r", "")
    if not text.strip():
        raise PipelineError(f"Input file is empty: {path}")
    if not text.lstrip().startswith(">"):
        raise PipelineError(f"Input file is not valid FASTA (no '>' header first): {path}")

    allowed = IUPAC_NT | {"-"} if allow_gaps else IUPAC_NT
    records: list[Record] = []
    with path.open() as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            seq = "".join(str(rec.seq).split()).upper().replace("U", "T")
            if not seq:
                raise PipelineError(f"Sequence '{rec.id}' is empty")
            for i, ch in enumerate(seq, start=1):
                if ch not in allowed:
                    raise PipelineError(
                        f"Sequence '{rec.id}' has an invalid character at position {i}: '{ch}'"
                    )
            records.append(Record(rec.id, seq))

    if not records:
        raise PipelineError(f"Input file is empty: {path}")
    dupes = sorted(rid for rid, n in Counter(r.id for r in records).items() if n > 1)
    if dupes:
        raise PipelineError(f"Duplicate IDs: {', '.join(dupes)}")
    log.info("Loaded %d sequences from %s", len(records), path)
    return records


def alignment_from_records(records: list[Record]) -> Alignment:
    """Build an Alignment from equal-length gapped records, dropping all-gap columns."""
    if not records:
        raise PipelineError("Alignment is empty")
    width = len(records[0].seq)
    for r in records:
        if len(r.seq) != width:
            raise PipelineError(
                f"Aligned sequences differ in length: '{r.id}' has {len(r.seq)} columns, "
                f"'{records[0].id}' has {width}"
            )
    keep = [c for c in range(width) if any(r.seq[c] != "-" for r in records)]
    rows = ["".join(r.seq[c] for c in keep) for r in records]
    return Alignment(ids=[r.id for r in records], rows=rows)


def load_alignment(path: Path) -> Alignment:
    """Load and validate a pre-aligned FASTA (gap character '-')."""
    return alignment_from_records(load_fasta(path, allow_gaps=True))


def write_alignment(aln: Alignment, path: Path) -> None:
    Path(path).write_text("".join(f">{i}\n{row}\n" for i, row in zip(aln.ids, aln.rows)))


TSV_COLUMNS = [
    "region", "rank", "forward_seq", "reverse_seq", "forward_tm", "reverse_tm", "forward_gc",
    "reverse_gc", "forward_len", "reverse_len", "amplicon_start", "amplicon_end", "amplicon_len",
    "coverage_pct", "n_covered", "n_total", "degenerate", "forward_tm_range", "reverse_tm_range",
    "penalty", "below_threshold", "not_covered",
]


def _fmt_range(r: tuple[float, float] | None) -> str:
    return f"{r[0]:.1f}-{r[1]:.1f}" if r else ""


def write_tsv(results: list[RegionResult], cmap: CoordMap, path: Path) -> None:
    """Write ranked pairs of every region; rank restarts at 1 per region."""
    with Path(path).open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(TSV_COLUMNS)
        for res in results:
            for rank, sp in enumerate(res.rank.pairs, start=1):
                start = cmap.col_to_ref(sp.fwd_cols[0], "start")
                end = cmap.col_to_ref(sp.rev_cols[-1], "end")
                cov = sp.coverage
                writer.writerow([
                    res.region.name, rank, sp.fwd_seq, sp.rev_seq,
                    f"{sp.pair.forward.tm:.1f}", f"{sp.pair.reverse.tm:.1f}",
                    f"{sp.pair.forward.gc:.1f}", f"{sp.pair.reverse.gc:.1f}",
                    sp.pair.forward.length, sp.pair.reverse.length,
                    start, end, end - start + 1,
                    f"{cov.coverage_pct:.1f}", len(cov.covered), len(cov.covered) + len(cov.not_covered),
                    sp.degenerate, _fmt_range(sp.fwd_tm_range), _fmt_range(sp.rev_tm_range),
                    f"{sp.pair.penalty:.3f}", res.rank.below_threshold,
                    ",".join(f"{sid}({why})" for sid, why in cov.not_covered.items()),
                ])


def _region_section(res: RegionResult, min_coverage: float) -> list[str]:
    r = res.region
    title = f"Region {r.name} ({r.start}-{r.end})"
    start, length = res.included
    lines = [title, "-" * len(title),
             f"Consensus positions (primer3 included region): {start + 1}-{start + length}",
             f"primer3 candidate pairs: {res.n_candidates}",
             f"Pairs meeting --min-coverage ({min_coverage:g}%): {res.n_passing}"]
    if res.n_candidates == 0:
        lines.append("primer3 found no candidate pairs in this region. primer3 explanation:")
        lines += [f"  {k}: {v}" for k, v in res.explain.items()]
        return lines
    best = res.rank.pairs[0].coverage.coverage_pct
    lines.append(f"Best coverage: {best:.1f}%")
    if res.rank.below_threshold:
        lines.append(
            f"WARNING: no pair reached {min_coverage:g}% coverage; the best available pairs are "
            "reported with below_threshold=True."
        )
        if res.failures:
            lines.append("Why the best pair misses sequences (count of uncovered sequences):")
            lines += [f"  {reason}: {n}" for reason, n in res.failures]
    return lines


def write_report(report: RunReport, path: Path) -> None:
    lf = report.length
    out = ["viraprime report", "=" * 24,
           f"Date: {datetime.now().isoformat(timespec='seconds')}",
           f"Command: {report.command}", "",
           "Parameters", "----------"]
    out += [f"{k}: {v}" for k, v in report.params.items()]
    out += ["", "Input", "-----",
            f"Sequences loaded: {report.n_input}",
            f"Length median/min/max: {report.median_len:.0f} / {report.min_len} / {report.max_len} bp"]
    if lf.outliers:
        verb = "dropped" if lf.action == "drop" else "kept (warn)"
        out.append(f"Length outliers ({verb}):")
        out += [f"  {sid}: {n} bp" for sid, n in lf.outliers]
    else:
        out.append("Length outliers: none")
    out += [f"Sequences used: {len(lf.kept)}", f"Reference: {report.reference_id}", "",
            "Alignment", "---------",
            f"Alignment columns: {report.aln_length}",
            f"Consensus length: {report.cons_length}",
            f"Fully conserved columns: {report.n_conserved}",
            f"Mean majority-base frequency: {report.mean_conservation:.3f}", ""]
    for res in report.regions:
        out += _region_section(res, report.min_coverage) + [""]
    failed = [res.region.name for res in report.regions if res.rank.below_threshold]
    out += ["Summary", "-------"]
    if failed:
        out.append(
            f"Regions without a pair meeting --min-coverage ({report.min_coverage:g}%): "
            + ", ".join(failed)
        )
    else:
        out.append("Every region has at least one pair meeting --min-coverage.")
    out += ["", "Log", "---", *report.log_lines]
    Path(path).write_text("\n".join(out) + "\n")
