"""File input/output: FASTA loading/validation and output writers."""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path

from Bio import SeqIO

from consensus_primers.models import Alignment, PipelineError, Record

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
