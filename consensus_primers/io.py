"""File input/output: FASTA loading/validation and output writers."""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path

from Bio import SeqIO

from consensus_primers.models import PipelineError, Record

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
