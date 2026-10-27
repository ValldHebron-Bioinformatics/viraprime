"""Shared data types for every pipeline stage."""

from __future__ import annotations

from dataclasses import dataclass


class PipelineError(Exception):
    """A user-facing error: bad input, bad parameters or a failed external tool."""


@dataclass(frozen=True)
class Record:
    id: str
    seq: str


@dataclass
class Alignment:
    ids: list[str]
    rows: list[str]

    @property
    def length(self) -> int:
        return len(self.rows[0]) if self.rows else 0

    def row(self, id: str) -> str:
        try:
            return self.rows[self.ids.index(id)]
        except ValueError:
            raise PipelineError(f"Sequence '{id}' not found in alignment") from None


@dataclass(frozen=True)
class Region:
    name: str
    start: int
    end: int
