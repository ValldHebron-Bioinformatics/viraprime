"""Tiling the whole genome with N overlapping amplicons in two pools."""

from __future__ import annotations

import math

from viraprime.models import Geometry

TOLERANCES = (0.20, 0.35, 0.50)  # amplicon length tolerances tried in order
LONG_AMPLICON = 3000  # warn above this target length


def target_length(genome_len: int, n: int, overlap: int) -> int:
    """Amplicon length that tiles genome_len with n amplicons overlapping by overlap."""
    return round((genome_len + (n - 1) * overlap) / n)


def geometry(
    genome_len: int, n: int, overlap: int, tolerance: float, primer_len_max: int
) -> Geometry | None:
    """Size limits for n amplicons, or None when amplicons would be too short.

    Neighbours may overlap (outer primer edges) by at most half the shortest
    amplicon, which keeps amplicons k and k+2 (same pool) apart. That overlap must
    still hold two primers plus the required insert overlap.
    """
    target = target_length(genome_len, n, overlap)
    len_min = math.floor(target * (1 - tolerance))
    max_outer = (len_min - 1) // 2
    if max_outer < overlap + 2 * primer_len_max:
        return None
    return Geometry(
        target_len=target,
        len_min=len_min,
        len_max=math.ceil(target * (1 + tolerance)),
        tolerance=tolerance,
        end_slack=target // 4,
        max_outer_overlap=max_outer,
    )


def max_amplicons(genome_len: int, overlap: int, primer_len_max: int) -> int:
    """Largest n whose amplicons are still long enough at the default tolerance (0 if none)."""
    best = 0
    for n in range(1, genome_len + 1):
        if geometry(genome_len, n, overlap, TOLERANCES[0], primer_len_max) is None:
            break
        best = n
    return best
