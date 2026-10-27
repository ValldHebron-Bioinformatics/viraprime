"""Ranking of scored primer pairs: coverage first, then primer3 penalty."""

from __future__ import annotations

from collections import Counter

from consensus_primers.models import RankResult, ScoredPair


def _sort_key(p: ScoredPair) -> tuple[float, float]:
    return (-p.coverage.coverage_pct, p.pair.penalty)


def rank(scored: list[ScoredPair], min_coverage: float, num_pairs: int) -> RankResult:
    """Top num_pairs pairs meeting min_coverage; if none do, the best pairs flagged."""
    ordered = sorted(scored, key=_sort_key)
    passing = [p for p in ordered if p.coverage.coverage_pct >= min_coverage]
    if passing:
        return RankResult(passing[:num_pairs], below_threshold=False)
    return RankResult(ordered[:num_pairs], below_threshold=True)


def failure_summary(scored: list[ScoredPair]) -> list[tuple[str, int]]:
    """Failure reasons of the best pair's uncovered sequences, most common first."""
    if not scored:
        return []
    best = min(scored, key=_sort_key)
    counts: Counter[str] = Counter()
    for reason in best.coverage.not_covered.values():
        counts.update(part.strip() for part in reason.split(";"))
    return counts.most_common()
