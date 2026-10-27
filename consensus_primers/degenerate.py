"""Greedy introduction of IUPAC degenerate bases to rescue uncovered sequences."""

from __future__ import annotations

import itertools
import math

import primer3

from consensus_primers.models import Alignment, PairCoverage, ScoringParams
from consensus_primers.scoring import CODE_FOR, IUPAC, score_pair

BASES = "ACGT"


def fold(site: str) -> int:
    """Number of distinct oligos a degenerate primer stands for."""
    return math.prod(len(IUPAC[b]) for b in site)


def tm_range(primer: str) -> tuple[float, float]:
    """Min and max primer3 Tm over every expansion of a (possibly degenerate) primer."""
    options = [sorted(IUPAC[b]) for b in primer]
    tms = [primer3.calc_tm("".join(v)) for v in itertools.product(*options)]
    return min(tms), max(tms)


def degenerate_pair(
    fwd_site: str,
    fwd_cols: list[int],
    rev_site: str,
    rev_cols: list[int],
    aln: Alignment,
    params: ScoringParams,
    max_fold: int,
) -> tuple[str, str, PairCoverage]:
    """Add one IUPAC base at a time, choosing the change that rescues most sequences.

    Sites are in forward-strand orientation. The 3'-terminal base of each primer (last
    base of the forward site, first base of the reverse site) is never changed.
    """
    sites = {"fwd": fwd_site, "rev": rev_site}
    cols = {"fwd": fwd_cols, "rev": rev_cols}
    rows = dict(zip(aln.ids, aln.rows))
    cov = score_pair(sites["fwd"], fwd_cols, sites["rev"], rev_cols, aln, params)

    while cov.not_covered:
        best: tuple[tuple[int, int, int, int], dict[str, str], PairCoverage] | None = None
        for which in ("fwd", "rev"):
            site = sites[which]
            n = len(site)
            terminal = n - 1 if which == "fwd" else 0
            for i in range(n):
                if i == terminal:
                    continue
                current = IUPAC[site[i]]
                seen = {rows[sid][cols[which][i]] for sid in cov.not_covered}
                for b in sorted((seen & set(BASES)) - current, key=BASES.index):
                    trial_site = site[:i] + CODE_FOR[current | {b}] + site[i + 1:]
                    if fold(trial_site) > max_fold:
                        continue
                    trial = dict(sites, **{which: trial_site})
                    trial_cov = score_pair(trial["fwd"], fwd_cols, trial["rev"], rev_cols, aln, params)
                    gain = len(trial_cov.covered) - len(cov.covered)
                    dist_from_3p = n - 1 - i if which == "fwd" else i
                    key = (gain, dist_from_3p, which == "fwd", -BASES.index(b))
                    if best is None or key > best[0]:
                        best = (key, trial, trial_cov)
        if best is None or best[0][0] <= 0:
            break
        _, sites, cov = best
    return sites["fwd"], sites["rev"], cov
