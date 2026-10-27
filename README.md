# viraprime

Design PCR primer pairs that work across a family of related sequences (variants of a gene,
strains of a virus, …). The tool aligns the sequences, builds a consensus, asks
[primer3](https://primer3.org) for candidate pairs on the consensus inside one or more target
regions, and then **re-checks every candidate against every individual sequence**. Pairs are
ranked by how many sequences they actually amplify.

## Installation

```bash
conda create -n primers -c conda-forge -c bioconda python=3.11 biopython "primer3-py>=2" pytest mafft muscle
conda activate primers
pip install -e .
```

Python ≥ 3.10 is required. Python dependencies are listed in `requirements.txt`; MAFFT (default)
or MUSCLE 5 must be on `PATH` unless you pass a pre-aligned file with `--aligned`.

## Example

```bash
python examples/make_example.py          # writes examples/family.fasta (10 synthetic sequences)
viraprime -i examples/family.fasta -o out \
    --region cons:61-540 --region hv:701-1100 \
    --amplicon-min 100 --amplicon-max 300 --allow-degenerate
```

The example family has a conserved block (61–540) and a hypervariable block (701–1100). The run
finds pairs with 100 % coverage in `cons`, reports the best available (30 %) pairs for `hv` with
`below_threshold=True`, explains why in `out/report.txt`, and exits with code 2.

`python -m viraprime …` works the same as `viraprime …`.

## Checking existing primers

`viraprime check` tests primer pairs you already have (a published assay, pairs from an earlier
run, pairs from another tool) against your sequences, without designing new ones. Each primer is
placed where it fits the consensus best, then the pair is scored on every sequence with exactly
the same matching rules as design mode.

List the pairs in a tab-separated file, written 5′→3′ as you would order them (IUPAC codes
allowed; a header line, blank lines and `#` comments are optional):

```
name	forward	reverse
cons	CAGTCCCCTACGTCGCATATC	TAGTTGGGCGTCCTGAGAAAC
hv	CTCCAGTGGMGGTCTTGGTKG	GCAYKCATCCMAGTAAGCACAC
```

```bash
viraprime check -i examples/family.fasta --primers assays.tsv -o out_check
```

It writes `alignment.fasta`, `report.txt` (one section per pair, with why it misses sequences)
and `check.tsv`, one row per pair: `name, forward_seq, reverse_seq, forward_tm, reverse_tm,
forward_gc, reverse_gc, forward_len, reverse_len, degenerate, forward_tm_range,
reverse_tm_range, located, forward_site_mismatches, reverse_site_mismatches, amplicon_start,
amplicon_end, amplicon_len, coverage_pct, n_covered, n_total, below_threshold, not_covered,
note`. On the example family the two pairs above give 100 % and 30 % coverage, as in the design
run.

A pair is reported as not located (0 % coverage, reason in `note`) when a primer fits two places
equally well or the reverse primer does not bind downstream of the forward one. Options shared
with design mode (`-i`/`--aligned`, `--aligner`, `--reference-id`, `--length-tolerance`,
`--outlier-action`, `--max-mismatches`, `--three-prime-window`, `--min-coverage`, `-v`) work
the same; design-only options such as `--region` are rejected. Exit codes: `0` every pair meets
`--min-coverage`, `2` at least one does not (outputs are still written), `1` error.

## Pipeline

1. **Load & validate** the FASTA: errors on missing/empty files, empty records, duplicate IDs
   and non-IUPAC characters. Lowercase and `U` are accepted (converted to upper-case / `T`).
2. **Length check**: sequences whose length differs from the median by more than
   `--length-tolerance` × median are dropped, kept with a warning, or abort the run
   (`--outlier-action`).
3. **Align** with MAFFT or MUSCLE, or use `--aligned` to skip this step.
4. **Consensus**: majority base per column; columns that are gaps in ≥ 50 % of sequences are left
   out of the consensus.
5. **Design** candidate pairs with primer3 on the consensus, with the whole amplicon inside each
   region.
6. **Re-score** each pair against every sequence (see *Matching rules*), optionally adding IUPAC
   degenerate bases (`--allow-degenerate`).
7. **Rank** by coverage, then primer3 penalty, and report the top `--num-pairs` per region.

### Matching rules

A primer matches a sequence when it has at most `--max-mismatches` mismatches **and** none in the
last `--three-prime-window` bases of its 3′ end. A pair covers a sequence when both primers match.

- A gap, `N` or any ambiguous base in the sequence counts as a mismatch.
- Each base the sequence has in addition to the primer (an insertion) counts as one mismatch.
- A degenerate primer base matches any base it stands for (`R` matches `A` and `G`).

With `--allow-degenerate`, IUPAC bases are added one at a time, each time choosing the change
that rescues the most sequences (then the one removing the most mismatches, then the one
farthest from the 3′ end), as long as the primer's fold (number of distinct oligos) stays at or
below `--max-degeneracy`. The 3′-terminal base is never made degenerate, and changes that do not
end up rescuing a sequence are undone.

## Coordinates and regions

All coordinates are 1-based, inclusive, on the reference sequence (`--reference-id`, default the
first sequence). They are converted to alignment columns internally.

- `--region START-END` or `--region NAME:START-END`, repeatable. Unnamed regions are called
  `region1`, `region2`, … in command-line order.
- `--regions-file FILE`: tab-separated `name`, `start`, `end`; an optional header line, blank lines
  and `#` comments are allowed.
- Without either, the whole reference is one region called `full`.

Regions are designed independently and may overlap.

## Options

| Option | Default | Description |
|---|---|---|
| `-i/--input` | — | multi-FASTA (exactly one of `-i` / `--aligned` is required) |
| `--aligned` | — | pre-aligned FASTA; skips alignment |
| `-o/--outdir` | required | output directory (created if missing) |
| `--aligner` | `mafft` | `mafft` or `muscle` |
| `--reference-id` | first sequence | reference for coordinates |
| `--region` / `--regions-file` | whole reference | target regions (see above) |
| `--amplicon-min`, `--amplicon-max` | **required** | amplicon length range (bp) |
| `--tm-min/--tm-opt/--tm-max` | 57 / 60 / 63 | primer Tm (°C) |
| `--max-tm-diff` | 2 | max Tm difference within a pair (°C) |
| `--primer-len-min/--primer-len-max` | 18 / 25 | primer length |
| `--gc-min/--gc-max` | 40 / 60 | primer GC % |
| `--length-tolerance` | 0.10 | allowed deviation from median length (fraction) |
| `--outlier-action` | `drop` | `warn`, `drop` or `abort` |
| `--max-mismatches` | 2 | max mismatches per primer |
| `--three-prime-window` | 5 | 3′ bases that must match exactly |
| `--min-coverage` | 95 | minimum coverage (%) |
| `--allow-degenerate` | off | allow IUPAC bases at variable positions |
| `--max-degeneracy` | 8 | max fold per degenerate primer |
| `--candidates` | 200 | primer3 pairs re-scored per region |
| `--num-pairs` | 5 | pairs reported per region |
| `-v/--verbose` | off | debug logging and tracebacks |

## Outputs

- `alignment.fasta` — the alignment used (after outlier removal).
- `primers.tsv` — one row per reported pair: `region, rank, forward_seq, reverse_seq, forward_tm,
  reverse_tm, forward_gc, reverse_gc, forward_len, reverse_len, amplicon_start, amplicon_end,
  amplicon_len, coverage_pct, n_covered, n_total, degenerate, forward_tm_range, reverse_tm_range,
  penalty, below_threshold, not_covered`. Amplicon coordinates are on the reference; `rank`
  restarts at 1 for each region; `not_covered` lists `ID(reason)` for every missed sequence. Tm
  ranges are given for degenerate primers (min–max over all expansions).
- `report.txt` — parameters, input statistics, dropped outliers, alignment statistics, one section
  per region (candidates, coverage, and why the best pair misses sequences when no pair reaches
  `--min-coverage`), a summary, and the full log.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | every region has at least one pair meeting `--min-coverage` |
| 1 | error (bad input or parameters, aligner missing or failed) |
| 2 | finished, but at least one region has no pair meeting `--min-coverage` (or no primer3 candidates); outputs are still written |

## Tests

```bash
pytest
```

Tests use a synthetic family (`tests/synthetic.py`); those that run MAFFT are skipped if it is not
installed.
