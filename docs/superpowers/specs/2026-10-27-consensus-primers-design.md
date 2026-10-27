# consensus-primers — Design

Date: 2026-10-27
Status: approved in brainstorming, pending written-spec review

## Goal

A Python CLI that designs consensus PCR primer pairs from a multi-FASTA of
related sequences (variants of one gene/virus), and reports for each pair how
many of the input sequences it actually amplifies under strict matching rules.

Success: on a synthetic dataset the tool runs end to end, finds pairs in
conserved windows, rejects 3′-end mismatches, and explains clearly when no pair
reaches the coverage threshold.

## Environment

- Python ≥ 3.10 (developed in a conda env `primers`, Python 3.11).
- Dependencies: `biopython`, `primer3-py` ≥ 2.0 (snake_case API), `pytest` (dev).
- External binaries: `mafft` (default) or `muscle` (v5 CLI), from bioconda.
- `requirements.txt` lists Python deps; README documents the conda setup.

## Pipeline

1. **Load & validate** (`io.load_fasta`). Biopython `SeqIO`. Errors on: missing
   or empty file, unparseable FASTA, duplicate IDs (all listed), any character
   outside IUPAC nucleotides `ACGTURYSWKMBDHVN` (case-insensitive; `-` allowed
   only with `--aligned`). Error gives ID, 1-based position and character.
   `U` is converted to `T`. Sequences are upper-cased.
2. **Length check** (`lengths.filter_by_length`). Median of ungapped lengths;
   tolerance `--length-tolerance` (fraction, default 0.10) → keep sequences with
   `|len − median| ≤ tol × median`. `--outlier-action`:
   - `drop` (default): remove, log and report IDs with their lengths.
   - `warn`: keep, log and report.
   - `abort`: error listing the outliers.
   If the reference sequence is an outlier and would be dropped → error.
   If fewer than 2 sequences remain → error.
3. **Align** (`align.run_aligner`). `--aligner {mafft,muscle}` (default mafft).
   `shutil.which` check; missing binary → error with hint
   `conda install -c bioconda <tool>`. MAFFT: `mafft --auto --quiet in.fa`;
   MUSCLE 5: `muscle -align in.fa -output out.fa`. Temp files in a
   `TemporaryDirectory`. Non-zero exit → error with stderr tail.
   `--aligned FILE` skips alignment: the file is validated (equal lengths, gap
   char `-`) and the length filter is applied to ungapped lengths, dropping
   rows from the alignment (all-gap columns then removed).
   Output: `alignment.fasta` (upper-case, original order).
4. **Conservation & consensus** (`conservation`). Per column: counts of
   A/C/G/T, gap count, ambiguous count; `majority_base`, `majority_freq`
   (fraction of all sequences), `gap_frac`. Consensus = majority base of every
   column with `gap_frac < 0.5`; gap-majority columns are omitted. `Consensus`
   holds the string plus `cons_to_col: list[int]`. Ties are broken in A,C,G,T
   order. Column stats are written to the report as a summary (mean
   conservation, number of fully conserved columns).
5. **Coordinates** (`coords.CoordMap`). Built from the reference row
   (`--reference-id`, default first record after filtering):
   - `ref_to_col(pos_1based)` → alignment column.
   - `col_to_ref(col, side)` → 1-based ref position; if the column is a gap in
     the reference, take the nearest ref base to the right (`side="start"`)
     or left (`side="end"`).
   - `col_to_cons(col, side)` analogous for consensus positions.
   Each region `[start, end]` (1-based, inclusive, on the reference) →
   columns → consensus interval → primer3 `SEQUENCE_INCLUDED_REGION` (0-based
   start, length). Adjustments onto gaps are logged. All regions are
   validated before any design work; an invalid region (outside the
   reference, start ≥ end, shorter than `--amplicon-min`) → error.

Steps 1–5 run once per run. Steps 6–9 run **independently for each region**
(see *Regions* below); the amplicon of a pair lies entirely inside its region.

6. **Design** (`design.run_primer3`). `primer3.bindings.design_primers` on the
   consensus with: `PRIMER_MIN/OPT/MAX_TM`, `PRIMER_PAIR_MAX_DIFF_TM`,
   `PRIMER_MIN/MAX_SIZE` (opt = midpoint), `PRIMER_MIN/MAX_GC`,
   `PRIMER_PRODUCT_SIZE_RANGE=[[amin, amax]]`, `PRIMER_NUM_RETURN=--candidates`
   (default 200), `PRIMER_MAX_NS_ACCEPTED=0`. Returns `CandidatePair`s
   (sequences, consensus positions, Tm, GC, penalty). Zero candidates → keep
   primer3's `PRIMER_*_EXPLAIN` strings for the report.
7. **Re-score** (`scoring`). Each primer is mapped to its span of alignment
   columns `[c0, c1]` (from `cons_to_col` of its first and last base).
   For each sequence row and each primer, walk the columns in the span:
   - consensus column (in `cons_to_col`): compare primer base `p` to row char
     `s`. Match iff `s ∈ IUPAC[p]` and `s ∈ {A,C,G,T}`. Gap or any ambiguous
     `s` (including N) = mismatch.
   - non-consensus column (gap-majority) where the row has a base: insertion
     → 1 mismatch.
   The reverse primer is compared as its reverse complement on the forward
   strand; its 3′ end is therefore the leftmost bases of the span.
   A primer matches a sequence iff total mismatches ≤ `--max-mismatches`
   (default 2) **and** 0 mismatches inside the `--three-prime-window`
   (default 5) primer bases at the 3′ end (insertions between those bases
   count as inside the window).
   A pair covers a sequence iff both primers match. Coverage = covered /
   sequences after filtering × 100. Each `ScoredPair` also records, per
   uncovered sequence, the reason (`fwd 3′ mismatch`, `rev too many
   mismatches`, …).
8. **Degeneracy** (`degenerate`, only with `--allow-degenerate`). For each
   primer of a pair with coverage < 100 %: candidate positions are primer
   positions (excluding the 3′-terminal base) where some uncovered sequence
   has a different A/C/G/T base. Greedy loop: pick the position+added base
   that rescues the most sequences (ties → further from the 3′ end), update
   the IUPAC code, stop when no improvement or when the fold (product of
   bases per position) would exceed `--max-degeneracy` (default 8). Re-score
   the pair. Tm reported as primer3 Tm of the consensus primer plus the
   min–max Tm over expanded variants (`primer3.calc_tm`).
9. **Rank** (`ranking.rank`, per region). Keep pairs with coverage ≥
   `--min-coverage` (default 95); sort by coverage desc, then primer3 penalty
   asc; output top `--num-pairs` (default 5). If none pass: output the top
   `--num-pairs` of all candidates with `below_threshold=True`, and the report
   explains the best coverage reached and the most common failure reasons.

## Regions

- `--region SPEC`, repeatable. `SPEC` = `START-END` or `NAME:START-END`
  (1-based inclusive reference coordinates). Unnamed regions are named
  `region1`, `region2`, … in command-line order.
- `--regions-file FILE`: tab-separated `name`, `start`, `end`; an optional
  header line (detected when `start` is not an integer); blank lines and
  `#` comments ignored.
- `--region` and `--regions-file` are mutually exclusive. Neither given →
  one region `full` = the whole reference.
- Names must be unique and match `[A-Za-z0-9_.-]+`. Regions may overlap.
- Parse errors, duplicate names and invalid ranges → error before alignment
  where possible (range-vs-reference checks happen once the reference is
  known, still before any primer design).
- Regions are independent: no cross-region Tm matching, dimer checks or
  overlap constraints (multiplexing is out of scope).
- Parsing lives in `regions.py` (`parse_region_spec`, `load_regions_file`
  → `list[Region]`).

## CLI

| Option | Default | Notes |
|---|---|---|
| `-i/--input` | — | multi-FASTA; exactly one of `-i` / `--aligned` required |
| `--aligned` | — | pre-aligned FASTA; skips the alignment step |
| `-o/--outdir` | required | created if missing |
| `--aligner` | mafft | `mafft`, `muscle` |
| `--reference-id` | first sequence | |
| `--region` | whole reference | repeatable, `[NAME:]START-END`, 1-based inclusive ref coords |
| `--regions-file` | — | TSV `name start end`; exclusive with `--region` |
| `--amplicon-min`, `--amplicon-max` | **required** | bp |
| `--tm-min/--tm-opt/--tm-max` | 57/60/63 | °C |
| `--max-tm-diff` | 2 | °C |
| `--primer-len-min/--primer-len-max` | 18/25 | |
| `--gc-min/--gc-max` | 40/60 | % |
| `--length-tolerance` | 0.10 | fraction of median |
| `--outlier-action` | drop | `warn`, `drop`, `abort` |
| `--max-mismatches` | 2 | per primer |
| `--three-prime-window` | 5 | bases |
| `--min-coverage` | 95 | % |
| `--allow-degenerate` | off | |
| `--max-degeneracy` | 8 | max fold per primer |
| `--candidates` | 200 | primer3 pairs to re-score |
| `--num-pairs` | 5 | |
| `-v/--verbose` | off | DEBUG logging, tracebacks |

Parameter sanity checks (`tm-min ≤ tm-opt ≤ tm-max`, len/gc/amplicon min ≤ max,
tolerance > 0, `three-prime-window ≤ primer-len-min`, counts ≥ 1, coverage in
0–100) → error before any work.

## Outputs (in `--outdir`)

- `alignment.fasta`
- `primers.tsv` columns: `region, rank, forward_seq, reverse_seq, forward_tm,
  reverse_tm, forward_gc, reverse_gc, forward_len, reverse_len,
  amplicon_start, amplicon_end, amplicon_len, coverage_pct, n_covered,
  n_total, degenerate, forward_tm_range, reverse_tm_range, penalty,
  below_threshold, not_covered`. Amplicon start/end are 1-based reference
  coordinates (outer primer edges); `amplicon_len` is measured on the
  reference. `not_covered` = comma-separated `ID(reason)` list, empty if
  none. `rank` restarts at 1 for each region; rows grouped by region in
  input order. A region where primer3 returned nothing contributes no rows
  (header only if that is true for all regions).
- `report.txt`: run date, command line, all parameter values, input stats
  (n sequences, median/min/max length), dropped/flagged outliers, alignment
  stats, warnings, then **one section per region**: coordinate mapping
  (incl. adjustments), number of primer3 candidates, number passing
  coverage, best coverage, and the no-valid-pair explanation when relevant.
  Ends with a summary line listing regions that failed.
  Logging goes to stderr and is also appended to `report.txt`'s log section.

## Exit codes

- `0` success: every region has ≥ 1 pair meeting `--min-coverage`.
- `1` error (`PipelineError` or argument error); message without traceback
  unless `--verbose`.
- `2` completed, but at least one region has no pair meeting
  `--min-coverage` (or primer3 found none there). Failing regions are named
  in the log and report. Outputs are still written for all regions.

## Code layout

```
consensus_primers/
  __init__.py, __main__.py
  cli.py          argparse, logging, orchestration, exit codes
  io.py           load_fasta, load_alignment, write_alignment, write_tsv, write_report
  lengths.py      filter_by_length
  align.py        check_binary, run_aligner
  conservation.py score_columns, build_consensus
  coords.py       CoordMap
  regions.py      Region, parse_region_spec, load_regions_file
  design.py       run_primer3
  scoring.py      match_primer, score_pair
  degenerate.py   degenerate_primer
  ranking.py      rank
  models.py       dataclasses + PipelineError
tests/
  conftest.py     synthetic dataset builder (fixed seed)
  test_io.py, test_lengths.py, test_coords.py, test_regions.py, test_scoring.py,
  test_degenerate.py, test_no_valid_pairs.py, test_cli_e2e.py
pyproject.toml    console script `consensus-primers`
requirements.txt, README.md
```

Stages are pure functions on dataclasses; only `cli.py`/`io.py` touch files
and only `align.py` runs subprocesses.

## Testing

Synthetic family: random 600 bp base (seed fixed) + ~10 variants with SNPs,
small indels and planted mismatches; per-test fixtures add specific features.

- io: empty file, duplicate IDs, invalid character → `PipelineError`.
- lengths: outlier dropped / kept (warn) / abort; IDs reported; reference-as-outlier error.
- coords: ref↔column↔consensus round-trip with gaps in reference and others;
  region edge on a gap; 1-based off-by-one.
- regions: `START-END` and `NAME:START-END` parsing, auto-naming, regions
  file with/without header and comments; errors for bad syntax, start ≥ end,
  duplicate names, `--region` + `--regions-file` together.
- scoring: mismatches outside window accepted up to the limit; one mismatch in
  the 3′ window rejected (forward and reverse); gap/N = mismatch; insertion counted.
- degenerate: SNP → correct IUPAC code; fold cap respected; 3′-terminal base fixed.
- no valid pairs: highly variable dataset → exit 2, `below_threshold`, reason in report.
- CLI e2e: full run with MAFFT on synthetic data (skip if MAFFT absent);
  all outputs exist, TSV header correct, top pair coverage 100 %.
- CLI multi-region: two regions → `region` column populated, ranks restart
  per region, one report section each; one conserved + one hypervariable
  region → exit 2 with the failing region named in the report.

Done = `pytest` green + README example command run on the test data with
outputs inspected.

## Out of scope

Off-target / genome-wide specificity checks, probe (TaqMan) design, primer
dimer screening beyond primer3's own checks, multiplexing, local-search
re-scoring (possible later fallback).
