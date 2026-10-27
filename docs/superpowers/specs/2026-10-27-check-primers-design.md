# viraprime check — Design

Date: 2026-10-27
Status: approved in brainstorming, pending written-spec review

## Goal

Let a user test primer pairs they already have (published assays, pairs from
an earlier design run, pairs from another tool) against a set of related
sequences, without designing new ones. For each pair the tool reports how
many sequences it amplifies and why it misses the others, using **exactly the
same matching rules as design mode**.

Success: checking the rank-1 pair of a design run reproduces that run's
coverage and `not_covered` list exactly; a multi-pair file is checked in one
run; bad input fails early with a clear message.

Decisions taken in brainstorming:

- Primers come from a **file of pairs** (no single-pair command-line flags).
- **No probe** column for now. The file format leaves room to add one later.
- Binding sites are found on the **consensus** (approach A), then every
  sequence is scored with the existing `scoring.score_pair`. Alignment-free
  per-sequence search (approach B) was rejected because its indel handling
  would not match design mode.

## Usage

```
viraprime check (-i FASTA | --aligned ALN) --primers ASSAYS.tsv -o OUTDIR
    [--aligner mafft|muscle] [--reference-id ID]
    [--length-tolerance 0.10] [--outlier-action drop|warn|abort]
    [--max-mismatches 2] [--three-prime-window 5] [--min-coverage 95]
    [-v]
```

`check` is a subcommand. The existing design invocation
(`viraprime -i … --region … --amplicon-min … --amplicon-max …`) is unchanged:
`main()` dispatches to the check parser only when the first argument is
`check`. The design parser's `--help` gains one epilog line pointing to
`viraprime check --help`. `python -m viraprime check …` works the same way.

Options shared with design mode keep their names, defaults and meaning.
Design-only options (`--region`, `--amplicon-*`, Tm/GC/length limits,
`--allow-degenerate`, `--candidates`, `--num-pairs`, …) are not accepted by
`check`; argparse rejects them with exit 1.

Parameter checks before any work (exit 1): `--max-mismatches ≥ 0`,
`--three-prime-window ≥ 1`, `0 ≤ --min-coverage ≤ 100`,
`--length-tolerance > 0`.

## Primer file (`--primers`)

Tab-separated, three columns: `name`, `forward`, `reverse`.

- `forward` and `reverse` are written 5′→3′, as ordered from a supplier.
  IUPAC codes `ACGTRYSWKMBDHVN` allowed; case-insensitive; `U` → `T`;
  stored upper-case.
- Blank lines and lines starting with `#` are ignored.
- **Header:** the first non-comment line is treated as a header when its
  `forward` or `reverse` field is not a valid oligo (e.g. `name forward reverse`).
- Names follow the region-name rule: letters, digits, `_ . -`.

Errors (exit 1), each naming the file and line:

- file missing, or no pairs in it;
- a line without exactly 3 tab-separated fields;
- empty or invalid name; duplicate names (all duplicates listed);
- a character outside the IUPAC set (field, 1-based position, character);
- a primer shorter than `--three-prime-window`.

## Pipeline

1. **Prepare** (shared with design mode): load and validate the FASTA or
   `--aligned` file, length filter, choose the reference (`--reference-id` or
   first kept sequence), align, write `alignment.fasta`, column stats,
   consensus, `CoordMap`. The code currently inline in `cli.run()` moves into
   one helper used by both modes; design behaviour and outputs do not change.
2. **Load assays** (`check.load_assays`). Validation as above.
3. **Locate** (`check.locate(site, consensus) -> Location | None`). Slide
   `site` along `consensus.seq`; a site base matches a consensus base when the
   consensus base is in the site base's IUPAC set. Return the start with the
   fewest mismatches and that count.
   - The forward primer is located as written; the reverse primer as its
     reverse complement (its binding site on the forward strand).
   - **Ambiguous:** two or more starts share the minimum → not located.
   - Site longer than the consensus → not located.
4. **Pair geometry.** Both sites located, and the reverse site starts after
   the forward site ends (no overlap) → located. Otherwise not located, with
   the reason `reverse site is not downstream of forward site`.
5. **Score.** For a located pair, site columns are
   `cons_to_col[start : start + len]`, then
   `score_pair(fwd_site, fwd_cols, rev_site, rev_cols, aln, ScoringParams(...))`,
   the same call design mode makes. A pair is below threshold when its
   coverage is under `--min-coverage`.
   - A site whose consensus mismatch count exceeds `--max-mismatches` is
     still scored (it will usually cover nothing) and logs a warning that the
     primer may not belong to this family.
   - A pair that is not located gets coverage 0 %, every sequence in
     `not_covered` with reason `not located`, and the reason in `note`.
6. **Primer properties.** Tm with `primer3.calc_tm` and GC % for primers
   without degenerate bases. For degenerate primers, `degenerate.tm_range`
   gives the Tm range, the Tm column holds its midpoint, and GC % is the mean
   over all expansions.

## Outputs (in `--outdir`)

- `alignment.fasta`, as in design mode.
- `check.tsv`, one row per pair in file order. Columns:
  `name, forward_seq, reverse_seq, forward_tm, reverse_tm, forward_gc,
  reverse_gc, forward_len, reverse_len, degenerate, forward_tm_range,
  reverse_tm_range, located, forward_site_mismatches, reverse_site_mismatches,
  amplicon_start, amplicon_end, amplicon_len, coverage_pct, n_covered,
  n_total, below_threshold, not_covered, note`.
  Amplicon coordinates are 1-based on the reference, outer primer edges,
  computed as in `write_tsv`; empty when not located. `not_covered` uses the
  design-mode format `ID(reason)`. Site-mismatch columns count mismatches
  against the consensus.
- `report.txt`: date, command, parameters, input and alignment statistics
  (same header blocks as design mode), then one section per pair: sequences,
  located site (consensus positions and reference coordinates) or the reason
  it was not located, consensus mismatches, coverage, and, when below
  threshold, counts of why sequences are missed (same wording as
  `failure_summary`). Ends with a summary listing pairs below threshold, then
  the log.

## Exit codes

- `0`: every pair reaches `--min-coverage`.
- `1`: error (bad arguments, bad input, aligner missing or failed).
- `2`: finished, but at least one pair is below threshold or not located.
  All outputs are still written.

## Code layout

- `viraprime/check.py` (new): `load_assays`, `locate`, `check_assay`.
- `viraprime/models.py`: `Assay(name, forward, reverse)`,
  `Location(start, mismatches)` and `CheckResult` (assay, located flag,
  note, fwd/rev locations and columns, coverage, Tm/GC values and ranges,
  below_threshold).
- `viraprime/io.py`: `write_check_tsv`, `write_check_report`; the shared
  report header blocks are factored out of `write_report` rather than copied.
- `viraprime/cli.py`: shared `_prepare()`, `build_check_parser()`,
  `run_check()`, `main()` dispatch, design-parser epilog.
- `README.md`: a "Checking existing primers" section with an example file
  and command.

## Testing (`tests/test_check.py`, written test-first)

- `load_assays`: header detected; comments and blank lines skipped;
  lowercase and `U` normalised; each error case above.
- `locate`: exact match; best site with mismatches; degenerate base matches
  either base; ambiguous (repeated site) → not located; site longer than the
  consensus → not located.
- Pair geometry: reverse upstream of forward, and overlapping sites → not
  located with the note.
- **Agreement with design mode:** run design on the synthetic family with
  `--aligned` (no MAFFT needed), once without and once with
  `--allow-degenerate`, feed its rank-1 pair of each region to
  `check`, and assert identical `coverage_pct`, `n_covered` and
  `not_covered`.
- CLI: exit 0 when all pairs pass; exit 2 with one failing pair (named in the
  report); exit 1 for a bad primer file and for a design-only option passed
  to `check`; outputs exist and the TSV header matches.
- The existing suite stays green (design mode unchanged by the refactor).

Done = `pytest` green, plus `viraprime check` run on `examples/family.fasta`
with the README example's pairs and the outputs inspected.

## Out of scope

Probes (TaqMan) and multiplex panels; alignment-free search; per-sequence
amplicon lengths; primer-dimer or hairpin checks; off-target specificity;
checking primers in the web page.
