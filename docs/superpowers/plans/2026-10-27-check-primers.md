# viraprime check Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `viraprime check`, which scores primer pairs the user already has against every input sequence with the same rules as design mode.

**Architecture:** A new `viraprime/check.py` loads a TSV of pairs, locates each primer on the consensus (fewest IUPAC-aware mismatches) and calls the existing `scoring.score_pair`. Input loading and alignment move out of `cli.run()` into two shared helpers used by both modes. `main()` dispatches to a separate check parser when `argv[0] == "check"`.

**Tech Stack:** Python ≥ 3.10, Biopython, primer3-py ≥ 2, pytest. Run everything inside the `primers` conda env (`conda activate primers`).

**Spec:** `docs/superpowers/specs/2026-10-27-check-primers-design.md`

## Global Constraints

- Design mode's CLI, outputs and exit codes do not change; the existing test suite stays green after every task.
- Matching rules are exactly `scoring.score_pair` / `match_primer`; check mode never re-implements them.
- Exit codes: `0` all pairs ≥ `--min-coverage`; `1` error (`PipelineError` or argparse); `2` any pair below threshold or not located (outputs still written).
- Primer alphabet `ACGTRYSWKMBDHVN`, case-insensitive, `U` → `T`, stored upper-case. Names match `regions.NAME_RE`.
- `check.tsv` column order exactly as in the spec's Outputs section.
- `not_covered` format `ID(reason)` joined by `,`; unlocated pairs use reason `not located`.
- Note texts: geometry `reverse site is not downstream of forward site`; ambiguity `ambiguous binding site for <fwd|rev> primer`; length `<fwd|rev> primer is longer than the consensus`. Notes join with `"; "`.
- Tm and Tm range are left empty when a primer's fold exceeds `TM_FOLD_LIMIT = 256`, with note `Tm not computed: <fwd|rev> fold <n> > 256`. This guards `degenerate.tm_range` against exponential expansion (extends the spec).
- In check mode the primer file is validated **before** alignment, so a bad file never starts MAFFT/MUSCLE.

## Review Focus

1. Primer files exported from Excel: UTF-8 BOM before the header, `\r\n` line ends, trailing empty columns (`name\tfwd\trev\t\t`). Expect them to parse like a clean file. → Task 1 test.
2. Leading/trailing spaces around fields (copied from a paper). Expect them stripped; internal spaces are an invalid character reported with position. → Task 1 tests.
3. Highly degenerate primers (e.g. 12 × `N`). Expect the run to finish in seconds with Tm empty and the note above, not to hang in `tm_range`. → Task 2 test.
4. A reverse primer pasted as its reverse complement (wrong orientation). Expect the site-mismatch warning to say `… may not belong to this family or may be written in the wrong orientation`. → Task 2 test.
5. A pair whose forward site starts where the reference has a gap. Expect amplicon coordinates snapped like design mode (`CoordMap.col_to_ref`), no crash. → Task 4 test.

---

## File map

- Create `viraprime/check.py`: `load_assays`, `locate`, `check_assay`, `TM_FOLD_LIMIT`.
- Modify `viraprime/models.py`: add `Assay`, `Location`, `CheckResult`, `CheckReport`.
- Modify `viraprime/ranking.py`: extract `count_reasons`.
- Modify `viraprime/io.py`: extract `_header_lines`; add `CHECK_TSV_COLUMNS`, `write_check_tsv`, `write_check_report`.
- Modify `viraprime/cli.py`: `_load_inputs`, `_align_inputs`, `build_check_parser`, `check_check_params`, `run_check`, dispatch in `main`, design epilog.
- Modify `tests/synthetic.py`: `make_aligned_family`.
- Create `tests/test_check.py` (units), `tests/test_check_cli.py` (CLI/e2e).
- Modify `README.md`.

---

### Task 1: Data types and primer-file loader

**Files:**
- Modify: `viraprime/models.py` (append)
- Create: `viraprime/check.py`
- Test: `tests/test_check.py`

**Interfaces:**
- Produces:
  - `@dataclass(frozen=True) class Assay: name: str; forward: str; reverse: str` (5′→3′, upper-case, U→T)
  - `load_assays(path: Path, window: int) -> list[Assay]` in `viraprime/check.py`; raises `PipelineError`.

- [ ] **Step 1: Write the failing tests** in `tests/test_check.py` (helper `_write(tmp_path, data: bytes | str) -> Path`):
  - `test_load_assays_header_comments_blank`: `"name\tforward\treverse\n# c\n\nA1\tacgtacgtacgtacgtac\tTTGGCCAATTGGCCAATT\n"` → `[Assay("A1", "ACGTACGTACGTACGTAC", "TTGGCCAATTGGCCAATT")]`.
  - `test_load_assays_u_to_t`: forward `"ACGUACGUACGUACGUAC"` → `"ACGTACGTACGTACGTAC"`.
  - `test_load_assays_excel_export`: bytes `b"\xef\xbb\xbfname\tforward\treverse\r\nA1\tACGTACGTACGTACGTAC\tTTGGCCAATTGGCCAATT\t\t\r\n"` → one assay `A1`.
  - `test_load_assays_strips_field_spaces`: `" A1 \t ACGTACGTACGTACGTAC \tTTGGCCAATTGGCCAATT\n"` → name `"A1"`, forward without spaces.
  - Error tests, each `pytest.raises(PipelineError, match=...)`:
    - missing file → `"not found"`;
    - header only → `"no primer pairs"`;
    - 2 fields → `"line 1: expected 3 tab-separated columns"`;
    - name `"a b"` → `"invalid name"`;
    - names `A1`, `A1` → `"duplicate.*A1"`;
    - forward `"ACGTXCGTACGTACGTAC"` on a data line after the header → `"forward.*position 5.*'X'"`;
    - internal space `"ACGT ACGTACGTACGTAC"` → `"position 5"`;
    - forward of length 4 with `window=5` → `"shorter than --three-prime-window"`.

- [ ] **Step 2: Run to verify they fail**
  Run: `pytest tests/test_check.py -v` → FAIL (ImportError: `viraprime.check`).

- [ ] **Step 3: Add `Assay` to `models.py` and implement `load_assays`**
  Read with `encoding="utf-8-sig"`; `splitlines()`; strip each field; drop trailing empty fields before counting. Header rule: only the first non-comment line, and only when its forward or reverse field (after U→T, upper-case) contains a character outside `scoring.IUPAC`, is skipped; later bad lines are errors. Names validated with `regions.NAME_RE`. Messages start `"{path.name} line {n}: "`.

- [ ] **Step 4: Run to verify they pass**
  Run: `pytest tests/test_check.py -v` → PASS.

- [ ] **Step 5: Commit**
  `git add viraprime/models.py viraprime/check.py tests/test_check.py && git commit -m "feat(check): load primer pairs from a TSV file"`

---

### Task 2: Locating primers and checking one assay

**Files:**
- Modify: `viraprime/models.py`, `viraprime/check.py`
- Test: `tests/test_check.py`

**Interfaces:**
- Consumes: `Assay`; `Alignment`, `Consensus`, `ScoringParams`, `PairCoverage` (models); `scoring.score_pair`, `scoring.IUPAC`; `degenerate.fold`, `degenerate.tm_range`; `Bio.Seq.reverse_complement`.
- Produces:
  - `@dataclass(frozen=True) class Location: start: int; mismatches: int` (0-based consensus start)
  - `locate(site: str, consensus_seq: str) -> Location | None` (None when ambiguous or the site is longer than the consensus)
  - `@dataclass class CheckResult: assay: Assay; located: bool; note: str; fwd: Location | None; rev: Location | None; fwd_cols: list[int]; rev_cols: list[int]; coverage: PairCoverage; fwd_tm: float | None; rev_tm: float | None; fwd_tm_range: tuple[float, float] | None; rev_tm_range: tuple[float, float] | None; fwd_gc: float; rev_gc: float; degenerate: bool; below_threshold: bool`
  - `check_assay(assay: Assay, aln: Alignment, cons: Consensus, params: ScoringParams, min_coverage: float) -> CheckResult`
  - `TM_FOLD_LIMIT = 256`

Fixture for this task: a hand-built alignment (6 rows × 140 columns from a seeded random base, a few planted SNPs outside the test sites, one all-but-one gap column), consensus via `build_consensus(score_columns(aln))`.

- [ ] **Step 1: Write the failing tests**
  - `test_locate_exact`: site = `cons.seq[30:50]` → `Location(30, 0)`.
  - `test_locate_best_with_mismatches`: same site with one base changed → `Location(30, 1)`.
  - `test_locate_degenerate`: site with `R` where the consensus has `A` → mismatches `0`.
  - `test_locate_ambiguous`: `locate("ACGTAC", "TTACGTACGGACGTACTT")` → `None`.
  - `test_locate_too_long`: `locate("A" * 30, "ACGT")` → `None`.
  - `test_check_assay_matches_score_pair`: forward = `cons.seq[10:30]`, reverse = reverse complement of `cons.seq[80:100]` → `located is True`, `coverage == score_pair(fwd, cons_to_col[10:30], cons.seq[80:100], cons_to_col[80:100], aln, params)`.
  - `test_check_assay_wrong_order`: forward from 80–99, reverse from 10–29 → `located is False`, `note == "reverse site is not downstream of forward site"`, `coverage.coverage_pct == 0.0`, `not_covered == {id: "not located" for id in aln.ids}`.
  - `test_check_assay_overlap`: sites 10–29 and 29–48 → not located, same note.
  - `test_check_assay_degenerate_tm`: forward with one `R` → `degenerate is True`, `fwd_tm_range` set, `fwd_tm == sum(fwd_tm_range) / 2`, `fwd_gc` equals the mean GC of the two expansions.
  - `test_check_assay_high_fold_skips_tm`: forward with 12 `N` → completes in < 5 s, `fwd_tm is None`, `fwd_tm_range is None`, `"Tm not computed: fwd fold 16777216 > 256" in note`.
  - `test_check_assay_wrong_orientation_warns(caplog)`: reverse given as `cons.seq[80:100]` (not reverse-complemented) → a WARNING containing `"wrong orientation"`.
  - `test_check_assay_below_threshold`: one row mutated at the forward primer's 3′-terminal column, `min_coverage=100` → `below_threshold is True`.

- [ ] **Step 2: Run to verify they fail**
  Run: `pytest tests/test_check.py -v -k "locate or check_assay"` → FAIL (ImportError).

- [ ] **Step 3: Implement `Location`, `CheckResult`, `locate`, `check_assay`**
  - `locate`: scan every start; count positions where the consensus base is not in `IUPAC[site_base]`; unique minimum or `None`.
  - Reverse primer located as `reverse_complement(reverse)`; columns `cons.cons_to_col[start:start + len]`.
  - GC % of any primer: `100 * mean over bases of |IUPAC[b] ∩ {G, C}| / |IUPAC[b]|`.
  - Tm: `primer3.calc_tm` when `fold == 1`; `tm_range` and its midpoint when `1 < fold <= TM_FOLD_LIMIT`; otherwise `None` plus note.
  - Warning (logger `"viraprime"`) when a located site has `mismatches > params.max_mismatches`: `"Pair {name}: {fwd|rev} primer matches the consensus with {n} mismatches; it may not belong to this family or may be written in the wrong orientation"`.
  - `below_threshold = coverage.coverage_pct < min_coverage` (always `True` when not located).

- [ ] **Step 4: Run to verify they pass**
  Run: `pytest tests/test_check.py -v` → PASS.

- [ ] **Step 5: Commit**
  `git add viraprime/models.py viraprime/check.py tests/test_check.py && git commit -m "feat(check): locate primers on the consensus and score a pair"`

---

### Task 3: Shared input loading and alignment in the CLI

**Files:**
- Modify: `viraprime/cli.py` (`run`, lines ~202–227)

**Interfaces:**
- Produces, in `cli.py`:
  - `@dataclass class _Inputs: records: list[Record]; lengths: list[int]; lf: LengthFilterResult; reference_id: str; ref_len: int`
  - `_load_inputs(a: argparse.Namespace) -> _Inputs`: load `-i`/`--aligned`, validate `--reference-id`, length filter, default reference = first kept record, `ref_len`. Uses `a.input, a.aligned, a.reference_id, a.length_tolerance, a.outlier_action`.
  - `@dataclass class _Aligned: aln: Alignment; stats: list[ColumnStats]; cons: Consensus; cmap: CoordMap`
  - `_align_inputs(a: argparse.Namespace, inp: _Inputs) -> _Aligned`: align (or build from `--aligned`), create `a.outdir`, write `alignment.fasta`, stats, consensus, `CoordMap`.
- `run()` becomes: `check_params` → `_load_inputs` → `resolve_regions(..., inp.ref_len)` → `_align_inputs` → region validation and design as today.

- [ ] **Step 1: Baseline**
  Run: `pytest -q` → note the pass/skip counts (all green).

- [ ] **Step 2: Refactor `run()` onto the two helpers.** No behaviour change; `RunReport` fields come from `_Inputs` / `_Aligned`.

- [ ] **Step 3: Verify**
  Run: `pytest -q` → identical counts, no failures.

- [ ] **Step 4: Commit**
  `git add viraprime/cli.py && git commit -m "refactor(cli): share input loading and alignment between modes"`

---

### Task 4: Check outputs (TSV and report)

**Files:**
- Modify: `viraprime/models.py`, `viraprime/ranking.py`, `viraprime/io.py`
- Test: `tests/test_check.py`

**Interfaces:**
- Consumes: `CheckResult`, `CoordMap`, `RunReport`.
- Produces:
  - `count_reasons(not_covered: dict[str, str]) -> list[tuple[str, int]]` in `ranking.py`; `failure_summary` delegates to it (behaviour unchanged).
  - `@dataclass class CheckReport` in `models.py`: `command, params, n_input, median_len, min_len, max_len, length, reference_id, aln_length, cons_length, n_conserved, mean_conservation, min_coverage, results: list[CheckResult], log_lines` (same types as `RunReport`).
  - `_header_lines(report: RunReport | CheckReport, title: str) -> list[str]` in `io.py` (title, date, command, parameters, input, alignment blocks); `write_report` uses it with `"viraprime report"`, unchanged output.
  - `CHECK_TSV_COLUMNS: list[str]`, `write_check_tsv(results: list[CheckResult], cmap: CoordMap, path: Path) -> None`, `write_check_report(report: CheckReport, path: Path) -> None` (title `"viraprime check report"`).

- [ ] **Step 1: Write the failing tests**
  - `test_count_reasons`: `{"s1": "fwd 3' mismatch; rev too many mismatches", "s2": "fwd 3' mismatch"}` → `[("fwd 3' mismatch", 2), ("rev too many mismatches", 1)]`.
  - `test_write_check_tsv_columns`: header equals `CHECK_TSV_COLUMNS`; a located result has `amplicon_start == cmap.col_to_ref(fwd_cols[0], "start")`, `amplicon_end == cmap.col_to_ref(rev_cols[-1], "end")`, `amplicon_len == end - start + 1`; an unlocated result has empty amplicon and site-mismatch fields, `located == "False"`, `note` filled.
  - `test_write_check_tsv_reference_gap`: alignment whose reference row has `-` at the forward site's first column → no exception, `amplicon_start` is the next reference base.
  - `test_write_check_report_sections`: contains `"viraprime check report"`, `"Pair A1"`, the consensus positions of a located pair, `"Why it misses sequences"` with counts for a below-threshold pair, `"not located"` for an unlocated pair, and a summary line naming every below-threshold pair (or `"Every pair meets --min-coverage."` when none), followed by the `Log` section with `log_lines`.

- [ ] **Step 2: Run to verify they fail**
  Run: `pytest tests/test_check.py -v -k "count_reasons or write_check"` → FAIL.

- [ ] **Step 3: Implement.** Number formats as `write_tsv`: `:.1f` for Tm, GC and coverage; `_fmt_range` for ranges; booleans as Python `str()`.

- [ ] **Step 4: Verify**
  Run: `pytest -q` → all green (existing report tests included).

- [ ] **Step 5: Commit**
  `git add viraprime tests/test_check.py && git commit -m "feat(check): write check.tsv and report"`

---

### Task 5: `viraprime check` command

**Files:**
- Modify: `viraprime/cli.py`, `tests/synthetic.py`
- Create: `tests/test_check_cli.py`

**Interfaces:**
- Consumes: `_load_inputs`, `_align_inputs`, `load_assays`, `check_assay`, `CheckReport`, `write_check_tsv`, `write_check_report`.
- Produces:
  - `build_check_parser() -> argparse.ArgumentParser` (`_Parser`, `prog="viraprime check"`). Options: `-i/--input` | `--aligned` (exclusive, required), `--primers` (Path, required), `-o/--outdir` (required), `--aligner`, `--reference-id`, `--length-tolerance`, `--outlier-action`, `--max-mismatches`, `--three-prime-window`, `--min-coverage`, `-v/--verbose`; defaults as in design mode.
  - `check_check_params(a: argparse.Namespace) -> None`: the spec's four checks, one `PipelineError` listing all problems.
  - `run_check(a: argparse.Namespace, command: str, log_lines: list[str]) -> int`: params → `_load_inputs` → `load_assays(a.primers, a.three_prime_window)` → `_align_inputs` → `check_assay` per assay → outputs → exit code.
  - `main()`: `argv[:1] == ["check"]` → `build_check_parser().parse_args(argv[1:])`, command `"viraprime check " + shlex.join(argv[1:])`, run `run_check`; logging setup and error handling shared with design mode. Design parser gains `epilog="To test primer pairs you already have: viraprime check --help"`.
  - `make_aligned_family(seed: int = 42) -> list[tuple[str, str]]` in `tests/synthetic.py`: `make_family` as a gapped alignment (3 gap columns in seq03 at base positions 30–32; a 2-column gap after base position 570 in every row except seq05).

- [ ] **Step 1: Write the failing tests** in `tests/test_check_cli.py` (MAFFT-free, all via `--aligned`):
  - `test_make_aligned_family`: equal row lengths; each row with `-` removed equals its `make_family` record.
  - `test_check_agrees_with_design` parametrized over `[]` and `["--allow-degenerate"]`: design `main(["--aligned", aln, "-o", d, "--region", "cons:61-540", "--region", "hv:701-1100", "--amplicon-min", "100", "--amplicon-max", "300", *extra])`; write the rank-1 row of each region as `name\tforward\treverse` (names `cons`, `hv`); `main(["check", "--aligned", aln, "--primers", p, "-o", c])`; per name, `coverage_pct`, `n_covered` and `not_covered` equal design's `primers.tsv`.
  - `test_check_exit_codes`: file with only the `cons` pair → `0`; with `cons` and `hv` → `2` and the report's summary line names `hv`.
  - `test_check_bad_primer_file_exit_1`: invalid character → `1`, `check.tsv` and `alignment.fasta` not written (validated before alignment).
  - `test_check_rejects_design_option`: `main(["check", "--aligned", aln, "--primers", p, "-o", c, "--region", "1-100"])` → `SystemExit` with code `1`.
  - `test_check_bad_params_exit_1`: `--min-coverage 101` → `1`; `--three-prime-window 0` → `1`.
  - `test_check_outputs`: `alignment.fasta`, `check.tsv` (header == `CHECK_TSV_COLUMNS`) and `report.txt` exist.
  - `test_design_help_mentions_check(capsys)`: `main(["--help"])` → `SystemExit(0)`, stdout contains `"viraprime check --help"`.

- [ ] **Step 2: Run to verify they fail**
  Run: `pytest tests/test_check_cli.py -v` → FAIL.

- [ ] **Step 3: Implement** parser, params check, `run_check`, `main` dispatch and epilog; `make_aligned_family`.

- [ ] **Step 4: Verify**
  Run: `pytest -q` → all green.

- [ ] **Step 5: Commit**
  `git add viraprime/cli.py tests/synthetic.py tests/test_check_cli.py && git commit -m "feat(check): viraprime check subcommand"`

---

### Task 6: Documentation and a real run

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Add a "Checking existing primers" section** after "Example": one paragraph on what it does and that it uses the same matching rules; an example `assays.tsv` (header plus the README example's two rank-1 pairs); the command `viraprime check -i examples/family.fasta --primers assays.tsv -o out_check`; the `check.tsv` columns; exit codes 0/1/2; design-only options are rejected.

- [ ] **Step 2: Real run**
  Run: `python examples/make_example.py && viraprime -i examples/family.fasta -o out --region cons:61-540 --region hv:701-1100 --amplicon-min 100 --amplicon-max 300 --allow-degenerate`; write the two rank-1 rows to `assays.tsv`; `viraprime check -i examples/family.fasta --primers assays.tsv -o out_check; echo $?`.
  Expected: exit `2`; `check.tsv` shows `cons` 100.0 and `hv` 30.0, the same `not_covered` as `out/primers.tsv`; `report.txt` has a section per pair.

- [ ] **Step 3: Commit**
  `git add README.md && git commit -m "docs: viraprime check usage"`
