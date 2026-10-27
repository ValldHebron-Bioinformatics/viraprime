# consensus-primers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A Python CLI that designs consensus PCR primer pairs for one or more reference-coordinate regions of a multi-FASTA, and scores each pair's coverage across all individual sequences.

**Architecture:** Package `consensus_primers/` with one module per pipeline stage, operating on dataclasses from `models.py`. Load/filter/align/consensus/coordinate-map run once; design → re-score → degenerate → rank run per region. Only `io.py`/`cli.py` touch files, only `align.py` runs subprocesses.

**Tech Stack:** Python 3.11 (conda env `primers`), Biopython, primer3-py ≥ 2.0, MAFFT / MUSCLE 5, pytest.

**Spec:** `docs/superpowers/specs/2026-10-27-consensus-primers-design.md` — read it alongside this plan; values cited below come from it.

## Global Constraints

- Python ≥ 3.10; full type hints; `logging` only, never `print` (except argparse's own output).
- primer3-py ≥ 2.0 snake_case API (`primer3.bindings.design_primers`, `primer3.calc_tm`).
- All user-facing coordinates are 1-based inclusive on the reference sequence.
- Defaults (verbatim): tm 57/60/63, max-tm-diff 2, primer-len 18/25, gc 40/60, length-tolerance 0.10, outlier-action drop, aligner mafft, max-mismatches 2, three-prime-window 5, min-coverage 95, max-degeneracy 8, candidates 200, num-pairs 5. `--amplicon-min/--amplicon-max` required, no default.
- Exit codes: 0 all regions OK, 1 error, 2 at least one region below threshold / no candidates.
- User errors raise `PipelineError`; `cli.py` logs the message and exits 1 (traceback only with `-v`).
- Failure-reason strings use ASCII `3'` (e.g. `fwd 3' mismatch`) so TSVs stay plain ASCII.
- Run tests with `conda run -n primers pytest` (or inside the activated env).

## Review Focus

1. Lowercase / soft-masked input, `U` bases, and Windows `\r\n` line endings → accepted and normalised (test in Task 1).
2. A record with an empty sequence among valid ones → `PipelineError` naming the ID, not a silent zero-length row (Task 1).
3. FASTA headers with descriptions (`>seqA some text`) → `--reference-id seqA` matches Biopython's `record.id` (Task 10).
4. Reference has a gap exactly at a primer's outer edge → amplicon start/end snap to the nearest reference base and stay inside the region (Task 4).
5. One region gets zero primer3 candidates while another succeeds → TSV still holds the good region's rows, report explains the empty one, exit 2 (Task 10).

---

### Task 1: Environment, scaffold, models, FASTA loading

**Files:**
- Create: `pyproject.toml`, `requirements.txt`, `.gitignore`, `consensus_primers/__init__.py`, `consensus_primers/models.py`, `consensus_primers/io.py`, `tests/synthetic.py`
- Test: `tests/conftest.py`, `tests/test_io.py`

**Interfaces:**
- Produces (models.py):
  - `class PipelineError(Exception)`
  - `@dataclass(frozen=True) class Record: id: str; seq: str`
  - `@dataclass class Alignment: ids: list[str]; rows: list[str]` with `length -> int` property and `row(id: str) -> str`
  - `@dataclass(frozen=True) class Region: name: str; start: int; end: int`
- Produces (io.py): `load_fasta(path: Path, allow_gaps: bool = False) -> list[Record]`
- Produces (tests/synthetic.py): `make_family(seed: int = 42) -> tuple[str, list[tuple[str, str]]]` (base sequence, records) and `write_fasta(path: Path, records: list[tuple[str, str]]) -> Path`
- Produces (conftest.py): fixtures `family(tmp_path) -> Path` (writes `family.fasta`) and `family_base -> str`; marker `needs_mafft = pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft not installed")`.

- [ ] **Step 1: Create env and scaffold**

```bash
conda create -n primers -y -c conda-forge -c bioconda python=3.11 biopython "primer3-py>=2" pytest mafft muscle
```
`requirements.txt`: `biopython>=1.81`, `primer3-py>=2.0`. `pyproject.toml`: setuptools, `requires-python = ">=3.10"`, console script `consensus-primers = "consensus_primers.cli:main"`, `[tool.pytest.ini_options] testpaths = ["tests"]`, `pythonpath = ["."]`. Then `conda run -n primers pip install -e .`. `.gitignore`: `__pycache__/`, `*.egg-info/`, `.pytest_cache/`, `out*/`, `examples/family.fasta`.

`make_family` uses `random.Random(seed)`: a 1200 bp random ACGT base (spec says ~600; 1200 gives room for one conserved and one hypervariable region); 10 sequences `seq01`..`seq10` (`seq01` = base). Positions 1–600 (1-based): each variant gets 3 SNPs, all inside 1–60 or 541–600 (keeps 61–540 fully conserved). Positions 701–1100: each variant independently substitutes every base with probability 0.25 (hypervariable). `seq03` has a 3 bp deletion at 30–32 and `seq05` a 2 bp insertion after 570.

- [ ] **Step 2: Write failing tests in `tests/test_io.py`**

```python
def test_empty_file_raises(tmp_path): (tmp_path/"e.fa").write_text(""); with pytest.raises(PipelineError, match="empty"): load_fasta(tmp_path/"e.fa")
def test_missing_file_raises(tmp_path): ... match="not found"
def test_duplicate_ids_listed(tmp_path): records a, b, a, b → match="Duplicate IDs: a, b"
def test_invalid_char_reports_id_and_position(tmp_path): ("s1", "ACGTXA") → match=r"s1.*position 5.*'X'"
def test_gap_rejected_unless_allowed(tmp_path): ("s1","AC-GT") raises; allow_gaps=True returns seq "AC-GT"
def test_empty_record_raises(tmp_path): ("s1","ACGT"),("s2","") → match="s2.*empty"
def test_normalises_case_u_and_crlf(tmp_path): write b">s1\r\nacgu\r\n" → [Record("s1","ACGT")]
def test_family_fixture_loads(family): len(load_fasta(family)) == 10
```

- [ ] **Step 3: Run** `conda run -n primers pytest tests/test_io.py -v` → FAIL (ImportError).

- [ ] **Step 4: Implement `models.py` dataclasses and `load_fasta`**

Parse with `SeqIO.parse(handle, "fasta")`; upper-case, `U→T`, strip whitespace/`\r`. Valid alphabet `ACGTRYSWKMBDHVN` (+ `-` when `allow_gaps`). Duplicate IDs listed sorted, comma-separated.

- [ ] **Step 5: Run tests** → PASS. **Step 6: Commit** `feat: scaffold package, models and FASTA loading`.

---

### Task 2: Length filter

**Files:** Create `consensus_primers/lengths.py`; Test `tests/test_lengths.py`

**Interfaces:**
- Consumes: `Record`, `PipelineError`
- Produces: `@dataclass class LengthFilterResult: kept: list[Record]; outliers: list[tuple[str, int]]; median: float; action: str` (in models.py) and `filter_by_length(records: list[Record], tolerance: float, action: Literal["warn","drop","abort"], reference_id: str) -> LengthFilterResult`. Length = count of non-`-` characters (so it works on aligned rows too).

- [ ] **Step 1: Failing tests**

```python
# 4 seqs of 100 bp + "long" of 150 bp, tolerance 0.10 → median 100
def test_drop_removes_outlier(): r = filter_by_length(recs, 0.10, "drop", "s1"); assert [x.id for x in r.kept] == ["s1","s2","s3","s4"]; assert r.outliers == [("long",150)]
def test_warn_keeps_outlier(): assert len(r.kept) == 5 and r.outliers == [("long",150)]
def test_abort_raises(): with pytest.raises(PipelineError, match="long"): ...
def test_boundary_is_inclusive(): a 110 bp seq is kept (|110-100| == 0.10*100)
def test_reference_outlier_dropped_raises(): reference_id="long", action="drop" → PipelineError match="reference"
def test_fewer_than_two_left_raises(): ...
```

- [ ] **Step 2: Run** → FAIL. **Step 3: Implement** (`statistics.median`). Log each outlier at WARNING. **Step 4: Run** → PASS. **Step 5: Commit** `feat: length outlier filter`.

---

### Task 3: Alignment

**Files:** Create `consensus_primers/align.py`; Modify `consensus_primers/io.py`; Test `tests/test_align.py`

**Interfaces:**
- Consumes: `Record`, `Alignment`, `load_fasta(allow_gaps=True)`
- Produces:
  - `check_binary(tool: Literal["mafft","muscle"]) -> str` (path; else `PipelineError` containing `conda install -c bioconda <tool>`)
  - `run_aligner(records: list[Record], tool: Literal["mafft","muscle"]) -> Alignment` (rows upper-case, input order restored by ID)
  - io: `load_alignment(path: Path) -> Alignment` (equal-length check → `PipelineError` naming the first mismatching ID)
  - io: `alignment_from_records(records: list[Record]) -> Alignment` (drops all-gap columns; used after length-filtering an `--aligned` file)
  - io: `write_alignment(aln: Alignment, path: Path) -> None`

Commands: `mafft --auto --quiet in.fa` (stdout = alignment); `muscle -align in.fa -output out.fa`. Use `subprocess.run(..., capture_output=True, text=True, check=False)` in a `TemporaryDirectory`; non-zero → `PipelineError` with the last 20 stderr lines.

- [ ] **Step 1: Failing tests**

```python
def test_missing_binary_hint(monkeypatch): monkeypatch.setattr(shutil,"which",lambda _: None); with pytest.raises(PipelineError, match="conda install -c bioconda mafft"): check_binary("mafft")
@needs_mafft
def test_mafft_aligns_indel(): recs = [Record("a","ACGTACGTTTGCA"), Record("b","ACGTACTTGCA")] → aln.ids == ["a","b"], equal row lengths, "-" in aln.row("b")
def test_load_alignment_unequal_lengths(tmp_path): match="b"
def test_alignment_from_records_drops_all_gap_columns(): rows "A-C","A-G" → "AC","AG"
```

- [ ] **Step 2: Run** → FAIL. **Step 3: Implement.** **Step 4: Run** → PASS. **Step 5: Commit** `feat: MAFFT/MUSCLE alignment and alignment IO`.

---

### Task 4: Conservation, consensus and coordinate map

**Files:** Create `consensus_primers/conservation.py`, `consensus_primers/coords.py`; Test `tests/test_conservation.py`, `tests/test_coords.py`

**Interfaces:**
- Consumes: `Alignment`, `Region`, `PipelineError`
- Produces (models.py): `@dataclass(frozen=True) class ColumnStats: majority_base: str; majority_freq: float; gap_frac: float; n_ambiguous: int` and `@dataclass class Consensus: seq: str; cons_to_col: list[int]`
- Produces: `score_columns(aln: Alignment) -> list[ColumnStats]`, `build_consensus(stats: list[ColumnStats]) -> Consensus`
- Produces: `class CoordMap` with `__init__(self, aln: Alignment, reference_id: str, consensus: Consensus)`, `ref_length: int`, `ref_to_col(pos: int) -> int`, `col_to_ref(col: int, side: Literal["start","end"]) -> int`, `col_to_cons(col: int, side: Literal["start","end"]) -> int`, `region_to_included(region: Region) -> tuple[int, int]` (primer3 0-based start, length; logs INFO when an edge was snapped), `validate_region(region: Region, amplicon_min: int) -> None`.

Consensus rule: column kept iff `gap_frac < 0.5`; majority among A/C/G/T counts, ties in A,C,G,T order; `majority_freq` = count / n_sequences. `col_to_*` with `side="start"` snaps right to the next non-gap/kept column, `"end"` snaps left; no such column → `PipelineError`.

- [ ] **Step 1: Failing tests**

```python
# alignment: ref "AC--GTA", s2 "ACTTGTA", s3 "AC--GAA"
def test_gap_majority_columns_dropped(): cons.seq == "ACGTA"; cons.cons_to_col == [0,1,4,5,6]
def test_majority_tie_order(): two rows "A","C" → majority_base "A"
def test_ref_to_col_skips_gaps(): ref_to_col(1)==0; ref_to_col(3)==4  # 1-based ref pos 3 = 'G' at column 4
def test_col_to_ref_snaps_on_gap(): col_to_ref(2,"start")==3; col_to_ref(3,"end")==2
def test_region_to_included(): Region("r",2,4) → (1, 3)
def test_region_outside_reference_raises(): Region("r",1,99) → match="outside"
def test_region_shorter_than_amplicon_min_raises(): ...
def test_ref_gap_at_primer_edge(): reference with "-" in the column of a primer's first base → col_to_ref(...,"start") returns the next ref base, which is within the region
```

- [ ] **Step 2: Run** → FAIL. **Step 3: Implement.** **Step 4: Run** → PASS. **Step 5: Commit** `feat: conservation scoring, consensus and coordinate mapping`.

---

### Task 5: Region parsing

**Files:** Create `consensus_primers/regions.py`; Test `tests/test_regions.py`

**Interfaces:**
- Consumes: `Region`, `PipelineError`
- Produces: `parse_region_spec(spec: str, index: int) -> Region` (index 1-based for auto-name `region{index}`), `load_regions_file(path: Path) -> list[Region]`, `resolve_regions(specs: list[str] | None, file: Path | None, ref_length: int) -> list[Region]` (neither → `[Region("full", 1, ref_length)]`; both → `PipelineError`; checks unique names, name regex `[A-Za-z0-9_.-]+`, `start < end`, `start >= 1`).

- [ ] **Step 1: Failing tests**

```python
def test_plain_spec(): parse_region_spec("120-480", 1) == Region("region1",120,480)
def test_named_spec(): parse_region_spec("E1:900-1350", 2) == Region("E1",900,1350)
@pytest.mark.parametrize("bad", ["abc","10-","5-5","9-3","x y:1-9","0-10"]) def test_bad_spec(bad): raises PipelineError
def test_file_with_header_and_comments(tmp_path): "name\tstart\tend\n# c\n\nA\t1\t50\nB\t40\t90\n" → [Region("A",1,50), Region("B",40,90)]
def test_file_without_header(tmp_path): ...
def test_duplicate_names_raise(): resolve_regions(["A:1-10","A:5-20"], None, 100) match="A"
def test_both_sources_raise(tmp_path): ...
def test_default_full_region(): resolve_regions(None, None, 600) == [Region("full",1,600)]
```

- [ ] **Step 2: Run** → FAIL. **Step 3: Implement.** **Step 4: Run** → PASS. **Step 5: Commit** `feat: region specs and regions file`.

---

### Task 6: Primer re-scoring

**Files:** Create `consensus_primers/scoring.py`; Test `tests/test_scoring.py`

**Interfaces:**
- Consumes: `Alignment`
- Produces (models.py): `@dataclass(frozen=True) class MatchResult: ok: bool; mismatches: int; three_prime_mismatches: int; reason: str | None` (reason ∈ `"3' mismatch"`, `"too many mismatches"`, None); `@dataclass class PairCoverage: covered: list[str]; not_covered: dict[str, str]; coverage_pct: float`; `@dataclass(frozen=True) class ScoringParams: max_mismatches: int; window: int`
- Produces: `match_primer(site: str, cols: list[int], row: str, three_prime: Literal["right","left"], params: ScoringParams) -> MatchResult` and `score_pair(fwd_site: str, fwd_cols: list[int], rev_site: str, rev_cols: list[int], aln: Alignment, params: ScoringParams) -> PairCoverage`

`site` = the oligo as it reads on the forward strand (forward primer as-is; reverse primer reverse-complemented with `Bio.Seq.reverse_complement`, which handles IUPAC). `cols[i]` = alignment column of site base i (consensus columns, ascending). Rules (spec §7): match iff row char ∈ {A,C,G,T} and ∈ IUPAC set of the site base; gap/ambiguous = mismatch; every column in `range(cols[0], cols[-1]+1)` not in `cols` where the row has a base = 1 insertion mismatch. 3′ window = last `window` site bases for `"right"`, first `window` for `"left"`; an insertion counts as in-window when it lies between two window columns. Prefer `"3' mismatch"` over `"too many mismatches"` when both apply. `not_covered[id]` = `"fwd <reason>"`, `"rev <reason>"` or both joined by `"; "`. Coverage = len(covered)/n×100.

- [ ] **Step 1: Failing tests** (use a 3-row alignment built by hand; `P = ScoringParams(2, 5)`)

```python
def test_exact_match(): ok, mismatches == 0
def test_two_mismatches_outside_window_ok(): mismatches at site positions 0 and 3 of a 20-mer → ok
def test_three_mismatches_rejected(): reason == "too many mismatches"
def test_single_3prime_mismatch_forward_rejected(): mismatch at position 19, three_prime="right" → ok False, reason "3' mismatch"
def test_single_3prime_mismatch_reverse_rejected(): mismatch at position 0, three_prime="left" → reason "3' mismatch"
def test_same_mismatch_outside_reverse_window_ok(): position 19 with three_prime="left" → ok
def test_gap_and_n_are_mismatches(): row with "-" at pos 2 and "N" at pos 6 → mismatches == 2
def test_insertion_counted(): cols skip one column where row has a base → mismatches == 1
def test_degenerate_site_base_matches(): site "R" at pos 2 matches both A and G rows
def test_score_pair_reasons(): seq with fwd 3' mismatch only → not_covered == {"s3": "fwd 3' mismatch"}; coverage_pct == pytest.approx(66.67, abs=0.01)
```

- [ ] **Step 2: Run** → FAIL. **Step 3: Implement.** **Step 4: Run** → PASS. **Step 5: Commit** `feat: per-sequence primer re-scoring with 3' window`.

---

### Task 7: primer3 design

**Files:** Create `consensus_primers/design.py`; Test `tests/test_design.py`

**Interfaces:**
- Consumes: `Consensus`
- Produces (models.py): `@dataclass(frozen=True) class Primer: seq: str; start: int; length: int; tm: float; gc: float` (`start` = 0-based leftmost consensus position, for both primers); `@dataclass(frozen=True) class CandidatePair: forward: Primer; reverse: Primer; penalty: float`; `@dataclass(frozen=True) class DesignParams` with fields tm_min, tm_opt, tm_max, max_tm_diff, len_min, len_max, gc_min, gc_max, amplicon_min, amplicon_max, candidates; `@dataclass class DesignResult: pairs: list[CandidatePair]; explain: dict[str, str]`
- Produces: `run_primer3(consensus: Consensus, included: tuple[int, int], params: DesignParams) -> DesignResult`

primer3 args per spec §6, `PRIMER_OPT_SIZE = (len_min+len_max)//2`. Right primer: primer3 reports `[three_prime_pos, len]` → `start = pos - len + 1`. `explain` = `PRIMER_LEFT_EXPLAIN`, `PRIMER_RIGHT_EXPLAIN`, `PRIMER_PAIR_EXPLAIN`.

- [ ] **Step 1: Failing tests**

```python
def test_pairs_inside_included_region(family_base): cons = Consensus(family_base, list(range(1200))); r = run_primer3(cons, (100, 400), params(100,300)); assert r.pairs; for p in r.pairs: assert p.forward.start >= 100 and p.reverse.start + p.reverse.length <= 500
def test_reverse_start_is_leftmost(family_base): consensus[rev.start: rev.start+rev.length] == reverse_complement(rev.seq)
def test_no_candidates_returns_explain(): Consensus("A"*300, ...) → pairs == [] and "PRIMER_LEFT_EXPLAIN" in explain
```

- [ ] **Step 2: Run** → FAIL. **Step 3: Implement.** **Step 4: Run** → PASS. **Step 5: Commit** `feat: primer3 candidate design on consensus`.

---

### Task 8: Degenerate primers

**Files:** Create `consensus_primers/degenerate.py`; Test `tests/test_degenerate.py`

**Interfaces:**
- Consumes: `score_pair`, `ScoringParams`, `PairCoverage`, `Alignment`
- Produces: `fold(site: str) -> int`, `degenerate_pair(fwd_site: str, fwd_cols: list[int], rev_site: str, rev_cols: list[int], aln: Alignment, params: ScoringParams, max_fold: int) -> tuple[str, str, PairCoverage]` (returns possibly-degenerate sites, forward-strand orientation), `tm_range(primer: str) -> tuple[float, float]` (min/max `primer3.calc_tm` over all expansions).

Greedy loop (spec §8): each iteration, for both sites, every position except the 3′-terminal one (last base of fwd site, first base of rev site), every A/C/G/T base seen at that column in a currently uncovered sequence and not yet in the code: trial-score the pair; keep the trial with the largest increase in covered count if `fold ≤ max_fold`. Ties: position farther from its 3′ end, then forward before reverse, then A,C,G,T. Stop when the best increase is 0. IUPAC code from the base set via a frozenset→code dict.

- [ ] **Step 1: Failing tests**

```python
def test_snp_becomes_iupac(): 4 rows, two with A, two with G at site pos 3 → fwd site has "R" at 3, coverage 100
def test_fold_cap_respected(): variation at 4 positions, max_fold=4 → fold(site) <= 4 and coverage < 100
def test_three_prime_terminal_never_degenerate(): variation only at last fwd base → site unchanged, coverage unchanged
def test_fold(): fold("ARYN") == 16
def test_tm_range_orders(): lo <= hi for "ACGTRACGTACGTACGTAC"
```

- [ ] **Step 2: Run** → FAIL. **Step 3: Implement.** **Step 4: Run** → PASS. **Step 5: Commit** `feat: greedy degenerate primers with fold cap`.

---

### Task 9: Ranking

**Files:** Create `consensus_primers/ranking.py`; Test `tests/test_ranking.py`

**Interfaces:**
- Consumes: `CandidatePair`, `PairCoverage`
- Produces (models.py): `@dataclass class ScoredPair: pair: CandidatePair; fwd_seq: str; rev_seq: str; fwd_cols: list[int]; rev_cols: list[int]; coverage: PairCoverage; degenerate: bool; fwd_tm_range: tuple[float, float] | None; rev_tm_range: tuple[float, float] | None` (`fwd_seq`/`rev_seq` in 5′→3′ primer orientation); `@dataclass class RankResult: pairs: list[ScoredPair]; below_threshold: bool`
- Produces: `rank(scored: list[ScoredPair], min_coverage: float, num_pairs: int) -> RankResult`; `failure_summary(scored: list[ScoredPair]) -> list[tuple[str, int]]` (reason counts across the best pair's not_covered, most common first).

- [ ] **Step 1: Failing tests**

```python
def test_coverage_then_penalty(): (cov 100, pen 2), (cov 100, pen 1), (cov 96, pen 0) → order pen1, pen2, cov96
def test_filters_below_min(): one at 90 dropped when another ≥ 95 exists
def test_none_pass_returns_best_flagged(): all < 95 → below_threshold True, len == min(num_pairs, n), sorted by coverage
def test_truncates_to_num_pairs(): ...
def test_failure_summary_counts(): best pair not_covered {"a":"fwd 3' mismatch","b":"fwd 3' mismatch","c":"rev too many mismatches"} → [("fwd 3' mismatch",2),("rev too many mismatches",1)]
```

- [ ] **Step 2: Run** → FAIL. **Step 3: Implement.** **Step 4: Run** → PASS. **Step 5: Commit** `feat: rank pairs by coverage then penalty`.

---

### Task 10: Outputs and CLI orchestration

**Files:** Create `consensus_primers/cli.py`, `consensus_primers/__main__.py`; Modify `consensus_primers/io.py`; Test `tests/test_cli.py`, `tests/test_no_valid_pairs.py`, `tests/test_cli_e2e.py`

**Interfaces:**
- Consumes: everything above.
- Produces (models.py): `@dataclass class RegionResult: region: Region; included: tuple[int, int]; n_candidates: int; n_passing: int; rank: RankResult; explain: dict[str, str]; failures: list[tuple[str, int]]` and `@dataclass class RunReport` (command, params dict, n_input, median/min/max length, length result, alignment length, consensus length, n fully conserved columns, mean majority_freq, `regions: list[RegionResult]`, `log_lines: list[str]`).
- Produces (io.py): `write_tsv(results: list[RegionResult], cmap: CoordMap, path: Path) -> None`, `write_report(report: RunReport, path: Path) -> None`.
- Produces (cli.py): `build_parser() -> argparse.ArgumentParser`, `main(argv: list[str] | None = None) -> int`.

TSV columns and order exactly as spec "Outputs" (`region, rank, …, not_covered`); Tm/GC 1 decimal, coverage 1 decimal, tm_range `lo-hi` or empty, `not_covered` = `ID(reason)` joined by `,`. Amplicon start = `cmap.col_to_ref(fwd_cols[0], "start")`, end = `cmap.col_to_ref(rev_cols[-1], "end")`. Logging: root handler to stderr plus an in-memory handler whose formatted lines go into `RunReport.log_lines`. `main` flow: parse → param sanity checks (spec CLI section) → load (`-i` → `load_fasta`, or `--aligned` → `load_fasta(allow_gaps=True)`) → `filter_by_length` → align or `alignment_from_records` → `write_alignment` → consensus → `CoordMap` → `resolve_regions` + `validate_region` for all → per region: `run_primer3`, score each candidate, degenerate if enabled and coverage < 100, `rank` → write TSV + report → return 0/2. `PipelineError` → log ERROR, return 1. The report's no-valid-pair text for a region includes `best coverage <x>%` and the `failures` list; for zero candidates it includes primer3's explain strings.

- [ ] **Step 1: Failing tests**

```python
# test_cli.py
def test_amplicon_args_required(): with pytest.raises(SystemExit): main(["-i","x.fa","-o","out"])
def test_input_and_aligned_exclusive(): ...SystemExit
def test_bad_tm_order_returns_1(family, tmp_path): main([... "--tm-min","62","--tm-opt","60"]) == 1
@needs_mafft
def test_reference_id_with_description(tmp_path): headers ">seqA some text" etc. → --reference-id seqA → exit != 1
# test_no_valid_pairs.py (needs_mafft)
def test_hypervariable_region_exit_2(family, tmp_path): main(["-i",family,"-o",out,"--region","hv:701-1100","--amplicon-min","100","--amplicon-max","300"]) == 2; every TSV row below_threshold == "True"; report contains "hv" and "best coverage"
def test_zero_candidates_region_explained(tmp_path): make_family records with positions 1-400 replaced by "A"*400 → --region polyA:1-400 --region ok:401-600 → exit 2, TSV rows only for "ok", report mentions primer3 explain text for polyA
# test_cli_e2e.py (needs_mafft)
def test_single_region_end_to_end(family, tmp_path): --region 100-500, amplicon 100-300 → exit 0; alignment.fasta, primers.tsv, report.txt exist; header == spec column list; first row coverage_pct == "100.0"
def test_two_regions(family, tmp_path): --region A:61-300 --region B:301-540 → region column {"A","B"}, ranks restart at 1, report has one section per region
def test_mixed_regions_exit_2(family, tmp_path): --region ok:100-500 --region hv:701-1100 → exit 2; report summary names "hv" only
def test_outlier_dropped_reported(family, tmp_path): append a 300 bp record "short" → report lists "short"; "short" not in alignment.fasta
def test_aligned_input_skips_aligner(monkeypatch, family, tmp_path): run once, then feed alignment.fasta via --aligned with shutil.which patched to return None → exit 0
```

For `test_zero_candidates_region_explained`, the "ok" region must lie in the conserved, non-poly-A part: use base positions 401–600 as conserved (build this dataset explicitly in the test with `make_family` then overwrite 1–400 with A in every record) and `--region ok:401-600`.

- [ ] **Step 2: Run** → FAIL. **Step 3: Implement writers, parser and `main`.** **Step 4: Run full suite** `conda run -n primers pytest -v` → all PASS. **Step 5: Commit** `feat: CLI orchestration, TSV and report outputs`.

---

### Task 11: README and end-to-end verification

**Files:** Create `README.md`, `examples/make_example.py` (writes `make_family()` output to `examples/family.fasta` via `tests/synthetic.py`)

- [ ] **Step 1: Write README** — purpose, conda install line from Task 1, example command, options table (from spec), output descriptions, exit codes, matching rules summary.

Example command:
```bash
consensus-primers -i examples/family.fasta -o out --region cons:61-540 --region hv:701-1100 --amplicon-min 100 --amplicon-max 300 --allow-degenerate
```

- [ ] **Step 2: Run** `conda run -n primers python examples/make_example.py` then the example command. Expected: exit code 2 (hv region fails by design), `out/primers.tsv` has 100 % rows for `cons` and `below_threshold` rows for `hv`, `out/report.txt` explains `hv`. Inspect all three files.
- [ ] **Step 3: Run** with only `--region cons:61-540` → exit 0.
- [ ] **Step 4: Run full test suite** → all PASS.
- [ ] **Step 5: Commit** `docs: README and example dataset`.
