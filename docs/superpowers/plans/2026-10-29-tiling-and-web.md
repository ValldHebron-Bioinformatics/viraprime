# Tiling mode and web app Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `viraprime -i genomes.fasta -n N` tiles the whole genome with N overlapping amplicons in 2 pools (global best fit), and the same runs from a hosted FastAPI web app.

**Architecture:** New engine modules `viraprime/candidates.py` (enumerate, filter, conservation-pick and score single primers) and `viraprime/tiling.py` (geometry, dynamic programming, pools, dimers). The CLI gains `-n`, a hidden `--progress-file`, and `viraprime web`. The web app (`viraprime/web/`) runs each job as `cli.main(argv)` in a child process and serves a static page. A Dockerfile and a GitHub Action deploy it to a Hugging Face Space.

**Tech Stack:** Python 3.11, primer3-py, biopython, FastAPI, uvicorn, python-multipart, httpx (tests), vanilla HTML/CSS/JS, Docker, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-10-29-tiling-and-web-design.md`

## Global Constraints

- Every command that works today behaves exactly as before; the existing suite (150 tests) stays green.
- Tiling defaults: overlap 50, degenerate **on**, max degeneracy 4, min coverage 95, Tm 57/60/63, max Tm diff 2, primer length 18–25, GC 40–60, max mismatches 2, 3′ window 5. Region-mode defaults unchanged (degenerate off, max 8).
- Same matching rules everywhere: final amplicon coverage is computed with `scoring.score_pair`.
- Exit codes: 0 all amplicons ≥ min coverage, 2 some below (outputs still written), 1 error.
- Primer names `viraprime_<k>_LEFT` / `viraprime_<k>_RIGHT`; odd amplicons pool 1, even pool 2.
- Web limits (env-configurable): upload ≤ 20 MB, ≤ 500 sequences, ≤ 35 000 bp each, N ≤ 200, 2 workers, 15 min timeout, files deleted 24 h after finishing.
- Web page follows `DESIGN.md` (untracked, in the repo root) and the landing page's look; taste lint and web-design-guidelines review must pass.
- Tests: `PATH=/home/vhir/anaconda3/envs/primers/bin:$PATH python -m pytest -q`. MAFFT tests skip if MAFFT is missing.

## Plan decisions (deviations from the spec's wording, same intent)

1. **Candidate primers are enumerated directly** (every 3′ anchor position, best length by penalty) with primer3's thermodynamic functions (`calc_tm`, `calc_hairpin`, `calc_homodimer`), instead of asking `design_primers` per window. Windowed `design_primers` ignores conservation and leaves gaps the DP needs filled.
2. **`degenerate_pair` is not rewritten**; a new `degenerate_site` follows the same greedy rule for one primer. Region-mode results therefore cannot change.
3. **"Same-pool amplicons don't touch" is enforced by a pairwise rule**: the outer overlap of neighbours is at most `(len_min − 1) // 2`, which guarantees amplicon k+2 starts after amplicon k ends. Consequence: the smallest workable amplicon is about 250 bp with defaults, so "N too large" comes from `geometry()` feasibility (not a flat 150 bp).
4. **The web app calls `cli.main(argv)` in a child process** rather than a separate library entry point: identical code path, isolated logging, killable on timeout. Results JSON is built from the output files.
5. Insert overlap (`overlap_next`) is counted in consensus positions.

## Review Focus

1. A stretch of the consensus with no acceptable primers (gap-majority columns, N, low complexity) → exit 1 naming the stretch in reference coordinates and suggesting N values that work, never a traceback.
2. `-n 1` → one amplicon spanning the genome, pool 1 only, no `overlap_next`, no dimer pairs across pools.
3. Reference with insertions or a partial reference → amplicon coordinates snap or are left empty with a note, as in check mode; BED rows skip unmappable primers with a warning.
4. Two web jobs at once and a job that times out → logs don't mix, the timed-out job is `failed` with a clear message, the other finishes.
5. Web download names and job ids from the URL → only whitelisted file names and uuid-hex ids are served (no path traversal).

---

### Task 1: Synthetic genome family, models and tiling geometry

**Files:**
- Modify: `tests/synthetic.py`, `viraprime/models.py`
- Create: `viraprime/tiling.py`, `tests/test_tiling_geometry.py`

**Interfaces:**
- Produces:
  - `tests.synthetic.make_genome_family(length=6000, n_seqs=12, seed=7, hyper=(3000, 3300)) -> list[tuple[str, str]]` — substitutions only (all rows the same length, usable directly with `--aligned`); scattered SNPs at rate 0.003 per base per sequence; every base in `hyper` (0-based half-open) substituted with probability 0.25 in each non-first sequence; first record is the unmodified base `g01`.
  - `models.TilingParams` (frozen dataclass): `n_amplicons: int, overlap: int, tm_min: float, tm_opt: float, tm_max: float, max_tm_diff: float, len_min: int, len_max: int, gc_min: float, gc_max: float, max_mismatches: int, window: int, degenerate: bool, max_degeneracy: int, min_coverage: float`; property `scoring -> ScoringParams`.
  - `models.Geometry` (frozen): `target_len: int, len_min: int, len_max: int, tolerance: float, end_slack: int, max_outer_overlap: int`.
  - `tiling.target_length(genome_len: int, n: int, overlap: int) -> int` = `round((genome_len + (n - 1) * overlap) / n)`.
  - `tiling.geometry(genome_len, n, overlap, tolerance, primer_len_max) -> Geometry | None` — `len_min = floor(L*(1-tol))`, `len_max = ceil(L*(1+tol))`, `end_slack = L // 4`, `max_outer_overlap = (len_min - 1) // 2`; returns `None` when `max_outer_overlap < overlap + 2 * primer_len_max` (infeasible).
  - `tiling.max_amplicons(genome_len, overlap, primer_len_max) -> int` — largest n ≥ 1 with `geometry(..., tolerance=0.20, ...)` not None (0 if none).
  - `tiling.TOLERANCES = (0.20, 0.35, 0.50)`; `tiling.LONG_AMPLICON = 3000`.

- [ ] **Step 1: Write the failing tests** in `tests/test_tiling_geometry.py`:

```python
def test_make_genome_family_shape():
    rows = make_genome_family()
    assert len(rows) == 12 and {len(s) for _, s in rows} == {6000}
    assert rows[0][0] == "g01" and set("".join(s for _, s in rows)) <= set("ACGT")

def test_target_length():
    assert target_length(6000, 5, 50) == 1240
    assert target_length(1000, 1, 50) == 1000

def test_geometry_values():
    g = geometry(6000, 5, 50, 0.20, 25)
    assert (g.target_len, g.len_min, g.len_max, g.end_slack) == (1240, 992, 1488, 310)
    assert g.max_outer_overlap == 495

def test_geometry_infeasible_when_amplicons_too_short():
    assert geometry(6000, 40, 50, 0.20, 25) is None  # L = 199

def test_max_amplicons():
    n = max_amplicons(6000, 50, 25)
    assert geometry(6000, n, 50, 0.20, 25) is not None
    assert geometry(6000, n + 1, 50, 0.20, 25) is None
```

- [ ] **Step 2: Run** `python -m pytest tests/test_tiling_geometry.py -q` — Expected: FAIL (ImportError).
- [ ] **Step 3: Implement** `make_genome_family`, the two dataclasses, and the geometry functions.
- [ ] **Step 4: Run** the same command — Expected: 5 passed. Then the full suite — Expected: all pass.
- [ ] **Step 5: Commit** `feat(tiling): genome test family and tiling geometry`.

---

### Task 2: Candidate primers (enumerate, conservation pick, score, degenerate)

**Files:**
- Create: `viraprime/candidates.py`, `tests/test_candidates.py`
- Modify: `viraprime/degenerate.py`, `viraprime/models.py`

**Interfaces:**
- Consumes: `TilingParams`, `scoring.match_primer`, `scoring.IUPAC`, `degenerate.fold`, `degenerate.tm_range`, `Alignment`, `Consensus`, `ColumnStats`.
- Produces:
  - `models.PrimerCandidate` (dataclass): `strand: Literal["+", "-"]`, `start: int` (0-based leftmost consensus position of the site), `length: int`, `site: str` (forward-strand orientation, possibly degenerate), `base_seq: str` (5′→3′ as ordered, non-degenerate), `seq: str` (5′→3′ as ordered, possibly degenerate), `tm: float` (midpoint if degenerate), `tm_range: tuple[float, float] | None`, `gc: float`, `penalty: float`, `mask: int` (bit i set ⇔ `aln.ids[i]` matched), `cols: list[int]`; property `end = start + length - 1`.
  - `candidates.enumerate_primers(cons_seq: str, params: TilingParams) -> tuple[list[PrimerCandidate], list[PrimerCandidate]]` — forward list anchored at each 3′ end e (site `cons[e-l+1 : e+1]`), reverse list anchored at each site start s (site `cons[s : s+l]`, seq = reverse complement). For each anchor, keep the length l in `[len_min, len_max]` with Tm and GC in range and the lowest penalty `|tm - tm_opt| + |l - (len_min+len_max)//2|`; skip sites with non-ACGT; reject if `calc_hairpin(seq).tm > 47.0`, `calc_homodimer(seq).tm > 47.0`, or more than 3 G/C in the last 5 bases of `seq`. `mask` 0 and `cols` empty at this stage.
  - `candidates.pick_conserved(cands, stats, cons, params, bin_size=40, per_bin=3) -> list[PrimerCandidate]` — conservation score = Σ over site columns of `(1 - majority_freq)`, weight 3 inside the 3′ window (right end for "+", left end for "−"); per bin of 3′-anchor positions keep the `per_bin` lowest `(score, penalty)`.
  - `candidates.score_primers(cands, aln, cons, params, progress=None) -> list[PrimerCandidate]` — sets `cols = cons.cons_to_col[start:start+length]`; scores each row with `match_primer` (three_prime "right" for "+", "left" for "−"), memoising by the row slice `row[cols[0]:cols[-1]+1]` (identical slices ⇒ identical result); if `params.degenerate` and some rows miss, replaces `site`/`seq`/`tm`/`tm_range` with `degenerate_site`'s result; sets `mask`. Calls `progress("scoring", i / len(cands))` every 50 primers when given.
  - `degenerate.degenerate_site(site, cols, aln, params: ScoringParams, max_fold, three_prime: ThreePrime) -> str` — same greedy rule as `degenerate_pair` for a single site: candidate changes at every position except the 3′-terminal base, bases seen at that column in currently unmatched rows; key `(rows rescued, mismatches removed over unmatched rows, distance from 3′ end, -BASES.index(b))`; stop when neither rescued nor removed is positive; return the state with most rows matched (earliest on ties); never exceed `max_fold`.

- [ ] **Step 1: Write the failing tests** in `tests/test_candidates.py`:
  - `test_enumerate_respects_limits`: on `make_genome_family()[0][1]` with default `TilingParams` (n=5), every candidate has `18 <= length <= 25`, `57 <= tm <= 63`, `40 <= gc <= 60`; reverse `seq == reverse_complement(site)`; forward `seq == site`.
  - `test_enumerate_skips_non_acgt`: a consensus with `N` at position 100 yields no candidate whose site covers 100.
  - `test_pick_conserved_avoids_variable_block`: with `make_genome_family()` as `alignment_from_records`, stats and consensus, picked candidates inside the hyper block (3000–3300) are fewer than 10 % of picked candidates elsewhere per kb; at most 3 per 40-bp bin per strand.
  - `test_score_matches_match_primer`: for 20 picked candidates, `mask` bits equal `match_primer(...).ok` per row computed directly.
  - `test_degenerate_site_rescues_and_keeps_3prime`: a site with one variable non-terminal column where half the rows carry G instead of A becomes `R` there, fold 2, all rows matched; the 3′-terminal base is never changed even when it is the variable column.
  - `test_degenerate_site_respects_max_fold`: with `max_fold=1` the site is returned unchanged.
  - `test_region_mode_unchanged` is implicit: the existing suite must stay green.
- [ ] **Step 2: Run** `python -m pytest tests/test_candidates.py -q` — Expected: FAIL (ImportError).
- [ ] **Step 3: Implement** the model, `degenerate_site`, and the three candidate functions.
- [ ] **Step 4: Run** the file — Expected: all pass; full suite — Expected: all pass.
- [ ] **Step 5: Commit** `feat(tiling): candidate primer enumeration, conservation pick and scoring`.

---

### Task 3: Best-fit tiling, pools, dimers, ends

**Files:**
- Modify: `viraprime/tiling.py`, `viraprime/models.py`
- Create: `tests/test_tiling.py`

**Interfaces:**
- Consumes: Task 1 geometry, Task 2 `PrimerCandidate` lists (scored), `score_pair`, `CoordMap`.
- Produces:
  - `models.Amplicon` (dataclass): `index: int` (1-based), `pool: int`, `fwd: PrimerCandidate`, `rev: PrimerCandidate`, `coverage: PairCoverage`, `below_threshold: bool`, `overlap_next: int | None`.
  - `models.TilingResult` (dataclass): `amplicons: list[Amplicon]`, `geometry: Geometry`, `n_fwd: int`, `n_rev: int`, `dimer_warnings: list[str]`, `uncovered_ends: tuple[int, int]` (reference bases outside the first/last insert).
  - `tiling.choose(fwds, revs, n, geom, params, genome_len) -> list[tuple[PrimerCandidate, PrimerCandidate]] | None`.
  - `tiling.tile(aln, cons, stats, cmap, params, progress=None) -> TilingResult` — raises `PipelineError` with the no-fit message.
  - `tiling.DIMER_DG = -9000.0` (cal/mol).

Algorithm for `choose` (not determined by the signature):

```
amplicon (f, r): r.end - f.start + 1 in [len_min, len_max], |f.tm - r.tm| <= max_tm_diff,
                 cov = popcount(f.mask & r.mask), pen = f.penalty + r.penalty
neighbours a -> b valid iff   a.r.start - b.f.end - 1 >= overlap          (insert overlap)
                          and a.r.end - b.f.start + 1 <= max_outer_overlap
first amplicon: f.start <= end_slack;  last: r.end >= genome_len - 1 - end_slack
feasible(t): DP over layers k = 1..n with state = reverse primer r ending amplicon k
   V_1(r)   = best over f (first-amplicon rule, cov >= t) of (cov, -pen)
   W_k(f)   = best V_k(r') over r' valid as predecessor of f (bisect revs by start)
   V_k+1(r) = best over amplicons (f, r) with cov >= t and W_k(f) set of W_k(f) + (cov, -pen)
   answer   = best V_n(r) with the last-amplicon rule; keep back-pointers
values compare lexicographically (sum of coverage, -sum of penalty)
binary-search t over the distinct cov values (feasibility is monotone); return the
path for the largest feasible t, or None
```

`tile` steps: progress `("finding primers", 0.1)` → enumerate, pick, score → `("choosing amplicons", 0.85)` → for tolerance in `TOLERANCES`: `geometry(...)`; skip None; `choose(...)`; stop at the first success and log `"used ±35 % amplicon length tolerance"` when not the first. On failure: find the longest consensus stretch with no candidate on either strand, convert to reference coordinates with `cmap` and raise `PipelineError(f"No tiling with {n} amplicons: no acceptable primers between {a:,} and {b:,}; try -n {list}")`, where the list holds the values in `n-3..n+3` (1..`max_amplicons`) for which `choose` succeeds at tolerance 0.50 with t=0 (or `"a different -n"` if none). If L > `LONG_AMPLICON`, log a warning. Final coverage of each amplicon = `score_pair(f.site, f.cols, r.site, r.cols, aln, params.scoring)`; pools `1 if index odd else 2`; dimers: for every unordered pair of primers in the same pool (fwd and rev of all amplicons in that pool), `primer3.calc_heterodimer(x.base_seq, y.base_seq).dg < DIMER_DG` → warning `"pool 1: viraprime_3_LEFT × viraprime_5_RIGHT ΔG -10.2 kcal/mol"`. `uncovered_ends` from the first insert start and the last insert end, via `cmap.col_to_ref`.

- [ ] **Step 1: Write the failing tests** in `tests/test_tiling.py` (fixture: `make_genome_family()` → alignment, stats, consensus, `CoordMap(aln, "g01", cons)`, default params with n=5):
  - `test_exact_count_pools_and_order`: 5 amplicons, indices 1..5, pools 1,2,1,2,1, starts strictly increasing.
  - `test_overlap_and_same_pool_spacing`: for each neighbour pair insert overlap ≥ 50 and `overlap_next` equals it; for k, k+2: `a[k+2].fwd.start > a[k].rev.end`.
  - `test_reaches_both_ends`: first `fwd.start <= end_slack`, last `rev.end >= 6000 - 1 - end_slack`.
  - `test_avoids_hypervariable_block`: no primer site overlaps 3000–3299 and every amplicon has coverage 100 %.
  - `test_mask_agrees_with_score_pair`: `popcount(f.mask & r.mask) == len(coverage.covered)` for each amplicon.
  - `test_n1_single_amplicon`: `-n 1` on a 1500-bp family: one amplicon, pool 1, `overlap_next is None`.
  - `test_widening_reported`: a geometry where ±20 % fails but ±35 % works (construct a family of length 3000 with a 500-bp stretch of `N`-free but primer-less low-complexity sequence `"AT" * 250` placed so the ±20 % windows cannot avoid it) → result `geometry.tolerance == 0.35` or `0.50`, and the log mentions it (`caplog`).
  - `test_no_fit_message`: a consensus whose middle 3000 bp is `"A" * 3000` with n=3 → `PipelineError` whose message contains `"no acceptable primers between"` and `"try -n"`.
  - `test_dimer_warning_format`: monkeypatch `primer3.calc_heterodimer` to return an object with `dg = -12000` → warnings name pool and primers.
- [ ] **Step 2: Run** `python -m pytest tests/test_tiling.py -q` — Expected: FAIL.
- [ ] **Step 3: Implement** `choose`, `tile`, the two models.
- [ ] **Step 4: Run** the file — Expected: all pass; full suite — Expected: all pass.
- [ ] **Step 5: Commit** `feat(tiling): best-fit dynamic programming, pools, dimer checks`.

---

### Task 4: Tiling outputs (TSV, BED, report)

**Files:**
- Modify: `viraprime/io.py`, `viraprime/models.py`
- Create: `tests/test_tiling_io.py`

**Interfaces:**
- Consumes: `TilingResult`, `CoordMap`, `_header_lines`, `count_reasons`, `_fmt_range`.
- Produces:
  - `io.TILING_TSV_COLUMNS` = `amplicon, pool, forward_name, forward_seq, reverse_name, reverse_seq, forward_tm, reverse_tm, forward_gc, reverse_gc, forward_len, reverse_len, degenerate, forward_tm_range, reverse_tm_range, amplicon_start, amplicon_end, amplicon_len, overlap_next, coverage_pct, n_covered, n_total, below_threshold, not_covered` (spec order).
  - `io.write_tiling_tsv(result, cmap, path)`; amplicon coordinates as `write_tsv` (outer primer edges via `col_to_ref`), empty when mapping raises `PipelineError` or start > end.
  - `io.write_scheme_bed(result, cmap, path)` — 7 tab-separated columns `chrom start end name pool strand sequence`, 0-based half-open on the reference, `chrom = cmap.reference_id`, LEFT rows `+`, RIGHT rows `-`, pool `1`/`2`; rows in amplicon order, LEFT before RIGHT; an unmappable primer is skipped with `log.warning`.
  - `models.TilingReport` — same header fields as `RunReport` plus `result: TilingResult`.
  - `io.write_tiling_report(report, cmap, path)` — `_header_lines(report, "viraprime tiling report")`, then a `Tiling` block (`Amplicons: N`, `Target length: L bp (allowed a–b, ±20 %)`, `Candidate primers: x forward, y reverse`, `Not covered at the ends: a bp (start), b bp (end)`), one line per amplicon (`  1 pool 1  101-1340  1240 bp  100.0 %`), for each amplicon below threshold its `count_reasons` lines, a `Primer dimers` block (`none` or the warnings), `Summary` (`All amplicons meet X % coverage` or `Amplicons below X % coverage: 3, 7`), then `Log`.

- [ ] **Step 1: Write the failing tests**: TSV header equals `TILING_TSV_COLUMNS`; names `viraprime_1_LEFT`; last row's `overlap_next` empty; BED has 2·N lines of 7 fields, `start < end`, `end - start == primer length`, strands alternate `+ -`; report contains `"Target length:"`, `"Primer dimers"`, and the below-threshold summary when a result has a low amplicon (build the `TilingResult` by hand with a 50 % amplicon).
- [ ] **Step 2: Run** `python -m pytest tests/test_tiling_io.py -q` — Expected: FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** — Expected: pass; full suite — Expected: pass.
- [ ] **Step 5: Commit** `feat(tiling): primers.tsv, scheme.bed and report`.

---

### Task 5: CLI `-n`, progress file, example data, README

**Files:**
- Modify: `viraprime/cli.py`, `examples/make_example.py`, `README.md`, `pyproject.toml` (version 0.2.0)
- Create: `tests/test_tiling_cli.py`, `examples/genomes.fasta` (generated), `viraprime/web/static/example.fasta` (generated; `viraprime/web/__init__.py` empty)

**Interfaces:**
- Consumes: `tiling.tile`, Task 4 writers, `_load_inputs`, `_align_inputs`.
- Produces:
  - Parser: `-n/--amplicons int`, `--overlap int (50)`, `--no-degenerate`, `--max-degeneracy` default `None` (→ 8 region mode, 4 tiling), `--amplicon-min/--amplicon-max` no longer `required`, hidden `--progress-file Path` (`help=argparse.SUPPRESS`) on both parsers.
  - `cli.parse_args(argv) -> tuple[argparse.Namespace, Callable[..., int], str]` (runner and command string; used by `main`).
  - `cli.validate_args(a: argparse.Namespace) -> None` — mode check + the mode's param checks; raises `PipelineError`. Rules: `-n` with any of `--region`, `--regions-file`, `--amplicon-min`, `--amplicon-max`, `--allow-degenerate` → `"-n cannot be combined with <flags>"`; `--no-degenerate` without `-n` → error; neither `-n` nor both amplicon bounds → `"use -n N to tile the whole genome with N amplicons, or --amplicon-min/--amplicon-max (and optionally --region) for region design"`; `-n < 1` → error; tiling checks reuse the Tm/length/GC/mismatch/window/coverage rules of `check_params` plus `--overlap >= 0`, `--max-degeneracy >= 1`.
  - `cli.make_progress(path: Path | None) -> Callable[[str, float], None]` — writes `{"stage": s, "progress": p}` atomically (temp file + `os.replace`); no-op for `None`.
  - `cli.run_tiling(a, command, log_lines) -> int` — progress `("aligning", 0.05)` before `_align_inputs`; `tile(...)`; writers; `("done", 1.0)`. `run` and `run_check` also emit `aligning`, `scoring`, `done`.
  - `cli.main(argv)` additionally dispatches `argv[:1] == ["web"]` to `viraprime.web.server.main(argv[1:])` (added in Task 6; until then `"viraprime web needs: pip install viraprime[web]"` on ImportError, exit 1).

- [ ] **Step 1: Write the failing tests** in `tests/test_tiling_cli.py` (write `make_genome_family()` to `aln.fasta`, use `--aligned`):
  - `test_tiling_run_outputs`: `main(["--aligned", aln, "-n", "5", "-o", out]) == 0`; `primers.tsv` has 5 rows; `scheme.bed` 10 lines; `report.txt` and `alignment.fasta` exist.
  - `test_tiling_agrees_with_check`: write each amplicon's `forward_seq`/`reverse_seq` to a primer file, run `check`; per amplicon `coverage_pct`, `n_covered`, `not_covered` identical.
  - `test_n_with_region_rejected` (exit 1, message names `--region`), `test_no_mode_message` (exit 1, message contains `use -n N`), `test_no_degenerate_without_n_rejected`.
  - `test_region_mode_still_full_reference`: the existing `--amplicon-min 100 --amplicon-max 300` run without `--region` still writes region `full`.
  - `test_exit_2_below_threshold`: `--min-coverage 100` on `make_genome_family(hyper=(0, 6000))` (everything variable) with `--no-degenerate` → exit 2, outputs written.
  - `test_progress_file`: with `--progress-file p.json`, the file ends as `{"stage": "done", "progress": 1.0}`.
  - `test_n_too_large`: `-n 40` on the 6000-bp family → exit 1, stderr contains `use at most`.
- [ ] **Step 2: Run** `python -m pytest tests/test_tiling_cli.py -q` — Expected: FAIL.
- [ ] **Step 3: Implement** the CLI changes. In `tile`/`run_tiling`, the `N too large` error is `f"N too large: amplicons would be ~{L} bp; use at most {max_amplicons(...)}"` when no tolerance gives a geometry.
- [ ] **Step 4: Run** — Expected: pass; full suite — Expected: 150 old + new all pass.
- [ ] **Step 5: Example data.** `examples/make_example.py` also writes `examples/genomes.fasta` and `viraprime/web/static/example.fasta` from `make_genome_family(length=9000, n_seqs=24, seed=11, hyper=(4200, 4500))`. Run it, then `viraprime -i examples/genomes.fasta -n 8 -o <scratchpad>/out` (MAFFT) and read `report.txt` and `primers.tsv`. Expected: exit 0, 8 amplicons, coverage 100 %.
- [ ] **Step 6: Performance check.** Generate `make_genome_family(length=30000, n_seqs=500, seed=3, hyper=(15000, 15400))` in the scratchpad and time `viraprime --aligned … -n 30`. Expected ≤ 5 min. If slower, profile (`python -m cProfile -s cumtime`) and fix the hot spot (scoring memo, bisect windows) before continuing; record the time in the ledger.
- [ ] **Step 7: README.** New top section "Quick start" with `viraprime -i genomes.fasta -n 10 -o out`, what it writes, and a pointer to "Run it from the web" (filled in Task 8); a "Tiling mode" section with the options table for tiling; existing region-design and check sections move under "Advanced: target regions" and keep their content.
- [ ] **Step 8: Commit** `feat(tiling): viraprime -n, progress file, examples, docs`.

---

### Task 6: Web backend (jobs, API, `viraprime web`)

**Files:**
- Create: `viraprime/web/app.py`, `viraprime/web/jobs.py`, `viraprime/web/params.py`, `viraprime/web/server.py`, `tests/test_web_params.py`, `tests/test_web_jobs.py`, `tests/test_web_api.py`
- Modify: `pyproject.toml` (extras `web = ["fastapi>=0.110", "uvicorn>=0.29", "python-multipart>=0.0.9"]`, `dev` adds `httpx`; packages `["viraprime", "viraprime.web"]`; package-data `viraprime.web = ["static/*"]`)

Install first: `PATH=… pip install -e ".[web,dev]"`.

**Interfaces:**
- Consumes: `cli.parse_args`, `cli.validate_args`, `cli.main`, `io.load_fasta`.
- Produces:
  - `params.Limits` (frozen dataclass) with `from_env()`: `max_upload_mb=20, max_seqs=500, max_seq_len=35000, max_amplicons=200, workers=2, timeout_s=900, retention_s=86400` from `VIRAPRIME_*` env vars.
  - `params.build_argv(mode: str, params: dict) -> list[str]` — whitelist per mode (raises `ValueError("unknown setting: x")`):
    - tile: `amplicons→-n`, `overlap`, `tm_min`, `tm_opt`, `tm_max`, `max_tm_diff`, `primer_len_min`, `primer_len_max`, `gc_min`, `gc_max`, `max_mismatches`, `three_prime_window`, `degenerate` (bool; false → `--no-degenerate`), `max_degeneracy`, `min_coverage`, `aligner`, `reference_id`;
    - regions: the shared Tm/len/GC/mismatch/window/coverage/aligner/reference keys plus `regions` (list of `[NAME:]START-END` → repeated `--region`), `amplicon_min`, `amplicon_max`, `allow_degenerate` (bool), `max_degeneracy`;
    - check: `max_mismatches`, `three_prime_window`, `min_coverage`, `aligner`, `reference_id`.
    Values are converted with the declared type (int/float/bool/str); `aligner` must be `mafft|muscle`.
  - `params.validate_upload(path: Path, limits) -> tuple[int, int]` (n sequences, max length) — `load_fasta` errors and limit breaches raise `ValueError` with a user-facing message.
  - `jobs.run_job(job_dir: Path, argv: list[str]) -> None` (child-process target): redirects stdout/stderr to `job_dir/log.txt`, runs `cli.main(argv + ["-o", out, "--progress-file", job_dir/"progress.json"])`, writes `job_dir/result.json` = `{"exit_code", "mode", "rows": [TSV rows as dicts], "reference_id", "reference_length", "error": last "ERROR: " log line or null}`; `rows` from `primers.tsv` or `check.tsv`; reference from the `Reference:` line of `report.txt` and the ungapped length of that row in `alignment.fasta`.
  - `jobs.JobManager(root: Path, limits: Limits, target=run_job)`: `submit(files: dict[str, bytes], argv: list[str]) -> str` (uuid4 hex; writes files into the job dir, queues); `status(job_id) -> dict` (`status` queued|running|done|failed, `stage`, `progress`, `queue_position`, `result`, `error`; `KeyError` if unknown); `tick()` (start queued jobs up to `workers` with `multiprocessing.get_context("spawn").Process`; mark finished; terminate and mark failed after `timeout_s` with `"took longer than 15 minutes"`); `cleanup(now: float | None = None) -> int` (delete job dirs whose finish time is older than `retention_s`); `start()`/`stop()` run `tick` every 0.5 s and `cleanup` every 10 min in a daemon thread.
  - `app.create_app(root: Path | None = None, limits: Limits | None = None, target=run_job) -> FastAPI` and module-level `app = create_app()`. Routes: `GET /` (static `index.html`), static mount `/static`, `POST /api/jobs` (form: `mode`, `params` JSON string, `fasta` file or `fasta_text`, `primers_text` for check; validates with `validate_upload`, `build_argv`, then `cli.parse_args` + `cli.validate_args` with stderr captured — any failure → 400 `{"detail": message}`; check mode writes `primers.tsv` from `primers_text`), `GET /api/jobs/{id}`, `GET /api/jobs/{id}/files/{name}` (`name` in `{"primers.tsv", "scheme.bed", "report.txt", "alignment.fasta", "check.tsv"}`; id must match `^[0-9a-f]{32}$`; else 404), `GET /api/jobs/{id}/zip`, `GET /api/example` (`static/example.fasta`), `GET /api/limits`.
  - `server.main(argv) -> int`: `--port 8000`, `--host 127.0.0.1`, `--no-browser`; runs uvicorn, opens the browser unless told not to.

- [ ] **Step 1: Write the failing tests:**
  - `test_web_params.py`: tile argv for `{"amplicons": 8, "degenerate": False}` equals `["-n", "8", "--no-degenerate"]`; unknown key → `ValueError`; bad aligner → `ValueError`; regions list → repeated `--region`; `validate_upload` rejects >`max_seqs`, a sequence over `max_seq_len`, and non-FASTA text, each with a clear message.
  - `test_web_jobs.py` (inject a fake `target` that sleeps/writes `result.json`, small `timeout_s`): queued → running → done; two workers run two jobs while a third reports `queue_position == 1`; a job exceeding the timeout becomes `failed` with the message; `cleanup` removes an old finished job dir and keeps a fresh one.
  - `test_web_api.py` (FastAPI `TestClient`, real `run_job`; mark MAFFT-dependent tests with the existing skip): POST the genome family as `fasta_text` with `{"amplicons": 4}` → 200 with id; poll `tick()` + `GET` until done (≤ 120 s) → `result.rows` has 4 entries and `reference_length` 6000; file download works; zip contains the 4 files; `GET /api/jobs/<id>/files/../../etc/passwd` and an unknown id → 404; 400 for `-n` too large, for unknown setting, for >500 sequences; check mode with two pairs returns 2 rows.
- [ ] **Step 2: Run** `python -m pytest tests/test_web_*.py -q` — Expected: FAIL.
- [ ] **Step 3: Implement** params, jobs, app, server and the `web` dispatch in `cli.main`. Static `index.html` can be a stub in this task (Task 7 builds it).
- [ ] **Step 4: Run** the web tests — Expected: pass; full suite — Expected: pass. Manual: `viraprime web --no-browser`, `curl -F mode=tile -F 'params={"amplicons":8}' -F fasta=@examples/genomes.fasta localhost:8000/api/jobs`, poll until done.
- [ ] **Step 5: Commit** `feat(web): FastAPI job API and viraprime web`.

---

### Task 7: Web page

**Files:**
- Create/replace: `viraprime/web/static/index.html`, `viraprime/web/static/app.css`, `viraprime/web/static/app.js`, `viraprime/web/static/favicon.svg` (copy of `site/favicon.svg`)

**Interfaces:**
- Consumes: the Task 6 API.
- Produces: the page described in the spec's "Page" and "Results" sections, with this exact behaviour:
  - Tabs `Design primers` (default) / `Check my primers` (ARIA tablist).
  - Sequences: drop zone + file input + paste textarea + "Try the example" (fetches `/api/example`). JS parses FASTA immediately and shows `N sequences · median M bp`; shows limit errors from `/api/limits` before upload.
  - Amplicons: number input (default 10, min 1, max from limits) with hint `≈ L bp each · O bp overlap · 2 pools` using `round((median + (n-1)*overlap)/n)`.
  - Advanced (`<details>`): every tile key from `build_argv`, with defaults from the Global Constraints; reference dropdown filled from parsed IDs; "Target specific regions instead" switch revealing regions textarea and amplicon min/max.
  - Check tab: sequences block (shared) + primers textarea (`name forward reverse` per line).
  - Run: POST, then poll `GET /api/jobs/{id}` every 1 s; progress bar with stage label and queue position; `?job=<id>` in the URL; on load with `?job=` resume polling.
  - Results: verdict line; SVG genome map (two lanes, bars coloured by coverage: ≥ min coverage one token colour, below another; `<title>` per bar and focusable bars showing primers and coverage); table with copy buttons per primer; download links for each file and the zip; error state shows the server's message.
  - Footer note: `Files are deleted 24 hours after your job finishes. For private or larger data, use the command line.`

- [ ] **Step 1: Read** `DESIGN.md`, `site/style.css` (tokens) and `site/index.html` (builder dialog pattern). Invoke `taste-skill:design-taste-frontend` (app-screen run) for the design pass.
- [ ] **Step 2: Build** the three files. Reuse the landing page's colour, type and spacing tokens; no external scripts.
- [ ] **Step 3: Run it end to end** with `viraprime web --no-browser` and playwright-cli: click "Try the example", set 8 amplicons, run, wait for results; verify verdict, 8 bars in alternating lanes, table rows, downloads respond 200; reload with `?job=` restores results; check tab with two pairs. Screenshots desktop (1440×900) and mobile (390×844) into the scratchpad; read them.
- [ ] **Step 4: Review** with the taste lint/inspect scripts and the `web-design-guidelines` skill on the three files; fix every finding (accepted exceptions recorded in the ledger).
- [ ] **Step 5: Run** the full suite — Expected: pass.
- [ ] **Step 6: Commit** `feat(web): design-and-check page with genome map`.

---

### Task 8: Docker image, Space deploy, landing page

**Files:**
- Create: `Dockerfile`, `.dockerignore`, `.github/workflows/space.yml`
- Modify: `site/index.html`, `site/app.js`, `site/data.js`, `site/style.css` (as needed), `README.md`

**Interfaces:**
- Produces:
  - `Dockerfile`: `python:3.11-slim`; `apt-get install -y --no-install-recommends mafft`; non-root user uid 1000 (HF requirement); copy `pyproject.toml`, `README.md`, `viraprime/`; `pip install --no-cache-dir ".[web]"`; `EXPOSE 7860`; `CMD ["uvicorn", "viraprime.web.app:app", "--host", "0.0.0.0", "--port", "7860"]`.
  - `.github/workflows/space.yml`: on push to master with paths `viraprime/**`, `Dockerfile`, `pyproject.toml`; job skipped unless `secrets.HF_TOKEN` and `vars.HF_SPACE` are set (check in a first step, `if:` on later steps); builds `deploy/` with Dockerfile, pyproject, viraprime/, and a `README.md` whose front matter is `title: viraprime`, `sdk: docker`, `app_port: 7860`; uploads with `huggingface_hub` `HfApi(token).upload_folder(repo_id=HF_SPACE, repo_type="space", folder_path="deploy")`.
  - Landing page: hero primary button `Open the web app` → `APP_URL` (in `site/data.js`), secondary `Install the CLI` → GitHub; "How it works" reduced to three steps (Drop your sequences / Choose how many amplicons / Get pooled primers); options list and command builder lead with `viraprime -i genomes.fasta -n 10 -o out`, builder default mode "Tile" with an amplicon-count field, region design and check under an "Advanced" heading; FAQ gains "Do I need to install anything?" (no: use the web app; install the CLI for private or larger data).
  - README: "Run it from the web" section (link) and "Deploy your own copy" (create a Docker Space, add repo secret `HF_TOKEN` with write access and repo variable `HF_SPACE=<owner>/<space>`).

- [ ] **Step 1: Docker.** Write the Dockerfile; `docker build -t viraprime-web .`; `docker run --rm -p 7860:7860 viraprime-web`; `curl localhost:7860/` returns the page and a tile job on `examples/genomes.fasta` finishes. Expected: both succeed.
- [ ] **Step 2: Workflow.** Write `space.yml`; validate with `python -c "import yaml; yaml.safe_load(open('.github/workflows/space.yml'))"`.
- [ ] **Step 3: Ask the user** for the Space id (`<owner>/<space>`). This is the one external dependency: if it is not available yet, set `APP_URL` to the README section "Run it from the web" and record a ruling.
- [ ] **Step 4: Landing page** changes; dev server (`python3 -m http.server 5173 --bind 127.0.0.1 -d site`), playwright desktop + mobile screenshots, taste lint, web-design-guidelines; fix findings.
- [ ] **Step 5: Full suite** — Expected: pass.
- [ ] **Step 6: Commit** `feat: Docker image, Hugging Face deploy, simpler landing page`.
