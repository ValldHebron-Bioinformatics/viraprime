# Tiling mode and web app — Design

Date: 2026-10-29
Status: approved in brainstorming, pending written-spec review

## Goal

Make viraprime simple enough that a user only has to say **how many amplicons**
they want. The tool tiles the whole genome with that many overlapping amplicons,
split into two pools, choosing primer positions that work on as many sequences
as possible. Everything else has sensible defaults and sits under "Advanced".
The same thing runs from a hosted web page, so nothing has to be installed.

Success:

- `viraprime -i genomes.fasta -n 10 -o out` produces 10 overlapping amplicons
  covering the genome end to end, in 2 pools, with `primers.tsv`,
  `scheme.bed` and `report.txt`.
- A user with only a browser drops a FASTA on the web page, picks a number,
  and gets the same result, a genome map and downloads.
- Existing commands (region design, `viraprime check`) keep working unchanged.

Decisions taken in brainstorming:

- Hosted Python web app (not in-browser WebAssembly, not local-only).
- Use case: whole-genome tiling for sequencing, **2 alternating pools**
  (ARTIC-style).
- Inputs up to ~30 kb genomes and ~500 sequences on the hosted app.
- Old modes kept but out of the way: tiling is the default everywhere.
- **Global best fit** (candidate primers + dynamic programming), not
  left-to-right greedy.
- **FastAPI + own page** in the landing page's design (not Gradio), hosted on a
  Hugging Face Docker Space.
- Defaults: overlap ≥ 50 bp, 2 pools, degenerate bases on (max 4 variants per
  primer), coverage target 95 %.

## 1. Tiling mode (engine and CLI)

### Usage

```
viraprime (-i FASTA | --aligned ALN) -n N -o OUTDIR [advanced options]
```

- `-n/--amplicons N` selects tiling mode. It cannot be combined with
  `--region`/`--regions-file`/`--amplicon-min`/`--amplicon-max`
  (argparse error, exit 1).
- Without `-n` and without `--amplicon-min/--amplicon-max`, the error message
  says `use -n N to tile the whole genome with N amplicons, or --amplicon-min/
  --amplicon-max (and optionally --region) for region design`. Every command
  that works today (including region design with no `--region`, which uses
  the whole reference as region `full`) behaves exactly as before.
- `N ≥ 1` (exit 1 otherwise).

Advanced options for tiling (defaults in brackets):
`--overlap` [50], `--tm-min/--tm-opt/--tm-max` [57/60/63], `--max-tm-diff` [2],
`--primer-len-min/--primer-len-max` [18/25], `--gc-min/--gc-max` [40/60],
`--max-mismatches` [2], `--three-prime-window` [5], `--no-degenerate`
(degenerate is **on** in tiling mode), `--max-degeneracy` [4 in tiling mode],
`--min-coverage` [95], `--aligner`, `--reference-id`, `--length-tolerance`,
`--outlier-action`, `-v`. Region-mode defaults (degenerate off, max 8) are
unchanged.

### Pipeline

1. **Prepare** — shared `_load_inputs` / `_align_inputs` (unchanged).
2. **Target size.** G = consensus length, O = `--overlap`.
   L = (G + (N − 1)·O) / N, rounded. Allowed amplicon length (outer primer
   edges) is L ± 20 %.
   - L < 150 bp → error (exit 1): `N too large: amplicons would be ~L bp;
     use at most K` where K is the largest N giving L ≥ 150.
   - L > 3000 bp → warning (long PCRs are harder), run continues.
3. **Candidate primers.** primer3 is asked for left-only and right-only primers
   in overlapping windows along the consensus (Tm, length, GC, hairpin,
   self-dimer and 3′-stability limits from the options; `PRIMER_MAX_NS_ACCEPTED`
   0). Duplicates across windows are merged.
   - Conservation pre-filter: a candidate whose consensus columns have mean
     majority frequency far below the best nearby candidates is dropped before
     scoring (cheap, keeps runtime down on large inputs).
   - Each remaining primer is scored **once** against every sequence with
     `match_primer` (same rules as design and check). Its result is the set of
     sequence IDs it matches, with the reason for each miss.
   - With degenerate on, a per-primer version of the existing greedy rule adds
     IUPAC bases (most sequences rescued first; never the 3′-terminal base;
     changes that rescue nothing are undone; fold ≤ `--max-degeneracy`).
     `degenerate_pair` is refactored so both use one per-site routine; region
     mode results do not change.
   - Pair coverage = sequences matched by both primers (the same as
     `score_pair`).
4. **Best fit (dynamic programming).** Choose exactly N amplicons
   a₁…a_N, in genome order, each a (forward, reverse) candidate pair, such that:
   - length of each amplicon is within the allowed range;
   - a₁'s forward primer starts within the first L/4 bases of the consensus and
     a_N's reverse primer ends within the last L/4 bases;
   - neighbours overlap: the region between a_{k+1}'s primers starts at least
     O bases before the region between a_k's primers ends
     (insert overlap ≥ O);
   - same-pool amplicons do not touch: a_{k+2}'s forward primer starts after
     a_k's reverse primer ends;
   - primer3 pair limits hold (max Tm difference).

   Objective, lexicographic: maximise the minimum amplicon coverage, then the
   sum of coverages, then minimise the summed primer3 penalties.
   (Implementation: binary search on the minimum-coverage threshold, DP for
   the best sum at each threshold.)
   - If no tiling exists at ±20 %, retry at ±35 %, then ±50 %; the report and
     log say which tolerance was used.
   - If still none, exit 1 with the reason (e.g. `no acceptable primers between
     12,400 and 13,900`) and a suggestion of N values worth trying.
5. **Pools.** Odd amplicons → pool 1, even → pool 2. For every pair of primers
   within a pool, `primer3.calc_heterodimer` ΔG is computed; ΔG below
   −9 kcal/mol is listed as a warning in the report (not a blocker).
6. **Ends.** The report states how many bases at each end of the reference are
   outside the region between primers (cannot be sequenced from these amplicons).

### Outputs (in `--outdir`)

- `alignment.fasta` — as today.
- `primers.tsv` — one row per amplicon: `amplicon, pool, forward_name,
  forward_seq, reverse_name, reverse_seq, forward_tm, reverse_tm, forward_gc,
  reverse_gc, forward_len, reverse_len, degenerate, forward_tm_range,
  reverse_tm_range, amplicon_start, amplicon_end, amplicon_len,
  overlap_next, coverage_pct, n_covered, n_total, below_threshold,
  not_covered`. Coordinates are 1-based on the reference (outer primer edges);
  `overlap_next` is the insert overlap with the next amplicon (empty for the
  last). Names are `viraprime_<k>_LEFT` / `viraprime_<k>_RIGHT`.
- `scheme.bed` — ARTIC primer-scheme BED: `chrom start end name pool strand
  sequence`, 0-based half-open, chrom = reference ID, pool `1`/`2`.
- `report.txt` — shared header blocks, tiling summary (N, target length,
  tolerance used, ends not covered), one line per amplicon, failure counts for
  amplicons below threshold (same wording as `failure_summary`), dimer
  warnings, summary, log.

### Exit codes

`0` every amplicon ≥ `--min-coverage`; `2` at least one below (outputs still
written); `1` error.

### Code layout

- `viraprime/tiling.py` (new): target-size maths, candidate collection,
  per-primer scoring, DP, pooling, dimer checks. Pure functions; a
  `progress(stage: str, fraction: float)` callback parameter.
- `viraprime/degenerate.py`: per-site routine shared by `degenerate_pair`.
- `viraprime/models.py`: `TilingParams`, `Amplicon`, `TilingResult`.
- `viraprime/io.py`: `write_tiling_tsv`, `write_scheme_bed`,
  `write_tiling_report`.
- `viraprime/cli.py`: `-n` flag, mode dispatch, `run_tiling`, a library entry
  point `design_tiling(records_or_path, params, progress) -> TilingResult`
  used by both CLI and web.

## 2. Web app

### Package

`viraprime/web/` (FastAPI app + static page). Installed with
`pip install viraprime[web]` (extra: `fastapi`, `uvicorn`,
`python-multipart`). `viraprime web [--port 8000]` serves it locally and opens
the browser. The hosted Space runs the same app.

### Page (follows DESIGN.md, consistent with the landing page)

**Tab "Design primers" (default)**

1. Sequences: drag-and-drop FASTA, paste, or "Try the example". Parsed in the
   browser immediately: shows `42 sequences · median 9,718 bp`, and flags bad
   format or over-limit input before upload.
2. Number of amplicons: number input with a live hint
   `≈ 1,020 bp each · 50 bp overlap · 2 pools`.
3. "Design primers" button.

**Advanced (collapsed):** overlap, Tm min/opt/max, primer length min/max,
GC min/max, max mismatches, 3′ window, degenerate on/off + max variants,
coverage target, aligner, reference sequence (dropdown of IDs). A switch
"Target specific regions instead" replaces the amplicon count with a regions
list and amplicon min/max (region mode).

**Tab "Check my primers":** sequences box + a box for `name forward reverse`
rows; runs `check`.

**Running:** progress bar with named stages (aligning, finding primers,
scoring, choosing amplicons), queue position when waiting. The job id is in the
URL (`?job=…`), so reload resumes.

**Results:** verdict line (`10 amplicons · worst coverage 97 % · all above
95 %`, or which amplicons fall short and the main reason); genome map (SVG, two
lanes for pool 1/2, bars coloured by coverage, hover/tap shows primers and
coverage); primer table with copy buttons; downloads (`primers.tsv`,
`scheme.bed`, `report.txt`, all as zip). Check results: one row per pair with
coverage and note.

### API

- `POST /api/jobs` — multipart: `fasta` file or text, `mode`
  (`tile`/`regions`/`check`), `params` JSON, `primers` text for check.
  Validates limits and parameters synchronously (HTTP 400 with the same
  messages as the CLI). Returns `{id}`.
- `GET /api/jobs/{id}` — `{status: queued|running|done|failed, stage,
  progress, queue_position, result?, error?}`; `result` is JSON of the
  amplicons/pairs, coverage and messages.
- `GET /api/jobs/{id}/files/{name}` — one output file (whitelisted names).
- `GET /api/jobs/{id}/zip` — all outputs.
- `GET /api/example` — the example FASTA.

Jobs run in a worker pool (2 workers, config by env var) inside the app
process; a FIFO queue; per-job working directory under a temp root; job IDs are
random (uuid4). Timeout 15 min per job (job marked failed). The engine is the
same library code as the CLI, with the progress callback.

### Limits and privacy

Upload ≤ 20 MB, ≤ 500 sequences, ≤ 35 kb per sequence, N ≤ 200 (env-var
configurable). Job files are deleted 24 h after completion (periodic
cleanup). No accounts, no database. The page states: files are deleted after
24 hours; for private or larger data, use the command line.

## 3. Hosting, landing page, docs

- `Dockerfile` at repo root: `python:3.11-slim`, `apt-get install mafft`,
  `pip install .[web]`, `uvicorn viraprime.web.app:app --host 0.0.0.0 --port
  7860`. Space metadata in a `README` header is generated by the deploy job
  (the repo README stays as is).
- `.github/workflows/space.yml`: on push to master touching `viraprime/**`,
  `Dockerfile` or `pyproject.toml`, push to the Space with `HF_TOKEN`; skipped
  when the secret is absent. One-time user steps (create Space, add secret)
  are documented in the README.
- Landing page (`site/`): primary hero CTA "Open the web app" (Space URL),
  secondary "Install the CLI"; "How it works" becomes three steps (drop your
  sequences → choose how many amplicons → get pooled primers); options list
  and command builder lead with `viraprime -i genomes.fasta -n 10`; region
  design and check shown lower, as Advanced.
- `README.md`: opens with the simple tiling command and the web app link;
  full option reference moves below.

## Testing (test-first)

Engine:
- target-length maths and N limits (too large → error with suggested max).
- synthetic genome: exactly N amplicons; each overlap ≥ O; same-pool
  amplicons don't touch; first/last amplicon reach the ends; pools alternate.
- family with a variable block: tiling avoids primers in it and beats an evenly
  spaced tiling's minimum coverage.
- agreement: each amplicon's primers fed to `check` give identical coverage.
- tolerance widening used and reported; no-fit error message.
- degenerate refactor: region-mode outputs on the example are unchanged.
- outputs: TSV header, BED format, exit codes, `-n` with `--region` rejected.
- existing suite stays green.

Web:
- FastAPI TestClient: job lifecycle on the synthetic family, limit
  rejections (size, count, length, bad FASTA, bad params), file whitelist,
  zip, cleanup of expired jobs.
- playwright-cli: example run end to end; desktop and mobile screenshots;
  web-design-guidelines review; taste lint.

Done = `pytest` green; `viraprime -i examples/family.fasta -n 4 -o out` run
and outputs inspected; web app run locally end to end with screenshots;
Docker image builds and serves the page.

## Out of scope

Accounts and saved history; probes; more than 2 pools; genomes > 35 kb on the
hosted app (CLI only); in-browser execution; per-pool primer concentration
balancing; off-target specificity.
