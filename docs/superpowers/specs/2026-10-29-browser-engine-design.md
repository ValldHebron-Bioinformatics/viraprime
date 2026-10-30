# Browser engine (GitHub Pages) — Design

> **Update 2026-10-30:** the browser now runs MAFFT 7.526 `disttbfast`, `tbfast`, `dndpre` and
> `dvtditr` with the method choice and arguments of `mafft --auto` (L-INS-i, FFT-NS-i or
> FFT-NS-2), so its alignment equals the command line's. Mentions of FFT-NS-2 only and
> MAFFT 7.520 below describe the first version.

Date: 2026-10-29
Status: approved in brainstorming, pending written-spec review
Builds on: `2026-10-29-tiling-and-web-design.md` (tiling mode and the FastAPI web app)

## Goal

Run the viraprime web app for free on GitHub Pages, with the analysis done inside the
visitor's browser: no server, and sequences never leave the user's computer. The same
page keeps working with the existing server engine (`viraprime web`, Docker) for large
inputs or a lab server.

Success:

- `https://valldhebron-bioinformatics.github.io/viraprime/app/` designs a tiling of
  `examples/genomes.fasta` with `-n 8` entirely in the browser, and its `primers.tsv`
  is identical to the command line's for the same input and settings.
- The same page, served by `viraprime web`, still runs jobs on the server.
- The landing page's "Open the web app" opens the browser app.

Evidence from the feasibility test (2026-10-29):

- primer3-py 2.2.0 builds for Pyodide 0.27 (Emscripten 3.1.58) and returns the same
  values as the native build; viraprime's tiling under Pyodide gives the same result as
  native, about 2.4× slower.
- MAFFT 7.520 `disttbfast` compiled with Emscripten runs FFT-NS-2 at native speed
  (160 × 13 kb HMPV: 19.0 s vs 18.1 s native) when the input is read with `-i <file>`
  (the program rewinds its input, which a stream cannot do), with `-D` for DNA, and
  built with the multithread code compiled in and run with `-C 0-0`.

Decisions taken in brainstorming:

- **One page, two engines.** Browser engine on GitHub Pages; server engine when the
  page is served by `viraprime web` / Docker.
- **The browser runs viraprime's own Python code** (Pyodide), not a JavaScript rewrite.
- **Hugging Face deploy is removed** (hosting there is now paid). Docker and
  `viraprime web` stay, presented as "Run it on your own server".

## 1. Engine design

### Hooks in viraprime (no change to CLI behaviour)

- `viraprime/align.py`: `register_aligner(name: str, fn: Callable[[list[Record]], Alignment] | None)`.
  `run_aligner(records, tool)` uses a registered function for `tool` when one exists,
  otherwise the external program as today. `check_binary` is skipped for registered
  aligners.
- `viraprime/cli.py`: `set_progress_reporter(fn: Callable[[str, float], None] | None)`.
  `make_progress(path)` also calls the reporter when one is set (the file, when given,
  is still written).

### Browser worker (`engine-worker.js`, a Web Worker)

1. Loads Pyodide and Biopython from the jsDelivr CDN, pinned to the Pyodide release the
   primer3 wheel was built for (0.27.x). Loads from the site itself: the primer3-py
   wheel, the viraprime wheel and the MAFFT module (`disttbfast.js` + `.wasm`).
2. Registers `mafft` as an aligner: writes the kept records to the worker's virtual
   file system, runs `disttbfast -D -q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0
   -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 1000 -i infile` (the arguments `mafft` uses for
   FFT-NS-2), parses the output into an `Alignment` in the input order. MUSCLE is not
   offered in the browser.
3. Registers a progress reporter that posts `{stage, progress}` to the page.
4. Writes the FASTA (and the primer table for check mode) to the virtual disk and calls
   `cli.main(argv)` with the same arguments the server builds (`params.build_argv`).
5. Posts back the same result object the server returns
   (`exit_code, mode, rows, reference_id, reference_length, warnings, error`), the
   output files as text, and a zip of them made by Python's `zipfile`.

### Page (`viraprime/web/static/`, one source for both engines)

- **Engine choice:** the server engine when `GET /api/limits` answers; otherwise the
  browser engine. The browser engine's limits are built in (same values as the server:
  500 sequences, 35 kb each, 200 amplicons) and its aligner list is `["mafft"]`.
- Both engines feed the existing results view (verdict, genome map, primer table,
  downloads). Browser downloads are Blob URLs; there is no `?job=` link.
- The analysis runs in the worker, so the page stays responsive. A **Cancel** button
  terminates the worker (a fresh one is started for the next run).
- Before leaving a page with a running job or unsaved results, the browser asks for
  confirmation (`beforeunload`).
- Browser mode shows: "Runs on your computer — your sequences are not uploaded." in
  place of the server's 24-hour deletion note.

### Files and build

- The page has one source, `viraprime/web/static/`. The Pages workflow copies it to
  `site/app/` at deploy time and adds the engine files; the web app lives at
  `…github.io/viraprime/app/`.
- The viraprime wheel is built at every deploy (`pip wheel . --no-deps`, pure Python).
- The primer3-py wheel and the MAFFT module change rarely: they are committed under
  `site/app/engine/` together with `tools/build_wasm.sh`, which rebuilds them
  (pyodide-build + Emscripten; MAFFT 7.520 from MAFFT's GitLab, BSD licence).
- `viraprime web` serves the same page; its engine files are not needed there, since the
  server engine is used.

## 2. Limits, errors and cleanup

- **First visit:** the page shows "Preparing the analysis tools…" with progress while
  about 15–20 MB load; the browser caches them afterwards.
- **Large inputs:** above about 3 million bases in total (e.g. 100 × 30 kb), the page
  asks for confirmation first: "This may take several minutes and use a lot of memory
  on your computer. For large sets, the command line is faster."
- **Out of memory / worker crash:** the run stops with a message in the results panel
  ("Your browser ran out of memory. Try fewer sequences or the command line." or what
  failed) and a Retry button.
- **viraprime's own errors** (bad FASTA, too many amplicons, no tiling) appear exactly as
  on the server, since the same code runs.
- **Cleanup:** remove `deploy/hf-space/`, `.github/workflows/space.yml` and the README's
  Hugging Face section; keep the Dockerfile and `viraprime web` under "Run it on your
  own server".
- **Landing page:** "Open the web app" points to `app/`; the FAQ "Do I need to install
  anything?" says the analysis runs on your computer and nothing is uploaded.

## 3. Testing

- **Python:** unit tests for `register_aligner` (a registered function is used instead of
  the external program, and `check_binary` is not called) and for the progress reporter
  (receives the stages, the file is still written). Existing suite stays green.
- **Browser engine, end to end in Node** (same Pyodide and WebAssembly as browsers):
  run `examples/genomes.fasta` through the worker's engine code with `-n 8`; its
  `primers.tsv` equals the command line's run on the alignment it produced; the MAFFT
  module's alignment equals the native MAFFT 7.520 `disttbfast` output for the same
  input.
- **Real browser** (playwright on the built site served locally): Try the example,
  design 8 amplicons, results, genome map, downloads (zip included), Cancel during a
  run, the Check tab; desktop and phone screenshots; taste inspector and
  web-design-guidelines review.
- **Server engine:** existing API tests stay green; the page served by `viraprime web`
  uses the server engine.

## Out of scope

Saving browser results between visits; several analyses at once; multithreaded MAFFT
(SharedArrayBuffer needs headers GitHub Pages cannot set); MUSCLE in the browser.
