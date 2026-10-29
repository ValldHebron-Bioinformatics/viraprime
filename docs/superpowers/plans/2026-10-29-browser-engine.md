# Browser engine Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The viraprime web page runs the analysis inside the browser on GitHub Pages (Pyodide + primer3 + MAFFT compiled to WebAssembly), and still uses the server engine when served by `viraprime web`.

**Architecture:** Two hooks in viraprime (a registrable aligner and a progress reporter) let the unchanged CLI run inside Pyodide. A module Web Worker loads Pyodide, the primer3-py wheel, the viraprime wheel and a MAFFT `disttbfast` WebAssembly module, registers MAFFT as the `mafft` aligner, and calls `cli.main(argv)`. The page picks the server engine when `/api/limits` answers and the browser engine otherwise; both feed the existing results view. The Pages workflow copies the page to `site/app/` with the engine files.

**Tech Stack:** Python 3.11/3.12, Pyodide 0.27.x (jsDelivr), primer3-py 2.2.0 built with pyodide-build, MAFFT 7.520 built with Emscripten 3.1.58, vanilla JS (module worker), Node 22+ for engine tests, playwright-cli, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-10-29-browser-engine-design.md`

## Global Constraints

- CLI, region design, check and the server web app behave exactly as before; the existing suite (230 tests) stays green.
- Pyodide version pinned to the 0.27 series whose ABI the primer3 wheel targets (`pyodide_2024_0_wasm32`); Pyodide and Biopython from `https://cdn.jsdelivr.net/pyodide/v0.27.7/full/`.
- MAFFT in the browser runs `disttbfast -D -q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 1000 -i infile` (FFT-NS-2), output on stdout.
- Browser limits = server defaults: 500 sequences, 35 000 bp each, 200 amplicons; aligners `["mafft"]`.
- Large-input confirmation above 3 000 000 bases in total.
- Browser-mode note: "Runs on your computer — your sequences are not uploaded."
- Page source stays single: `viraprime/web/static/`; deployed copy at `site/app/`.
- Tests: `PATH=/home/vhir/anaconda3/envs/primers/bin:$PATH python -m pytest -q -p no:warnings`.

## Plan decisions

1. **CLI keeps `mafft --auto`; the browser runs FFT-NS-2.** `--auto` chooses FFT-NS-2 for genome-sized inputs (the tiling use case) but may choose L-INS-i for small inputs. The engine test therefore compares the browser's `primers.tsv` with the CLI run on the alignment the browser produced (as the spec says), and separately compares the browser alignment with native `disttbfast` (same arguments). The README notes the difference.
2. **MAFFT is called synchronously from Python.** The worker compiles the `.wasm` once at start-up (async); each alignment creates a fresh module instance synchronously (Emscripten `-sWASM_ASYNC_COMPILATION=0`, `wasmBinary` provided), writes `infile`, `callMain`, collects stdout. If synchronous instantiation proves impossible inside Pyodide's call, fall back to a two-phase run (Python prepares and filters, JS aligns, Python continues with `--aligned`) and record the ruling.
3. **Built engine binaries are committed** (`site/app/engine/`: the primer3 wheel, `mafft-disttbfast.js`, `mafft-disttbfast.wasm`) with `tools/build_wasm.sh`. Building runs downloaded build scripts, which this session's permission check blocks: **the user runs `tools/build_wasm.sh` once** (Task 2, step 3).

## Review Focus

1. A FASTA with Windows line endings, lowercase, or `U` uploaded to the browser engine → same normalisation as the server (load_fasta), same result.
2. Cancel during the MAFFT step or during scoring → worker terminated, page usable, next run works with a fresh worker.
3. Browser out of memory (very large input) → message in the results panel with Retry; the page does not hang.
4. The page served by `viraprime web` → server engine, `?job=` links still work; on Pages, a `?job=` URL is ignored with no error.
5. The worker failing to load (offline, CDN blocked) → clear message: "Could not load the analysis tools. Check your connection and reload."

---

### Task 1: Hooks in viraprime

**Files:** Modify `viraprime/align.py`, `viraprime/cli.py`. Create `tests/test_hooks.py`.

**Interfaces — Produces:**
- `align.register_aligner(name: str, fn: Callable[[list[Record]], Alignment] | None) -> None` (None unregisters); `run_aligner` uses a registered function when present (no `check_binary`, logs `Aligning N sequences with <name>`), and checks that every input ID comes back (same error as today).
- `cli.set_progress_reporter(fn: Callable[[str, float], None] | None) -> None`; `make_progress(path)` returns an emitter that calls the reporter (if set) and writes the file (if a path is given).

- [ ] **Step 1: Write the failing tests** (`tests/test_hooks.py`):
  - `test_registered_aligner_is_used`: register `"mafft"` → a function returning `alignment_from_records` of the input; monkeypatch `shutil.which` to return None; `main(["-i", fasta, "-n", "2", "-o", out])` on `make_genome_family(length=3000)` written unaligned exits 0 and the function was called once with 12 records; unregister in `finally`.
  - `test_registered_aligner_missing_ids_error`: a function that drops one record → exit 1, stderr names the missing id.
  - `test_progress_reporter_receives_stages`: reporter collects stages for an `--aligned` tiling run → includes `aligning`, `finding primers`, `scoring`, `choosing amplicons`, `done` in that order; with `--progress-file` the file is still written.
- [ ] **Step 2: Run** `python -m pytest tests/test_hooks.py -q` — Expected: FAIL (ImportError).
- [ ] **Step 3: Implement** the two hooks (module-level dicts/variables; `make_progress` without path and without reporter stays a no-op).
- [ ] **Step 4: Run** the file, then the full suite — Expected: all pass.
- [ ] **Step 5: Commit** `feat: aligner and progress hooks for embedding viraprime`.

---

### Task 2: Engine build (primer3 wheel, MAFFT module)

**Files:** Create `tools/build_wasm.sh`, `site/app/engine/` (built files), `site/app/engine/README.md` (what each file is, versions, licences, how to rebuild).

**Interfaces — Produces:** `site/app/engine/primer3_py-2.2.0-cp312-cp312-pyodide_2024_0_wasm32.whl`; `site/app/engine/mafft-disttbfast.js` (classic script defining `createMafft`; CommonJS export in Node), `site/app/engine/mafft-disttbfast.wasm`.

- [ ] **Step 1: Write** `tools/build_wasm.sh` (bash, `set -euo pipefail`, work dir under `${TMPDIR:-/tmp}/viraprime-wasm`): creates a Python 3.12 venv with `pyodide-build`; `pyodide xbuildenv install 0.27.7` and `install-emscripten`; downloads the primer3-py 2.2.0 sdist, `pyodide build .`, copies the wheel renamed to the `pyodide_2024_0_wasm32` tag; downloads `https://gitlab.com/sysimm/mafft/-/archive/v7.520/mafft-v7.520.tar.gz`, `make clean`, then `emmake make disttbfast CC=emcc LIBS=-lm CFLAGS="-O3 -Wno-error=implicit-function-declaration" LDFLAGS="-O3 -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createMafft -sENVIRONMENT=web,worker,node -sINVOKE_RUN=0 -sWASM_ASYNC_COMPILATION=0 -sFORCE_FILESYSTEM=1 -sEXPORTED_RUNTIME_METHODS=FS,callMain -sALLOW_MEMORY_GROWTH=1 -sMAXIMUM_MEMORY=4GB"`, and copies `disttbfast` / `disttbfast.wasm` to `site/app/engine/mafft-disttbfast.js` / `.wasm`. Uses `EMSDK_PYTHON` from the venv (base conda's Python 3.8 is too old).
- [ ] **Step 2: Write** `site/app/engine/README.md`.
- [ ] **Step 3: Ask the user to run** `bash tools/build_wasm.sh` in a terminal (permission check blocks it here). Expected: the three files exist; the wasm is about 0.4 MB, the wheel about 1 MB.
- [ ] **Step 4: Node smoke test** (`tests/browser/mafft.test.mjs`, run with `node --test`): instantiate the module synchronously from the `.wasm` bytes, align `examples/genomes.fasta`, assert 24 rows of 9 000 columns and each row ungapped equals the input; run twice in one process (fresh instance each time).
- [ ] **Step 5: Commit** `build: primer3 and MAFFT for the browser engine`.

---

### Task 3: Worker and browser engine

**Files:** Create `tools/build_site.py`, `tests/test_build_site.py`, `viraprime/web/static/engine-worker.js`, `viraprime/web/static/engine.js`, `tests/browser/package.json` (pinned `pyodide@0.27.7`), `tests/browser/engine.test.mjs`, `tests/test_browser_engine.py` (runs the Node tests when `node` and `tests/browser/node_modules` exist, else skips).

**Interfaces — Produces:**
- `engine-worker.js` messages: in `{type: "run", mode, argv, fasta, primers}`; out `{type: "loading", stage, progress}`, `{type: "progress", stage, progress}`, `{type: "done", result, files: {name: text}, zip: ArrayBuffer}`, `{type: "error", message}`.
- `engine.js` (ES module): `serverEngine` and `browserEngine`, each `{name, limits(): Promise<Limits>, run({mode, settings, fasta, primers}, onProgress): {promise, cancel()}}`; `pickEngine(): Promise<engine>` (server if `GET /api/limits` answers within 2 s). `browserEngine.run` builds argv with a JS copy of `params.build_argv` rules (same whitelist and flags; covered by a test that compares with the Python output for 6 settings sets).
- Worker internals: loads Pyodide from the pinned CDN, `loadPackage(["biopython", "micropip"])`, installs the primer3 and viraprime wheels from `./engine/`; registers `mafft` via a Python shim `viraprime_browser.py` (written into the worker FS) that calls the JS `runMafft(text) -> text` and parses with `load_fasta` + `alignment_from_records`; sets the progress reporter to post messages; runs `cli.main(argv + ["-o", "out"])` with stderr captured; builds the result object exactly as `jobs.run_job` does (reuse its helper by importing `viraprime.web.jobs` functions that do not need FastAPI — move `_read_rows`/`_reference` result building into `viraprime/web/result.py` so both use it).

- [ ] **Step 0: Site builder (test first).** `tools/build_site.py OUT_DIR`: copies `site/` to OUT_DIR, copies `viraprime/web/static/*` to `OUT_DIR/app/`, builds the viraprime wheel (`pip wheel . --no-deps -w OUT_DIR/app/engine/`) next to the committed engine files, and writes `OUT_DIR/app/engine/manifest.json` (`{"viraprime": "<wheel name>", "primer3": "<wheel name>", "pyodide": "0.27.7"}`), which the worker reads. Failing test first (`tests/test_build_site.py`: `app/index.html`, `app/engine/manifest.json`, the three engine files and a `viraprime-*.whl` exist; manifest names match files), then implement, then pass.
- [ ] **Step 1: Refactor** `jobs.run_job`'s result building into `viraprime/web/result.py: build_result(job_dir: Path, argv, code) -> dict`; existing web tests stay green (run them).
- [ ] **Step 2: Write the failing Node test** `tests/browser/engine.test.mjs`: loads the worker code in Node (worker logic exported as `createEngineCore({fetchFile})` from `engine-worker.js` so it runs outside a Worker); runs mode `tile`, `-n 8`, on `examples/genomes.fasta`; asserts `result.exit_code === 0`, 8 rows, `reference_length === 9000`; writes the browser's `alignment.fasta` and `primers.tsv` to a temp dir; then runs `viraprime --aligned <that alignment> -n 8` natively and asserts identical `primers.tsv`. Second case: check mode with two pairs returns 2 rows. Third: argv builder equals Python `build_argv` for 6 settings sets (Python output pre-computed by the pytest wrapper into a JSON fixture).
- [ ] **Step 3: Run** `cd tests/browser && npm install && node --test` — Expected: FAIL (module not found).
- [ ] **Step 4: Implement** `engine-worker.js` (core + `self.onmessage` wrapper) and `engine.js`.
- [ ] **Step 5: Run** the Node tests — Expected: pass; then `python -m pytest -q` — Expected: all pass (the wrapper runs or skips the Node tests).
- [ ] **Step 6: Commit** `feat(web): browser engine (Pyodide worker with primer3 and MAFFT)`.

---

### Task 4: Page uses either engine

**Files:** Modify `viraprime/web/static/app.js`, `index.html`, `app.css`.

**Interfaces — Consumes:** `engine.js` from Task 3. `app.js` becomes a module (`<script type="module">`).

Behaviour (spec §1–2): engine picked at load; limits and aligner options from the engine; browser mode shows the privacy note instead of the 24-hour note and hides nothing else; "Preparing the analysis tools…" progress during the first load; confirmation above 3 000 000 total bases; a Cancel button beside the progress bar (both engines: server engine stops polling and leaves the job; browser engine terminates the worker); `beforeunload` while running or holding unsaved browser results; downloads as Blob URLs in browser mode (files + zip); errors (load failure, out of memory, viraprime errors) in the results panel with Retry. `?job=` only used by the server engine.

- [ ] **Step 1:** Implement; keep the server-engine path identical (the existing API tests and the page's server flow must not change).
- [ ] **Step 2: Verify server engine** with `viraprime web --no-browser` + playwright: example, `-n 8`, results, reload with `?job=` — Expected: as before.
- [ ] **Step 3: Verify browser engine**: build a local Pages tree (`python tools/build_site.py`, Task 3) and serve it with `python3 -m http.server`; playwright: example, `-n 8`, results identical to the server run, downloads and zip, Cancel during MAFFT then a new run, Check tab, page with no network for the CDN → load error message. Screenshots desktop + phone.
- [ ] **Step 4:** taste lint/inspect (`--app`, run from the project for the app sheet) and web-design-guidelines on the three files; fix findings.
- [ ] **Step 5: Commit** `feat(web): one page, server or browser engine`.

---

### Task 5: Pages build, landing page, cleanup

**Files:** Modify `.github/workflows/pages.yml`, `site/index.html`, `README.md`. Delete `deploy/hf-space/`, `.github/workflows/space.yml`.

**Interfaces — Consumes:** `tools/build_site.py` (Task 3).

- [ ] **Step 1: Workflow:** `pages.yml` triggers on `site/**`, `viraprime/**`, `tools/build_site.py`, the workflow; adds `actions/setup-python@v5` (3.11), `python tools/build_site.py _site`, uploads `_site`. Validate YAML.
- [ ] **Step 2: Landing page:** "Open the web app" → `app/`; FAQ a9 says the analysis runs on your computer and nothing is uploaded; README: "Run it from the web" links the Pages app, notes FFT-NS-2 in the browser vs `--auto` on the CLI, and "Run it on your own server" (Docker, `viraprime web`); remove the Hugging Face section and files.
- [ ] **Step 3:** Full suite; site inspector (page mode, from the scratchpad) — Expected: 0 FAIL.
- [ ] **Step 4: Commit** `feat: deploy the browser app on GitHub Pages; remove Hugging Face deploy`.
