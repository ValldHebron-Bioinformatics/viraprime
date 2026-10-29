# Browser engine files

Compiled parts of the viraprime browser engine, loaded by `app/engine-worker.js`.
Rebuild them with `bash tools/build_wasm.sh` (from the repository root).

| File | What it is | Source and licence |
|---|---|---|
| `primer3_py-2.2.0-cp312-cp312-pyodide_2024_0_wasm32.whl` | primer3-py 2.2.0 built for Pyodide 0.27 (Emscripten 3.1.58) | [primer3-py](https://github.com/libnano/primer3-py), GPL-2.0 |
| `mafft-disttbfast.js`, `mafft-disttbfast.wasm` | MAFFT 7.520 `disttbfast` (FFT-NS-2) as a WebAssembly module; factory `createMafft` | [MAFFT](https://mafft.cbrc.jp/alignment/software/), BSD |

The viraprime wheel and `manifest.json` are added at deploy time by `tools/build_site.py`.
Pyodide and Biopython are loaded from `https://cdn.jsdelivr.net/pyodide/v0.27.7/full/`.
