#!/usr/bin/env bash
# Rebuild the compiled parts of the browser engine into site/app/engine/:
#   * primer3-py 2.2.0 as a Pyodide wheel (pyodide-build, Pyodide 0.27 / Emscripten 3.1.58)
#   * MAFFT 7.520 `disttbfast` (FFT-NS-2) as a WebAssembly module
# Needs: uv (https://docs.astral.sh/uv/), curl, make, a few GB of disk for the toolchain.
# Takes about 10 minutes the first time (it downloads Emscripten); later runs are faster.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT="$ROOT/site/app/engine"
WORK="${TMPDIR:-/tmp}/viraprime-wasm"
XBUILDENV=0.27.8            # same ABI (pyodide_2024_0) as the Pyodide 0.27.7 runtime the page loads
PRIMER3=2.2.0
MAFFT=7.520
mkdir -p "$WORK" "$OUT"
cd "$WORK"

echo "== build tools"
command -v uv >/dev/null || { echo "uv is required: https://docs.astral.sh/uv/" >&2; exit 1; }
[ -x venv/bin/pyodide ] || { uv venv -q --seed -p 3.12 venv; uv pip install -q -p venv/bin/python "pyodide-build==0.39.3"; }
export PATH="$WORK/venv/bin:$PATH"
export EMSDK_PYTHON="$WORK/venv/bin/python"   # emsdk needs Python >= 3.10
pyodide xbuildenv install "$XBUILDENV"
pyodide xbuildenv install-emscripten
EMSDK_DIR=$(ls -td "$HOME"/.cache/pyodide-build/.pyodide-xbuildenv-*/"$XBUILDENV"/emsdk | head -n1)

echo "== primer3-py $PRIMER3"
python -m pip download -q --no-deps --no-binary :all: "primer3-py==$PRIMER3" -d src
mkdir -p build-primer3 && tar xzf "src/primer3_py-$PRIMER3.tar.gz" -C build-primer3
(cd "build-primer3/primer3_py-$PRIMER3" && pyodide build . -o "$WORK/dist")
WHEEL=$(ls "$WORK"/dist/primer3_py-"$PRIMER3"-*wasm32.whl | head -n1)
# newer pyodide-build names the platform pyemscripten_2024_0; Pyodide 0.27 expects pyodide_2024_0
cp "$WHEEL" "$OUT/primer3_py-$PRIMER3-cp312-cp312-pyodide_2024_0_wasm32.whl"

echo "== MAFFT $MAFFT disttbfast"
curl -fsSL -o "src/mafft-v$MAFFT.tar.gz" "https://gitlab.com/sysimm/mafft/-/archive/v$MAFFT/mafft-v$MAFFT.tar.gz"
mkdir -p build-mafft && tar xzf "src/mafft-v$MAFFT.tar.gz" -C build-mafft
# shellcheck disable=SC1091
source "$EMSDK_DIR/emsdk_env.sh" >/dev/null 2>&1
cd "build-mafft/mafft-v$MAFFT/core"
make clean >/dev/null
# The multithread code stays compiled in (disttbfast needs parts of it) but runs with -C 0-0.
# disttbfast must read its input with -i <file>: it rewinds the input, which a stream cannot do.
emmake make disttbfast CC=emcc LIBS=-lm \
  CFLAGS="-O3 -Wno-error=implicit-function-declaration" \
  LDFLAGS="-O3 -sMODULARIZE=1 -sEXPORT_NAME=createMafft -sENVIRONMENT=web,worker,node \
-sINVOKE_RUN=0 -sWASM_ASYNC_COMPILATION=0 -sFORCE_FILESYSTEM=1 -sEXPORTED_RUNTIME_METHODS=FS,callMain \
-sALLOW_MEMORY_GROWTH=1 -sMAXIMUM_MEMORY=4GB -sEXIT_RUNTIME=0" >build.log 2>&1 \
  || { tail -20 build.log >&2; exit 1; }
cp disttbfast "$OUT/mafft-disttbfast.js"
cp disttbfast.wasm "$OUT/mafft-disttbfast.wasm"
cp ../license "$OUT/MAFFT-LICENSE.txt"   # BSD: the notice travels with the binaries

echo "== done"
ls -la "$OUT"
