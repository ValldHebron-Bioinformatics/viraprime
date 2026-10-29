// MAFFT's disttbfast (FFT-NS-2, the method `mafft` uses for genome-sized input) compiled to
// WebAssembly, run synchronously so Python code can call it like a function.
// Classic script: importScripts() in a worker defines self.ViraMafft; require() in Node.
(function (root) {
'use strict';

// The arguments `mafft` passes to disttbfast for FFT-NS-2 on nucleotides (MAFFT 7.520).
const MAFFT_ARGS = ['-D', '-q', '0', '-E', '2', '-V', '-1.53', '-s', '0.0', '-W', '6', '-O',
  '-C', '0-0', '-b', '62', '-g', '0', '-f', '-1.53', '-Q', '100.0', '-h', '0', '-F', '-X', '0.1',
  '-x', '1000'];

function parseFasta(text) {
  const records = [];
  let current = null;
  for (const raw of text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line) continue;
    if (line.startsWith('>') || line.startsWith('=')) {
      current = { id: line.slice(1).split(/\s+/)[0], seq: '' };
      records.push(current);
    } else if (current) {
      current.seq += line;
    }
  }
  return records;
}

// createMafft: the module factory (mafft-disttbfast.js); wasmBinary: the .wasm bytes.
// Returns runMafft(fastaText) -> aligned FASTA text, in the input order.
function makeMafftRunner(createMafft, wasmBinary) {
  return function runMafft(fastaText) {
    let mod = null;
    const out = [];
    const err = [];
    // Built with -sWASM_ASYNC_COMPILATION=0: with the bytes given, the instance starts
    // synchronously and the pre-run hook hands us the module before createMafft returns.
    const started = createMafft({
      wasmBinary,
      noInitialRun: true,
      print: (line) => out.push(line),
      printErr: (line) => err.push(line),
      preRun: [(m) => { mod = m; }],
    });
    if (started && typeof started.catch === 'function') started.catch(() => {});
    if (!mod) throw new Error('MAFFT could not start');
    mod.FS.writeFile('/infile', fastaText);
    let code = 0;
    try {
      // -i: disttbfast rewinds its input, which a stream cannot do
      code = mod.callMain([...MAFFT_ARGS, '-i', '/infile']);
    } catch (e) {
      if (e && e.name === 'ExitStatus') code = e.status;
      else throw new Error(`MAFFT failed: ${e && e.message ? e.message : e}`);
    }
    const text = out.join('\n');
    if (code !== 0 || !text.trim()) {
      throw new Error(`MAFFT failed: ${err.filter(Boolean).slice(-3).join(' ') || `exit code ${code}`}`);
    }
    return text.replace(/^=/gm, '>');
  };
}

const api = { MAFFT_ARGS, parseFasta, makeMafftRunner };
if (typeof module === 'object' && module.exports) module.exports = api;
else root.ViraMafft = api;
})(typeof self !== 'undefined' ? self : globalThis);
