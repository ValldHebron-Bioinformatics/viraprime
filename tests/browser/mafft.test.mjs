// node --test tests/browser/mafft.test.mjs
import test from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { createRequire } from "node:module";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "../..");
const require = createRequire(import.meta.url);
const createMafft = require(path.join(ROOT, "site/app/engine/mafft-disttbfast.js"));
const { makeMafftRunner, parseFasta } = require(path.join(ROOT, "viraprime/web/static/mafft-runner.js"));

const wasm = fs.readFileSync(path.join(ROOT, "site/app/engine/mafft-disttbfast.wasm"));
const input = fs.readFileSync(path.join(ROOT, "examples/genomes.fasta"), "utf8");

test("aligns the example genomes synchronously, twice in one process", () => {
  const runMafft = makeMafftRunner(createMafft, wasm);
  for (let round = 0; round < 2; round += 1) {
    const aligned = runMafft(input);          // a plain value, not a Promise
    assert.equal(typeof aligned, "string");
    const rows = parseFasta(aligned);
    const source = parseFasta(input);
    assert.equal(rows.length, 24);
    assert.deepEqual(rows.map((r) => r.id), source.map((r) => r.id));
    assert.deepEqual(new Set(rows.map((r) => r.seq.length)), new Set([9000]));
    rows.forEach((r, i) => assert.equal(r.seq.replace(/-/g, "").toUpperCase(), source[i].seq.toUpperCase()));
  }
});

test("reports MAFFT failures as errors", () => {
  const runMafft = makeMafftRunner(createMafft, wasm);
  assert.throws(() => runMafft("no fasta here\n"), /MAFFT failed.*0 sequence/);
});
