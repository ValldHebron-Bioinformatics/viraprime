// Reproducible, diverse test genomes: a random ancestor, a few clades, then point mutations and
// short insertions/deletions per sequence (so alignment methods give different answers).
const BASES = 'ACGT';

function rng(seed) {
  let s = seed >>> 0;
  return () => {
    s = (s + 0x6d2b79f5) >>> 0;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function mutate(seq, rand, subRate, indels) {
  const out = seq.split('');
  for (let i = 0; i < out.length; i += 1) {
    if (rand() < subRate) out[i] = BASES[Math.floor(rand() * 4)];
  }
  for (let k = 0; k < indels; k += 1) {
    const at = Math.floor(rand() * out.length);
    const len = 1 + Math.floor(rand() * 12);
    if (rand() < 0.5) out.splice(at, len);
    else out.splice(at, 0, ...Array.from({ length: len }, () => BASES[Math.floor(rand() * 4)]));
  }
  return out.join('');
}

// nseq sequences about `length` bases long, as FASTA text.
export function synthGenomes(nseq, length, seed = 1) {
  const rand = rng(seed);
  const root = Array.from({ length }, () => BASES[Math.floor(rand() * 4)]).join('');
  const clades = Array.from({ length: 4 }, () => mutate(root, rand, 0.12, Math.ceil(length / 300)));
  const lines = [];
  for (let i = 0; i < nseq; i += 1) {
    const seq = mutate(clades[i % clades.length], rand, 0.04, Math.ceil(length / 600));
    lines.push(`>s${i + 1}`, seq);
  }
  return `${lines.join('\n')}\n`;
}
