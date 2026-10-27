"""Write the synthetic test family to examples/family.fasta (used by the README example)."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from tests.synthetic import make_family, write_fasta  # noqa: E402

if __name__ == "__main__":
    out = write_fasta(ROOT / "examples" / "family.fasta", make_family()[1])
    print(f"Wrote {out}")
