#!/usr/bin/env python3
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re
import sys

GSM_RE = re.compile(r"(GSM\d{7})")

def infer_gsm(path: Path) -> str:
    m = GSM_RE.search(path.name)
    if not m:
        raise ValueError(f"Cannot infer GSM from filename: {path.name}")
    return m.group(1)

def main() -> int:
    here = Path(__file__).resolve().parent          # .../scripts
    repo = here.parent                              # repo root (adjust if needed)

    # Your known output directory:
    h5ad_dir = (repo / "fastq" / "data" / "h5ad").resolve()

    if not h5ad_dir.exists():
        print(f"ERROR: h5ad_dir does not exist: {h5ad_dir}", file=sys.stderr)
        return 2

    files = sorted(h5ad_dir.glob("*.h5ad"))
    print(f"h5ad_dir: {h5ad_dir}")
    print(f"Found {len(files)} .h5ad files")

    by_gsm: dict[str, list[Path]] = defaultdict(list)
    for f in files:
        try:
            gsm = infer_gsm(f)
        except ValueError:
            continue
        by_gsm[gsm].append(f)

    # Keep exactly one file per GSM: newest by mtime (you can change policy)
    kept: list[Path] = []
    dupes: list[tuple[str, int]] = []

    for gsm, fs in sorted(by_gsm.items()):
        if len(fs) > 1:
            dupes.append((gsm, len(fs)))
        fs_sorted = sorted(fs, key=lambda p: p.stat().st_mtime)
        kept.append(fs_sorted[-1])

    print(f"Unique GSM: {len(kept)}")
    if dupes:
        print(f"Duplicate GSM count: {len(dupes)} (showing first 10)")
        for gsm, n in dupes[:10]:
            print(f"  {gsm}: {n}")

    # OPTIONAL: write a manifest of what to use for merge
    manifest = h5ad_dir / "dedup_manifest.txt"
    manifest.write_text("\n".join(str(p) for p in kept) + "\n")
    print(f"Wrote manifest: {manifest}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())

