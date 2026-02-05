#!/usr/bin/env python3
"""
load_into_anndata.py

Fast, safe ingestion of GEO "RAW" tarballs where each member is a gzipped
gene x cell count matrix with:

- Header row: cell IDs (no leading "gene" column name)
- First column: gene name
- Remaining columns: integer counts

This script:
- Streams from the .tar member (no extraction to disk)
- Decompresses in-memory streaming
- Builds a sparse matrix WITHOUT ever materializing a dense DataFrame
- Writes one .h5ad per member (or a selected subset)
- Optionally writes a shared genes.tsv once and verifies consistency

Designed for large batches (hundreds of GSMs) and high sparsity.
"""

from __future__ import annotations

import argparse
import gzip
import io
import sys
import tarfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix


@dataclass(frozen=True)
class MemberSpec:
    member_name: str                # e.g., GSM5836979_AB10128.txt.gz
    sample_id: str                  # e.g., GSM5836979
    out_h5ad: Path                  # output path


def _parse_sample_id(member_name: str) -> str:
    # Typical member naming: GSMxxxxxxx_*.txt.gz
    # We extract everything up to first underscore, else strip extensions.
    base = Path(member_name).name
    if "_" in base:
        return base.split("_", 1)[0]
    # Fallback
    for suf in (".txt.gz", ".gz", ".txt"):
        if base.endswith(suf):
            return base[: -len(suf)]
    return base


def _iter_members(tf: tarfile.TarFile) -> Iterable[str]:
    # Only regular files ending in .txt.gz
    for m in tf.getmembers():
        if m.isfile() and m.name.endswith(".txt.gz"):
            yield m.name


def _make_out_path(out_dir: Path, member_name: str) -> Path:
    # Output file name: member basename with .h5ad
    base = Path(member_name).name
    # Strip .txt.gz -> .h5ad (or .gz -> .h5ad)
    if base.endswith(".txt.gz"):
        stem = base[:-7]
    elif base.endswith(".gz"):
        stem = base[:-3]
    else:
        stem = base
    return out_dir / f"{stem}.h5ad"


def _read_gz_member_bytes(tf: tarfile.TarFile, member_name: str) -> bytes:
    f = tf.extractfile(member_name)
    if f is None:
        raise FileNotFoundError(f"Member not found or not extractable: {member_name}")
    return f.read()


def _split_gene_and_counts(line: str) -> Tuple[str, str]:
    """
    Split a data line into gene name and the rest (counts string).
    Uses first whitespace boundary; robust to tabs + variable spaces.
    """
    # Find first whitespace (space/tab). Faster than .split() for long lines.
    n = len(line)
    i = 0
    while i < n and not line[i].isspace():
        i += 1
    gene = line[:i]
    # Skip whitespace run
    while i < n and line[i].isspace():
        i += 1
    rest = line[i:]
    return gene, rest


def _build_sparse_from_gz_bytes(
    gz_bytes: bytes,
    *,
    dtype: np.dtype = np.int32,
    verify_ncols: Optional[int] = None,
) -> Tuple[csr_matrix, np.ndarray, np.ndarray]:
    """
    Returns: (X_cells_by_genes, cell_ids, gene_names)

    - Streams through decompressed text.
    - Parses counts using np.fromstring (C-speed).
    - Builds CSR via COO triplets (rows, cols, data) then converts.
    """
    # Decompress stream
    with gzip.GzipFile(fileobj=io.BytesIO(gz_bytes), mode="rb") as gz:
        header_bytes = gz.readline()
        if not header_bytes:
            raise ValueError("Empty gz member (no header).")
        header = header_bytes.decode("utf-8", errors="replace").strip()
        if not header:
            raise ValueError("Header row is blank.")
        cell_ids = np.array(header.split(), dtype=object)
        n_cells = cell_ids.size
        if n_cells == 0:
            raise ValueError("No cell IDs parsed from header.")

        # COO triplets (as lists of ndarrays for amortized efficiency)
        rows_chunks = []
        cols_chunks = []
        data_chunks = []
        gene_names = []

        gene_idx = 0
        for raw in gz:
            line = raw.decode("utf-8", errors="replace").rstrip("\n")
            if not line:
                continue

            gene, counts_str = _split_gene_and_counts(line)
            if not gene:
                continue

            # Parse counts quickly; np.fromstring ignores repeated whitespace
            counts = np.fromstring(counts_str, sep=" ", dtype=dtype)
            if verify_ncols is not None and counts.size != verify_ncols:
                raise ValueError(
                    f"Bad column count at gene_idx={gene_idx}, gene={gene!r}: "
                    f"got {counts.size}, expected {verify_ncols}"
                )
            if counts.size != n_cells:
                # Common issue: malformed row; fail fast with context.
                raise ValueError(
                    f"Row length mismatch at gene_idx={gene_idx}, gene={gene!r}: "
                    f"got {counts.size} counts but header has {n_cells} cells"
                )

            nz = counts.nonzero()[0]
            if nz.size:
                rows_chunks.append(nz.astype(np.int32, copy=False))
                cols_chunks.append(np.full(nz.size, gene_idx, dtype=np.int32))
                data_chunks.append(counts[nz])

            gene_names.append(gene)
            gene_idx += 1

        n_genes = gene_idx
        if n_genes == 0:
            raise ValueError("No gene rows parsed from member.")

        if rows_chunks:
            rows = np.concatenate(rows_chunks)
            cols = np.concatenate(cols_chunks)
            data = np.concatenate(data_chunks).astype(dtype, copy=False)
        else:
            # All zeros matrix (unlikely but handle)
            rows = np.array([], dtype=np.int32)
            cols = np.array([], dtype=np.int32)
            data = np.array([], dtype=dtype)

        X = csr_matrix((data, (rows, cols)), shape=(n_cells, n_genes), dtype=dtype)
        gene_names = np.array(gene_names, dtype=object)
        return X, cell_ids, gene_names


def _write_genes_tsv(path: Path, gene_names: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text("\n".join(map(str, gene_names)) + "\n", encoding="utf-8")
    tmp.replace(path)


def _load_existing_genes(path: Path) -> Optional[np.ndarray]:
    if not path.exists():
        return None
    genes = path.read_text(encoding="utf-8").splitlines()
    # Drop any empty trailing line artifacts
    genes = [g for g in genes if g != ""]
    return np.array(genes, dtype=object)


def _member_specs(
    tar_path: Path,
    out_dir: Path,
    *,
    only_prefix: Optional[str] = None,
    limit: Optional[int] = None,
) -> list[MemberSpec]:
    specs: list[MemberSpec] = []
    with tarfile.open(tar_path, "r") as tf:
        for name in _iter_members(tf):
            sample_id = _parse_sample_id(name)
            if only_prefix and not sample_id.startswith(only_prefix):
                continue
            specs.append(MemberSpec(name, sample_id, _make_out_path(out_dir, name)))
            if limit is not None and len(specs) >= limit:
                break
    return specs


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Stream GEO RAW tar members into sparse AnnData .h5ad files (fast, no dense pandas)."
    )
    ap.add_argument("--tar", dest="tar_path", required=True, type=Path, help="Path to GSE*_RAW.tar")
    ap.add_argument("--out-dir", dest="out_dir", required=True, type=Path, help="Directory for .h5ad outputs")
    ap.add_argument("--genes-tsv", dest="genes_tsv", type=Path, default=None,
                    help="Write/verify shared genes list at this path (recommended).")
    ap.add_argument("--compression", choices=["gzip", "none"], default="none",
                    help="H5AD compression. Use 'none' for speed; compress later if needed.")
    ap.add_argument("--dtype", choices=["int32", "int64"], default="int32", help="Count matrix dtype.")
    ap.add_argument("--only-prefix", default=None,
                    help="Only process members whose sample_id starts with this prefix (e.g. GSM5836979).")
    ap.add_argument("--limit", type=int, default=None, help="Process only first N matching members.")
    ap.add_argument("--verify-genes", action="store_true",
                    help="If --genes-tsv exists, verify gene order matches; else write it from first sample.")
    ap.add_argument("--progress-every", type=int, default=10, help="Print progress every N samples.")
    args = ap.parse_args()

    tar_path: Path = args.tar_path
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    dtype = np.int32 if args.dtype == "int32" else np.int64
    compression = None if args.compression == "none" else "gzip"

    specs = _member_specs(tar_path, out_dir, only_prefix=args.only_prefix, limit=args.limit)
    if not specs:
        print("No matching .txt.gz members found.", file=sys.stderr)
        return 2

    existing_genes = _load_existing_genes(args.genes_tsv) if args.genes_tsv else None
    wrote_genes = False

    t0 = time.time()
    with tarfile.open(tar_path, "r") as tf:
        for idx, spec in enumerate(specs, start=1):
            s0 = time.time()
            try:
                gz_bytes = _read_gz_member_bytes(tf, spec.member_name)
                X, cell_ids, gene_names = _build_sparse_from_gz_bytes(gz_bytes, dtype=dtype)
            except Exception as e:
                print(f"[{idx}/{len(specs)}] ERROR {spec.member_name}: {e}", file=sys.stderr)
                return 1

            # genes.tsv handling
            if args.genes_tsv and args.verify_genes:
                if existing_genes is None and not wrote_genes:
                    _write_genes_tsv(args.genes_tsv, gene_names)
                    existing_genes = gene_names
                    wrote_genes = True
                else:
                    # Verify identical gene order
                    if existing_genes is None:
                        existing_genes = gene_names
                    if gene_names.shape != existing_genes.shape or not np.array_equal(gene_names, existing_genes):
                        print(
                            f"[{idx}/{len(specs)}] ERROR gene order mismatch in {spec.member_name}",
                            file=sys.stderr,
                        )
                        return 3
            elif args.genes_tsv and (not args.verify_genes) and (existing_genes is None) and (not wrote_genes):
                # Write once from first sample, no verification
                _write_genes_tsv(args.genes_tsv, gene_names)
                wrote_genes = True

            # Build AnnData
            adata = ad.AnnData(X=X)
            adata.obs_names = cell_ids.astype(str)
            adata.var_names = gene_names.astype(str)
            adata.obs["sample_id"] = spec.sample_id

            # Write .h5ad
            spec.out_h5ad.parent.mkdir(parents=True, exist_ok=True)
            adata.write_h5ad(spec.out_h5ad.as_posix(), compression=compression)

            elapsed = time.time() - s0
            if (idx % args.progress_every) == 0 or idx == 1 or idx == len(specs):
                nnz = int(X.nnz)
                print(
                    f"[{idx}/{len(specs)}] {spec.sample_id} -> {spec.out_h5ad.name} "
                    f"shape={X.shape} nnz={nnz} time={elapsed:.1f}s",
                    flush=True,
                )

    total = time.time() - t0
    print(f"Done. Wrote {len(specs)} samples to {out_dir.resolve()} in {total/60:.1f} min.")
    if args.genes_tsv:
        print(f"genes.tsv: {Path(args.genes_tsv).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

