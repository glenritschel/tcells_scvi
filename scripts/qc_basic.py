#!/usr/bin/env python3
"""
qc_basic.py — conservative QC for merged GSE195452 scRNA-seq counts

Reads:
  fastq/data/combined_raw.h5ad

Writes:
  fastq/data/GSE195452_qc_basic.h5ad
  fastq/data/qc_reports/qc_summary.tsv
  fastq/data/qc_reports/qc_thresholds.txt

Notes:
- Conservative filters (good first pass for low-depth data):
    min_counts >= 100
    min_genes  >= 100
    pct_counts_mt <= 20  (only if MT genes exist)
- Does NOT do HVG selection or normalization/log1p.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp


def detect_mt(var_names: pd.Index) -> np.ndarray:
    """
    Return boolean mask for mitochondrial genes.
    Handles common conventions: MT-*, mt-*, and some Ensembl-style cases.
    """
    v = var_names.astype(str)
    upper = pd.Index([x.upper() for x in v])
    return upper.str.startswith("MT-")


def ensure_categorical_obs(adata: ad.AnnData, key: str) -> None:
    if key in adata.obs.columns:
        adata.obs[key] = adata.obs[key].astype(str).astype("category")


def sum_axis1(X) -> np.ndarray:
    """Fast per-row sums for dense or sparse."""
    if sp.issparse(X):
        return np.asarray(X.sum(axis=1)).ravel()
    return X.sum(axis=1)


def nnz_axis1(X) -> np.ndarray:
    """Fast per-row nonzero counts for dense or sparse."""
    if sp.issparse(X):
        return np.diff(X.indptr)
    return (X != 0).sum(axis=1)


def compute_qc(adata: ad.AnnData) -> None:
    """
    Adds basic QC fields to adata.obs:
      - total_counts
      - n_genes_by_counts
      - pct_counts_mt (if MT genes detected)
    """
    X = adata.X

    adata.obs["total_counts"] = sum_axis1(X).astype(np.float64)
    adata.obs["n_genes_by_counts"] = nnz_axis1(X).astype(np.int32)

    mt_mask = detect_mt(adata.var_names)
    if mt_mask.any():
        if sp.issparse(X):
            mt_counts = np.asarray(X[:, mt_mask].sum(axis=1)).ravel()
        else:
            mt_counts = X[:, mt_mask].sum(axis=1)
        mt_counts = mt_counts.astype(np.float64)
        adata.obs["pct_counts_mt"] = (mt_counts / np.maximum(adata.obs["total_counts"].to_numpy(), 1.0)) * 100.0
    else:
        # Keep column absent if no MT genes; caller can branch on that.
        pass


def summarize(adata: ad.AnnData, out_tsv: Path) -> None:
    cols = ["total_counts", "n_genes_by_counts"]
    if "pct_counts_mt" in adata.obs.columns:
        cols.append("pct_counts_mt")

    def q(s: pd.Series) -> dict:
        return {
            "min": float(s.min()),
            "p01": float(s.quantile(0.01)),
            "p05": float(s.quantile(0.05)),
            "median": float(s.median()),
            "p95": float(s.quantile(0.95)),
            "p99": float(s.quantile(0.99)),
            "max": float(s.max()),
        }

    rows = []
    for c in cols:
        d = q(adata.obs[c])
        d["metric"] = c
        rows.append(d)

    df = pd.DataFrame(rows).set_index("metric")
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--input",
        default="fastq/data/combined_raw.h5ad",
        help="Input merged raw h5ad",
    )
    ap.add_argument(
        "--output",
        default="fastq/data/GSE195452_qc_basic.h5ad",
        help="Output QC-filtered h5ad",
    )
    ap.add_argument("--min-counts", type=int, default=100)
    ap.add_argument("--min-genes", type=int, default=100)
    ap.add_argument("--max-pct-mt", type=float, default=20.0)
    args = ap.parse_args()

    in_path = Path(args.input).resolve()
    out_path = Path(args.output).resolve()

    if not in_path.exists():
        raise FileNotFoundError(in_path)

    print(f"Loading: {in_path}")
    a = ad.read_h5ad(in_path)
    print(f"Loaded: {a.shape}  (cells, genes)")

    # Ensure expected obs fields are categorical if present
    ensure_categorical_obs(a, "sample_id")
    ensure_categorical_obs(a, "batch")

    # Basic QC metrics
    compute_qc(a)

    # Build filter mask (conservative)
    keep = (a.obs["total_counts"].to_numpy() >= args.min_counts) & (
        a.obs["n_genes_by_counts"].to_numpy() >= args.min_genes
    )

    used_mt = False
    if "pct_counts_mt" in a.obs.columns:
        used_mt = True
        keep &= a.obs["pct_counts_mt"].to_numpy() <= args.max_pct_mt

    n0 = a.n_obs
    n_keep = int(keep.sum())
    print(f"Keeping {n_keep}/{n0} cells ({n_keep/n0:.3%})")

    a_qc = a[keep].copy()

    # Preserve raw counts as-is; do NOT normalize here.
    # Make obs_names unique as a safety net (should already be unique from merge step).
    if not a_qc.obs_names.is_unique:
        a_qc.obs_names_make_unique(join="-")

    # Write outputs
    out_path.parent.mkdir(parents=True, exist_ok=True)
    a_qc.write_h5ad(out_path, compression="gzip")
    print(f"Wrote: {out_path}  shape={a_qc.shape}")

    report_dir = out_path.parent / "qc_reports"
    summarize(a_qc, report_dir / "qc_summary.tsv")

    thresholds_txt = report_dir / "qc_thresholds.txt"
    thresholds_txt.write_text(
        "\n".join(
            [
                f"input={in_path}",
                f"output={out_path}",
                f"min_counts={args.min_counts}",
                f"min_genes={args.min_genes}",
                f"max_pct_mt={args.max_pct_mt} (applied={used_mt})",
                f"cells_before={n0}",
                f"cells_after={n_keep}",
            ]
        )
        + "\n"
    )
    print(f"Wrote QC reports to: {report_dir}")


if __name__ == "__main__":
    main()

