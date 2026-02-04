"""
paper3_geo_pbmc_pipeline.py

- Download scRNA-seq 10x matrices from GEO Series supplementary files
- Load into Scanpy (AnnData)
- QC + initial broad cell type annotation (marker-based + optional celltypist)
- Verify T cell representation
"""

from __future__ import annotations

import os
import re
import tarfile
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

# Optional: install if you want auto-annotation
# pip install celltypist
try:
    import celltypist
    HAS_CELLTYPIST = True
except Exception:
    HAS_CELLTYPIST = False

# Optional: GEO downloader
# pip install GEOparse
try:
    import GEOparse
    HAS_GEOPARSE = True
except Exception:
    HAS_GEOPARSE = False


# ---------------------------
# Filesystem helpers
# ---------------------------

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def extract_archive(archive_path: Path, out_dir: Path) -> None:
    out_dir = ensure_dir(out_dir)

    if archive_path.suffix in [".zip"]:
        with zipfile.ZipFile(archive_path, "r") as zf:
            zf.extractall(out_dir)
        return

    # Handle .tar, .tar.gz, .tgz
    if archive_path.suffix in [".tar"] or archive_path.name.endswith((".tar.gz", ".tgz")):
        mode = "r:gz" if archive_path.name.endswith((".tar.gz", ".tgz")) else "r"
        with tarfile.open(archive_path, mode) as tf:
            tf.extractall(out_dir)
        return

    raise ValueError(f"Unsupported archive type: {archive_path}")

def find_10x_matrix_dirs(root: Path) -> List[Path]:
    """
    Search for directories that look like 10x matrices (mtx + barcodes + features/genes).
    """
    candidates = []
    for d in root.rglob("*"):
        if not d.is_dir():
            continue
        has_mtx = (d / "matrix.mtx").exists() or (d / "matrix.mtx.gz").exists()
        has_barcodes = (d / "barcodes.tsv").exists() or (d / "barcodes.tsv.gz").exists()
        has_features = (
            (d / "features.tsv").exists() or (d / "features.tsv.gz").exists() or
            (d / "genes.tsv").exists() or (d / "genes.tsv.gz").exists()
        )
        if has_mtx and has_barcodes and has_features:
            candidates.append(d)
    return sorted(set(candidates))


# ---------------------------
# GEO download
# ---------------------------

def geo_download_supplementary(gse_id: str, out_dir: Path) -> List[Path]:
    """
    Downloads supplementary files for a GEO Series using GEOparse.
    Returns local file paths.
    """
    if not HAS_GEOPARSE:
        raise RuntimeError("GEOparse not installed. `pip install GEOparse`")

    out_dir = ensure_dir(out_dir)
    gse = GEOparse.get_GEO(geo=gse_id, destdir=str(out_dir), how="full")

    # This downloads supplementary files into destdir; GEOparse names vary.
    supp_paths = []
    if hasattr(gse, "supplementary_files") and gse.supplementary_files:
        for url in gse.supplementary_files:
            # GEOparse will download during get_GEO(how="full") for many series.
            # But in some cases you need to call get_GEO + manually fetch.
            pass

    # Collect actual files present after download
    for f in out_dir.rglob("*"):
        if f.is_file() and f.suffix.lower() in [".gz", ".zip", ".tar"] or f.name.endswith((".tar.gz", ".tgz")):
            supp_paths.append(f)

    if not supp_paths:
        raise RuntimeError(
            f"No archives found under {out_dir}. "
            f"This GSE may not provide 10x matrices as supplementary files."
        )

    return sorted(supp_paths)


# ---------------------------
# Load + merge
# ---------------------------

def read_10x_dir(matrix_dir: Path, sample_id: str) -> sc.AnnData:
    """
    Reads a 10x matrix directory into AnnData and tags sample_id.
    """
    ad = sc.read_10x_mtx(
        str(matrix_dir),
        var_names="gene_symbols",
        make_unique=True
    )
    ad.obs["sample_id"] = sample_id
    return ad

def load_all_10x_from_extracted(extracted_root: Path) -> sc.AnnData:
    """
    Finds all 10x matrix dirs under extracted_root and concatenates into one AnnData.
    """
    dirs = find_10x_matrix_dirs(extracted_root)
    if not dirs:
        raise RuntimeError(f"No 10x matrix directories found under: {extracted_root}")

    adatas = []
    for i, d in enumerate(dirs, start=1):
        # heuristic sample_id from path
        sample_id = d.parts[-2] if len(d.parts) >= 2 else f"sample_{i}"
        adatas.append(read_10x_dir(d, sample_id=sample_id))

    ad = sc.concat(adatas, join="outer", label="sample_id", keys=[a.obs["sample_id"][0] for a in adatas])
    return ad


# ---------------------------
# QC
# ---------------------------

def run_qc(
    ad: sc.AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
    max_pct_mito: float = 20.0,
) -> sc.AnnData:
    ad = ad.copy()

    # Basic filtering
    sc.pp.filter_genes(ad, min_cells=min_cells)
    sc.pp.filter_cells(ad, min_genes=min_genes)

    # Mito fraction
    ad.var["mt"] = ad.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    ad = ad[ad.obs["pct_counts_mt"] <= max_pct_mito].copy()
    return ad


# ---------------------------
# Initial broad annotation
# ---------------------------

BROAD_MARKERS: Dict[str, List[str]] = {
    "T_cells": ["CD3D", "CD3E", "TRAC"],
    "CD4_T": ["IL7R", "CCR7", "LTB"],
    "CD8_T": ["CD8A", "CD8B", "NKG7", "GZMB"],
    "B_cells": ["MS4A1", "CD79A", "CD74"],
    "Plasma": ["MZB1", "XBP1", "JCHAIN"],
    "NK": ["NKG7", "GNLY", "FCGR3A"],
    "Monocytes": ["LYZ", "S100A8", "S100A9", "FCN1"],
    "DC": ["FCER1A", "CST3"],
}

def score_markers(ad: sc.AnnData, markers: Dict[str, List[str]]) -> sc.AnnData:
    ad = ad.copy()
    for label, genes in markers.items():
        genes_present = [g for g in genes if g in ad.var_names]
        if not genes_present:
            ad.obs[f"score_{label}"] = 0.0
            continue
        sc.tl.score_genes(ad, gene_list=genes_present, score_name=f"score_{label}", use_raw=False)
    return ad

def annotate_broad(ad: sc.AnnData) -> sc.AnnData:
    """
    Very simple baseline: assign the max-scoring broad class.
    You can replace this with celltypist or reference mapping later.
    """
    ad = ad.copy()
    ad = score_markers(ad, BROAD_MARKERS)
    score_cols = [c for c in ad.obs.columns if c.startswith("score_")]
    scores = ad.obs[score_cols].to_numpy()
    labels = [c.replace("score_", "") for c in score_cols]
    ad.obs["broad_type"] = [labels[i] for i in scores.argmax(axis=1)]
    return ad

def annotate_with_celltypist(ad: sc.AnnData, model: str = "Immune_All_Low.pkl") -> sc.AnnData:
    if not HAS_CELLTYPIST:
        raise RuntimeError("celltypist not installed. `pip install celltypist`")
    ad = ad.copy()

    # CellTypist expects normalized/log1p expression
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    pred = celltypist.annotate(ad, model=model, majority_voting=True)
    ad.obs["celltypist_label"] = pred.predicted_labels["majority_voting"]
    return ad


# ---------------------------
# Standard preprocessing for quick UMAP sanity-check
# ---------------------------

def quick_umap(ad: sc.AnnData, color: List[str], out_png: Optional[Path] = None) -> None:
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat_v3")
    ad = ad[:, ad.var["highly_variable"]].copy()
    sc.pp.scale(ad, max_value=10)
    sc.tl.pca(ad, n_comps=50)
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30)
    sc.tl.umap(ad)
    sc.tl.leiden(ad, resolution=0.8)

    sc.pl.umap(ad, color=color, show=(out_png is None))
    if out_png is not None:
        sc.pl.umap(ad, color=color, show=False, save=f"_{out_png.name}")  # scanpy saves under ./figures/


def summarize_t_cells(ad: sc.AnnData, label_col: str = "broad_type") -> pd.DataFrame:
    """
    Quick quantification of T cell representation.
    """
    df = (
        ad.obs.groupby(["sample_id", label_col])
        .size()
        .reset_index(name="n_cells")
        .sort_values(["sample_id", "n_cells"], ascending=[True, False])
    )
    # add sample totals + percent
    totals = ad.obs.groupby("sample_id").size().rename("sample_total").reset_index()
    df = df.merge(totals, on="sample_id", how="left")
    df["pct"] = 100.0 * df["n_cells"] / df["sample_total"]
    return df


# ---------------------------
# Main runner
# ---------------------------

def run_geo_series(gse_id: str, workdir: Path) -> Tuple[sc.AnnData, pd.DataFrame]:
    workdir = ensure_dir(workdir)
    raw_dir = ensure_dir(workdir / "raw_download")
    extracted_dir = ensure_dir(workdir / "extracted")

    # 1) download archives
    archives = geo_download_supplementary(gse_id, raw_dir)

    # 2) extract
    for a in archives:
        # Extract each into its own folder to reduce collisions
        target = extracted_dir / a.stem.replace(".tar", "").replace(".gz", "")
        ensure_dir(target)
        try:
            extract_archive(a, target)
        except ValueError:
            # ignore non-archives
            continue

    # 3) load all 10x matrices found
    ad = load_all_10x_from_extracted(extracted_dir)

    # 4) qc
    ad = run_qc(ad, min_genes=200, max_pct_mito=20.0)

    # 5) baseline broad annotation
    ad = annotate_broad(ad)

    # 6) optional celltypist (better first-pass labels)
    # ad = annotate_with_celltypist(ad)

    # 7) quick sanity-check embedding
    quick_umap(ad, color=["leiden", "broad_type", "sample_id"])

    # 8) quantify T cells
    summary = summarize_t_cells(ad, label_col="broad_type")

    # 9) persist intermediate
    out_h5ad = workdir / f"{gse_id}.qc.annotated.h5ad"
    ad.write(out_h5ad)

    return ad, summary


if __name__ == "__main__":
    # Example: start with GSE195452 as the PBMC/blood anchor
    gse = "GSE195452"
    adata, tcell_summary = run_geo_series(gse, Path(f"./work_{gse}"))
    print(tcell_summary.head(30))

