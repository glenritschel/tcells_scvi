#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys
import numpy as np
import anndata as ad


H5AD_DIR = Path("/home/glen/repos/tcells_scvi/fastq/data/h5ad")
OUT = Path("/home/glen/repos/tcells_scvi/fastq/data/combined_raw.h5ad")
CHUNK_SIZE = 50  # GSMs per concat chunk


def add_basic_qc(a: ad.AnnData) -> None:
    # Works with sparse matrices efficiently
    a.obs["n_counts"] = np.asarray(a.X.sum(axis=1)).ravel()
    a.obs["n_genes"] = np.asarray((a.X > 0).sum(axis=1)).ravel()


def infer_sample_id_from_path(p: Path) -> str:
    # Handles "GSM5836979.h5ad" or "GSM5836979_AB10128.h5ad"
    stem = p.stem
    return stem.split("_")[0]


def read_one(p: Path) -> ad.AnnData:
    a = ad.read_h5ad(p)

    gsm = infer_sample_id_from_path(p)   # "GSM5836979" (repeats)
    file_id = p.stem                    # unique across 1446, e.g. "GSM5836979_AB10128" or similar

    # enforce sample_id as GSM (string)
    a.obs["sample_id"] = gsm

    # add unique file_id
    a.obs["file_id"] = file_id

    # make obs_names globally unique via file_id
    a.obs_names = a.obs["file_id"].iloc[0] + ":" + a.obs_names.astype(str)

    # choose batch definition
    a.obs["batch"] = a.obs["file_id"]   # recommended for now

    add_basic_qc(a)

    for col in ["sample_id", "file_id", "batch"]:
        a.obs[col] = a.obs[col].astype(str)

    return a


def concat_list(lst: list[ad.AnnData]) -> ad.AnnData:
    # index_unique=None preserves existing obs_names as-is; fine if already unique
    return ad.concat(
        lst,
        axis=0,
        join="outer",
        merge="same",
        index_unique=None,
    )


def assert_unique_obs(a: ad.AnnData, where: str):
    if not a.obs_names.is_unique:
        raise ValueError(f"{where}: obs_names not unique")


def main() -> int:
    paths = sorted(H5AD_DIR.glob("GSM*.h5ad"))
    if not paths:
        raise FileNotFoundError(f"No .h5ad files found under {H5AD_DIR}")

    chunks: list[ad.AnnData] = []
    current: list[ad.AnnData] = []

    for i, p in enumerate(paths, 1):
        a = read_one(p)
        assert_unique_obs(a, p.name)
        current.append(a)

        if len(current) == CHUNK_SIZE or i == len(paths):
            chunk = concat_list(current)
            chunks.append(chunk)
            current = []
            print(f"Chunked concat: {len(chunks)} chunk(s), processed {i}/{len(paths)} files")

    adata = concat_list(chunks)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT, compression="gzip")

    print("Done.")
    print(f"Samples (files): {len(paths)}")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    print(f"Wrote: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

