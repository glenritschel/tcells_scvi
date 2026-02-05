#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import scanpy as sc

def mean_expr(adata, genes):
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        return None, []
    x = adata[:, genes].X
    # Works for sparse/dense
    m = np.asarray(x.mean(axis=1)).ravel()
    return m, genes

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adata", required=True, help="Full-gene .h5ad with embeddings (X_scVI/X_umap)")
    ap.add_argument("--out", required=True, help="Output T-cell subset .h5ad")
    ap.add_argument("--t-quantile", type=float, default=0.85, help="Quantile threshold for tcell_score")
    ap.add_argument("--require-embeddings", action="store_true", help="Fail if embeddings missing")
    args = ap.parse_args()

    a = sc.read_h5ad(args.adata)

    if args.require_embeddings:
        for k in ["X_scVI","X_umap"]:
            if k not in a.obsm:
                raise SystemExit(f"Missing {k} in input object.")

    t_genes = ["CD3D","CD3E","CD247","TRAC","TRBC1"]
    b_genes = ["MS4A1","CD79A","CD74","HLA-DRA"]

    tscore, used_t = mean_expr(a, t_genes)
    if tscore is None:
        raise SystemExit("No T-cell marker genes found in var_names; cannot gate.")
    a.obs["tcell_score"] = tscore

    bscore, used_b = mean_expr(a, b_genes)
    if bscore is not None:
        a.obs["bcell_score"] = bscore

    thr = float(np.quantile(a.obs["tcell_score"], args.t_quantile))
    mask = a.obs["tcell_score"] >= thr

    t = a[mask].copy()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    t.write_h5ad(args.out, compression="gzip")
    print(f"Wrote: {args.out}")
    print(f"T-cell gate: quantile={args.t_quantile} thr={thr:.6g} n={t.n_obs}/{a.n_obs} used_t={used_t} used_b={used_b}")

if __name__ == "__main__":
    main()

