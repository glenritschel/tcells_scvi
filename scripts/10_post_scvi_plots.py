#!/usr/bin/env python3
import argparse
from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adata", required=True, help="HVG scVI+UMAP .h5ad (contains X_umap)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--batch-key", default="sample_id")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.adata)

    # Batch mixing plot (no legend)
    if args.batch_key in adata.obs.columns:
        sc.pl.umap(adata, color=args.batch_key, legend_loc=None, show=False)
        plt.savefig(outdir / f"umap_by_{args.batch_key}_nolegend.png", dpi=200, bbox_inches="tight")
        plt.close()

    # If any markers exist in HVG space, plot them (usually they won't)
    for g in ["CD3D","CD3E","TRAC","TRBC1","MS4A1","CD79A","LYZ","NKG7","GNLY"]:
        if g in adata.var_names:
            sc.pl.umap(adata, color=g, show=False)
            plt.savefig(outdir / f"umap_{g}.png", dpi=200, bbox_inches="tight")
            plt.close()

    print(f"Wrote plots to: {outdir}")

if __name__ == "__main__":
    main()

