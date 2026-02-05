#!/usr/bin/env python3
import argparse
import scanpy as sc
import numpy as np

def main():
    ap = argparse.ArgumentParser(description="Create HVG-subset AnnData for scVI training.")
    ap.add_argument("--input", required=True, help="Input .h5ad")
    ap.add_argument("--output", required=True, help="Output .h5ad (HVG subset)")
    ap.add_argument("--n-top-genes", type=int, default=3000, help="Number of HVGs")
    ap.add_argument("--min-cells-gene", type=int, default=10,
                    help="Drop genes expressed in fewer than this many cells (stability).")
    args = ap.parse_args()

    adata = sc.read_h5ad(args.input)

    # Preserve raw counts
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # 1) Drop ultra-rare genes to avoid mean-bin degeneracy
    # This does NOT change cell QC; it removes genes with ~no information.
    sc.pp.filter_genes(adata, min_cells=args.min_cells_gene)

    # 2) HVG selection on normalized/log data (no LOESS)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.n_top_genes,
        flavor="seurat",
        subset=True,
    )

    # Restore counts for scVI
    adata.X = adata.layers["counts"]

    adata.write_h5ad(args.output)
    print(f"Wrote HVG subset: {args.output} | cells={adata.n_obs} genes={adata.n_vars}")

if __name__ == "__main__":
    main()

