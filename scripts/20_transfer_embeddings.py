#!/usr/bin/env python3
import argparse
from pathlib import Path
import scanpy as sc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", required=True, help="Full-gene .h5ad (counts, all genes)")
    ap.add_argument("--hvg", required=True, help="HVG scVI+UMAP .h5ad (contains X_scVI and X_umap)")
    ap.add_argument("--out", required=True, help="Output .h5ad with embeddings attached")
    args = ap.parse_args()

    full = sc.read_h5ad(args.full)
    hvg  = sc.read_h5ad(args.hvg)

    if full.n_obs != hvg.n_obs:
        raise SystemExit(f"n_obs mismatch: full={full.n_obs} hvg={hvg.n_obs}")
    if not (full.obs_names == hvg.obs_names).all():
        raise SystemExit("Cell order mismatch between full and HVG objects. Refuse to transfer embeddings.")

    for k in ["X_scVI","X_umap"]:
        if k not in hvg.obsm:
            raise SystemExit(f"Missing {k} in HVG object.")
        full.obsm[k] = hvg.obsm[k]

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    full.write_h5ad(args.out, compression="gzip")
    print(f"Wrote: {args.out}")

if __name__ == "__main__":
    main()

