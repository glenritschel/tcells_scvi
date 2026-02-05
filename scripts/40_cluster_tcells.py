#!/usr/bin/env python3
import argparse
from pathlib import Path
import scanpy as sc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adata", required=True, help="T-cell subset .h5ad (must contain X_scVI)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-neighbors", type=int, default=15)
    ap.add_argument("--leiden", type=float, nargs="+", default=[0.5, 0.8], help="Leiden resolutions")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    t = sc.read_h5ad(args.adata)
    if "X_scVI" not in t.obsm:
        raise SystemExit("Missing X_scVI in obsm.")

    sc.pp.neighbors(t, use_rep="X_scVI", n_neighbors=args.n_neighbors)
    sc.tl.umap(t, random_state=args.seed)

    for r in args.leiden:
        key = f"leiden_r{str(r).replace('.','')}"
        sc.tl.leiden(t, resolution=r, key_added=key, random_state=args.seed)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    t.write_h5ad(args.out, compression="gzip")
    print(f"Wrote: {args.out}")
    for r in args.leiden:
        key = f"leiden_r{str(r).replace('.','')}"
        print(key, "n_clusters=", t.obs[key].nunique())

if __name__ == "__main__":
    main()

