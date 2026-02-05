#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import scanpy as sc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adata", required=True)
    ap.add_argument("--sample-key", default="sample_id")
    ap.add_argument("--state-key", default="functional_state")
    ap.add_argument("--out-counts", required=True)
    ap.add_argument("--out-frac", required=True)
    args = ap.parse_args()

    t = sc.read_h5ad(args.adata)
    sk, fk = args.sample_key, args.state_key
    if sk not in t.obs.columns or fk not in t.obs.columns:
        raise SystemExit(f"Missing required obs columns: {sk} or {fk}")

    ct = pd.crosstab(t.obs[sk], t.obs[fk])
    frac = ct.div(ct.sum(axis=1), axis=0)

    Path(args.out_counts).parent.mkdir(parents=True, exist_ok=True)
    ct.to_csv(args.out_counts)
    frac.to_csv(args.out_frac)

    print(f"Wrote: {args.out_counts}")
    print(f"Wrote: {args.out_frac}")
    print("Samples:", ct.shape[0], "States:", ct.shape[1])

if __name__ == "__main__":
    main()

