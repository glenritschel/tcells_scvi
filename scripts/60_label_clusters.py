#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import scanpy as sc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adata", required=True)
    ap.add_argument("--cluster-key", default="leiden_r05")
    ap.add_argument("--out-h5ad", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--label-col", default="functional_state")
    args = ap.parse_args()

    t = sc.read_h5ad(args.adata)
    ck = args.cluster_key
    if ck not in t.obs.columns:
        raise SystemExit(f"Missing cluster key in obs: {ck}")

    score_cols = [c for c in t.obs.columns if c.startswith("score_")]
    if not score_cols:
        raise SystemExit("No score_* columns found. Run scoring first.")

    df = t.obs.groupby(ck)[score_cols].mean()
    df["n_cells"] = t.obs[ck].value_counts().loc[df.index].values

    # best program per cluster
    best = df[score_cols].idxmax(axis=1).str.replace("score_", "", regex=False)
    best_score = df[score_cols].max(axis=1)

    # map to labels (your current decisions; adjust anytime)
    # - Cytotoxic → CD8/Cytotoxic
    # - Naive program → Naive/CM
    # - otherwise label by program name (or "Other/Unclear")
    label = []
    for prog in best:
        if prog == "Cytotoxic":
            label.append("Cytotoxic")
        elif prog == "Naive":
            label.append("Naive")
        elif prog in ("Th1","Th2","Th17","Tfh","Treg"):
            label.append(prog)
        else:
            label.append("Unpolarized")
    df["best_program"] = best
    df["best_score"] = best_score
    df[args.label_col] = label
    df = df.sort_values("n_cells", ascending=False)

    # apply labels to cells
    label_map = df[args.label_col].to_dict()
    t.obs[args.label_col] = t.obs[ck].map(label_map).astype("category")

    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv)
    print(f"Wrote: {args.out_csv}")

    Path(args.out_h5ad).parent.mkdir(parents=True, exist_ok=True)
    t.write_h5ad(args.out_h5ad, compression="gzip")
    print(f"Wrote: {args.out_h5ad}")

if __name__ == "__main__":
    main()

