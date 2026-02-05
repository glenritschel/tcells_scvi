#!/usr/bin/env python3
import argparse
from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt

PROGRAMS = {
    "Treg":      ["FOXP3","IL2RA","CTLA4","IKZF2"],
    "Tfh":       ["CXCR5","PDCD1","ICOS","BCL6"],
    "Th1":       ["TBX21","IFNG","CXCR3"],
    "Th2":       ["GATA3","IL4","IL13"],
    "Th17":      ["RORC","IL17A","IL17F","CCR6"],
    "Naive":     ["CCR7","LTB","IL7R"],
    "Cytotoxic": ["NKG7","GNLY","GZMB","PRF1","CTSW"],
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adata", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--outdir-plots", required=True)
    ap.add_argument("--use-raw", action="store_true", help="If set, score using adata.raw")
    args = ap.parse_args()

    t = sc.read_h5ad(args.adata)
    outdir = Path(args.outdir_plots)
    outdir.mkdir(parents=True, exist_ok=True)

    for name, genes in PROGRAMS.items():
        present = [g for g in genes if g in t.var_names]
        if not present:
            print(f"Skip {name}: no genes found")
            continue
        sc.tl.score_genes(t, gene_list=present, score_name=f"score_{name}", use_raw=args.use_raw)
        sc.pl.umap(t, color=f"score_{name}", show=False)
        plt.savefig(outdir / f"tcells_umap_score_{name}.png", dpi=200, bbox_inches="tight")
        plt.close()
        print(f"{name} genes_found={present}")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    t.write_h5ad(args.out, compression="gzip")
    print(f"Wrote: {args.out}")
    print(f"Wrote score PNGs to: {outdir}")

if __name__ == "__main__":
    main()

