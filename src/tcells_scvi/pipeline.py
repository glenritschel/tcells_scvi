import os
import scanpy as sc

def run_qc(cfg: dict, outdir: str) -> None:
    # Placeholder: you will adapt to your dataset loader
    # For now, assume an h5ad exists at data/processed/input.h5ad
    in_path = os.path.join(cfg["paths"]["processed_dir"], "input.h5ad")
    ad = sc.read_h5ad(in_path)

    # Basic QC example (adapt as needed)
    ad.var["mt"] = ad.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True)

    qc = cfg["qc"]
    ad = ad[ad.obs.n_genes_by_counts >= qc["min_genes"], :].copy()
    ad = ad[ad.obs.pct_counts_mt <= qc["max_mito_pct"], :].copy()

    out_path = os.path.join(outdir, "qc_filtered.h5ad")
    ad.write(out_path)

