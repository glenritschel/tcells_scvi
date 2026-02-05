from pathlib import Path
import anndata as ad

REPO = Path(__file__).resolve().parents[1]
H5AD = REPO / "fastq" / "data" / "combined_raw.h5ad"

a = ad.read_h5ad(H5AD)
print("loaded:", H5AD, a.shape)

