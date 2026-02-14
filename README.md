[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18644684.svg)](https://doi.org/10.5281/zenodo.18644684)

# Paper 3 – T Cell Functional States in Systemic Sclerosis

This repository contains the full, reproducible computational pipeline for **Paper 3** of a multi-paper PhD thesis investigating immune scarring and fibrosis in systemic sclerosis (SSc).

## Paper 3 Focus

Functional state characterization of peripheral blood T cells in systemic sclerosis using single-cell RNA sequencing and probabilistic latent variable modeling.

## Key Features

- **Batch correction**: scVI integration across 578 samples (67,592 cells)
- **T cell extraction**: Marker-based identification (14,617 T cells)
- **Functional states**: Program-based scoring for 7 states (Naive, Th1, Th2, Th17, Tfh, Treg, Cytotoxic)
- **Validation**: Leiden clustering at multiple resolutions
- **Reproducibility**: Fixed random seeds, complete documentation

## Quick Links

- **Manuscript**: [Frontiers in Immunology - Submitted]
- **Input Data**: [GEO GSE195452](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195452)
- **Archived Code**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18644684.svg)](https://doi.org/10.5281/zenodo.18644684)

## Dataset

- **GEO accession:** GSE195452
- **Sample type:** PBMCs (systemic sclerosis and controls)
- **Data format:** Processed GSM gene × cell count tables (not FASTQ / 10x matrices)

---

## Overview of the Pipeline

The pipeline performs the following steps:

1. Quality control and filtering of PBMC single-cell data  
2. Selection of highly variable genes (HVGs)  
3. scVI model training for batch-corrected latent representation  
4. UMAP embedding and batch-mixing diagnostics  
5. Extraction of T cells using conservative marker-based gating  
6. Clustering of T cells in scVI latent space  
7. Functional program scoring (Naive, Th1, Th2, Th17, Tfh, Treg, Cytotoxic)  
8. Assignment of functional state labels  
9. Export of per-sample functional state counts and fractions  
10. Generation of manuscript-ready figures  

All major steps are implemented as standalone scripts and orchestrated using a Makefile.

---

## Reproducibility Guarantees

- CPU-only execution (no GPU required)
- Fixed random seeds for scVI and downstream analysis
- Threading restricted to avoid nondeterministic floating-point behavior
- All parameters centralized in `configs/paper3.yaml`
- Automatic capture of run provenance and software versions

Re-running the pipeline reproduces:
- The same biological conclusions
- Near-identical embeddings and cluster assignments
- Identical tables and figures (up to floating-point noise)

---

## Requirements

- Linux or WSL2
- Python ≥ 3.10
- Conda or Mamba environment recommended

### Key Python Packages

- `scanpy`
- `anndata`
- `scvi-tools`
- `torch`
- `lightning`
- `numpy`
- `scipy`
- `pandas`

Exact package versions used in any run are recorded automatically.

### Setup
```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/paper3-tcell-functional-states.git
cd paper3-tcell-functional-states

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Dependencies

scanpy==1.9.3
scvi-tools==1.0.0
numpy==1.24.3
pandas==2.0.2
scipy==1.10.1
matplotlib==3.7.1
seaborn==0.12.2
anndata==0.9.1

---

## How to Run the Pipeline

### View Available Targets

```bash
make
```

### Expected Runtime
- **Total**: ~2-4 hours on standard laptop
- **scVI integration**: ~30-60 minutes (most time-consuming)

## Repository Structure
```bash
├── code/                # Analysis scripts
├── data/                # Gene programs and reference files
├── docs/                # Documentation
├── README.md            # This file
├── requirements.txt     # Python dependencies
└── LICENSE.txt          # MIT License
```

## Citation

If you use this code, please cite:

```bash
Ritschel, G.C. (2026). Single-cell transcriptomic analysis reveals program-based T cell functional states in systemic sclerosis peripheral blood.
Submitted to Frontiers in Immunology.

Code: https://doi.org/10.5281/zenodo.18644684
```

## License

MIT License - See LICENSE

## Contact

Glen C. Ritschel  
Email: glen.ritschel@gmail.com

## Acknowledgments

Analysis uses data from GEO accession GSE195452.

## Related Work

This work is part of a research program examining the immune-fibrosis pathway in systemic sclerosis:

- **Paper 1**: Fibroblast activation states (Submitted to Frontiers in Medicine)
- **Paper 2**: EBV immune scarring in B cells (Submitted to PLOS One)
- **Paper 3**: T cell functional states (This work)

