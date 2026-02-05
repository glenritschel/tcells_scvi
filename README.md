# Paper 3 – T Cell Functional States in Systemic Sclerosis

This repository contains the full, reproducible computational pipeline for **Paper 3** of a multi-paper PhD thesis investigating immune scarring and fibrosis in systemic sclerosis (SSc).

## Paper 3 Focus

Functional state characterization of peripheral blood T cells in systemic sclerosis using single-cell RNA sequencing and probabilistic latent variable modeling.

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

---

## How to Run the Pipeline

### View Available Targets

```bash
make

