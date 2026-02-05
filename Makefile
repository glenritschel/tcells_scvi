# ===============================
# Reproducibility hardening
# ===============================

export PYTHONHASHSEED ?= 0
export TZ ?= UTC
export LC_ALL ?= C
export LANG ?= C

# Prevent nondeterminism from threaded math libs (CPU)
export OMP_NUM_THREADS ?= 1
export MKL_NUM_THREADS ?= 1
export OPENBLAS_NUM_THREADS ?= 1
export NUMEXPR_NUM_THREADS ?= 1


# ---- Reproducibility hardening ----
export PYTHONHASHSEED ?= 0
export TZ ?= UTC
export LC_ALL ?= C
export LANG ?= C

# Reduce nondeterminism from threaded math libs (CPU)
export OMP_NUM_THREADS ?= 1
export MKL_NUM_THREADS ?= 1
export OPENBLAS_NUM_THREADS ?= 1
export NUMEXPR_NUM_THREADS ?= 1

SHELL := /bin/bash

.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c

PY ?= python
SCRIPTS := scripts

DATA := fastq/data

# ---------- Inputs you may change ----------
# Full pre-scVI QC'd dataset (all genes)
QC_FULL := $(DATA)/GSE195452_qc_75_75_min25cells.h5ad

# HVG subset output (used for scVI speed/stability)
HVG_N := 3000
HVG := $(DATA)/GSE195452_qc_75_75_min25cells_hvg$(HVG_N).h5ad

# scVI run directory + model artifact
SCVI_RUN := $(DATA)/scvi_hvg$(HVG_N)_run1
SCVI_MODEL := $(SCVI_RUN)/scvi_model/model.pt

# scVI latent+UMAP h5ad output (produced by post_scvi_umap.py)
SCVI_LATENT_UMAP := $(DATA)/GSE195452_scvi_hvg$(HVG_N)_latent_umap.h5ad

# Full-gene object with embeddings transferred from HVG latent/UMAP
FULL_WITH_EMB := $(DATA)/GSE195452_qc_75_75_min25cells_with_scvi_umap.h5ad

# T-cell subset and downstream objects
TCELLS := $(DATA)/GSE195452_Tcells_scvi.h5ad
TCELLS_CLUSTERED := $(DATA)/GSE195452_Tcells_scvi_clustered.h5ad
TCELLS_SCORED := $(DATA)/GSE195452_Tcells_scvi_clustered_scored.h5ad
TCELLS_LABELED := $(DATA)/GSE195452_Tcells_scvi_clustered_labeled.h5ad

# Tables/figures
R05_LABELS := $(DATA)/leiden_r05_cluster_labels.csv
STATE_COUNTS := $(DATA)/tcells_counts_by_sample_functional_state.csv
STATE_FRAC := $(DATA)/tcells_frac_by_sample_functional_state.csv

PLOTS_DIR := $(DATA)
SCORE_PLOTS_DIR := $(DATA)

# Parameters used in scripts
BATCH_KEY := sample_id
N_LATENT := 30
MAX_EPOCHS := 100
BATCH_SIZE := 512
NUM_WORKERS := 0
SEED := 0

TCELL_QUANTILE := 0.85
N_NEIGHBORS := 15
LEIDEN_RES := 0.5 0.8
CLUSTER_KEY := leiden_r05

# ---------- Convenience targets ----------
.PHONY: all help results

# Default target when running `make`
all: help

# What a reviewer / CI / you should run to reproduce results
results: provenance tables figures

help:
	@echo ""
	@echo "Paper 3 – T cell functional states pipeline"
	@echo ""
	@echo "Primary targets:"
	@echo "  make results       Reproduce all tables and figures (recommended)"
	@echo ""
	@echo "Stage targets:"
	@echo "  make qc            Run QC summary/filtering (qc_basic.py)"
	@echo "  make hvg           Create HVG subset: $(HVG)"
	@echo "  make scvi          Train scVI model: $(SCVI_MODEL)"
	@echo "  make postscvi      Compute latent space + UMAP embeddings"
	@echo "  make tcells        Extract, cluster, score, and label T cells"
	@echo "  make tables        Export per-sample functional state tables"
	@echo "  make figures       Generate manuscript figures (PNGs)"
	@echo ""
	@echo "Maintenance targets:"
	@echo "  make clean         Remove derived figures/CSVs (keep models and h5ad)"
	@echo "  make distclean     Remove all derived outputs including models"
	@echo ""

# ---------- Stage 0: QC ----------
# NOTE: If qc_basic.py creates a *new* file in your repo, wire it here.
# Right now your QC file is already created, so this is a placeholder target.
qc:
	@echo "QC stage: using existing QC input: $(QC_FULL)"
	@ls -lah "$(QC_FULL)" >/dev/null

# ---------- Stage 1: HVG ----------
$(HVG): configs/paper3.yaml $(QC_FULL) $(SCRIPTS)/make_hvg_input.py
	$(PY) $(SCRIPTS)/make_hvg_input.py \
	  --input "$(QC_FULL)" \
	  --output "$(HVG)" \
	  --n-top-genes "$(HVG_N)" \
	  --batch-key "$(BATCH_KEY)"

hvg: $(HVG)

# ---------- Stage 2: scVI training ----------
$(SCVI_MODEL): configs/paper3.yaml $(HVG) $(SCRIPTS)/train_scvi_cpu.py
	mkdir -p "$(SCVI_RUN)"
	$(PY) $(SCRIPTS)/train_scvi_cpu.py \
	  --input "$(HVG)" \
	  --outdir "$(SCVI_RUN)" \
	  --batch-key "$(BATCH_KEY)" \
	  --n-latent "$(N_LATENT)" \
	  --max-epochs "$(MAX_EPOCHS)" \
	  --batch-size "$(BATCH_SIZE)" \
	  --num-workers "$(NUM_WORKERS)" \
	  --seed "$(SEED)"

scvi: $(SCVI_MODEL)

# ---------- Stage 3: post-scVI latent + UMAP ----------
$(SCVI_LATENT_UMAP): configs/paper3.yaml $(HVG) $(SCVI_MODEL) $(SCRIPTS)/post_scvi_umap.py
	$(PY) $(SCRIPTS)/post_scvi_umap.py \
	  --input "$(HVG)" \
	  --model "$(SCVI_MODEL)" \
	  --out "$(SCVI_LATENT_UMAP)" \
	  --batch-key "$(BATCH_KEY)"

postscvi: $(SCVI_LATENT_UMAP)

# ---------- Stage 4: analysis steps (formerly here-docs) ----------
$(FULL_WITH_EMB): configs/paper3.yaml $(QC_FULL) $(SCVI_LATENT_UMAP) $(SCRIPTS)/20_transfer_embeddings.py
	$(PY) $(SCRIPTS)/20_transfer_embeddings.py \
	  --full "$(QC_FULL)" \
	  --hvg "$(SCVI_LATENT_UMAP)" \
	  --out "$(FULL_WITH_EMB)"

$(TCELLS): configs/paper3.yaml $(FULL_WITH_EMB) $(SCRIPTS)/30_extract_tcells.py
	$(PY) $(SCRIPTS)/30_extract_tcells.py \
	  --adata "$(FULL_WITH_EMB)" \
	  --out "$(TCELLS)" \
	  --t-quantile "$(TCELL_QUANTILE)" \
	  --require-embeddings

$(TCELLS_CLUSTERED): configs/paper3.yaml $(TCELLS) $(SCRIPTS)/40_cluster_tcells.py
	$(PY) $(SCRIPTS)/40_cluster_tcells.py \
	  --adata "$(TCELLS)" \
	  --out "$(TCELLS_CLUSTERED)" \
	  --n-neighbors "$(N_NEIGHBORS)" \
	  --leiden $(LEIDEN_RES)

$(TCELLS_SCORED): configs/paper3.yaml $(TCELLS_CLUSTERED) $(SCRIPTS)/50_score_programs.py
	$(PY) $(SCRIPTS)/50_score_programs.py \
	  --adata "$(TCELLS_CLUSTERED)" \
	  --out "$(TCELLS_SCORED)" \
	  --outdir-plots "$(SCORE_PLOTS_DIR)"

$(TCELLS_LABELED): configs/paper3.yaml $(TCELLS_SCORED) $(SCRIPTS)/60_label_clusters.py
	$(PY) $(SCRIPTS)/60_label_clusters.py \
	  --adata "$(TCELLS_SCORED)" \
	  --cluster-key "$(CLUSTER_KEY)" \
	  --out-h5ad "$(TCELLS_LABELED)" \
	  --out-csv "$(R05_LABELS)" \
	  --label-col functional_state

tcells: $(TCELLS_LABELED)

# ---------- Tables ----------
$(STATE_COUNTS) $(STATE_FRAC): configs/paper3.yaml $(TCELLS_LABELED) $(SCRIPTS)/70_export_state_tables.py
	$(PY) $(SCRIPTS)/70_export_state_tables.py \
	  --adata "$(TCELLS_LABELED)" \
	  --sample-key "$(BATCH_KEY)" \
	  --state-key functional_state \
	  --out-counts "$(STATE_COUNTS)" \
	  --out-frac "$(STATE_FRAC)"

tables: $(STATE_COUNTS) $(STATE_FRAC)

RUN_MANIFEST := $(DATA)/run_manifest.txt
VERSIONS := $(DATA)/versions.txt

.PHONY: provenance
provenance: $(RUN_MANIFEST) $(VERSIONS)

$(RUN_MANIFEST): configs/paper3.yaml Makefile
	@mkdir -p "$(DATA)"
	@echo "date_utc=$$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$(RUN_MANIFEST)"
	@echo "git_commit=$$(git rev-parse HEAD 2>/dev/null || echo 'NA')" >> "$(RUN_MANIFEST)"
	@echo "git_status=$$(git status --porcelain 2>/dev/null | wc -l | tr -d ' ')" >> "$(RUN_MANIFEST)"
	@echo "config_sha256=$$(sha256sum configs/paper3.yaml | awk '{print $$1}')" >> "$(RUN_MANIFEST)"
	@echo "makefile_sha256=$$(sha256sum Makefile | awk '{print $$1}')" >> "$(RUN_MANIFEST)"
	@echo "PYTHONHASHSEED=$(PYTHONHASHSEED)" >> "$(RUN_MANIFEST)"
	@echo "OMP_NUM_THREADS=$(OMP_NUM_THREADS)" >> "$(RUN_MANIFEST)"
	@echo "MKL_NUM_THREADS=$(MKL_NUM_THREADS)" >> "$(RUN_MANIFEST)"
	@echo "OPENBLAS_NUM_THREADS=$(OPENBLAS_NUM_THREADS)" >> "$(RUN_MANIFEST)"
	@echo "NUMEXPR_NUM_THREADS=$(NUMEXPR_NUM_THREADS)" >> "$(RUN_MANIFEST)"
	@echo "Wrote: $(RUN_MANIFEST)"

$(VERSIONS):
	@mkdir -p "$(DATA)"
	@echo "python=$$($(PY) -V 2>&1)" > "$(VERSIONS)"
	@echo "pip_freeze:" >> "$(VERSIONS)"
	@$(PY) -m pip freeze >> "$(VERSIONS)" || true
	@echo "" >> "$(VERSIONS)"
	@echo "selected_versions:" >> "$(VERSIONS)"
	@$(PY) - << 'EOF' >> "$(VERSIONS)" || true
import sys
pkgs = ["scanpy","anndata","scvi-tools","torch","lightning","numpy","scipy","pandas","umap-learn","leidenalg","igraph"]
try:
    import importlib.metadata as md
except Exception:
    import importlib_metadata as md
print("python_executable=", sys.executable)
for p in pkgs:
    try:
        print(p, md.version(p))
    except Exception:
        pass
EOF
	@echo "Wrote: $(VERSIONS)"

RUN_MANIFEST := $(DATA)/run_manifest.txt
VERSIONS := $(DATA)/versions.txt

.PHONY: provenance
provenance: $(RUN_MANIFEST) $(VERSIONS)

$(RUN_MANIFEST): configs/paper3.yaml Makefile
	@mkdir -p "$(DATA)"
	@echo "date_utc=$$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$(RUN_MANIFEST)"
	@echo "git_commit=$$(git rev-parse HEAD 2>/dev/null || echo NA)" >> "$(RUN_MANIFEST)"
	@echo "config_sha256=$$(sha256sum configs/paper3.yaml | awk '{print $$1}')" >> "$(RUN_MANIFEST)"
	@echo "makefile_sha256=$$(sha256sum Makefile | awk '{print $$1}')" >> "$(RUN_MANIFEST)"
	@echo "PYTHONHASHSEED=$(PYTHONHASHSEED)" >> "$(RUN_MANIFEST)"
	@echo "OMP_NUM_THREADS=$(OMP_NUM_THREADS)" >> "$(RUN_MANIFEST)"
	@echo "OPENBLAS_NUM_THREADS=$(OPENBLAS_NUM_THREADS)" >> "$(RUN_MANIFEST)"
	@echo "Wrote: $(RUN_MANIFEST)"

$(VERSIONS):
	@mkdir -p "$(DATA)"
	@echo "python=$$($(PY) -V 2>&1)" > "$(VERSIONS)"
	@$(PY) -m pip freeze >> "$(VERSIONS)" || true
	@echo "Wrote: $(VERSIONS)"


# ---------- Figures ----------
# Optional: wire more figure scripts here if you split by figure number.
figures: provenance postscvi tcells tables
	@echo "Figures are emitted by scripts into $(DATA) (UMAPs, score plots, etc.)."
	@echo "Key outputs:"
	@ls -lah "$(DATA)"/umap_by_*_nolegend.png 2>/dev/null || true
	@ls -lah "$(DATA)"/tcells_umap_score_*.png 2>/dev/null || true
	@ls -lah "$(DATA)"/tcells_functional_state_cytokine_heatmap.png 2>/dev/null || true

# ---------- Cleanup ----------
clean:
	rm -f "$(DATA)"/umap_by_*png
	rm -f "$(DATA)"/umap_*png
	rm -f "$(DATA)"/tcells_umap_*.png
	rm -f "$(DATA)"/tcells_functional_state_*.png
	rm -f "$(DATA)"/leiden_r05_program_means.csv "$(DATA)"/leiden_r08_program_means.csv
	rm -f "$(R05_LABELS)" "$(STATE_COUNTS)" "$(STATE_FRAC)"
	@echo "Cleaned figures/tables. Kept .h5ad and scVI model outputs."

distclean: clean
	rm -f "$(HVG)" "$(SCVI_LATENT_UMAP)" "$(FULL_WITH_EMB)" "$(TCELLS)" "$(TCELLS_CLUSTERED)" "$(TCELLS_SCORED)" "$(TCELLS_LABELED)"
	rm -rf "$(SCVI_RUN)"
	@echo "Removed derived .h5ad outputs and scVI run dir."

