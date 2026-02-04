SHELL := /bin/bash
PY := python

RUN_ID ?= $(shell date +%Y%m%d_%H%M%S)

.PHONY: help setup qc integrate annotate signatures de figures export all

help:
	@echo "Targets:"
	@echo "  make setup"
	@echo "  make qc"
	@echo "  make integrate"
	@echo "  make annotate"
	@echo "  make signatures"
	@echo "  make de"
	@echo "  make figures"
	@echo "  make export"
	@echo "  make all"

setup:
	$(PY) -m pip install -r requirements-pip.txt

qc:
	$(PY) scripts/01_qc.py --config configs/pipeline.yaml --run-id $(RUN_ID)

integrate:
	$(PY) scripts/02_integrate_scvi.py --config configs/pipeline.yaml --run-id $(RUN_ID)

annotate:
	$(PY) scripts/03_annotate_tcells.py --config configs/pipeline.yaml --run-id $(RUN_ID)

signatures:
	$(PY) scripts/04_signatures.py --config configs/pipeline.yaml --run-id $(RUN_ID)

de:
	$(PY) scripts/05_de_tests.py --config configs/pipeline.yaml --run-id $(RUN_ID)

figures:
	$(PY) scripts/06_make_figures.py --config configs/pipeline.yaml --run-id $(RUN_ID)

export:
	$(PY) scripts/07_export_tables.py --config configs/pipeline.yaml --run-id $(RUN_ID)

all: qc integrate annotate signatures de figures export

