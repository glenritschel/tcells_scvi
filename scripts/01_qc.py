#!/usr/bin/env python
import argparse, os, yaml
from tcells_scvi.pipeline import run_qc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    outdir = os.path.join(cfg["project"]["run_output_root"], args.run_id, "qc")
    os.makedirs(outdir, exist_ok=True)

    run_qc(cfg, outdir)

if __name__ == "__main__":
    main()

