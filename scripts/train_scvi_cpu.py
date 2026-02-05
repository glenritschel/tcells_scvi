#!/usr/bin/env python3
"""
Train scVI on an AnnData (.h5ad) using CPU (WSL-friendly).

Assumptions:
- adata.X contains raw counts (integers)
- adata.obs has:
    - 'sample_id' (or your batch_key)
    - optionally other covariates
- adata.var_names are aligned across samples (you already ensured this)

Outputs:
- model_dir/: scvi model
- output_h5ad: input adata + latent embedding in .obsm["X_scVI"]
- run_manifest.json: key run metadata

Usage:
  python scripts/train_scvi_cpu.py \
    --input fastq/data/GSE195452_qc_75_75_min25cells.h5ad \
    --outdir fastq/data/scvi_run1 \
    --batch-key sample_id \
    --latent-key X_scVI \
    --max-epochs 200 \
    --seed 0
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
import scvi


@dataclass
class RunManifest:
    started_utc: str
    duration_sec: float
    input_path: str
    output_h5ad: str
    model_dir: str
    n_obs: int
    n_vars: int
    batch_key: str
    latent_key: str
    n_latent: int
    max_epochs: int
    plan: dict
    versions: dict
    device: str
    seed: int


def utc_now_iso() -> str:
    # no timezone libs; just a coarse ISO-ish string
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def set_reproducibility(seed: int) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)
    # CPU determinism (best-effort)
    torch.use_deterministic_algorithms(False)
    os.environ["PYTHONHASHSEED"] = str(seed)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Input .h5ad (counts in .X)")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--batch-key", default="sample_id", help="Batch key in adata.obs")
    p.add_argument("--layer", default=None, help="Counts layer to use (default: .X)")
    p.add_argument("--latent-key", default="X_scVI", help='obsm key for latent, e.g. "X_scVI"')

    # core model/training knobs
    p.add_argument("--n-latent", type=int, default=30)
    p.add_argument("--max-epochs", type=int, default=200)
    p.add_argument("--early-stop", action="store_true", help="Enable early stopping")
    p.add_argument("--patience", type=int, default=20, help="Early stopping patience")
    p.add_argument("--min-delta", type=float, default=0.0, help="Early stopping min_delta")
    p.add_argument("--lr", type=float, default=1e-3)

    # performance for CPU
    p.add_argument("--num-workers", type=int, default=0, help="DataLoader workers (WSL: keep 0)")
    p.add_argument("--batch-size", type=int, default=512, help="Minibatch size (CPU-safe)")
    p.add_argument("--seed", type=int, default=0)

    # optional exports
    p.add_argument("--write-normalized", action="store_true",
                  help="Write scVI normalized expression to .layers['scvi_norm']")
    p.add_argument("--norm-library-size", type=float, default=1e4,
                  help="Library size for normalized expression if enabled")

    return p.parse_args()


def main() -> int:
    args = parse_args()
    t0 = time.time()
    started = utc_now_iso()

    in_path = Path(args.input).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    output_h5ad = outdir / "adata_scvi_latent.h5ad"
    model_dir = outdir / "scvi_model"
    manifest_path = outdir / "run_manifest.json"

    set_reproducibility(args.seed)

    # CPU setup
    scvi.settings.seed = args.seed
    scvi.settings.num_threads = max(1, os.cpu_count() or 1)

    # DataLoader worker control (Trainer doesn't accept num_workers)
    scvi.settings.dl_num_workers = args.num_workers
    # optional, safe defaults on WSL/CPU:
    scvi.settings.dl_pin_memory = False

    device = "cpu"
    torch.set_num_threads(scvi.settings.num_threads)

    print(f"Loading: {in_path}")
    adata = ad.read_h5ad(in_path)

    # basic checks
    if args.batch_key not in adata.obs.columns:
        raise KeyError(f"batch-key '{args.batch_key}' not in adata.obs.columns")

    # Ensure batch is a string/categorical that won't explode later
    adata.obs[args.batch_key] = adata.obs[args.batch_key].astype(str)

    # Make obs_names unique defensively (should already be unique after your prefixing)
    if not adata.obs_names.is_unique:
        # keep stable behavior: make unique rather than failing
        adata.obs_names_make_unique()

    # scVI expects raw counts; ensure sparse or ndarray ok
    # If you use a layer, scVI will read counts from that layer.
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key=args.batch_key,
        layer=args.layer,
    )

    # Instantiate model
    model = scvi.model.SCVI(
        adata,
        n_latent=args.n_latent,
    )

    # Training
    print(
        f"Training scVI (device={device}) "
        f"cells={adata.n_obs} genes={adata.n_vars} "
        f"batches={adata.obs[args.batch_key].nunique()} "
        f"max_epochs={args.max_epochs} batch_size={args.batch_size}"
    )

    train_kwargs = dict(
        max_epochs=args.max_epochs,
        accelerator=device,
        devices=1,
        batch_size=args.batch_size,
        plan_kwargs={"lr": args.lr},
        check_val_every_n_epoch=1,
        enable_progress_bar=True,
    )

    if args.early_stop:
        train_kwargs.update(
            dict(
                early_stopping=True,
                early_stopping_patience=args.patience,
                early_stopping_min_delta=args.min_delta,
                early_stopping_monitor="elbo_validation",
                # rule-of-thumb: keep some warmup
            )
        )

    model.train(**train_kwargs)

    # Save model
    model.save(model_dir, overwrite=True)
    print(f"Saved model: {model_dir}")

    # Latent embedding
    latent = model.get_latent_representation()
    adata.obsm[args.latent_key] = latent
    print(f"Wrote latent to adata.obsm['{args.latent_key}'] shape={latent.shape}")

    # Optional normalized expression (can be big in RAM; OK at ~68k x 58k only if you store sparse? It's dense.)
    # For your gene count (57,874), writing dense normalized expression will be HUGE.
    # So default is off; enable only if you later subset genes/cells.
    if args.write_normalized:
        print("Computing normalized expression (this can be large)...")
        norm = model.get_normalized_expression(
            library_size=args.norm_library_size,
            return_numpy=True
        )
        # WARNING: norm is dense float; storing it will balloon disk
        adata.layers["scvi_norm"] = norm
        print("Stored normalized expression in adata.layers['scvi_norm']")

    # Write output h5ad
    adata.write_h5ad(output_h5ad, compression="gzip")
    print(f"Wrote: {output_h5ad}")

    # Manifest
    duration = time.time() - t0
    plan = {
        "batch_size": args.batch_size,
        "lr": args.lr,
        "early_stop": args.early_stop,
        "patience": args.patience,
        "min_delta": args.min_delta,
        "layer": args.layer,
    }
    versions = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "torch": torch.__version__,
        "scvi-tools": scvi.__version__,
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }

    manifest = RunManifest(
        started_utc=started,
        duration_sec=duration,
        input_path=str(in_path),
        output_h5ad=str(output_h5ad),
        model_dir=str(model_dir),
        n_obs=int(adata.n_obs),
        n_vars=int(adata.n_vars),
        batch_key=args.batch_key,
        latent_key=args.latent_key,
        n_latent=args.n_latent,
        max_epochs=args.max_epochs,
        plan=plan,
        versions=versions,
        device=device,
        seed=args.seed,
    )

    manifest_path.write_text(json.dumps(asdict(manifest), indent=2) + "\n")
    print(f"Wrote manifest: {manifest_path}")
    print(f"Done in {duration/60:.2f} min.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

