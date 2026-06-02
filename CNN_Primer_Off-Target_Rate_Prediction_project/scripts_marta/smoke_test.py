#!/usr/bin/env python
# To run the srcipt, activate the uv environment and change th end line from CRLF to LF (in vscode, bottom right corner) and then run:
# smoke_test.py
from dataclasses import replace
import torch
import mlflow

from src.project_package.config import Config
from src.project_package.data import (
    load_npz_files, AlignmentDataset, split_by_primer
)
from src.project_package.train import train_one_run


def main():
    # 1) Override a few config values for a *fast* run.
    cfg = Config()
    cfg = replace(
        cfg,
        epochs=10,          # ceiling — early stopping should fire well before
        patience=4,         # don't wait 10 epochs to early-stop on a smoke test
        lr_patience=2,      # let the scheduler kick in faster too
        batch_size=64,      # smaller batches = more gradient updates per epoch
        run_name="smoke_test",
    )
    print("Config:", cfg)

    # 2) Pick 5 npz files. Just take the first 5 alphabetically for reproducibility.
    npz_paths = sorted(cfg.data_dir.glob("*.npz"))[:5]
    assert len(npz_paths) == 5, f"Expected 5 npz files, found {len(npz_paths)}"
    print(f"Loading from: {[p.name for p in npz_paths]}")

    data = load_npz_files(npz_paths)
    print(f"Loaded {len(data['images'])} sites across "
          f"{len(set(data['primer_ids']))} unique primers")
    print(f"Image shape: {data['images'].shape}, dtype: {data['images'].dtype}")
    print(f"Label range: [{data['labels'].min():.2f}, {data['labels'].max():.2f}], "
          f"mean: {data['labels'].mean():.2f}")

    # 3) Single primer-level split — no K-fold for a smoke test.
    train_idx, val_idx = split_by_primer(
        data['primer_ids'], val_fraction=cfg.val_fraction, seed=cfg.seed
    )
    print(f"Train sites: {len(train_idx)}, Val sites: {len(val_idx)}")

    train_ds = AlignmentDataset(data['images'][train_idx], data['labels'][train_idx])
    val_ds   = AlignmentDataset(data['images'][val_idx],   data['labels'][val_idx])

    # 4) Sanity-check one sample before training. This catches shape/dtype bugs
    #    that would otherwise only surface mid-epoch.
    x_sample, y_sample = train_ds[0]
    print(f"Sample x: shape={tuple(x_sample.shape)}, dtype={x_sample.dtype}, "
          f"range=[{x_sample.min():.3f}, {x_sample.max():.3f}]")
    print(f"Sample y: {y_sample.item():.4f}")

    # 5) Train with MLflow tracking so you can inspect curves later.
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    mlflow.set_experiment(cfg.experiment_name)
    with mlflow.start_run(run_name=cfg.run_name) as run:
        mlflow.log_params({
            k: str(v) for k, v in cfg.__dict__.items()
        })
        model, best_val = train_one_run(
            cfg, train_ds, val_ds, device=device, mlflow_run=run
        )
        print(f"\n✅ Smoke test complete. Best val MSE: {best_val:.4f}")


if __name__ == "__main__":
    main()