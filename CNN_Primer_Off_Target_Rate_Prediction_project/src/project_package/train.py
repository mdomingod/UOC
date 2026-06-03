# src/project_package/train.py

"""
Training loop with full MLflow instrumentation.
 
Design:
  - train_one_run(...) is the inner trainer; it can run inside any MLflow
    run started by the caller (single-run main, CV fold, hyperparam sweep).
  - Reproducibility is handled by `utils.set_global_seed` + a controlled
    DataLoader generator (PyTorch reproducibility docs).
  - All scalar metrics are logged per epoch with `step=epoch` so the MLflow UI
    can compare learning curves across runs (this is what powers `Compare`).
"""

from __future__ import annotations
 
import math
import time
from pathlib import Path
import json
 
from flask import json
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import torch.nn.utils as nn_utils
 
from .config import Config
from .model import PordleCNN
from .utils import (
    MLFLOW_AVAILABLE,
    get_device,
    log_figure,
    log_model_summary,
    make_generator,
    plot_loss_curves,
    plot_lr_schedule,
    seed_worker,
    set_global_seed,
)
 
if MLFLOW_AVAILABLE:
    import mlflow
    import mlflow.pytorch


# ============================================================================
# 1. Per-epoch helpers
# ============================================================================
 
def _train_one_epoch(model, loader, optimizer, loss_fn, device) -> tuple[float, float]:
    """Returns (mean train loss, total grad norm before clipping)."""
    model.train()
    total_loss = 0.0
    total_grad_norm = 0.0

    for batch_idx, (x, y) in enumerate(loader, 1):
        x = x.to(device, non_blocking=True)
        y = y.to(device, non_blocking=True)

        optimizer.zero_grad(set_to_none=True)

        loss = loss_fn(model(x), y)
        loss.backward()

        # Use PyTorch's optimized gradient norm computation (no clipping)
        grad_norm = nn_utils.clip_grad_norm_(model.parameters(), float('inf'))
 
        optimizer.step()

        # Accumulate with detach to prevent graph tracking
        total_loss += loss.detach()
        total_grad_norm += grad_norm


    num_batches = batch_idx
    return float(total_loss / num_batches), float(total_grad_norm / num_batches)


 
 
@torch.no_grad()
def _validate_one_epoch(model, loader, loss_fn, device) -> float:
    model.eval()
    total_loss = 0.0
    
    for batch_idx, (x, y) in enumerate(loader, 1):
        x = x.to(device, non_blocking=True)
        y = y.to(device, non_blocking=True)
        
        loss = loss_fn(model(x), y)
        total_loss += loss.detach()
    
    num_batches = batch_idx
    return float(total_loss / num_batches)


# ============================================================================
# 2. Public training function
# ============================================================================
 
def train_one_run(
    cfg: Config,
    train_ds,
    val_ds,
    device: str | torch.device | None = None,
    use_mlflow: bool = False,
):
    """
    Train a single model.
 
    Returns:
        model:        the trained model with best-validation weights restored
        best_val:     best (lowest) validation MSE seen during training
        stopped_at:   last epoch actually executed
        best_epoch:   epoch where best_val was reached
        history:      dict of per-epoch metrics (train_loss, val_loss, lr, grad_norm)
 
    If `use_mlflow=True`, per-epoch scalars and final artifacts are written to
    the currently active MLflow run. The CALLER starts/ends the run.
    """
    # ---- 0. Reproducibility ------------------------------------------------
    set_global_seed(cfg.seed, deterministic=cfg.deterministic,
                    cudnn_benchmark=cfg.cudnn_benchmark)
 
    device = torch.device(device) if device is not None else get_device(cfg.device)
 
    # ---- 1. DataLoaders with controlled randomness -------------------------
    generator = make_generator(cfg.seed)
 
    train_loader = DataLoader(
        train_ds,
        batch_size=cfg.batch_size,
        shuffle=True,
        num_workers=cfg.n_workers,
        worker_init_fn=seed_worker,   # seeds NumPy/random inside workers
        generator=generator,          # controls shuffle order
        pin_memory=cfg.pin_memory and device.type == "cuda",
        persistent_workers=cfg.n_workers > 0,
    )
    val_loader = DataLoader(
        val_ds,
        batch_size=cfg.batch_size,
        shuffle=False,
        num_workers=cfg.n_workers,
        worker_init_fn=seed_worker,
        pin_memory=cfg.pin_memory and device.type == "cuda",
        persistent_workers=cfg.n_workers > 0,
    )
 
    # ---- 2. Model / loss / optimizer ---------------------------------------
    model = PordleCNN(cfg).to(device)
    #model = torch.compile(model)  # optional, requires PyTorch 2.0+ and compatible hardware
    loss_fn = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg.learning_rate)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.1, patience=cfg.lr_patience
    )
 
    if use_mlflow and MLFLOW_AVAILABLE:
        log_model_summary(model)
        mlflow.log_param("device", str(device))
 
    # Sanity check: what does an untrained model predict?
    model.eval()
    with torch.no_grad():
        sample_x, _ = next(iter(train_loader))
        raw = model(sample_x.to(device))
        print(
            f"Pre-training preds: min={raw.min():.4f}, max={raw.max():.4f}, "
            f"Mean ± Std:    {raw.mean():.4f} ± {raw.std():.4f}"
            f"%zero={((raw == 0).float().mean() * 100):.1f}%"
        )
 
    # ---- 3. Training loop ---------------------------------------------------
    history = {"train_loss": [], "val_loss": [], "lr": [], "grad_norm": []}
    best_val, best_state, best_epoch = float("inf"), None, 0
    epochs_no_improve, stopped_at = 0, cfg.epochs
    t0 = time.time()
 
    for epoch in range(1, cfg.epochs + 1):
        train_loss, grad_norm = _train_one_epoch(model, train_loader, optimizer, loss_fn, device)
        val_loss              = _validate_one_epoch(model, val_loader,   loss_fn, device)

        scheduler.step(val_loss)
        current_lr = optimizer.param_groups[0]["lr"]
 
        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)
        history["lr"].append(current_lr)
        history["grad_norm"].append(grad_norm)
 
        if use_mlflow and MLFLOW_AVAILABLE:
            
            mlflow.log_metrics(
                {
                    "train_mse": train_loss,
                    "val_mse":   val_loss,
                    "val_rmse":  float(np.sqrt(val_loss)),
                    "lr":        current_lr,
                    "grad_norm": grad_norm,
                },
                step=epoch,
            )
 
        print(f"Epoch {epoch:3d} | train loss MSE {train_loss:9.2f} "
              f"| val loss MSE {val_loss:9.2f} | lr {current_lr:.2e}")
 
        # ---- 4. Early stopping --------------------------------------------
        if val_loss < best_val:
            best_val   = val_loss
            best_state = {k: v.detach().cpu().clone() 
                          for k, v in model.state_dict().items()}
            best_epoch = epoch
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1
            if epochs_no_improve >= cfg.patience:
                stopped_at = epoch
                print(f"Early stopping at epoch {epoch}: no improvement for "
                      f"{epochs_no_improve} epochs")
                break
 
        stopped_at = epoch
 
    elapsed = time.time() - t0
 
    # ---- 5. Restore best weights -------------------------------------------
    if best_state is not None:
        model.load_state_dict(best_state)
 
    history["best_epoch"] = best_epoch
 
    # ---- 6. Log final summary, plots, and model artifact -------------------
    if use_mlflow and MLFLOW_AVAILABLE:

        

        mlflow.log_metrics({
            "best_val_mse":  float(best_val),
            "best_val_rmse": float(np.sqrt(best_val)),
            "best_epoch":    float(best_epoch),
            "stopped_at":    float(stopped_at),
            "training_time_s": float(elapsed),
        })
        log_figure(plot_loss_curves(history,
                                    title=f"Loss curves — {cfg.run_name}"),
                   "plots/loss_curves.png")
        log_figure(plot_lr_schedule(history),
                   "plots/lr_schedule.png")
 
        if cfg.log_model:
            # Provide an input example so MLflow can infer the model signature.
            # This makes the logged model reusable via `mlflow.pytorch.load_model`.
            example_x = next(iter(val_loader))[0][:2].cpu().numpy()
            mlflow.pytorch.log_model(
                pytorch_model=model,
                artifact_path="model",
                input_example=example_x,
            )
 
    return model, best_val, stopped_at, best_epoch, history