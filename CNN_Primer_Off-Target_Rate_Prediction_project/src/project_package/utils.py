# src/project_package/utils.py
"""
Shared utilities: reproducibility, device detection, MLflow logging helpers,
and plotting. Keep this file framework-agnostic where possible.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import random
import subprocess
import sys
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import torch
import shutil 

try:
    import mlflow
    MLFLOW_AVAILABLE = True
except ImportError:
    MLFLOW_AVAILABLE = False


# ============================================================================
# 1. REPRODUCIBILITY  (PyTorch reproducibility docs)
# ============================================================================

def set_global_seed(seed: int, deterministic: bool = True, cudnn_benchmark: bool = False) -> None:
    """
    Seed every RNG used in the training stack.

    Andrew Ng's ML specialization stresses that reproducible experiments are
    the foundation of any iterative model-development loop. PyTorch's own
    docs (https://pytorch.org/docs/stable/notes/randomness.html) require:
      - seed Python's `random`, NumPy, and torch on both CPU and CUDA
      - disable cudnn benchmarking
      - enable deterministic algorithms
      - set CUBLAS_WORKSPACE_CONFIG for full CUDA determinism
    """
    os.environ["PYTHONHASHSEED"] = str(seed)
    # Required by torch for deterministic CUDA matmul; harmless on CPU/MPS.
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    torch.backends.cudnn.benchmark     = cudnn_benchmark
    torch.backends.cudnn.deterministic = deterministic

    if deterministic:
        # warn_only=True: some ops (e.g. some pooling backward) have no det
        # implementation; we warn rather than crash so the CNN still trains.
        torch.use_deterministic_algorithms(True, warn_only=True)


def seed_worker(worker_id: int) -> None:
    """
    Re-seed each DataLoader worker.
    PyTorch only seeds its own RNG in workers, so NumPy and `random` must
    be reseeded manually here.
    """
    worker_seed = torch.initial_seed() % 2**32
    np.random.seed(worker_seed)
    random.seed(worker_seed)


def make_generator(seed: int) -> torch.Generator:
    """Generator that controls DataLoader shuffle order."""
    g = torch.Generator()
    g.manual_seed(seed)
    return g


# ============================================================================
# 2. DEVICE
# ============================================================================

def get_device(preferred: str = "auto") -> torch.device:
    """Resolve 'auto' / 'cpu' / 'cuda' / 'mps' to a torch.device."""
    if preferred == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    return torch.device(preferred)


# ============================================================================
# 3. MLFLOW: setup and parameter logging
# ============================================================================

def init_mlflow(cfg, experiment_name: str | None = None) -> None:
    """Set tracking URI + experiment. Safe to call multiple times."""
    if not MLFLOW_AVAILABLE:
        return
    if cfg.mlflow_tracking_uri:
        mlflow.set_tracking_uri(cfg.mlflow_tracking_uri)
    mlflow.set_experiment(experiment_name or cfg.experiment_name)


def log_params_from_cfg(cfg) -> None:
    """Log every config field as an MLflow param (Paths -> strings)."""
    if not MLFLOW_AVAILABLE:
        return
    flat = cfg.to_dict()
    # `tags` and `npz_files` are special — log them differently
    tags = flat.pop("tags", {}) or {}
    npz_files = flat.pop("npz_files", []) or []

    # MLflow params must be strings/numbers. Cast complex values.
    safe = {}
    for k, v in flat.items():
        if isinstance(v, (list, dict)):
            safe[k] = json.dumps(v, default=str)
        else:
            safe[k] = v
    mlflow.log_params(safe)

    if tags:
        mlflow.set_tags({str(k): str(v) for k, v in tags.items()})

    if npz_files:
        mlflow.log_text(
            "\n".join(str(p) for p in npz_files),
            "inputs/npz_files.txt",
        )


# ============================================================================
# 4. MLFLOW: environment, git, dataset provenance
# ============================================================================

def log_environment() -> None:
    """
    Log everything needed to reproduce the Python environment.

    For a uv-managed project, `uv.lock` is the authoritative artifact: it
    pins every transitive dependency with hashes and can be replayed
    exactly with `uv sync`. We also log pyproject.toml (the declared
    constraints), `uv pip freeze` (current resolved state), and basic
    Python/Torch info.
    """
    if not MLFLOW_AVAILABLE:
        return

    # ---- Python / Torch info ------------------------------------------
    info = {
        "python_version": sys.version,
        "platform":       platform.platform(),
        "torch_version":  torch.__version__,
        "cuda_available": torch.cuda.is_available(),
        "cuda_version":   torch.version.cuda,
    }
    mlflow.log_dict(info, "env/python_info.json")

    # Find project root: src/project_package/utils.py -> root
    project_root = Path(__file__).resolve().parents[2]

    # ---- 1. uv.lock (authoritative reproducibility artifact) ----------
    uv_lock = project_root / "uv.lock"
    if uv_lock.exists():
        mlflow.log_artifact(str(uv_lock), artifact_path="env")

    # ---- 2. pyproject.toml (declared dependency constraints) ----------
    pyproject = project_root / "pyproject.toml"
    if pyproject.exists():
        mlflow.log_artifact(str(pyproject), artifact_path="env")

    # ---- 3. Resolved package list -------------------------------------
    # Prefer `uv pip freeze`; fall back to `python -m pip freeze` if uv
    # isn't on PATH (e.g. when running a packaged wheel outside the repo).
    freeze_text: str | None = None
    if shutil.which("uv"):
        try:
            freeze_text = subprocess.check_output(
                ["uv", "pip", "freeze"], text=True, timeout=30
            )
        except (subprocess.SubprocessError, OSError):
            pass
    if freeze_text is None:
        try:
            freeze_text = subprocess.check_output(
                [sys.executable, "-m", "pip", "freeze"], text=True, timeout=30
            )
        except (subprocess.SubprocessError, OSError):
            pass
    if freeze_text:
        mlflow.log_text(freeze_text, "env/pip_freeze.txt")

    # ---- 4. uv version itself (so the lock-file producer is recorded) -
    if shutil.which("uv"):
        try:
            uv_version = subprocess.check_output(
                ["uv", "--version"], text=True, timeout=5
            ).strip()
            mlflow.set_tag("uv.version", uv_version)
        except (subprocess.SubprocessError, OSError):
            pass


def log_git_state() -> None:
    """Log current commit SHA + dirty-flag. Skip silently if not in a git repo."""
    if not MLFLOW_AVAILABLE:
        return
    try:
        sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, timeout=5
        ).strip()
        dirty = subprocess.check_output(
            ["git", "status", "--porcelain"], text=True, timeout=5
        ).strip()
        mlflow.set_tag("git.commit", sha)
        mlflow.set_tag("git.dirty", str(bool(dirty)))
        if dirty:
            mlflow.log_text(dirty, "env/git_status.txt")
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        pass


def hash_array(arr: np.ndarray) -> str:
    """Stable hash of a NumPy array (first 12 hex chars of SHA256)."""
    return hashlib.sha256(arr.tobytes()).hexdigest()[:12]


def log_dataset_info(data: dict, split_indices: dict | None = None) -> None:
    """
    Log shapes, label stats, primer counts, and a content hash of each array.
    `data` is the dict returned by `load_npz_files`.
    `split_indices` (optional) is {'train': idx, 'val': idx, 'test': idx}.
    """
    if not MLFLOW_AVAILABLE:
        return

    images, labels = data["images"], data["labels"]
    info: dict[str, Any] = {
        "n_sites":         int(len(images)),
        "image_shape":     list(images.shape[1:]),
        "n_unique_primers": int(len(np.unique(data["primer_seqs"]))),
        "label_min":       float(labels.min()),
        "label_max":       float(labels.max()),
        "label_mean":      float(labels.mean()),
        "label_zero_frac": float((labels == 0).mean()),
        "hash_images":     hash_array(images),
        "hash_labels":     hash_array(labels),
        "hash_site_ids": hash_array(data["site_ids"]),
    }
    if "is_target" in data:
        info["frac_target"] = float(np.asarray(data["is_target"]).mean())

    mlflow.log_dict(info, "data/dataset_info.json")

    if split_indices:
        split_info = {
            name: {
                "n_sites":   int(len(idx)),
                "n_primers": int(len(np.unique(data["primer_seqs"][idx]))),
            }
            for name, idx in split_indices.items()
        }
        mlflow.log_dict(split_info, "data/split_info.json")


def log_model_summary(model: torch.nn.Module) -> None:
    """Log architecture as text + parameter count."""
    if not MLFLOW_AVAILABLE:
        return
    mlflow.log_text(str(model), "model/architecture.txt")
    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    mlflow.log_param("n_parameters", n_params)


# ============================================================================
# 5. PLOTS  (return fig so caller can also save/display)
# ============================================================================
def save(fig, run_dir, name, subdir=None):
    target = run_dir / subdir if subdir else run_dir
    #run_dir.mkdir(exist_ok=True)
    target.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(target / name, dpi=150)
    plt.close(fig)
    print("saved", target/ name)

def plot_loss_curves(history: dict, title: str = "Training history") -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7, 4))
    epochs = range(1, len(history["train_loss"]) + 1)
    ax.plot(epochs, history["train_loss"], label="train", linewidth=1.5)
    ax.plot(epochs, history["val_loss"],   label="val",   linewidth=1.5)
    if "best_epoch" in history:
        ax.axvline(history["best_epoch"], color="green", linestyle="--",
                   alpha=0.6, label=f"best (epoch {history['best_epoch']})")
    ax.set_xlabel("epoch"); ax.set_ylabel("MSE loss")
    ax.set_title(title); ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def plot_lr_schedule(history: dict) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7, 3))
    epochs = range(1, len(history["lr"]) + 1)
    ax.plot(epochs, history["lr"], linewidth=1.5, color="purple")
    ax.set_yscale("log")
    ax.set_xlabel("epoch"); ax.set_ylabel("learning rate")
    ax.set_title("Learning rate schedule"); ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def plot_pred_vs_true(y_true: np.ndarray, y_pred: np.ndarray,
                      title: str = "Predicted vs true",
                      xlabel: str = "true", ylabel: str = "predicted") -> plt.Figure:
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(y_true, y_pred, s=4, alpha=0.3, edgecolors="none")
    lo = float(min(y_true.min(), y_pred.min()))
    hi = float(max(y_true.max(), y_pred.max()))
    ax.plot([lo, hi], [lo, hi], "k--", linewidth=1, label="y = x")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title)
    ax.legend(); ax.grid(alpha=0.3); ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    return fig


def plot_pred_vs_true_log(y_true: np.ndarray, y_pred: np.ndarray, 
                          title: str = "Predicted vs true (log scale)") -> plt.Figure:
    eps = 1e-6
    y_true = y_true + eps
    y_pred = y_pred + eps

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(y_true, y_pred, s=5, alpha=0.3)

    ax.set_xscale("log")
    ax.set_yscale("log")

    lo = min(y_true.min(), y_pred.min())
    hi = max(y_true.max(), y_pred.max())
    ax.plot([lo, hi], [lo, hi], "k--", linewidth=1)

    ax.set_xlabel("true")
    ax.set_ylabel("predicted")
    ax.set_title(title)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    return fig


def plot_residuals(y_true: np.ndarray, y_pred: np.ndarray,
                   title: str = "Residuals") -> plt.Figure:
    residuals = y_pred - y_true
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].scatter(y_true, residuals, s=4, alpha=0.3, edgecolors="none")
    axes[0].axhline(0, color="k", linewidth=1)
    axes[0].set_xlabel("true"); axes[0].set_ylabel("residual (pred - true)")
    axes[0].set_title("Residuals vs true"); axes[0].grid(alpha=0.3)
    axes[1].hist(residuals, bins=60, edgecolor="black", linewidth=0.5)
    axes[1].axvline(0, color="k", linewidth=1)
    axes[1].set_xlabel("residual"); axes[1].set_ylabel("count")
    axes[1].set_title("Residual distribution"); axes[1].grid(alpha=0.3)
    fig.suptitle(title); fig.tight_layout()
    return fig


def plot_residuals_clip(y_true: np.ndarray, y_pred: np.ndarray,
                        title: str = "Residuals (clipped)") -> plt.Figure:
    residuals = (y_pred - y_true) / (y_true + 1e-6)
    residuals = np.clip(residuals, -10, 10)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].scatter(y_true, residuals, s=4, alpha=0.3, edgecolors="none")
    axes[0].axhline(0, color="k", linewidth=1)
    axes[0].set_xlabel("true"); axes[0].set_ylabel("residual (pred - true)")
    axes[0].set_title("Residuals vs true"); axes[0].grid(alpha=0.3)
    axes[1].hist(residuals, bins=60, edgecolor="black", linewidth=0.5)
    axes[1].set_xscale("symlog")
    axes[1].axvline(0, color="k", linewidth=1)
    axes[1].set_xlabel("residual"); axes[1].set_ylabel("count")
    axes[1].set_title("Residual distribution"); axes[1].grid(alpha=0.3)
    fig.suptitle(title); fig.tight_layout()
    return fig


def plot_panel_offtarget_comparison(panel_df) -> plt.Figure:
    """Bar plot of predicted vs true off-target rate per panel."""
    import pandas as pd  # local: only needed here
    valid = panel_df.dropna(subset=["pred_offtarget_rate", "true_offtarget_rate"])
    x = np.arange(len(valid))
    fig, ax = plt.subplots(figsize=(max(7, len(valid) * 0.6), 4))
    width = 0.4
    ax.bar(x - width/2, valid["true_offtarget_rate"], width, label="true")
    ax.bar(x + width/2, valid["pred_offtarget_rate"], width, label="pred")
    ax.set_xticks(x); ax.set_xticklabels(valid["panel"].astype(str), rotation=45, ha="right")
    ax.set_ylabel("off-target rate"); ax.set_title("Panel-level off-target rate")
    ax.legend(); ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    return fig

def _save(fig, run_dir, name, subdir=None):
    target = run_dir / subdir if subdir else run_dir
    #run_dir.mkdir(exist_ok=True)
    target.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(target / name, dpi=150)
    plt.close(fig)
    print("saved", target/ name)


def log_figure(fig: plt.Figure, artifact_path: str) -> None:
    """Log a figure to MLflow then close it."""
    if MLFLOW_AVAILABLE:
        mlflow.log_figure(fig, artifact_path)
    plt.close(fig)