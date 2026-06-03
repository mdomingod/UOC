# src/project_package/test.py
"""
Held-out test evaluation.

Loads a trained model from an MLflow run (or directly from a state-dict) and
runs the full evaluation suite (site / primer / panel) on the test set.
The test set must be primer-disjoint from train and val — this is the
"evaluate ONCE" step in Andrew Ng's framework.

Usage:
    # Either pass a run_id to download the model from MLflow:
    run_test(cfg, test_data, source={"run_id": "abc123"})

    # Or pass an in-memory model (e.g. right after training):
    run_test(cfg, test_data, source={"model": trained_model})
"""

from __future__ import annotations
from zipfile import Path

from flask import json
import numpy as np
import torch
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import pandas as pd
import matplotlib.pyplot as plt

from .config import Config
from .data import AlignmentDataset
from .evaluate import (
    compute_offtarget_rate,
    level_bias,
    predict,
)
from .utils import (
    MLFLOW_AVAILABLE,
    get_device,
    init_mlflow,
    log_dataset_info,
    log_environment,
    log_git_state,
    log_params_from_cfg,
)

from src.project_package.utils            import (
    get_device,
    init_mlflow,
    log_dataset_info,
    log_environment,
    log_git_state,
    log_params_from_cfg,
    plot_pred_vs_true, 
    plot_residuals, 
    plot_residuals_clip,
    plot_pred_vs_true_log, 
    save,  
)

if MLFLOW_AVAILABLE:
    import mlflow
    import mlflow.pytorch


# ============================================================================
# 1. Model loading
# ============================================================================

def load_model_from_run(run_id: str, device: torch.device) -> torch.nn.Module:
    """Load a PyTorch model that was logged with mlflow.pytorch.log_model."""
    if not MLFLOW_AVAILABLE:
        raise RuntimeError("mlflow is not installed")
    model = mlflow.pytorch.load_model(f"runs:/{run_id}/model", map_location=device)
    model.to(device).eval()
    return model


# ============================================================================
# 2. Test runner
# ============================================================================

def run_test(
    cfg: Config,
    test_data: dict,
    source: dict,
    y_mean: float | None = None,
    y_std: float | None = None,
    run_dir_test: Path | None = None,
    run_name: str | None = None,
    device: str | torch.device | None = None,
) -> dict:
    """
    Evaluate a trained model on the hold-out test set.

    `source` must contain ONE of:
        {"run_id": "<mlflow_run_id>"}  – load model from MLflow
        {"model": <torch.nn.Module>}   – use a model already in memory

    `test_data` is the dict from `data.load_npz_files` restricted to the
    test split (primer-disjoint from train and val).

    Returns the evaluation `report_test` dict.
    """
    device = torch.device(device) if device is not None else get_device(cfg.device)
    init_mlflow(cfg, experiment_name=cfg.experiment_name + "_test")

    # ---- Resolve model -----------------------------------------------------
    if "model" in source:
        model = source["model"].to(device).eval()
        parent_run_id = source.get("run_id")  # optional
    elif "run_id" in source:
        model = load_model_from_run(source["run_id"], device)
        parent_run_id = source["run_id"]
    else:
        raise ValueError("source must contain 'model' or 'run_id'")

    # ---- Predict -----------------------------------------------------------
    images, labels, primer_seqs = test_data["images"], test_data["labels"], test_data["primer_seqs"]
    is_target = test_data.get("is_target")
    test_ds = AlignmentDataset(images, labels, y_mean, y_std)
    site_ids = test_data['site_ids']
    panel_ids = test_data.get("panel_ids")

    y_pred = predict(model, test_ds, batch_size=cfg.batch_size, device=str(device))
    y_true = (labels - test_ds.y_mean.numpy()) / test_ds.y_std.numpy()

    # ---- Compute metrics ---------------------------------------------------
    y_pred_norm = y_pred
    y_true_norm = y_true

    y_true_denorm = test_data['labels']
    y_pred_denorm = y_pred * test_ds.y_std.numpy() + test_ds.y_mean.numpy()
    
    offtarget_rate = compute_offtarget_rate(y_pred_denorm, test_data['is_target']) # return n_offtarget / total if total > 0 else 0.0
    primer_level_bias = level_bias(y_pred_denorm, y_true_denorm) # return a dict

    df_test = pd.DataFrame({
        'site_id': test_data['site_ids'],
        'y_true_norm': y_true_norm,
        'y_pred_norm': y_pred_norm,
        'y_true_denorm': y_true_denorm,
        'y_pred_denorm': y_pred_denorm,
        'is_target': test_data['is_target'],
    })

    df_test['offtarget_rate'] = offtarget_rate
    df_test['mean_bias'] = primer_level_bias['mean_bias']
    df_test['mean_abs_bias'] = primer_level_bias['mean_abs_bias']

    df_test.to_csv(run_dir_test / f"{run_name}_test_validation_predictions.csv", index=False)

    test_metrics = ({
        'r2': r2_score(y_true_denorm,  y_pred_denorm),
        'mse' : mean_squared_error(y_true_denorm,  y_pred_denorm),
        'rmse': np.sqrt(mean_squared_error(y_true_denorm,  y_pred_denorm)),
        'mae': mean_absolute_error(y_true_denorm,  y_pred_denorm),
    })

    with open(run_dir_test / f"{run_name}test_metrics.json", "w") as f:
        json.dump(test_metrics, f, indent=2)

    # denormalised plots
    plot_pred_vs_true_fig = plot_pred_vs_true(y_true_denorm, y_pred_denorm, "Predicted vs true")
    save(plot_pred_vs_true_fig, run_dir_test, f"{run_name}_test_plot_pred_vs_true.png")

    plot_pred_vs_true_log_fig = plot_pred_vs_true_log(y_true_denorm, y_pred_denorm, "Predicted vs true (log scale)")
    save(plot_pred_vs_true_log_fig, run_dir_test, f"{run_name}_test_plot_pred_vs_true_log.png")
    
    plot_residuals_log_fig = plot_residuals_clip(y_true_denorm, y_pred_denorm, "Residuals (log scale)")
    save(plot_residuals_log_fig, run_dir_test, f"{run_name}_test_plot_residuals_log.png")

    plot_residuals_fig = plot_residuals(y_true_denorm, y_pred_denorm, "Residuals")
    save(plot_residuals_fig, run_dir_test, f"{run_name}_test_plot_residuals.png")

    # normalised plots
    plot_pred_vs_true_fig = plot_pred_vs_true(y_true_norm, y_pred_norm, "Predicted vs true (normalised)")
    save(plot_pred_vs_true_fig, run_dir_test, f"{run_name}_test_plot_pred_vs_true_norm.png")

    plot_residuals_fig = plot_residuals(y_true_norm, y_pred_norm, "Residuals (normalised)")
    save(plot_residuals_fig, run_dir_test, f"{run_name}_test_plot_residuals_norm.png")


    print("\n── Hold-out test results ──")
    print(f"  Site RMSE     = {test_metrics['rmse']:.4f}  ")
    print(f"  Site MAE  = {test_metrics['mae']:.4f}")
    print(f"  Site R²   = {test_metrics['r2']:.4f}")
          #f"[{lo_rmse:.4f}, {hi_rmse:.4f}]  ({cfg.bootstrap_confidence:.0%} CI)")
    #print(f"  Panel MAPB    = {report_test['panel_level']['mean_abs_bias']:.4f}")

    # ---- Log everything to a fresh MLflow run ------------------------------
    with mlflow.start_run(run_name=f"TEST_{run_name}") as run:
        mlflow.set_tag("run_type", "test_eval")
        if parent_run_id:
            mlflow.set_tag("source_run_id", parent_run_id)

        log_params_from_cfg(cfg)
        if cfg.log_environment:
            log_environment()
        if cfg.log_git_state:
            log_git_state()
        if cfg.log_dataset_info:
            log_dataset_info(test_data)

        #log_evaluation_to_mlflow(report_test, prefix="test")
        #make_evaluation_plots(report_test, prefix="test")

        
    