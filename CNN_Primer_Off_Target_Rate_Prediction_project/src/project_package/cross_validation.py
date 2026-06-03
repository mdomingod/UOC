# src/project_package/cross_validation.py
"""
K-fold cross-validation with nested MLflow runs.

Layout in the MLflow UI:
  experiment: <cfg.experiment_name>
    └── parent run "CV_<cfg.run_name>"   (tag: run_type=cv_parent)
          ├── child run "fold_0"          (tag: run_type=cv_child, fold=0)
          ├── child run "fold_1"
          └── ...

The parent run holds aggregated metrics (mean / std across folds).
Each child run holds per-fold training history and evaluation.

This is the layout Andrew Ng's specialization recommends for honest
hyperparameter selection: tune on CV folds, then evaluate ONCE on a
held-out test set (handled in test.py).
"""

from __future__ import annotations

import numpy as np
import torch
import time

from .config import Config
from .data import AlignmentDataset, kfold_by_primer_seq, kfold_by_panel
from .evaluate_complex import (
    evaluate_predictions,
    log_evaluation_to_mlflow,
    make_evaluation_plots,
    predict,
)
from .train import train_one_run
from .utils import (
    MLFLOW_AVAILABLE,
    init_mlflow,
    log_dataset_info,
    log_environment,
    log_git_state,
    log_params_from_cfg,
)

if MLFLOW_AVAILABLE:
    import mlflow


def run_cv(cfg: Config, data: dict, device: str | torch.device | None = None) -> dict:
    """
    Run K-fold CV at the primer level. The held-out test set should NOT be in
    `data` — use split_train_val_test upstream and pass only the train+val portion here.

    Returns a dict with per-fold and aggregated metrics.
    """
    init_mlflow(cfg, experiment_name=cfg.experiment_name + "_cv")

    images = data["images"]
    labels = data["labels"]
    primer_seqs = data["primer_seqs"]
    panel_ids = data["panel_ids"]
    is_target = data.get("is_target")

    fold_results: list[dict] = []

    # ---- Parent run --------------------------------------------------------
    with mlflow.start_run(run_name=f"CV_{cfg.run_name}") as parent:
        mlflow.set_tag("run_type", "cv_parent")
        log_params_from_cfg(cfg)

        if cfg.log_environment:
            log_environment()
        if cfg.log_git_state:
            log_git_state()
        if cfg.log_dataset_info:
            log_dataset_info(data)

        print(f"Parent CV run started: {parent.info.run_id}")

        # ---- Per-fold child runs ------------------------------------------
        print("Splitting at kfold by panel...")
        for fold_idx, train_idx, val_idx in kfold_by_panel(
            panel_ids, n_folds=cfg.n_folds, seed=cfg.seed
        ):
            with mlflow.start_run(run_name=f"fold_{fold_idx}", nested=True) as child:
                mlflow.set_tags({"run_type": "cv_child", "fold": str(fold_idx), "fold_run_id": child.info.run_id})
                log_params_from_cfg(cfg)
                mlflow.log_param("fold", fold_idx)

                if cfg.log_dataset_info:
                    log_dataset_info(
                        data,
                        split_indices={"train": train_idx, "val": val_idx},
                    )

                # ---- Verify primer or panel_id disjointness — mandatory check ----

                train_panel = set(panel_ids[train_idx])
                val_panel = set(panel_ids[val_idx])
                overlap = train_panel & val_panel

                assert len(overlap) == 0, (
                    f"Panel ID leakage detected: {len(overlap)} IDs "
                    f"appear in both train and val. Split is broken."
                )

                print(f"  Unique primers in train: {len(train_panel):,}")
                print(f"  Unique primers in val:   {len(val_panel):,}")
                print(f"  Panel overlap:          0 ✓")

                # Compute normalisation stats from training labels only
                labels = data["labels"].astype(np.float32)
                train_labels = labels[train_idx]

                val_labels   = data["labels"][val_idx].astype(np.float32)

                print(f"Label Train — mean: {train_labels.mean():.2f}, "
                    f"std: {train_labels.std():.2f}, "
                    f"zeros: {(train_labels==0).mean():.1%}, "
                    f"max: {train_labels.max():.0f}, "
                    f"p99: {np.percentile(train_labels, 99):.0f} ")

                print(f"Label Val   — mean: {val_labels.mean():.2f}, "
                    f"std: {val_labels.std():.2f}, "
                    f"zeros: {(val_labels==0).mean():.1%}, "
                    f"max: {val_labels.max():.0f}, "
                    f"p99: {np.percentile(val_labels, 99):.0f} ")
                
                # Data normalisation
                y_mean = float(train_labels.mean())
                y_std  = float(train_labels.std())

                val_norm = (val_labels - y_mean) / (y_std + 1e-8)
                zero_pred_val_mse = float((val_norm ** 2).mean())
                print(f"\nExpected val MSE for zero-predictor: {zero_pred_val_mse:.4f}")

                train_ds = AlignmentDataset(images[train_idx], train_labels, y_mean, y_std)
                val_ds   = AlignmentDataset(images[val_idx],   labels[val_idx], y_mean, y_std)

                print(f"\n── Fold {fold_idx} | "
                      f"train: {len(train_idx)} sites / "
                      f"val: {len(val_idx)} sites ──")
                
                print("[3/5] Training one run...")
                init_mlflow(cfg)
                mlflow.set_experiment(cfg.experiment_name)
                with mlflow.start_run(run_name=cfg.run_name) as run:
                    with mlflow.start_run(nested=True, run_name='metadata'):
                        mlflow.set_tag("run_type", "single_run")
                        log_params_from_cfg(cfg)
                        if cfg.log_environment:
                            log_environment()
                        if cfg.log_git_state:
                            log_git_state()
                        if cfg.log_dataset_info:
                            log_dataset_info(
                                data,
                                split_indices={"train": train_idx, "val": val_idx}) #, "test": test_idx},
                            #)

                t0 = time.time()
                model, best_val, stopped_at, best_epoch, _ = train_one_run(
                    cfg, train_ds, val_ds, device=device, use_mlflow=True
                )
                duration = time.time() - t0
                print(f"  Training finished in {duration:.1f}s")
                print(f"  Stopped at epoch {stopped_at}, best val MSE {best_val:.2f} at epoch {best_epoch}")
                print()

                # Log training summary as metrics
                mlflow.log_metrics({
                    "training_duration_sec": duration,
                    "best_val_mse":          best_val,
                    "stopped_at_epoch":      stopped_at,
                    "best_epoch":            best_epoch,
                })

                # ---- Evaluate this fold on its validation set ----------
                print("[4/5] Evaluating on validation set...")
                val_preds = predict(model, val_ds, batch_size=cfg.batch_size, device=device)

                y_true_norm = data['labels'][val_idx] - val_ds.y_mean.numpy() / val_ds.y_std.numpy()
                y_pred_norm = val_preds

                y_true_denorm = data['labels'][val_idx]
                y_pred_denorm = val_preds * val_ds.y_std.numpy() + val_ds.y_mean.numpy()
    

                report = evaluate_predictions(
                    y_pred_denorm, 
                    y_true_denorm, 
                    data['primer_seqs'][val_idx], 
                    data['is_target'][val_idx], 
                    data['panel_ids'][val_idx],
                )
                log_evaluation_to_mlflow(report, prefix=f"fold{fold_idx}_val")
                make_evaluation_plots(report, prefix=f"fold{fold_idx}_val")

                fold_results.append({
                    "fold":       fold_idx,
                    "best_val":   best_val,
                    "best_epoch": best_epoch,
                    "stopped_at": stopped_at,
                    "site_level":   report["site_level"],
                    "primer_level": report["primer_level"],
                })

        # ---- Aggregate to the parent run ---------------------------------
        agg = _aggregate_folds(fold_results)
        mlflow.log_metrics(agg)
        mlflow.log_dict({"folds": fold_results, "aggregate": agg},
                        "cv_summary.json")

        print("\n── CV summary ──")
        for k, v in agg.items():
            print(f"  {k:30s} = {v:.4f}")

    return {"folds": fold_results, "aggregate": agg}


def _aggregate_folds(fold_results: list[dict]) -> dict:
    """Compute mean ± std across folds for every numeric metric."""
    agg: dict[str, float] = {}

    # best_val across folds
    best_vals = [f["best_val"] for f in fold_results]
    agg["cv_best_val_mean"] = float(np.mean(best_vals))
    agg["cv_best_val_std"]  = float(np.std(best_vals))

    # site-level metrics
    for key in fold_results[0]["site_level"]:
        vals = [f["site_level"][key] for f in fold_results
                if isinstance(f["site_level"][key], (int, float))
                and not (isinstance(f["site_level"][key], float) and np.isnan(f["site_level"][key]))]
        if not vals:
            continue
        agg[f"cv_site_{key}_mean"] = float(np.mean(vals))
        agg[f"cv_site_{key}_std"]  = float(np.std(vals))

    # primer-level metrics
    for key in fold_results[0]["primer_level"]:
        vals = [f["primer_level"][key] for f in fold_results
                if isinstance(f["primer_level"][key], (int, float))
                and not (isinstance(f["primer_level"][key], float) and np.isnan(f["primer_level"][key]))]
        if not vals:
            continue
        agg[f"cv_primer_{key}_mean"] = float(np.mean(vals))
        agg[f"cv_primer_{key}_std"]  = float(np.std(vals))

    return agg