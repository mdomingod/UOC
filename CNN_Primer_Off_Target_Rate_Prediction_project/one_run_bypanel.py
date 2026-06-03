#!/usr/bin/env python3
# To run the srcipt, activate the uv environment and change th end line from CRLF to LF (in vscode, bottom right corner) and then run:
# one_run_main.py

"""
Single-run pipeline:
    1. Load .npz files
    2. Primer-disjoint split into train / val / test
    3. Train one model, logging everything to MLflow
    4. Evaluate on the held-out test set in a separate MLflow run
 
This is the discipline Andrew Ng's ML specialization recommends:
TUNE on validation, EVALUATE ONCE on test.
"""

import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import torch
import mlflow
import random
import pandas as pd
import matplotlib.pyplot as plt


from torchgen import model

from src.project_package.config           import Config
from src.project_package.data             import (
    AlignmentDataset,
    load_npz_files,
    split_by_panel,
    split_by_site_id,
    split_by_primer_seq,
    kfold_by_primer_seq
)
from src.project_package.train            import train_one_run
from src.project_package.test             import run_test
from src.project_package.utils            import (
    get_device,
    init_mlflow,
    log_dataset_info,
    log_environment,
    log_git_state,
    log_params_from_cfg,
    plot_pred_vs_true, 
    plot_residuals, 
    plot_panel_offtarget_comparison,
    plot_loss_curves, 
    save, 
    plot_pred_vs_true_log, 
    plot_residuals_clip, 
)

from src.project_package.evaluate          import (
    compute_offtarget_rate,
    level_bias,
    predict
)

from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from scipy.stats import pearsonr, spearmanr



# ============================================================================
# 1. SETUP
# ============================================================================

def setup_run_directory(cfg: Config) -> Path:
    """Create a timestamped folder for this run's outputs."""
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    run_dir = cfg.reports_dir / f"{timestamp}_{cfg.run_name}"
    run_dir.mkdir(exist_ok=True, parents=True)
    run_dir_test = cfg.report_test_dir / f"{timestamp}_{cfg.run_name}"
    run_dir_test.mkdir(exist_ok=True, parents=True)
    return run_dir, run_dir_test

def set_global_seed(seed: int):
    """Make every source of randomness deterministic."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def setup_mlflow(cfg: Config):
    """Configure MLflow tracking location and experiment."""
    if cfg.mlflow_tracking_uri:
        mlflow.set_tracking_uri(cfg.mlflow_tracking_uri)
    mlflow.set_experiment(cfg.experiment_name)


class TeeOutput:
    """Write stdout to both terminal AND a log file."""
    def __init__(self, log_path: Path):
        self.log_path = log_path
        self.terminal = sys.stdout
        self.log = open(log_path, "w", buffering=1)  # line-buffered
    def write(self, msg):
        self.terminal.write(msg)
        self.log.write(msg)
    def flush(self):
        self.terminal.flush()
        self.log.flush()


# ============================================================================
# 2. REPORT WRITING
# ============================================================================

def render_markdown_report(cfg: Config, report: dict, train_info: dict) -> str:
    """Build a human-readable Markdown summary of the run."""
    md = []
    md.append(f"# Training Run Report — {cfg.run_name}")
    md.append(f"")
    md.append(f"_Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}_")
    md.append(f"")

    # Training summary
    md.append(f"## Training")
    md.append(f"")
    md.append(f"| Item | Value |")
    md.append(f"|---|---|")
    md.append(f"| Epochs run        | {train_info['epochs_run']} |")
    md.append(f"| Training time     | {train_info['duration_sec']:.1f} s |")
    md.append(f"| Best val MSE      | {train_info['best_val_mse']:,.2f} |")
    md.append(f"| Train samples     | {train_info['n_train']:,} |")
    md.append(f"| Val samples       | {train_info['n_val']:,} |")
    md.append(f"")

    # Site-level metrics
    s = report['site_level']
    md.append(f"## Site-level metrics")
    md.append(f"")
    md.append(f"How well does the model predict individual UMI counts?")
    md.append(f"")
    md.append(f"| Metric | Value |")
    md.append(f"|---|---|")
    md.append(f"| n_sites     | {s['n_sites']:,} |")
    md.append(f"| MSE         | {s['mse']:,.2f} |")
    md.append(f"| RMSE        | {s['rmse']:,.2f} |")
    md.append(f"| MAE         | {s['mae']:,.2f} |")
    #md.append(f"| Pearson r   | {s['pearson']:.4f} |")
    #md.append(f"| Spearman r  | {s['spearman']:.4f} |")
    md.append(f"")

    # Primer-level metrics
    p = report['primer_level']
    md.append(f"## Primer-level metrics")
    md.append(f"")
    md.append(f"Can the model rank primers by off-target risk?")
    md.append(f"")
    md.append(f"| Metric | Value |")
    md.append(f"|---|---|")
    md.append(f"| n_primers              | {p['n_primers']:,} |")
    md.append(f"| n_dropped (zero total) | {p['n_dropped_zero_total']:,} |")
    md.append(f"| Mean bias              | {p['mean_bias']:+.4f} |")
    md.append(f"| Mean abs bias          | {p['mean_abs_bias']:.4f} |")
    md.append(f"| RMSE on rates          | {p['rmse']:.4f} |")
    #md.append(f"| Pearson r              | {p['pearson']:.4f} |")
    #md.append(f"| Spearman r             | {p['spearman']:.4f} |")
    md.append(f"")

    # Panel-level metrics (optional)
    if 'panel_level' in report:
        pn = report['panel_level']
        md.append(f"## Panel-level metrics")
        md.append(f"")
        md.append(f"The headline metric from the QIAGEN paper.")
        md.append(f"")
        md.append(f"| Metric | Value |")
        md.append(f"|---|---|")
        md.append(f"| n_panels       | {pn['n_panels']} |")
        md.append(f"| Mean bias      | {pn['mean_bias']:+.4%} |")
        md.append(f"| MAPB           | {pn['mean_abs_bias']:.4%} |")
        md.append(f"| Max abs bias   | {pn['max_abs_bias']:.4%} |")
        md.append(f"")
        md.append(f"### Per-panel breakdown")
        md.append(f"")
        md.append(f"| Panel | Predicted rate | Observed rate | Bias |")
        md.append(f"|---|---:|---:|---:|")
        for row in pn['per_panel']:
            md.append(
                f"| {row['panel']} "
                f"| {row['pred_offtarget_rate']:.4%} "
                f"| {row['true_offtarget_rate']:.4%} "
                f"| {row['bias']:+.4%} |"
            )

    return "\n".join(md)


def save_reports(run_dir: Path, cfg: Config, train_info: dict): #report: dict
    """Persist everything: config, metrics, tables, and the markdown summary."""
    # 1. Config used for this run
    (run_dir / "config.json").write_text(
        json.dumps(cfg.to_dict(), indent=2, default=str)
    )

    # 2. All metrics as JSON (machine-readable). Strip non-serializable bits.
    metrics = {
        #'site_level':   report['site_level'],
        #'primer_level': report['primer_level'],
        'training':     train_info,
    }
    

    # 4. Human-readable Markdown report
    md = render_markdown_report(cfg, train_info) # report
    (run_dir / "report.md").write_text(md)

# ============================================================================
# 3. MAIN
# ============================================================================

def main():
    # ----- Build config -----
    cfg = Config()
    # Override config here for this specific run (or load from YAML/CLI):
    cfg.run_name      = "baseline_onerun_train_val_test_splitbyPanel_bs550_ep30_lr3e-3_v1"
    cfg.epochs        = 30
    cfg.batch_size    = 550
    cfg.learning_rate = 3e-3
    cfg.allow_negative_output = True  
    cfg.seed = 42

    # ----- Setup the run directory, MLflow, and tee stdout -----
    run_dir, run_dir_test = setup_run_directory(cfg)
    setup_mlflow(cfg)
    set_global_seed(cfg.seed)
    sys.stdout = TeeOutput(run_dir / "training_log.txt")
    print(f"=== Run started: {datetime.now().isoformat()} ===")
    print(f"Output directory: {run_dir}")
    print()

    # ----- Start the MLflow run as a context manager -----
    with mlflow.start_run(run_name=cfg.run_name) as mlflow_run:
        print(f"MLflow run ID: {mlflow_run.info.run_id}")
        print()

        # ===== LOG ALL CONFIG PARAMETERS =====
        with mlflow.start_run(nested=True):
            mlflow.log_params(cfg.to_dict())

        # Tag the run with useful metadata for filtering later
        with mlflow.start_run(nested=True):
            mlflow.set_tag("phase", "diagnostic")
        with mlflow.start_run(nested=True):    
            mlflow.set_tag("model_status", "investigating_dying_relu")
        with mlflow.start_run(nested=True):    
            mlflow.set_tag("dataset_files", str(len(sorted(cfg.data_dir.glob("*.npz")))))

    # ----- Load data -----
    print("[1/5] Loading data...")
    npz_files = sorted(cfg.data_dir.glob("*.npz"))
    print(f"  Found {len(npz_files)} npz files")

    data = load_npz_files(npz_files)
    
    n_total = len(data['labels'])
    n_unique_primers = len(np.unique(data['primer_seqs']))

    #panel_ids =  np.concatenate(data['panel_ids']) if isinstance(data['panel_ids'], (list, np.ndarray)) else data['panel_ids']
    n_unique_panels = len(np.unique(data['panel_ids']))

    print(f"  Total alignments: {n_total:,}")
    print(f"  Unique primers:   {n_unique_primers:,}")
    print(f"  Unique panels:    {n_unique_panels:,}")
    print(f"Zero-UMI:     {(data['labels'] == 0).sum():,}  ({100*(data['labels']==0).mean():.1f}%)")
    print(f"Positive-UMI: {(data['labels'] > 0).sum():,}  ({100*(data['labels']>0).mean():.1f}%)")
    print()


    # ----- Panel or Primer-level split -----
    print("[2/5] Splitting at Panel level...")
    train_idx, val_idx = split_by_panel(
        data["panel_ids"], cfg.val_fraction, cfg.seed
    )

    # ---- Verify primer or panel_id disjointness — mandatory check ----
    train_panel = set(data['panel_ids'][train_idx])
    val_panel = set(data['panel_ids'][val_idx])
    overlap = train_panel & val_panel

    assert len(overlap) == 0, (
        f"Panel ID leakage detected: {len(overlap)} IDs "
        f"appear in both train and val. Split is broken."
    )

    print(f"  Unique panels in train: {len(train_panel):,}")
    print(f"  Unique panels in val:   {len(val_panel):,}")
    print(f"  Primers overlap:          0 ✓")

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

    #print(f"  Label stats train: mean={y_mean:.2f}, std={y_std:.2f}")


    train_ds = AlignmentDataset(data["images"][train_idx], train_labels, y_mean, y_std)
    val_ds   = AlignmentDataset(data["images"][val_idx],   labels[val_idx], y_mean, y_std)
 
    device = get_device(cfg.device)
    print(f"Using device: {device}")
    print(f"Train: {len(train_idx)} sites | "
          f"Val: {len(val_idx)} sites | ")
         # f"Test: {len(test_idx)} sites")


     # ---- 2. Train one run --------------------------------------------------
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
        model, best_val, stopped_at, best_epoch, history = train_one_run(
            cfg, train_ds, val_ds, device=device, use_mlflow=True,
        )
        training_run_id = run.info.run_id
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

    # ----- Save model and reports -----
    print("Saving the model...")
    model_path = run_dir / "model.pt"
    torch.save(model.state_dict(), model_path)

    # ----- Evaluate on the validation set -----
    print("[4/5] Evaluating on validation set...")
    val_preds = predict(model, val_ds, batch_size=cfg.batch_size, device=device)

    y_true_norm = (data['labels'][val_idx] - val_ds.y_mean.numpy()) / val_ds.y_std.numpy()
    y_pred_norm = val_preds

    y_true_denorm = data['labels'][val_idx]
    y_pred_denorm = val_preds * val_ds.y_std.numpy() + val_ds.y_mean.numpy()
    
    df_val = pd.DataFrame({
        'site_id': data['site_ids'][val_idx],
        'y_true_norm': y_true_norm,
        'y_pred_norm': y_pred_norm,
        'y_true_denorm': y_true_denorm,
        'y_pred_denorm': y_pred_denorm,
        'is_target': data['is_target'][val_idx],
        'r2': r2_score(y_true_denorm,  y_pred_denorm),
        'mse' : mean_squared_error(y_true_denorm,  y_pred_denorm),
        'rmse': np.sqrt(mean_squared_error(y_true_denorm,  y_pred_denorm)),
        'mae': mean_absolute_error(y_true_denorm,  y_pred_denorm),
    })


    df_val.to_csv(run_dir / f"{cfg.run_name}_validation_predictions.csv", index=False)

    # Plots 
    plot_pred_vs_true_fig = plot_pred_vs_true(y_true_denorm, y_pred_denorm, "Predicted vs true")
    save(plot_pred_vs_true_fig, run_dir, f"{cfg.run_name}_plot_pred_vs_true.png")

    plot_residuals_fig = plot_residuals(y_true_denorm, y_pred_denorm, "Residuals")
    save(plot_residuals_fig, run_dir, f"{cfg.run_name}_plot_residuals.png")


    plot_pred_vs_true_fig = plot_pred_vs_true(y_true_norm, y_pred_norm, "Predicted vs true (normalised)")
    save(plot_pred_vs_true_fig, run_dir, f"{cfg.run_name}_plot_pred_vs_true_norm.png")

    plot_residuals_fig = plot_residuals(y_true_norm, y_pred_norm, "Residuals (normalised)")
    save(plot_residuals_fig, run_dir, f"{cfg.run_name}_plot_residuals_norm.png")


    plot_loss_curves_fig = plot_loss_curves(history)
    save(plot_loss_curves_fig, run_dir, f"{cfg.run_name}_plot_loss_curves.png")

    # Save history to disk
    with open(run_dir / f"{cfg.run_name}_training_history.json", "w") as f:
        json.dump(history, f, indent=2)


    # ---- 5. Evaluate ONCE on the test set ----------------------------------
    print("[5/5] Evaluating on test set...")
    npz_files_test = sorted(cfg.data_test_dir.glob("*.npz"))
    test_data = load_npz_files(npz_files_test)

    n_total          = len(test_data["labels"])
    n_unique_primers = len(np.unique(test_data["primer_seqs"]))
    n_unique_panels  = len(np.unique(test_data["panel_ids"]))
    
    print(f"  Total alignments: {n_total:,}")
    print(f"  Unique primers:   {n_unique_primers:,}")
    print(f"  Unique panels:    {n_unique_panels:,}")
    print(f"  Zero-UMI:     {(data['labels'] == 0).sum():,}  "
          f"({100 * (data['labels'] == 0).mean():.1f}%)")
    print(f"  Positive-UMI: {(data['labels'] > 0).sum():,}  "
          f"({100 * (data['labels'] > 0).mean():.1f}%)")
    print()
 
    run_test(
        cfg,
        test_data=test_data,
        source={"model": model, "run_id": training_run_id},
        device=device,
        y_mean = float(train_labels.mean()),
        y_std  = float(train_labels.std()),
        run_dir_test=run_dir_test,
        run_name=cfg.run_name

    )


    # Later, load:
    #model.load_state_dict(torch.load("model.pt"))
    print("Saving artifacts...")
    train_info = {
        'epochs_run':   stopped_at,    # default: ran all epochs
        'best_epoch':       best_epoch,
        'duration_sec': duration,
        'best_val_mse': best_val,
        'n_train':      len(train_idx),
        'n_val':        len(val_idx),
        'device':       device,
    }
    save_reports(run_dir, cfg, train_info) 
    print(f"  All train artifacts saved to: {run_dir}")

    #save_report_test(run_dir_test, cfg. report_test)

    # ----- Log artifacts to MLflow too -----
    mlflow.log_artifact(str(model_path))
    mlflow.log_artifact(str(run_dir / "report.md"))
    mlflow.log_artifact(str(run_dir / "training_log.txt"))
    mlflow.log_artifact(str(run_dir / "config.json"))
    mlflow.log_artifact(str(run_dir / "training_history.json"))

    print(f"  MLflow run ID: {mlflow_run.info.run_id}")
    print()
    print(f"=== Run finished: {datetime.now().isoformat()} ===")


if __name__ == "__main__":
    main()