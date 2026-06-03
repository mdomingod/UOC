#!/usr/bin/env python3
# test_main.py
"""
Held-out TEST entry-point.

Run from the project ROOT:
    python test_main.py

What it does (mirrors one_run.py, but evaluates instead of trains):
    1. Load the held-out TEST .npz files (must be primer-disjoint from train/val).
    2. Load a trained model — either from a saved state-dict on disk OR an MLflow run.
    3. Run the full evaluation suite on the test set via
       src/project_package/test.run_test (site / primer / panel + bootstrap CIs).
    4. Save the model, metrics, per-primer / per-panel tables and a Markdown
       report to a timestamped output folder, and let run_test log to MLflow.

This is the "evaluate ONCE on test" step (Andrew Ng's ML discipline):
tune on the validation set in one_run.py, then evaluate here exactly once.
"""

import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import torch
import mlflow

# NOTE: import path matches one_run.py ("src.project_package..."). main.py uses
# "project_package..." instead — adjust these to whichever layout actually runs.
from src.project_package.config import Config
from src.project_package.data   import load_npz_files
from src.project_package.model  import PordleCNN
from src.project_package.utils  import get_device
from src.project_package.test   import run_test, load_model_from_run


# ============================================================================
# 1. SETUP  (copied from one_run.py so this script is self-contained)
# ============================================================================

def setup_run_directory(cfg: Config) -> Path:
    """Create a timestamped folder for this test run's outputs."""
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    run_dir = cfg.output_dir_dir / f"{timestamp}_TEST_{cfg.run_name}"
    run_dir.mkdir(exist_ok=True, parents=True)
    return run_dir


def set_global_seed(seed: int):
    """Make every source of randomness deterministic (bootstrap CIs use it)."""
    import random
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
    # run_test sets its own experiment (cfg.experiment_name + "_test"); this is
    # just the tracking URI / fallback experiment.
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
# 2. REPORT WRITING  (test variant — no training section)
# ============================================================================

def render_test_markdown_report(cfg: Config, report: dict) -> str:
    """Human-readable Markdown summary of the hold-out test evaluation."""
    md = []
    md.append(f"# Hold-out TEST Report — {cfg.run_name}")
    md.append("")
    md.append(f"_Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}_")
    md.append("")

    # Site-level metrics
    s = report["site_level"]
    md.append("## Site-level metrics")
    md.append("")
    md.append("How well does the model predict individual UMI counts?")
    md.append("")
    md.append("| Metric | Value |")
    md.append("|---|---|")
    md.append(f"| n_sites     | {s['n_sites']:,} |")
    md.append(f"| MSE         | {s['mse']:,.2f} |")
    md.append(f"| RMSE        | {s['rmse']:,.2f} |")
    md.append(f"| MAE         | {s['mae']:,.2f} |")
    md.append("")

    # Primer-level metrics
    p = report["primer_level"]
    md.append("## Primer-level metrics")
    md.append("")
    md.append("Can the model rank primers by off-target risk?")
    md.append("")
    md.append("| Metric | Value |")
    md.append("|---|---|")
    md.append(f"| n_primers              | {p['n_primers']:,} |")
    md.append(f"| n_dropped (zero total) | {p['n_dropped_zero_total']:,} |")
    md.append(f"| Mean bias              | {p['mean_bias']:+.4f} |")
    md.append(f"| Mean abs bias          | {p['mean_abs_bias']:.4f} |")
    md.append(f"| RMSE on rates          | {p['rmse']:.4f} |")
    md.append("")

    # Panel-level metrics (optional)
    if "panel_level" in report:
        pn = report["panel_level"]
        md.append("## Panel-level metrics")
        md.append("")
        md.append("The headline metric from the QIAGEN paper.")
        md.append("")
        md.append("| Metric | Value |")
        md.append("|---|---|")
        md.append(f"| n_panels     | {pn['n_panels']} |")
        md.append(f"| Mean bias    | {pn['mean_bias']:+.4%} |")
        md.append(f"| MAPB         | {pn['mean_abs_bias']:.4%} |")
        md.append(f"| Max abs bias | {pn['max_abs_bias']:.4%} |")
        md.append("")
        md.append("### Per-panel breakdown")
        md.append("")
        md.append("| Panel | Predicted rate | Observed rate | Bias |")
        md.append("|---|---:|---:|---:|")
        for row in pn["per_panel"]:
            md.append(
                f"| {row['panel']} "
                f"| {row['pred_offtarget_rate']:.4%} "
                f"| {row['true_offtarget_rate']:.4%} "
                f"| {row['bias']:+.4%} |"
            )

    return "\n".join(md)


def save_test_reports(run_dir: Path, cfg: Config, report: dict):
    """Persist config, metrics, tables and the Markdown summary for the test run."""
    # 1. Config used for this run
    (run_dir / "config.json").write_text(
        json.dumps(cfg.to_dict(), indent=2, default=str)
    )

    # 2. All metrics as JSON (machine-readable)
    metrics = {
        "site_level":   report["site_level"],
        "primer_level": report["primer_level"],
    }
    if "panel_level" in report:
        metrics["panel_level"] = report["panel_level"]
    (run_dir / "test_metrics.json").write_text(
        json.dumps(metrics, indent=2, default=str)
    )

    # 3. Per-primer / per-panel tables (CSV for downstream analysis)
    if "primer_table" in report:
        report["primer_table"].to_csv(run_dir / "test_primer_table.csv", index=False)
    if "panel_table" in report:
        report["panel_table"].to_csv(run_dir / "test_panel_table.csv", index=False)

    # 4. Human-readable Markdown report
    (run_dir / "test_report.md").write_text(render_test_markdown_report(cfg, report))


# ============================================================================
# 3. MAIN
# ============================================================================

def main():
    # ----- Build config -----
    cfg = Config()
    cfg.run_name = "baseline_DHS001Z_SplitByPrimer_tv04_bs600_ep50_lr1e-4"
    cfg.seed = 42

    # ===================== EDIT THESE TWO BLOCKS ============================
    # (a) Where the held-out TEST .npz files live. This must be a DIFFERENT,
    #     primer-disjoint set from what one_run.py trained on. If your Config
    #     already exposes a test directory, use that instead of a literal path.
    TEST_DATA_DIR = getattr(cfg, "test_data_dir", cfg.data_dir)

    # (b) Where the trained model comes from. Pick ONE:
    #     - From disk: the model.pt that one_run.py saved (run_dir / "model.pt").
    #     - From MLflow: the run_id of the training run.
    USE_MLFLOW_MODEL = False
    MODEL_STATE_DICT_PATH = Path("reports/<your_training_run_dir>/model.pt")
    MODEL_RUN_ID = "<paste_training_mlflow_run_id_here>"
    # =======================================================================

    # ----- Setup the run directory, MLflow, seed, and tee stdout -----
    run_dir = setup_run_directory(cfg)
    setup_mlflow(cfg)
    set_global_seed(cfg.seed)
    sys.stdout = TeeOutput(run_dir / "test_log.txt")
    print(f"=== TEST run started: {datetime.now().isoformat()} ===")
    print(f"Output directory: {run_dir}")
    print()

    device = get_device(cfg.device)
    print(f"Using device: {device}")

    # ----- 1. Load the held-out TEST data ----------------------------------
    print("[1/4] Loading TEST data...")
    npz_files = sorted(cfg.data_test_dir.glob("*.npz"))
    #npz_files = sorted(Path(TEST_DATA_DIR).glob("*.npz"))
    if not npz_files:
        raise FileNotFoundError(f"No .npz files found in {TEST_DATA_DIR}")
    print(f"  Found {len(npz_files)} test npz files")

    data = load_npz_files(npz_files)

    n_total          = len(data["labels"])
    n_unique_primers = len(np.unique(data["primer_seqs"]))
    n_unique_panels  = len(np.unique(data["panel_ids"]))
    
    print(f"  Total alignments: {n_total:,}")
    print(f"  Unique primers:   {n_unique_primers:,}")
    print(f"  Unique panels:    {n_unique_panels:,}")
    print(f"  Zero-UMI:     {(data['labels'] == 0).sum():,}  "
          f"({100 * (data['labels'] == 0).mean():.1f}%)")
    print(f"  Positive-UMI: {(data['labels'] > 0).sum():,}  "
          f"({100 * (data['labels'] > 0).mean():.1f}%)")
    print()

    # The full loaded dict already carries every key run_test needs
    # (images, labels, primer_seqs, is_target, site_ids, panel_ids).
    test_data = data

    # ----- 2. Load the trained model ---------------------------------------
    print("[2/4] Loading trained model...")
    if USE_MLFLOW_MODEL:
        model = load_model_from_run(MODEL_RUN_ID, device)
        source = {"model": model, "run_id": MODEL_RUN_ID}
        print(f"  Loaded model from MLflow run {MODEL_RUN_ID}")
    else:
        model = PordleCNN(cfg).to(device)
        state = torch.load(MODEL_STATE_DICT_PATH, map_location=device)
        model.load_state_dict(state)
        model.eval()
        source = {"model": model}
        print(f"  Loaded model state-dict from {MODEL_STATE_DICT_PATH}")

    # Persist the exact weights used for this evaluation alongside the results
    model_out = run_dir / "model.pt"
    torch.save(model.state_dict(), model_out)
    print(f"  Saved model copy to {model_out}")
    print()

    # ----- 3. Evaluate ONCE on the test set --------------------------------
    print("[3/4] Evaluating on the hold-out test set...")
    report = run_test(
        cfg,
        test_data=test_data,
        source=source,
        device=device,
    )
    print()

    # ----- 4. Save evaluation + results to the output folder ---------------
    print("[4/4] Saving artifacts...")
    save_test_reports(run_dir, cfg, report)
    print(f"  All artifacts saved to: {run_dir}")

    # Echo headline numbers to the log
    s = report["site_level"]
    print()
    print(f"  Site-level RMSE: {s['rmse']:,.2f}")
    print(f"  Site-level MAE:  {s['mae']:,.2f}")
    if "primer_level" in report:
        print(f"  Primer-level MAPB: {report['primer_level']['mean_abs_bias']:.4f}")
    if "panel_level" in report:
        print(f"  Panel-level MAPB:  {report['panel_level']['mean_abs_bias']:.4%}")

    print()
    print(f"=== TEST run finished: {datetime.now().isoformat()} ===")


if __name__ == "__main__":
    main()
