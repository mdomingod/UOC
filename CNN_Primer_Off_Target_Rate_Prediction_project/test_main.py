import json
import sys
import time
from datetime import datetime
from pathlib import Path

from mlflow import data
from tensorboard import data
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


def main():
    # ----- Build config -----
    cfg = Config()
    # Override config here for this specific run (or load from YAML/CLI):
    cfg.run_name      = "baseline_onerun_train_val_test_splitbyPrimer_bs600_ep50_lr1e-3_v1" 
    cfg.epochs        = 50
    cfg.batch_size    = 600
    cfg.learning_rate = 1e-3
    cfg.allow_negative_output = True  
    cfg.seed = 42

    device = get_device(cfg.device)

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

    #  mean: 61.85, std: 186.56 for DHS001Z_SplitByPrimer_tv04_bs600_ep50_lr1e-4
    # mean: 27.74, std: 381.55 splitbyPanel_bs600_ep50_lr1e-4_v1


    run_test(
        cfg,
        test_data=test_data,
        source={"model": model},
        run_name=cfg.run_name,
        run_dir_test=cfg.run_dir_test,
        y_mean=cfg.y_mean,
        y_std=cfg.y_std,
        device=device,
        

    )