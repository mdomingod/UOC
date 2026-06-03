# src/project_package/evaluate.py
"""
Evaluation module — site, primer, and panel level metrics.
 
Computes the metrics needed for TFM-quality model evaluation:
  - Site-level:   MSE, MAE, Pearson, Spearman on raw UMI predictions
  - Primer-level: per-primer off-target rate and prediction error
  - Panel-level:  per-panel off-target rate bias and MAPB
  - Bootstrap CIs for any metric
 
Pure functions operate on numpy arrays / DataFrames so they're reusable for
baseline models (SVR, Random Forest) too. MLflow logging is kept in a separate
set of functions (`log_evaluation_to_mlflow`, `make_evaluation_plots`) so the
metric code remains side-effect free.
"""

from __future__ import annotations
 
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, mean_absolute_error
 
from .utils import (
    MLFLOW_AVAILABLE,
    log_figure,
    plot_panel_offtarget_comparison,
    plot_pred_vs_true,
    plot_residuals,
)
 
if MLFLOW_AVAILABLE:
    import mlflow


# ============================================================================
# 1. INFERENCE
# ============================================================================

@torch.no_grad()
def predict(model, dataset, batch_size: int = 128, device: str = 'cpu') -> np.ndarray:
    """Run the model on a dataset and return predictions as a 1D numpy array."""
    model.to(device).eval()
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
    preds = []
    for batch in loader:
        x = batch[0] if isinstance(batch, (tuple, list)) else batch
        preds.append(model(x.to(device)).cpu().numpy())
    return np.concatenate(preds)


# ============================================================================
# 2. SITE-LEVEL METRICS  (one prediction per alignment)
# ============================================================================

def site_level_metrics(y_pred: np.ndarray, y_true: np.ndarray) -> dict:
    """
    Site-level regression metrics. Answers: how well does the model predict
    individual UMI counts at each alignment?
    """
    y_pred = np.asarray(y_pred, dtype=np.float64)
    y_true = np.asarray(y_true, dtype=np.float64)

    # Clip negative predictions to zero — UMI counts can't be negative.
    # (Only matters if allow_negative_output=True or log-transform used.)
    #y_pred_clipped = np.clip(y_pred, 0, None)

    mse  = float(mean_squared_error(y_true, y_pred))
    mae  = float(mean_absolute_error(y_true, y_pred))
    rmse = float(np.sqrt(mse))

    """# Correlations — only meaningful if predictions and labels both vary
    if y_pred.std() > 0 and y_true.std() > 0:
        pearson_r, _ = pearsonr(y_pred, y_true)
        spearman_r, _ = spearmanr(y_pred, y_true)
    else:
        pearson_r, spearman_r = float('nan'), float('nan')"""

    return {
        'mse':      mse,
        'rmse':     rmse,
        'mae':      mae,
        #'pearson':  float(pearson_r),
        #'spearman': float(spearman_r),
        'n_sites':  len(y_true),
    }


# ============================================================================
# 3. PRIMER-LEVEL AGGREGATION
# ============================================================================

def aggregate_to_primer(
    y_pred: np.ndarray,
    y_true: np.ndarray,
    primer_seqs: np.ndarray,
    is_target: np.ndarray,
) -> pd.DataFrame:
    """
    Aggregate site-level predictions to one row per primer.

    Returns a DataFrame with columns:
        primer_seqs, pred_target_umi, pred_offtarget_umi, pred_total_umi,
        true_target_umi, true_offtarget_umi, true_total_umi,
        pred_offtarget_rate, true_offtarget_rate
    """
    #y_pred_clipped = np.clip(y_pred, 0, None)
    is_target = is_target.astype(bool)

    df = pd.DataFrame({
        'primer_seqs':     primer_seqs,
        'is_target':     is_target,
        'pred_umi':      y_pred,
        'true_umi':      y_true,
    })

    # Split each primer's sites into target and off-target totals
    grouped = df.groupby('primer_seqs').apply(lambda g: pd.Series({
        'pred_target_umi':    g.loc[g['is_target'],  'pred_umi'].sum(),
        'pred_offtarget_umi': g.loc[~g['is_target'], 'pred_umi'].sum(),
        'true_target_umi':    g.loc[g['is_target'],  'true_umi'].sum(),
        'true_offtarget_umi': g.loc[~g['is_target'], 'true_umi'].sum(),
        'n_sites':            len(g),
    })).reset_index()

    grouped['pred_total_umi'] = grouped['pred_target_umi'] + grouped['pred_offtarget_umi']
    grouped['true_total_umi'] = grouped['true_target_umi'] + grouped['true_offtarget_umi']

    grouped['pred_offtarget_rate'] = np.where(
        grouped['pred_total_umi'] > 0,
        grouped['pred_offtarget_umi'] / grouped['pred_total_umi'],
        np.nan
    )
    grouped['true_offtarget_rate'] = np.where(
        grouped['true_total_umi'] > 0,
        grouped['true_offtarget_umi'] / grouped['true_total_umi'],
        np.nan
    )

    return grouped


def primer_level_metrics(primer_df: pd.DataFrame) -> dict:
    """
    Primer-level metrics from the aggregated DataFrame.
    Answers: can the model rank primers by off-target risk?
    """
    valid = primer_df.dropna(subset=['pred_offtarget_rate', 'true_offtarget_rate'])
    n_dropped = len(primer_df) - len(valid)

    pred = valid['pred_offtarget_rate'].values
    true = valid['true_offtarget_rate'].values
    bias = (pred - true)

    """if pred.std() > 0 and true.std() > 0:
        pearson_r, _ = pearsonr(pred, true)
        spearman_r, _ = spearmanr(pred, true)
    else:
        pearson_r, spearman_r = float('nan'), float('nan')"""

    return {
        'n_primers':              int(len(valid)),
        'n_dropped_zero_total':   int(n_dropped),
        'mean_bias':              float(bias.mean()),
        'mean_abs_bias':          float(np.abs(bias).mean()),
        'rmse':                   float(np.sqrt((bias ** 2).mean())),
        #'pearson':                float(pearson_r),
        #'spearman':               float(spearman_r),
    }


# ============================================================================
# 4. PANEL-LEVEL AGGREGATION
# ============================================================================

def _to_str_array(arr: np.ndarray) -> np.ndarray:
    """
    Decode bytes -> str if needed. npz files often store strings as `bytes`
    (dtype kind 'S') or as object arrays — pandas groupby keys are cleaner
    with proper Python strings.
    """
    arr = np.asarray(arr)
    if arr.dtype.kind in ("S", "O"):
        return np.asarray([
            x.decode("utf-8") if isinstance(x, bytes) else str(x)
            for x in arr.ravel()
        ]).reshape(arr.shape)
    return arr.astype(str)


def aggregate_to_panel(
    y_pred: np.ndarray,
    y_true: np.ndarray,
    panel_ids: np.ndarray,
    is_target: np.ndarray,
    primer_seqs: np.ndarray | None = None,
) -> pd.DataFrame:
    """
    Aggregate site-level predictions to one row per panel.

    Goes directly from sites to panels (not via primers) so that primers
    appearing in multiple panels — explicitly noted in the QIAGEN paper —
    are correctly counted in each panel they belong to.

    `primer_seqs` is optional; if provided, the output gains an `n_primers`
    column (unique primers per panel).
    """
    #y_pred_clipped = np.clip(np.asarray(y_pred, dtype=np.float64), 0, None)
    y_pred         = np.asarray(y_pred, dtype=np.float64)
    y_true         = np.asarray(y_true, dtype=np.float64)
    is_target_b    = np.asarray(is_target).astype(bool)
    panel_ids    = _to_str_array(panel_ids)

    # Pre-compute conditional sums so we can use the fast .agg path
    df = pd.DataFrame({
        "panel":              panel_ids,
        "pred_target_umi":    np.where(is_target_b,  y_pred, 0.0),
        "pred_offtarget_umi": np.where(~is_target_b, y_pred, 0.0),
        "true_target_umi":    np.where(is_target_b,  y_true,         0.0),
        "true_offtarget_umi": np.where(~is_target_b, y_true,         0.0),
    })
    if primer_seqs is not None:
        df["primer_seqs"] = primer_seqs

    agg_spec = {
        "n_sites":            ("pred_target_umi", "count"),
        "pred_target_umi":    ("pred_target_umi", "sum"),
        "pred_offtarget_umi": ("pred_offtarget_umi", "sum"),
        "true_target_umi":    ("true_target_umi", "sum"),
        "true_offtarget_umi": ("true_offtarget_umi", "sum"),
    }
    if primer_seqs is not None:
        agg_spec["n_primers"] = ("primer_seq", "nunique")

    grouped = df.groupby("panel", sort=True).agg(**agg_spec).reset_index()

    grouped["pred_total_umi"] = grouped["pred_target_umi"] + grouped["pred_offtarget_umi"]
    grouped["true_total_umi"] = grouped["true_target_umi"] + grouped["true_offtarget_umi"]

    grouped["pred_offtarget_rate"] = np.where(
        grouped["pred_total_umi"] > 0,
        grouped["pred_offtarget_umi"] / grouped["pred_total_umi"],
        np.nan,
    )
    grouped["true_offtarget_rate"] = np.where(
        grouped["true_total_umi"] > 0,
        grouped["true_offtarget_umi"] / grouped["true_total_umi"],
        np.nan,
    )
    grouped["bias"]     = grouped["pred_offtarget_rate"] - grouped["true_offtarget_rate"]
    grouped["abs_bias"] = np.abs(grouped["bias"])

    return grouped


def panel_level_metrics(panel_df: pd.DataFrame) -> dict:
    """
    Panel-level metrics — the headline numbers from the QIAGEN paper.

    Reports MAPB (mean absolute prediction bias) across panels — the same
    metric used for model selection in the paper's cross-validation.
    """
    valid = panel_df.dropna(subset=["pred_offtarget_rate", "true_offtarget_rate"])
    n_dropped = len(panel_df) - len(valid)
    bias = valid["bias"].values

    out = {
        "n_panels":              int(len(valid)),
        "n_dropped_zero_total":  int(n_dropped),
        "mean_bias":             float(bias.mean())            if len(bias) else float("nan"),
        "mean_abs_bias":         float(np.abs(bias).mean())    if len(bias) else float("nan"),  # MAPB
        "max_abs_bias":          float(np.abs(bias).max())     if len(bias) else float("nan"),
        "rmse_offtarget_rate":   float(np.sqrt((bias ** 2).mean())) if len(bias) else float("nan"),
        "per_panel": valid[
            [c for c in ["panel", "n_sites", "n_primers",
                         "pred_offtarget_rate", "true_offtarget_rate", "bias"]
             if c in valid.columns]
        ].to_dict("records"),
    }
    return out


# ============================================================================
# 5. UNCERTAINTY — bootstrap CIs
# ============================================================================

def bootstrap_ci(
    y_pred: np.ndarray,
    y_true: np.ndarray,
    metric_fn,
    n_iterations: int = 1000,
    confidence: float = 0.95,
    seed: int = 42,
) -> tuple[float, float, float]:
    """
    Bootstrap confidence interval for any metric.
    `metric_fn` takes (y_pred, y_true) and returns a scalar.
    Returns (point_estimate, lower, upper).
    """
    rng = np.random.default_rng(seed)
    n = len(y_pred)
    point = metric_fn(y_pred, y_true)

    samples = np.empty(n_iterations)
    for i in range(n_iterations):
        idx = rng.integers(0, n, size=n)
        samples[i] = metric_fn(y_pred[idx], y_true[idx])

    alpha = (1 - confidence) / 2
    lower = float(np.quantile(samples, alpha))
    upper = float(np.quantile(samples, 1 - alpha))
    return float(point), lower, upper


# ============================================================================
# 6. TOP-LEVEL CONVENIENCE
# ============================================================================

def evaluate_predictions(
    y_pred: np.ndarray,
    y_true: np.ndarray,
    primer_seqs: np.ndarray,
    is_target: np.ndarray,
    panel_ids: np.ndarray | None = None,
) -> dict:
    """
    Run all evaluations in one call. 
    Returns a nested dict with site, primer, and (optionally) panel-level metrics.

    All inputs are per-site arrays of the same length.
    """
    report = {
        'site_level': site_level_metrics(y_pred, y_true),
    }

    primer_df = aggregate_to_primer(y_pred, y_true, primer_seqs, is_target)
    report['primer_level'] = primer_level_metrics(primer_df)
    report['primer_table'] = primer_df

    if panel_ids is not None:
        # Build a primer → panel mapping (one panel per primer)
        primer_to_panel = (
            pd.DataFrame({'primer_seqs': primer_seqs, 'panel': panel_ids})
            .drop_duplicates('primer_seqs')
            .set_index('primer_seqs')['panel']
        )
        panel_per_primer = primer_df['primer_seqs'].map(primer_to_panel).values

        panel_df = aggregate_to_panel(primer_df, panel_per_primer)
        report['panel_level'] = panel_level_metrics(panel_df)
        report['panel_table'] = panel_df

    # Raw predictions for later artifact logging
    report["_y_pred_norm"] = y_pred
    report["_y_true_norm"] = y_true

    return report
