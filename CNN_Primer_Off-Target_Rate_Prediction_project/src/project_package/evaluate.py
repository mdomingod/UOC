# src/project_package/evaluate.py
import numpy as np
import torch
from torch.utils.data import DataLoader

@torch.no_grad()
def predict(model, dataset, batch_size=128, device='cpu'):
    model.to(device).eval()
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
    preds = []
    for batch in loader:
        x = batch[0] if isinstance(batch, (tuple, list)) else batch
        preds.append(model(x.to(device)).cpu().numpy())
    return np.concatenate(preds)

def compute_offtarget_rate(predictions, is_ontarget_mask):
    predictions = np.clipclear(predictions, 0, None)
    n_target    = predictions[is_ontarget_mask].sum()
    n_offtarget = predictions[~is_ontarget_mask].sum()
    total = n_target + n_offtarget
    return n_offtarget / total if total > 0 else 0.0

def level_bias(predicted_rates, observed_rates):
    biases = np.array(predicted_rates) - np.array(observed_rates)
    return {
        'mean_bias': float(biases.mean()),
        'mean_abs_bias': float(np.abs(biases).mean()),
        'per_panel_bias': biases.tolist()
    }

def evaluate_model(model, dataset, true_labels, is_ontarget, device='cpu'):
    """Single function that produces all the metrics you want at evaluation time."""
    preds = predict(model, dataset, device=device)
    mse = float(((preds - true_labels) ** 2).mean())
    mae = float(np.abs(preds - true_labels).mean())
    offtarget_rate = compute_offtarget_rate(preds, is_ontarget)
    return {
        'mse': mse,
        'mae': mae,
        'predicted_offtarget_rate': offtarget_rate,
    }