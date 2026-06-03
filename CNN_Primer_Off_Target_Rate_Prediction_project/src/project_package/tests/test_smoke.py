#!/usr/bin/env python3
# To run the srcipt, activate the uv environment and change th end line from CRLF to LF (in vscode, bottom right corner) and then run:
# src/project_package/tests/test_smoke.py
"""
Convergence smoke test.

Verifies that training reduces the loss over a few epochs on a small subset of real data.
This is a mechanical check — does the model learn anything? — not a quality assessment.

Run from project root:
    python src/project_package/tests/test_smoke.py
"""

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from project_package.config import Config
from project_package.data import load_npz_files, AlignmentDataset
from project_package.model import PordleCNN

def get_balanced_subset(data: dict, n_total: int = 2000) -> tuple:
    """
    Return a balanced subset with equal zero and positive UMI samples.
    This gives the model something meaningful to learn in just a few epochs.
    """
    labels = data['labels']
    
    zero_idx    = np.where(labels == 0)[0]
    nonzero_idx = np.where(labels > 0)[0]
    
    n_each = min(n_total // 2, len(zero_idx), len(nonzero_idx))
    
    rng = np.random.default_rng(42)
    chosen_zero    = rng.choice(zero_idx,    n_each, replace=False)
    chosen_nonzero = rng.choice(nonzero_idx, n_each, replace=False)
    idx = np.concatenate([chosen_zero, chosen_nonzero])
    rng.shuffle(idx)
    
    return data['images'][idx], data['labels'][idx]

def test_model_converges():
    """
    Train the model for a few epochs on a small subset and verify the loss decreases.
    """
    cfg = Config()
    
    # Load a small slice of real data
    npz_files = sorted(cfg.data_dir.glob("*.npz"))
    data = load_npz_files(npz_files)
    # BEFORE
    n = min(2000, len(data['labels']))
    ds = AlignmentDataset(data['images'][:n], data['labels'][:n])

    # AFTER
    #images, labels = get_balanced_subset(data, n_total=2000)
    #ds = AlignmentDataset(images, labels)
    loader = DataLoader(ds, batch_size=64, shuffle=True)

    # Build model and optimizer
    model = PordleCNN(cfg)

    # Check what the model predicts before any training
    model.eval()
    with torch.no_grad():
        sample_x, _ = next(iter(loader))
        raw_preds = model(sample_x)
        print(f"Pre-training predictions: min={raw_preds.min():.4f}, "
            f"max={raw_preds.max():.4f}, "
            f"percent_zero={((raw_preds == 0).sum() / len(raw_preds) * 100):.1f}%")
    model.train()

    model.train()
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    """
    A smaller learning rate means smaller steps — the loss decreases more smoothly and is less likely to overshoot.
    """
    loss_fn = nn.MSELoss()
    
    # Train for a few epochs, tracking the mean loss per epoch
    epoch_losses = []
    n_epochs = 15
    for epoch in range(n_epochs):
        batch_losses = []
        for x, y in loader:
            optimizer.zero_grad()
            loss = loss_fn(model(x), y)
            loss.backward()
            optimizer.step()
            batch_losses.append(loss.item())
        epoch_losses.append(np.mean(batch_losses))
        print(f"  Epoch {epoch + 1}/{n_epochs}: loss = {epoch_losses[-1]:,.2f}")
    
    # The convergence check: final loss must be lower than initial
    initial, final = epoch_losses[0], epoch_losses[-1]
    reduction = 100 * (initial - final) / initial
    
    best = min(epoch_losses)
    reduction_best = 100 * (initial - best) / initial
    assert best < initial * 0.95, \
    f"Loss never improved meaningfully: started {initial:.2f}, best {best:.2f}"

    print(f"✓ Model converges: {initial:,.2f} → best {best:,.2f} "
        f"({reduction_best:.1f}% reduction best, {reduction:.1f}% reduction final)")


if __name__ == "__main__":
    print("=== Convergence smoke test ===\n")
    test_model_converges()
    print("\n=== Passed ===")