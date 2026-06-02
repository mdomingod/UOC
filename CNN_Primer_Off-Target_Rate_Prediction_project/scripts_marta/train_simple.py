# src/project_package/train.py
import numpy as np
import torch
import torch.nn as nn
import mlflow
from torch.utils.data import DataLoader
from .config import Config
from .model import PordleCNN
from .evaluate import evaluate_model



def train_one_run(cfg: Config, train_ds, val_ds, device='cpu', mlflow_run=None):
    """Train a single model. Returns the trained model and its best validation MSE."""
    train_loader = DataLoader(train_ds, batch_size=cfg.batch_size, shuffle=True,  num_workers=2)
    val_loader   = DataLoader(val_ds,   batch_size=cfg.batch_size, shuffle=False, num_workers=2)

    model = PordleCNN(cfg).to(device)
    loss_fn = nn.MSELoss()

    """ Types of lost functions:
- MSELoss: Mean Squared Error, good for regression tasks where you want to penalize larger errors more heavily.
- L1Loss: Mean Absolute Error, good for regression tasks where you want to treat all errors equally.
- CrossEntropyLoss: Good for multi-class classification tasks, combines LogSoftmax and NLLLoss in one single class.
- BCEWithLogitsLoss: Good for binary classification tasks, combines a Sigmoid layer and the BCELoss in one single class.
    """

    optimizer = torch.optim.Adam(model.parameters(), lr=cfg.learning_rate, momentum=cfg.momentum)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', factor=0.1, patience=cfg.lr_patience
    )

    # Log architecture info once
    if mlflow_run is not None:
        mlflow.log_param("n_parameters", model.count_parameters())

    best_val = float('inf')
    best_state = None
    epochs_no_improve = 0

    for epoch in range(1, cfg.epochs + 1):
        model.train()
        train_losses = []
        for x, y in train_loader:
            x, y = x.to(device), y.to(device)
            optimizer.zero_grad()
            loss = loss_fn(model(x), y)
            loss.backward()
            optimizer.step()
            train_losses.append(loss.item())

        model.eval()
        val_losses = []
        with torch.no_grad():
            for x, y in val_loader:
                x, y = x.to(device), y.to(device)
                val_losses.append(loss_fn(model(x), y).item())

        train_loss = float(np.mean(train_losses))
        val_loss   = float(np.mean(val_losses))
        scheduler.step(val_loss)

        # ---- MLflow logging happens per-epoch ----
        if mlflow_run is not None:
            mlflow.log_metrics(
                {"train_mse": train_loss, "val_mse": val_loss,
                 "lr": optimizer.param_groups[0]['lr']},
                step=epoch
            )

        print(f"Epoch {epoch:3d} | train MSE {train_loss:9.2f} | val MSE {val_loss:9.2f}")

        if val_loss < best_val:
            best_val = val_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1
            if epochs_no_improve >= cfg.patience:
                print(f"Early stopping at epoch {epoch}")
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    return model, best_val