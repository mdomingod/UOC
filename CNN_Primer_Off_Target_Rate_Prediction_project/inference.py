#!/usr/bin/env python3
"""Minimal inference with a locally saved PordleCNN model.

Usage:
    python infer.py  <model.pt>  <data.npz>  <y_mean>  <y_std>

- <model.pt>  : a saved state_dict (e.g. models/bypanel_bs600_ep50_lr1e4_model.pt)
- <data.npz>  : an encoded file with an 'images' array, shape (N, 5, 66, 2), uint8
- <y_mean>,<y_std> : the TRAINING label mean/std used when that model was trained
                     (printed by the one_run_*.py script; needed to denormalise)
"""
import sys
import numpy as np
import torch

from src.project_package.config import Config
from src.project_package.model import PordleCNN

model_path, npz_path = sys.argv[1], sys.argv[2]
y_mean, y_std = float(sys.argv[3]), float(sys.argv[4])

# 1. Rebuild the model and load the saved weights
cfg = Config()                                  # defaults must match how it was trained
model = PordleCNN(cfg)
model.load_state_dict(torch.load(model_path, map_location="cpu"))
model.eval()

# 2. Load and prepare the images exactly as the dataset does
images = np.load(npz_path, allow_pickle=True)["images"]   # (N, 5, 66, 2) uint8
x = torch.from_numpy(images.astype(np.float32) / 255.0)   # scale to [0, 1]
x = x.permute(0, 3, 1, 2)                                 # (N, 5, 66, 2) -> (N, 2, 5, 66)

# 3. Predict (output is in normalised space), then denormalise to real counts
with torch.no_grad():
    pred_norm = model(x)
pred_counts = pred_norm.numpy() * (y_std + 1e-8) + y_mean

# 4. Show results
for i, p in enumerate(pred_counts):
    print(f"site {i}: predicted UMI count = {p:.2f}")