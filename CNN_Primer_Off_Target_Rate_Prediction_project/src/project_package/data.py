# src/project_package/data.py
import numpy as np
import torch
from pathlib import Path
from torch.utils.data import Dataset
from collections import defaultdict

import torchvision
from torchvision.transforms import v2

# PyTorch TensorBoard support

from torch.utils.tensorboard import SummaryWriter
from datetime import datetime


def load_npz_files(npz_paths: list[Path]) -> dict[str, np.ndarray]:
    """
    Load one or more .npz files and concatenate their arrays.
    Expects each file to contain: 'images', 'labels', 'site_ids'.
    """
    images_list, labels_list, site_ids_list = [], [], []
    target_list = []
    """has_target_per_file = []"""
    primer_seq_list = []
    panel_id_list = []
    """has_target_per_file, has_panel_per_file = [], []"""

    for path in npz_paths:
        with np.load(path, allow_pickle=True) as npz:
            images_list.append(npz['images'])
            labels_list.append(npz['labels'])

            # Primer IDs
            raw_site_ids = npz["site_ids"]            # shape (n, 2) or object array of tuples
            # Collapse (alignChrom, read1_anchor) -> "alignChrom:read1_anchor"
            if raw_site_ids.ndim == 2:
                pid_1d = np.array(
                    [f"{row[0]}:{row[1]}" for row in raw_site_ids],
                    dtype=object,
                )
            else:
                # already 1-D (object array of tuples, or already-joined strings)
                pid_1d = np.array(
                    [f"{t[0]}:{t[1]}" if isinstance(t, (tuple, list, np.ndarray)) else str(t)
                    for t in raw_site_ids],
                    dtype=object,
                )
            site_ids_list.append(pid_1d)

            # is_target
            target_list.append(npz['is_target'])
            
            """if 'is_target' in npz.files:
                target_list.append(npz['is_target'])
                has_target_per_file.append(True)
            else:
                has_target_per_file.append(False)"""
            

            # primer_seq
            primer_seq_list.append(npz['primer_seq'])

            
            # panel_id scalar in npz -> broadcast to (n,) so it aligns site-wise
            pid = npz["panel_id"]
            pid_scalar = pid.item() if pid.ndim == 0 else pid[0]
            if isinstance(pid_scalar, bytes):
                pid_scalar = pid_scalar.decode("utf-8")
            panel_id_list.append(np.full(len(npz['images']), pid_scalar, dtype=object))


            """ if "panel_id" in npz.files:
                # scalar in npz -> broadcast to (n,) so it aligns site-wise
                pid = npz["panel_id"]
                pid_scalar = pid.item() if pid.ndim == 0 else pid[0]
                if isinstance(pid_scalar, bytes):
                    pid_scalar = pid_scalar.decode("utf-8")
                panel_id_list.append(np.full(len(npz['images']), pid_scalar, dtype=object))
                has_panel_per_file.append(True)
            else:
                has_panel_per_file.append(False)"""

    """# Consistency checks (same pattern you already use for is_target)
    for name, flags in [("is_target", has_target_per_file),
                        ("panel_id",  has_panel_per_file)]:
        if any(flags) and not all(flags):
            missing = [p.name for p, has in zip(npz_paths, flags) if not has]
            raise ValueError(
                f"Inconsistent dataset: '{name}' is missing from {len(missing)} files. "
                f"Missing in: {missing}."
            )"""

    out = {
        "images":     np.concatenate(images_list,     axis=0),
        "labels":     np.concatenate(labels_list,     axis=0),
        "site_ids": np.concatenate(site_ids_list, axis=0),
        "is_target": np.concatenate(target_list,      axis=0),
        "panel_ids": np.concatenate(panel_id_list,    axis=0),
        "primer_seqs": np.concatenate(primer_seq_list, axis=0)
    }
    """if all(has_target_per_file):
        out["is_target"] = np.concatenate(target_list, axis=0)
    if all(has_panel_per_file):
        out["panel_ids"] = np.concatenate(panel_id_list, axis=0)"""

    return out


class AlignmentDataset(Dataset):
    """
    Wraps pre-encoded images + labels.
    - Images: lazy uint8 → float32 /255 per sample (memory-efficient)
    - Labels: z-score normalised using externally provided stats
    - Channel layout: detected once at init, not re-checked per sample
    """
    """Wraps pre-encoded images + labels. Scales uint8 -> float32 in [0, 1]."""
    def __init__(self, images: np.ndarray, labels: np.ndarray, y_mean: float, y_std: float):
        self.images = images
        #self.labels = labels

        # ---- Detect channel layout once, not per sample ----
        # images shape is either (N, H, W, C) [channels-last] or (N, C, H, W) [channels-first]
        # The model expects (C, H, W), so we need to permute if channels-last
        sample_shape = images.shape[1:]   # (H, W, C) or (C, H, W)
        self._needs_permute = (
            len(sample_shape) == 3 and sample_shape[-1] in (1, 2, 3)
        )


        #y = torch.tensor(self.labels, dtype=torch.float32)
        # Normalise labels upfront — vectorised, done once
        y = torch.from_numpy(np.asarray(labels, dtype=np.float32))
        self.y_mean  = torch.tensor(y_mean,         dtype=torch.float32)
        self.y_std   = torch.tensor(y_std + 1e-8,   dtype=torch.float32)
        self.labels  = (y - self.y_mean) / self.y_std

        # 


    def __len__(self):
        """
        PyTorch's DataLoader calls this to know how many batches to produce per epoch. 
        Trivial but required — without __len__ the DataLoader falls back to streaming mode and you lose shuffling.
        """
        return len(self.images)
        

    def __getitem__(self, idx):
        x = torch.from_numpy(self.images[idx].astype(np.float32) / 255.0)
        
        # Make sure shape is (C, H, W) — adjust if your npz stores (H, W, C)
        """if x.ndim == 3 and x.shape[-1] in (1, 2, 3):
            x = x.permute(2, 0, 1)
        if self.labels is None:
            return x"""
        
        if self._needs_permute:        # branch eliminated at runtime by Python — bool lookup only
            x = x.permute(2, 0, 1)

        y_norm = self.labels[idx]

        return x, y_norm


def split_by_site_id(site_ids: np.ndarray, val_fraction: float = 0.2, seed: int = 42): # in fact it is split by site.
    """Return (train_indices, val_indices) sites (chr:read1pos) into two groups."""

    rng = np.random.default_rng(seed)
    unique_site_ids = np.unique(site_ids)
    rng.shuffle(unique_site_ids)

    n_val = int(len(unique_site_ids) * val_fraction)
    #val_primers = set(unique_site_ids[:n_val]) BEFORE
    val_site_ids = unique_site_ids[:n_val]          # numpy slice — no set, no list

    # BEFORE
    #train_idx = np.where(~np.isin(site_ids, list(val_primers)))[0] 
    #val_idx   = np.where( np.isin(site_ids, list(val_primers)))[0]

    # Compute the boolean mask ONCE
    val_mask   = np.isin(site_ids, val_site_ids)

    # Derive both index arrays from the same mask
    val_idx   = np.where( val_mask)[0]
    train_idx = np.where(~val_mask)[0]

    return train_idx, val_idx



def split_by_primer_seq(primer_seqs: np.ndarray, val_fraction: float = 0.2, seed: int = 42):
    """Return (train_indices, val_indices) splitting primers sequences into two groups."""

    rng = np.random.default_rng(seed)
    unique_primers_seqs = np.unique(primer_seqs)
    rng.shuffle(unique_primers_seqs)

    n_val = int(len(unique_primers_seqs) * val_fraction)
    #val_primers = set(unique_primers_seqs[:n_val]) BEFORE
    val_primers_seqs = unique_primers_seqs[:n_val]          # numpy slice — no set, no list

    # BEFORE
    #train_idx = np.where(~np.isin(site_ids, list(val_primers)))[0] 
    #val_idx   = np.where( np.isin(site_ids, list(val_primers)))[0]

    # Compute the boolean mask ONCE
    val_mask   = np.isin(primer_seqs, val_primers_seqs)

    # Derive both index arrays from the same mask
    val_idx   = np.where( val_mask)[0]
    train_idx = np.where(~val_mask)[0]

    return train_idx, val_idx


def kfold_by_primer_seq(primer_seqs: np.ndarray, n_folds: int = 5, seed: int = 42):
    """Yield (train_indices, val_indices) for K-fold splitting at the primer level."""
    rng = np.random.default_rng(seed)

    _, unique_primers = np.unique(primer_seqs, return_inverse=True)
    primer_seqs = primer_seqs.astype(np.int32)

    unique_primers_ids = np.arange(primer_seqs.max() + 1, dtype=np.int32)
    rng.shuffle(unique_primers_ids)

    fold_assignments = np.array_split(unique_primers_ids, n_folds)

    for fold_idx, val_primers in enumerate(fold_assignments):
        val_set = set(val_primers)
        val_mask = np.isin(primer_seqs, val_set)

        val_idx   = np.where(val_mask)[0]
        train_idx = np.where(~val_mask)[0]
        yield fold_idx, train_idx, val_idx

def kfold_by_panel(panel_ids: np.ndarray, n_folds: int = 5, seed: int = 42):
    """Yield (train_indices, val_indices) for K-fold splitting at the primer level."""
    rng = np.random.default_rng(seed)

    _, unique_panel = np.unique(panel_ids, return_inverse=True)
    panel_ids = panel_ids.astype(np.int32)

    unique_panel_ids = np.arange(panel_ids.max() + 1, dtype=np.int32)
    rng.shuffle(unique_panel_ids)

    fold_assignments = np.array_split(unique_panel_ids, n_folds)

    for fold_idx, val_panels in enumerate(fold_assignments):
        val_set = set(val_panels)
        val_mask = np.isin(panel_ids, val_set)

        val_idx   = np.where(val_mask)[0]
        train_idx = np.where(~val_mask)[0]
        yield fold_idx, train_idx, val_idx


# CHECK IF IT WORKS WHEN OTHER PANELS EXIST
def split_by_panel(panel_ids: np.ndarray, val_fraction: float = 0.2, seed: int = 42):
    """Return (train_indices, val_indices) splitting panels into two groups."""
    rng = np.random.default_rng(seed)

    #panel_ids =  np.concatenate(panel_ids) if isinstance(panel_ids, (list, np.ndarray)) else panel_ids
    
    unique_panels = np.unique(panel_ids)
    rng.shuffle(unique_panels)
    n_val = int(len(unique_panels) * val_fraction)
    val_panels = set(unique_panels[:n_val])

    # Compute the boolean mask ONCE
    val_mask   = np.isin(panel_ids, list(val_panels))

    # Derive both index arrays from the same mask
    val_idx   = np.where( val_mask)[0]
    train_idx = np.where(~val_mask)[0]

    return train_idx, val_idx