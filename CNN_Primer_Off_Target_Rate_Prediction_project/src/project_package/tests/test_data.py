#!/usr/bin/env python3
# To run the srcipt, activate the uv environment and change th end line from CRLF to LF (in vscode, bottom right corner) and then run:
# src/project_package/tests/test_data.py
"""
Tests for data loading and dataset behavior.
Run from project root: python src/project_package/tests/test_data.py
"""

import numpy as np
import torch
from pathlib import Path
from torch.utils.data import DataLoader

from project_package.config import Config
from project_package.data import (
    load_npz_files,
    AlignmentDataset,
    split_by_primer,
    kfold_by_primer,
    split_train_val_test
)


# ---------------------------------------------------------------------------
# Fixture: pick one .npz file to use for tests
# ---------------------------------------------------------------------------

cfg = Config()
ALL_NPZ_FILES = sorted(cfg.data_dir.glob("*.npz"))
assert len(ALL_NPZ_FILES) > 0, f"No .npz files found in {cfg.data_dir}"
SINGLE_FILE = [ALL_NPZ_FILES[0]]


# ---------------------------------------------------------------------------
# Test 1: load_npz_files reads the file correctly
# ---------------------------------------------------------------------------

def test_load_single_file():
    """Verify load_npz_files returns the expected keys and consistent array lengths."""
    data = load_npz_files(SINGLE_FILE)

    # Required keys
    assert 'images'     in data, "Missing 'images' key in loaded data"
    assert 'labels'     in data, "Missing 'labels' key in loaded data"
    assert 'primer_ids' in data, "Missing 'primer_ids' key in loaded data"

    # Optional but expected after regeneration
    assert 'is_target'  in data, "Missing 'is_target' — did you regenerate the .npz files?"

    # Length consistency: every array has one entry per alignment
    n = len(data['labels'])
    assert len(data['images'])     == n, f"images has {len(data['images'])} entries, expected {n}"
    assert len(data['primer_ids']) == n, f"primer_ids has {len(data['primer_ids'])} entries, expected {n}"
    assert len(data['is_target'])  == n, f"is_target has {len(data['is_target'])} entries, expected {n}"

    print(f"✓ Loaded {n:,} alignments from {SINGLE_FILE[0].name}")


# ---------------------------------------------------------------------------
# Test 2: image data has the expected shape and dtype
# ---------------------------------------------------------------------------

def test_image_shape_and_dtype():
    """Images should be (N, 5, 66, 2) or (N, 2, 5, 66), stored as uint8."""
    data = load_npz_files(SINGLE_FILE)
    images = data['images']

    # uint8 in storage — that's the whole point of the storage format
    assert images.dtype == np.uint8, \
        f"Expected images to be uint8 (storage format), got {images.dtype}"

    # 4D array: (N, height, width, channels) or (N, channels, height, width)
    assert images.ndim == 4, f"Expected 4D image array, got {images.ndim}D"

    # The non-batch dimensions should be (5, 66, 2) in some order
    other_dims = sorted(images.shape[1:])
    assert other_dims == [2, 5, 66], \
        f"Expected dimensions {{2, 5, 66}}, got {other_dims}"

    # uint8 means values must be in [0, 255]
    assert images.min() >= 0 and images.max() <= 255, \
        f"uint8 values out of range: [{images.min()}, {images.max()}]"

    print(f"✓ Image shape {images.shape}, dtype {images.dtype}, "
          f"value range [{images.min()}, {images.max()}]")


# ---------------------------------------------------------------------------
# Test 3: labels look like UMI counts
# ---------------------------------------------------------------------------

def test_labels_are_valid():
    """Labels should be non-negative integers (or floats), with a mix of zeros and positives."""
    data = load_npz_files(SINGLE_FILE)
    labels = data['labels']

    # All non-negative
    assert (labels >= 0).all(), \
        f"Found negative labels — UMI counts can't be negative. Min: {labels.min()}"

    # Most labels should be small or zero; some should be large
    n_zero = (labels == 0).sum()
    n_pos  = (labels > 0).sum()
    assert n_pos > 0, "No positive labels found — dataset has no signal"

    print(f"✓ Labels: min={labels.min()}, max={labels.max()}, "
          f"zero={n_zero:,}, positive={n_pos:,}, mean={labels.mean():.2f}")


# ---------------------------------------------------------------------------
# Test 4: AlignmentDataset returns correctly normalized tensors
# ---------------------------------------------------------------------------

def test_dataset_single_sample():
    """A single Dataset sample should be (2, 5, 66) float32 in [0, 1] paired with a scalar label."""
    data = load_npz_files(SINGLE_FILE)
    ds = AlignmentDataset(data['images'], data['labels'])

    x, y = ds[0]

    # Type check
    assert isinstance(x, torch.Tensor), f"Expected x to be a Tensor, got {type(x)}"
    assert isinstance(y, torch.Tensor), f"Expected y to be a Tensor, got {type(y)}"

    # Dtype check — model expects float32
    assert x.dtype == torch.float32, f"Expected float32 image, got {x.dtype}"
    assert y.dtype == torch.float32, f"Expected float32 label, got {y.dtype}"

    # Shape check — channels first, as the model expects
    assert x.shape == (2, 5, 66), f"Expected x shape (2, 5, 66), got {tuple(x.shape)}"
    assert y.shape == (), f"Expected y to be a scalar, got shape {tuple(y.shape)}"

    # Normalization check — values must be in [0, 1] after /255 scaling
    assert x.min() >= 0.0 and x.max() <= 1.0, \
        f"Image values out of [0, 1]: [{x.min():.3f}, {x.max():.3f}]"

    print(f"✓ Sample 0: x.shape={tuple(x.shape)}, y={y.item():.1f}, "
          f"x range [{x.min():.3f}, {x.max():.3f}]")


# ---------------------------------------------------------------------------
# Test 5: DataLoader produces correctly-shaped batches
# ---------------------------------------------------------------------------

def test_dataloader_batch():
    """A batch from the DataLoader should be (batch_size, 2, 5, 66) and (batch_size,)."""
    data = load_npz_files(SINGLE_FILE)
    ds = AlignmentDataset(data['images'], data['labels'])
    loader = DataLoader(ds, batch_size=32, shuffle=False)

    x_batch, y_batch = next(iter(loader))

    assert x_batch.shape == (32, 2, 5, 66), f"Got x batch shape {tuple(x_batch.shape)}"
    assert y_batch.shape == (32,),           f"Got y batch shape {tuple(y_batch.shape)}"
    assert x_batch.dtype == torch.float32, f"Expected float32 x batch, got {x_batch.dtype}"
    assert y_batch.dtype == torch.float32, f"Expected float32 y batch, got {y_batch.dtype}"

    print(f"✓ Batch shapes: x={tuple(x_batch.shape)}, y={tuple(y_batch.shape)}")


# ---------------------------------------------------------------------------
# Test 6: iterate over 10 samples without errors
# ---------------------------------------------------------------------------

def test_iterate_10_samples():
    """Verify __getitem__ works repeatedly without raising or producing inconsistent shapes."""
    data = load_npz_files(SINGLE_FILE)
    ds = AlignmentDataset(data['images'], data['labels'])

    for i in range(10):
        x, y = ds[i]
        assert x.shape == (2, 5, 66)
        assert isinstance(y.item(), float)

    print(f"✓ Iterated over 10 samples without error")


# ---------------------------------------------------------------------------
# Test 7: split_by_primer produces non-overlapping primer sets
# ---------------------------------------------------------------------------

def test_split_by_primer_no_leakage():
    """The most important test in this file. Train and val sets must share NO primers."""
    data = load_npz_files(SINGLE_FILE)
    primer_ids = data['primer_ids']

    train_idx, val_idx = split_by_primer(primer_ids, val_fraction=0.2, seed=42)

    # No site appears in both splits
    assert len(set(train_idx) & set(val_idx)) == 0, \
        "Site-level overlap between train and val — split is broken"

    # No primer appears in both splits — THIS is the critical anti-leakage check
    train_primers = set(primer_ids[train_idx])
    val_primers   = set(primer_ids[val_idx])
    overlap = train_primers & val_primers
    assert len(overlap) == 0, \
        f"Primer leakage: {len(overlap)} primers appear in both train and val"

    # Sanity: roughly the expected fraction (allow some tolerance because we split by primer, not site)
    val_frac_actual = len(val_idx) / len(primer_ids)
    assert 0.10 <= val_frac_actual <= 0.30, \
        f"Validation fraction {val_frac_actual:.2%} is far from expected 20%"

    print(f"✓ Primer-level split: {len(train_idx):,} train / {len(val_idx):,} val "
          f"({val_frac_actual:.1%}), no primer overlap")


# ---------------------------------------------------------------------------
# Test 8: kfold_by_primer produces clean disjoint folds
# ---------------------------------------------------------------------------

def test_kfold_no_leakage():
    """Every fold's val set must be primer-disjoint from its train set."""
    data = load_npz_files(SINGLE_FILE)
    primer_ids = data['primer_ids']

    seen_val_primers = set()
    n_folds = 5

    for fold_idx, train_idx, val_idx in kfold_by_primer(primer_ids, n_folds=n_folds, seed=42):
        train_primers = set(primer_ids[train_idx])
        val_primers   = set(primer_ids[val_idx])

        # Anti-leakage within each fold
        overlap = train_primers & val_primers
        assert len(overlap) == 0, \
            f"Fold {fold_idx}: primer leakage detected ({len(overlap)} primers in both splits)"

        # Across folds, each primer should appear in exactly one val set
        repeat = seen_val_primers & val_primers
        assert len(repeat) == 0, \
            f"Fold {fold_idx}: primers appear as val in multiple folds: {len(repeat)}"
        seen_val_primers |= val_primers

    print(f"✓ {n_folds}-fold CV: each fold primer-disjoint, all primers used exactly once for val")


# ---------------------------------------------------------------------------
# Test 9: end-to-end — feed real data into the model
# ---------------------------------------------------------------------------

def test_real_data_through_model():
    """Take real samples, batch them, and verify the model accepts them without errors."""
    from project_package.model import PordleCNN

    data = load_npz_files(SINGLE_FILE)
    ds = AlignmentDataset(data['images'][:100], data['labels'][:100])
    loader = DataLoader(ds, batch_size=32, shuffle=False)

    cfg = Config()
    model = PordleCNN(cfg)
    model.eval()

    x_batch, y_batch = next(iter(loader))
    with torch.no_grad():
        y_hat = model(x_batch)

    assert y_hat.shape == (32,), f"Model output shape wrong: {tuple(y_hat.shape)}"
    assert (y_hat >= 0).all(), "Model produced negative predictions"

    print(f"✓ Real data passes through model: predictions in [{y_hat.min():.2f}, {y_hat.max():.2f}]")


# ---------------------------------------------------------------------------
# Test 10: confirm the three sets are primer-disjoint
# ---------------------------------------------------------------------------

def test_three_way_split_no_leakage():
    """Train, val, and test sets must be pairwise primer-disjoint."""
    data = load_npz_files(SINGLE_FILE)
    primer_ids = data['primer_ids']

    train_idx, val_idx, test_idx = split_train_val_test(
        primer_ids, val_fraction=0.15, test_fraction=0.15, seed=42
    )

    train_primers = set(primer_ids[train_idx])
    val_primers   = set(primer_ids[val_idx])
    test_primers  = set(primer_ids[test_idx])

    assert len(train_primers & val_primers)  == 0, "train/val leakage"
    assert len(train_primers & test_primers) == 0, "train/test leakage"
    assert len(val_primers   & test_primers) == 0, "val/test leakage"

    print(f"✓ Three-way split: {len(train_primers)} train / "
          f"{len(val_primers)} val / {len(test_primers)} test primers, "
          f"all disjoint")
    

# ---------------------------------------------------------------------------
# Run all tests
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=== Testing data.py ===\n")
    test_load_single_file()
    test_image_shape_and_dtype()
    test_labels_are_valid()
    test_dataset_single_sample()
    test_dataloader_batch()
    test_iterate_10_samples()
    test_split_by_primer_no_leakage()
    test_kfold_no_leakage()
    test_real_data_through_model()
    test_three_way_split_no_leakage()
    print("\n=== All tests passed ===")