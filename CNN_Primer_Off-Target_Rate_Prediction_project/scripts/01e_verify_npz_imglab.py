#!/usr/bin/env python3

# 19/05/2026

import numpy as np
from pathlib import Path

from tqdm import tqdm

def verify_dataset(path, expected_pixel_values={0, 100, 250}):
    data = np.load(path, allow_pickle=True)
    imgs, lbls, site_ids, is_target, panel_id, primer_seq = data["images"], data["labels"], data["site_ids"], data["is_target"], data["panel_id"], data['primer_seq']

    assert imgs.shape[0] == lbls.shape[0], "N mismatch"
    assert imgs.dtype == np.uint8, f"images dtype is {imgs.dtype}, expected uint8"
    assert set(np.unique(imgs).tolist()).issubset(expected_pixel_values), \
        "unexpected pixel values"
    assert (lbls >= 0).all(), "negative labels"

    # full transform
    x = imgs.astype(np.float32) / 255.0
    assert 0.0 <= x.min() and x.max() <= 1.0
    assert not np.isnan(x).any() and not np.isinf(x).any()

    print(f"OK — {imgs.shape[0]} examples, "
          f"file size {Path(path).stat().st_size/1024**2:.1f} MB, "
          f"labels in [min: {lbls.min()}, max: {lbls.max()}]")
    
    primer_seq = [str(pid) for pid in primer_seq]
    site_ids = [str(pid) for pid in site_ids]  # convert to strings for easier checking
    # assert len(set(primer_ids)) == len(primer_ids), "duplicate primer_ids"
    assert set(is_target).issubset({0, 1}), "is_target should be binary"
    # assert len(set(panel_id)) == 1, "multiple panel_ids in one file"

    print(f"OK — Sites IDs: {len(set(site_ids))} unique,"
          f"OK — Primer: {len(set(primer_seq))} unique," 
          f"Is Target distribution: {np.bincount(is_target)},"
          f"Panel ID: {panel_id}")


csv_file_names = []
with open("/home/domingom/projects/SPE_Primer_Alignments/data/need_to_be_combined/images_name_npz.txt") as f:
    for line in f:
        line = line.strip()  # remove trailing newline and whitespace
        csv_file_names.append(line)


for name in tqdm(csv_file_names, desc='processing the file'):

    verify_dataset(f"/home/domingom/projects/SPE_Primer_Alignments/data/need_to_be_combined/processed/images_labels_CollapsedAlign_npz/{name}")