#!/usr/bin/env python3

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import islice
import json
import pandas as pd
import matplotlib.pyplot as plt
import json
import os
from tqdm import tqdm
from Bio import Align
import numpy as np

def prepare_for_encoding(primer_align, gen_seq_align, strand, gap="-"):
    """Strip extraction-buffer gaps and orient so the primer's 3' end is on the left."""
    # 1. Strip leading/trailing columns where the primer is gapped (template buffer)
    left  = len(primer_align) - len(primer_align.lstrip(gap))
    right = len(primer_align) - len(primer_align.rstrip(gap))
    end   = len(primer_align) - right if right else len(primer_align)
    p, t  = primer_align[left:end], gen_seq_align[left:end]

    # 2. Forward strand has 3' on the right -> reverse to put it on the left.
    #    Reverse strand already has 3' on the left -> leave alone. Updated on 4 May 26: 
    # Now the reverse primer/template have been reverted back before the alignment therefore now are 5'->3' oriented too.
    # Hence, the conditional is removed and the reverting is applied on both sense.
    # Ask Manfred's opinion since maybe it would not needed to do this reverse step since it is an agreement to CNN to have a persistance sense (orientation).
    
    p, t = p[::-1], t[::-1]
    
    return p, t

import numpy as np
 
# Row order: ACGT then deletion
ROWS = ['A', 'C', 'G', 'T', '-']
ROW_INDEX = {b: i for i, b in enumerate(ROWS)}
 
# Encoding values (as paper of reference described at Section 4.2)
CORRECT_VAL = 250 
OTHER_VAL   = 100
LEFT_PAD    = 4     # "padded 4 bases of zeros to the left of all images"
TARGET_W    = 66    # final image width (from Fig. 6: Input is 5x66x2)
 
 
def encode_sequence(seq: str) -> np.ndarray:
    """Turn a single aligned sequence (with '-' for gaps) into a 5xL matrix.
    The correct row gets 250, the others get 100."""
    L = len(seq)
    M = np.full((5, L), OTHER_VAL, dtype=np.int16)
    for j, ch in enumerate(seq):
        ch = ch.upper()
        if ch not in ROW_INDEX:
            raise ValueError(f"Unexpected character {ch!r} at position {j}")
        M[ROW_INDEX[ch], j] = CORRECT_VAL
    return M
 
 
def encode_alignment(primer_aln: str, template_aln: str,
                     target_width: int = TARGET_W,
                     left_pad: int = LEFT_PAD) -> np.ndarray:
    """Produce the (5, target_width, 2) image array as in the paper.
 
    Channel 0: primer
    Channel 1: template
 
    NOTE: both inputs must be the same length (the alignment length) and must
    be oriented so that the 3' end of the primer is on the LEFT of the string.
    See Fig. 1b in the paper - the right side is the 'less important' 5' end
    of the primer. If your ntthal output is the other way around, reverse both
    strings before calling this function.
    """
    if len(primer_aln) != len(template_aln):
        raise ValueError("Primer and template alignments must have equal length")
 
    # Step 1: encode each sequence into a 5xL matrix
    primer_mat   = encode_sequence(primer_aln)     # (5, L)
    template_mat = encode_sequence(template_aln)   # (5, L)
 
    # Step 2: stack as 2 channels -> (5, L, 2)
    img = np.stack([primer_mat, template_mat], axis=-1)
    L = img.shape[1]
 
    # Step 3: left-pad with 4 zero columns (3' side) and right-pad/trim to width 66
    needed_right = target_width - left_pad - L
    if needed_right >= 0:
        img = np.pad(img,
                     pad_width=((0, 0), (left_pad, needed_right), (0, 0)),
                     mode='constant', constant_values=0)
    else:
        # Trim from the right (5' end of primer) and then add only the left pad
        keep = target_width - left_pad
        img = img[:, :keep, :]
        img = np.pad(img,
                     pad_width=((0, 0), (left_pad, 0), (0, 0)),
                     mode='constant', constant_values=0)
 
    # Step 4: scale to [0, 1] by dividing by 255
    #img = img.astype(np.float32) / 255.0
    img_uint8 = img.astype(np.uint8)   # values 0, 100, 250

    # When loading for training:
    #X = data["images"].astype(np.float32) / 255.0

    return img_uint8

Alignfile_names = []
with open("/home/marta/SPE_alignments_project/data/alignments_json_files_names.txt") as f:
    for line in f:
        line = line.strip()  # remove trailing newline and whitespace
        Alignfile_names.append(line)
    
for name in Alignfile_names:
     # Load the Alignments
    with open(f"/home/marta/SPE_alignments_project/src/data/ResultFiles_marta/Align_collap_v3_1slip_RevToFw/{name}", "r") as f: 
        AlignAligner_collapse_v3_1slippage= json.load(f)
    print(name) # name example: DHS-001Z_20220525_HA_S1_R_AlignAligner_collapsed_v3_1slippage_RevertingRevToFw.json

    image_list = []
    labels_umiC = []
    meta_list = []
    for chr, strands in tqdm(AlignAligner_collapse_v3_1slippage.items(), desc="Processing chromosome"):

        for strand, locations in strands.items():
            
            for loc, site in locations.items():

                umi_count = site["umi_count"]

                for key, val in site.items():
                    if key == "umi_count":
                        continue

                    primer_seq = key 

                    for aln_key, (tpl, prm) in val.items(): # level to iterate to prepare the prm and tpl to later encode them into an image. 

                        primer, template = prepare_for_encoding(prm, tpl, strand)
                        image = encode_alignment(primer, template)      # (5, 66, 2), 3' on the left
                    # ... pair with umi_count and feed to the CNN

                        image_list.append(image) # an image of each alignment
                        labels_umiC.append(umi_count)
                        meta_list.append({"primer": primer_seq, "chr": chr,
                                    "loc": int(loc), "strand": strand})

        npz_name = name.replace("AlignAligner_collapsed_v3_1slippage_RevertingRevToFw.json", "img_and_labels_uint8.npz")  
        npz_folder = "/home/marta/SPE_alignments_project/src/data/img_labels_npz/"

        # Build the full path
        full_path_npz = os.path.join(npz_folder, npz_name)

        np.savez_compressed(
            full_path_npz,
            images=np.stack(image_list, axis=0),       # (N, 5, 66, 2) float32
            labels=np.array(labels_umiC, dtype=np.float32),
        )

        parquet_name = name.replace("AlignAligner_collapsed_v3_1slippage_RevertingRevToFw.json", "metadata.parquet")  
        parquet_folder = "/home/marta/SPE_alignments_project/src/data/metadata_parquet/"

        # Build the full path
        full_path_parquet = os.path.join(parquet_folder, parquet_name)

        pd.DataFrame(meta_list).to_parquet(full_path_parquet)

