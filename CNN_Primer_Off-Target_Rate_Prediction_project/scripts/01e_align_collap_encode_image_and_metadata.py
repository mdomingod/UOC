#!/usr/bin/env python3

# 24/05/2026 - 

# This script performs the alignment of the primers and genomic sequences, 
# collapses the alignments by summing the umi counts for each unique alignment, 
# encodes the collapsed alignments into images and saves them in npz format.

# Each npz file also contains the labels (umi counts) for each image, as well as the primer ids and 
# whether the alignment is on-target or off-target.
#
# It also generates a metadata json file for each npz file containing information about the encoding parameters 
# and the number of alignments. 
# 
# The script processes multiple csv files containing the alignments and epcr07 data, 
# concatenates them, and then applies the alignment, collapsing, 
# and encoding steps to produce the final images and labels for training a CNN model. 
# 
# The images are saved as uint8 to save space, and the labels are saved as float32.

# 15_CollapsingAlignments.ipynb contains the code for collapsing the alignments and summing the umi counts for each unique alignment, which is then used in this script to produce the final images and labels.


from datetime import datetime
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import islice
import os
import json
import pandas as pd
from collections import Counter
from collections import defaultdict
from tqdm import tqdm
from Bio import Align
from pathlib import Path
from data_quality_report import compute_mixed_labels_report, save_report, render_markdown_report

""" # When loading images for training:
    # X = data["images"].astype(np.float32) / 255.0
    # Above is the step to do when loading the images for training,
    # not here since we want to save them as uint8 to save space and only convert to float32 and scale when loading for training."""


# 1. Alignments function definitions
# Needleman Wunch Algorithm alignment for primers and genomic sequence function definition)
def calculate_cigar_nwa(primer, read):
    """
    Semi-global alignment between primer and read.
    The full primer is forced to align (no soft-clipping),
    while the read has free end gaps (primer can land anywhere).
    Models PCR annealing where the primer binds with its full length.
    Returns CIGAR string and error count (mismatches + gaps across the entire primer).

    Implementation: Needleman-Wunsch global alignment recursion modified
    to behave as semi-global by zeroing end-gap penalties on the read only.

    - mode = "global" applies the Needleman-Wunsch global alignment recursion,
      which forces both sequences to participate in the alignment end-to-end.
    - Setting free end-gap penalties on the read only overrides this for the read,
      making the read's terminal gaps unpenalized (semi-global behaviour).
    - The primer retains full-length forced alignment (it cannot be soft-clipped),
      because its end-gap penalties remain penalized.
    - The read may contain extra sequence at either end, modelling how a primer
      can anneal at any internal location of a read.
    - counts() reflects mismatches and gaps across the entire primer,
      because none of the primer bases can be clipped away (unlike Smith-Waterman,
      which may hide mismatches by localizing the alignment to the best-scoring
      subsequence and trimming misaligned ends).

    """

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    # Free end gaps on target (read) only
    aligner.open_left_insertion_score = 0
    aligner.extend_left_insertion_score = 0
    aligner.open_right_insertion_score = 0
    aligner.extend_right_insertion_score = 0

    alignments = aligner.align(read, primer)

    if not alignments:
        return "*", len(primer)

    alignment = alignments[0]

    # CIGAR from built-in SAM formatter (6th field)
    cigar = alignment.format("sam").split("\t")[5]

    # Error count from built-in counts()
    c = alignment.counts()
    err = c.mismatches + c.gaps

    return cigar, err, alignment


# The alignments funcion returns the aligned primer and genomic sequence, 
# which are then fed to the encoding function to produce the image.
def alignment(primer, genomic_seq):
    
    cigar, err, alignment = calculate_cigar_nwa(primer, genomic_seq)
    align_tuple = str(alignment).split("0 ")[1].split()[0], str(alignment).split("0 ")[3].split()[0]
    gsq_aligned = align_tuple[0] # genomic seq
    prm_aligned = align_tuple[1] # primer

    return prm_aligned, gsq_aligned



# 2. Encoding to image function definitions
def prepare_for_encoding(primer_align, gen_seq_align, gap="-"):
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


def encoding_to_image(prm_aligned, gsq_aligned):
    
    primer, gen_seq = prepare_for_encoding(prm_aligned, gsq_aligned)
    image = encode_alignment(primer, gen_seq)      # (5, 66, 2), 3' on the left
# ... pair with umi_count and feed to the CNN

    return image


def extract_panel_id(name):
    if name[7] == "_":
        return name[:7]
    if name[8] == "_" and name[:8] != "DHS5000Z":
        return name[:8]
    if name[9] == "_":
        return name[:9]
    if name[:8] == "DHS5000Z":
        return "DHS-5000Z"
    return None


csv_file_names = []
with open("/home/domingom/projects/SPE_Primer_Alignments/data/need_to_be_combined/csv_align_ONLYnames.txt") as f:
    for line in f:
        line = line.strip()  # remove trailing newline and whitespace
        csv_file_names.append(line)


for name in tqdm(csv_file_names, desc='processing the file'):

    #Load data
    alignments_file = f'/home/domingom/projects/SPE_Primer_Alignments/data/need_to_be_combined/raw/{name}_alignments_from_data.csv'
    data_to_align = pd.read_csv(alignments_file, sep=",")

    alignments_file = f'/home/domingom/projects/SPE_Primer_Alignments/data/need_to_be_combined/raw/{name}_EPCR07_not_in_data.csv'
    data_to_align_epcr07 = pd.read_csv(alignments_file, sep=",")


    # Concatenate alignments and epcr07 files.
    combined_data = pd.concat([data_to_align, data_to_align_epcr07], ignore_index=True)




    ALLOWED = set("ACGT")  # no '-' here because input sequences before alignment shouldn't have gaps

    # Identify rows where either primer or genomic_seq contains a non-ACGT character
    bad_primer  = ~combined_data["primer"].str.upper().str.fullmatch(r"[ACGT]+")
    bad_genomic = ~combined_data["genomic_seq"].str.upper().str.fullmatch(r"[ACGT]+")
    bad_mask    = bad_primer | bad_genomic

    n_dropped = bad_mask.sum()
    print(f"Dropping {n_dropped} rows with non-ACGT characters "
        f"({100 * n_dropped / len(combined_data):.2f}%)")
    print(f"  - {bad_primer.sum()} rows due to primer")
    print(f"  - {bad_genomic.sum()} rows due to genomic_seq")
    print(f"  - {(bad_primer & bad_genomic).sum()} rows due to both")

    combined_data = combined_data[~bad_mask].reset_index(drop=True)

    """ .str.fullmatch(r"[ACGT]+") returns True only if the entire string consists of A/C/G/T characters. 
    Anything else — N, lowercase letters, whitespace, IUPAC ambiguity codes like Y or R — fails the match. 
    The leading ~ inverts it, so bad_primer is True when the primer contains something unexpected. """

    print("Sample primers with non-ACGT:")
    print(combined_data.loc[bad_primer, "primer"].head(5).tolist())

    print("\nSample genomic_seqs with non-ACGT:")
    print(combined_data.loc[bad_genomic, "genomic_seq"].head(5).tolist())



    # Executing the alignment for all the data
    for index, row in tqdm(combined_data.iterrows(), total=combined_data.shape[0]):
        primer = row['primer']
        genomic_seq = row['genomic_seq']
        prm_aligned, gsq_aligned = alignment(primer, genomic_seq)
        combined_data.at[index, 'prm_aligned'] = prm_aligned
        combined_data.at[index, 'gsq_aligned'] = gsq_aligned


    # Collapsing the alignments by unique primer-genomic sequence combinations, summing the UMI group counts and counting how many original sites were collapsed into each unique combination. 
    # Also keeping track of whether any of the original sites in each group were on-target (isIntendedSite=1).
    collapsed_alignments = (
        combined_data
        .groupby(["primer", "prm_aligned", "gsq_aligned"], as_index=False)
        .agg(
            umi_groups_cnt=('umi_groups_cnt', 'sum'),
            read1_anchor=('read1_anchor', 'first'),
            alignChrom=('alignChrom', 'first'),
            is_target=('isIntendedSite', 'max'),
            n_genomic_occurrences=('umi_groups_cnt', 'count'),
        )
        .rename(columns={'primer': 'primer_seq'})   # ← normalize the name for downstream code
        .sort_values(by='umi_groups_cnt', ascending=False)
    )


    # Data quality report: identify how many unique (prm_aligned, gsq_aligned) groups have mixed on-target/off-target labels among their contributing rows, 
    # and provide a per-group breakdown of the label distribution and genomic positions.
    report = compute_mixed_labels_report(combined_data)
    save_report(
        report,
        output_dir=Path(f"/home/domingom/projects/SPE_Primer_Alignments/reports/{name}_data_quality_reports_01e"),
        dataset_name=f"{name} collapsed alignments",
    )


    images_labels_list = []
    for index, row in tqdm(collapsed_alignments.iterrows(), total=collapsed_alignments.shape[0]):
        prm_aligned = row['prm_aligned']
        gsq_aligned = row['gsq_aligned']

        images_labels_list.append((encoding_to_image(prm_aligned, gsq_aligned), row['umi_groups_cnt']))


    npz_name = f'{name}_images_and_labels_uint8_.npz'
    npz_folder = "/home/domingom/projects/SPE_Primer_Alignments/data/need_to_be_combined/processed/images_labels_CollapsedAlign_npz"

    # Build the full path
    full_path_npz = os.path.join(npz_folder, npz_name)

    np.savez_compressed(
        full_path_npz,
        images=np.stack([img for img, _ in images_labels_list], axis=0),       # (N, 5, 66, 2) float32
        labels=np.array([label for _, label in images_labels_list], dtype=np.float32), 
        site_ids=tuple(
            (str(a) if isinstance(a, np.ndarray) else a,
             str(b) if isinstance(b, np.ndarray) else b)
            for a, b in zip(
                collapsed_alignments['alignChrom'].values,
                collapsed_alignments['read1_anchor'].values)),  # list of tuples (alignChrom, read1_anchor)
        is_target=collapsed_alignments['is_target'].values.astype(np.uint8),
        primer_seq=collapsed_alignments['primer_seq'].values.astype(object),  # object dtype for variable-length strings
        panel_id=extract_panel_id(name)
)
    
    
    """
    line 132: CORRECT_VAL = 250 # "the correct row gets 250, the others get 100" as described in the paper at Section 4.2
    line 133: OTHER_VAL   = 100 # the others get 100" 
    line 134: LEFT_PAD     = 4     # "padded 4 bases of zeros to the left of all images"
    line 135: TARGET_W     = 66    # final image width, from paper is 66 (Fig. 6: Input is 5x66x2)

    """

    metadata = {
        "img_height": 5,
        "img_width": TARGET_W,
        "img_channels": 2,
        "left_pad": LEFT_PAD,
        "on_value": CORRECT_VAL,
        "off_value": OTHER_VAL,
        "encoding_script_version": "01e_align_collap_encode_image_and_metadata.py",
        "encoding_date": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        "n_alignments": len(collapsed_alignments)
    }


    metadata_json_name = f'{name}_metadata.json'
    full_path_metadata_json = os.path.join(npz_folder, metadata_json_name)
    with open(full_path_metadata_json, 'w') as f:
        json.dump(metadata, f, indent=4)
