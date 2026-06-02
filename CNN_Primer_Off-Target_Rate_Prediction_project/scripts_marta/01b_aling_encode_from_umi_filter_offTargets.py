#!/usr/bin/env python3

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

BASES = ["A", "C", "G", "T"]

ref_hg19 = pysam.FastaFile("/home/marta/SPE_alignments_project/data/genome/hg19.fa")


# Load alignment file
alignments_file = '/home/marta/SPE_alignments_project/data/umi_filter_alignments/DHS-001Z_20220525_HA_S1_R.umi_filter.alignments.txt'

aligns = pd.read_csv(alignments_file, sep="|")

aligns.columns = aligns.columns.str.replace('.', '_')

aligns.columns = ['primerChrom', 'loc5', 'primerStrand', 'primer', 'umiSeq', 'isIntendedSite', 'alignChrom', 'alignStrand', 'alignLocRand', 'alignLoc', 'readId', 'read1_pos', 'read1_aend', 'cigar1', 'read2_pos', 'read2_aend', 'cigar2']


# Load off-target
off_targets = aligns[aligns['isIntendedSite'] == 0]

off_target = off_targets.iloc[0]

off_target_primer_end = off_target['read1_pos']+8
off_target_primer_start = off_target_primer_end - len(off_target['primer'])

# genomic sequence generation
genomic_seq = ref_hg19.fetch(off_target['alignChrom'], off_target_primer_start, off_target_primer_end)


# UMI cluestering function definition
def hamming_distance(a, b):
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def cluster_umis(off_targets):
    """Greedy UMI clustering: merge UMIs within Hamming distance <= 1.

    Parameters
    ----------
    off_targets : pandas DataFrame

    Returns
    -------
    list of (centroid_umi: str, reads: list[dict])
        One entry per collapsed molecule.
    """
    
    print('preparing UMI clusters...')
    umi_counts = defaultdict(int)
    umi_to_reads = defaultdict(list)
    for read_data in off_targets.itertuples():
        umiSeq = read_data.umiSeq
        umi_counts[umiSeq] += 1
        umi_to_reads[umiSeq].append(read_data)

    # Process UMIs from most to least frequent (greedy)
    sorted_umis = sorted(umi_counts, key=lambda u: umi_counts[u], reverse=True)

    #print('assigning UMI clusters...')
    assigned = {}   # umi -> centroid umi
    for umi in sorted_umis:
        if umi in assigned:
            continue
        # This UMI becomes a centroid; absorb all unassigned neighbours
        assigned[umi] = umi
        for other in sorted_umis:
            if other in assigned:
                continue
            if hamming_distance(umi, other) <= 1:
                assigned[other] = umi

    # Collect reads per centroid
    #print('collecting reads per centroid...')
    centroid_reads = defaultdict(list)
    for umi, reads in umi_to_reads.items():
        if assigned[umi] in centroid_reads:
            centroid_reads[assigned[umi]]['umi_group_size'] += len(reads)
        else:
            centroid_reads[assigned[umi]] = pd.DataFrame([reads[0]._asdict()])
            centroid_reads[assigned[umi]]['umi_group_size'] = len(reads)

    umi_groups = pd.concat(centroid_reads.values(), ignore_index=True)
    
    return umi_groups

columns_umi_collapse = ['primer', 'umiSeq', 'alignChrom', 'alignStrand', 'read1_pos', 'cigar1']

#Performing UMI collapse
off_targets_per_position_list = [group for _, group in off_targets.groupby(['primer', 'alignChrom', 'alignStrand'])]
before_cluster_sizes =[x.shape[0] for x in off_targets_per_position_list]
largest_umi_group = off_targets_per_position_list[before_cluster_sizes.index(max(before_cluster_sizes))]

off_targets_per_position_list_clustered = [cluster_umis(group) for group in off_targets_per_position_list]

umi_groups = pd.concat(off_targets_per_position_list_clustered, ignore_index=True)

# Addition of the genomic sequence (with some padding)
def get_genomic_seq(row):
    alignChrom = row.alignChrom
    alignStrand = row.alignStrand
    primer = row.primer
    read1_pos = row.read1_pos
    read1_aend = row.read1_aend
    if alignStrand == 0:
        gen_seq = ref_hg19.fetch(alignChrom, read1_pos+3-len(primer), read1_pos+13) # extra padding +/- 5bp on each side
    elif alignStrand == 1:
        gen_seq = ref_hg19.fetch(alignChrom, read1_aend-13, read1_aend-3+len(primer)) # extra padding +/- 5bp on each side
        gen_seq = Seq(gen_seq).reverse_complement()  # Reverse complement for reverse strand
    gen_seq = str(gen_seq).upper()
    return gen_seq

umi_groups['genomic_seq'] = umi_groups.apply(lambda row: get_genomic_seq(row), axis=1)

# identical primer collapsing :genome pairs and count how many umi_groups support this.
umi_groups_final = umi_groups[['primer', 'genomic_seq']].groupby(['primer', 'genomic_seq']).aggregate(len).reset_index()
umi_groups_final.columns = ['primer', 'genomic_seq', 'umi_groups_cnt']


# # Needleman Wunch Algorithm alignment for primers and genomic sequence function definition)
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

# Encoding to image function definition
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



# Alignment and Encoding execution.
image_list = []
labels_umiC = []


po = 0
for row in umi_groups_final.itertuples():
    
    primer = row[1] 
    
    gen_seq = row[2] 

    umi_count = row[3]
    
    
    cigar, err, align_fw = calculate_cigar_nwa(primer, gen_seq)
    align_tuple_fw = str(align_fw).split("0 ")[1].split()[0], str(align_fw).split("0 ")[3].split()[0]

    print(align_tuple_fw)

    gsq = align_tuple_fw[0] # genomic seq
    prm = align_tuple_fw[1] # primer

    primer, gen_seq = prepare_for_encoding(prm, gsq)
    image = encode_alignment(primer, gen_seq)      # (5, 66, 2), 3' on the left
# ... pair with umi_count and feed to the CNN

    image_list.append(image) # an image of each alignment
    labels_umiC.append(umi_count)

    po +=  1
    if po > 5:
        break


    #npz_name = name.replace("AlignAligner_collapsed_v3_1slippage_RevertingRevToFw.json", "img_and_labels_uint8.npz")  
    npz_name = 'DHS-001Z_20220525_HA_S1_R_imgages_and_labels_uint8_firsts5.npz'
    npz_folder = "/home/marta/SPE_alignments_project/src/data/images_and_labels_npz/"

    # Build the full path
    full_path_npz = os.path.join(npz_folder, npz_name)

    np.savez_compressed(
        full_path_npz,
        images=np.stack(image_list, axis=0),       # (N, 5, 66, 2) float32
        labels=np.array(labels_umiC, dtype=np.float32),
    )

    