#!/usr/bin/env python3
"""
Script does 3 steps
. Load and extract read info from alignments file. 
. Collapse UMI groups.
. Extract genomic sequence for the primer match position.

Usage:
    python extract_from_alignments_file.py <input_alignments_file> [--output OUTPUT_FILE]

- <input_alignments_file>: Path to the umi_filter.alignments.txt file (pipe-separated columns)
- --output OUTPUT_FILE: Optional output CSV file (default: aligns_umi_grouped.csv)

Requires: pandas, pysam, biopython
"""
import argparse
import pandas as pd
import pysam
from Bio.Seq import Seq
from collections import defaultdict
from tqdm.auto import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

import multiprocessing
    
# ---- Argument parsing ----
def parse_args():
    parser = argparse.ArgumentParser(description="Collapse UMI and extract genomic sequences.")
    parser.add_argument("input", help="Input alignments .txt file (pipe-separated)")
    parser.add_argument("--output", default="aligns_umi_grouped.csv", help="Output CSV file (default: aligns_umi_grouped.csv)")
    parser.add_argument("--ref", default="../data/genome/hg19.fa", help="Reference genome FASTA (default: ../data/genome/hg19.fa)")
    return parser.parse_args()

# ---- UMI clustering ----
def hamming_distance(a, b):
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(c1 != c2 for c1, c2 in zip(a, b))

def cluster_umis(off_targets):
    """Greedy UMI clustering: merge UMIs within Hamming distance <= 1."""
    umi_counts = defaultdict(int)
    umi_to_reads = defaultdict(list)
    for read_data in off_targets.itertuples():
        umiSeq = read_data.umiSeq
        umi_counts[umiSeq] += 1
        umi_to_reads[umiSeq].append(read_data)
    sorted_umis = sorted(umi_counts, key=lambda u: umi_counts[u], reverse=True)
    assigned = {}   # umi -> centroid umi
    for umi in sorted_umis:
        if umi in assigned:
            continue
        assigned[umi] = umi
        for other in sorted_umis:
            if other in assigned:
                continue
            if hamming_distance(umi, other) <= 1:
                assigned[other] = umi
    centroid_reads = {}
    for umi, reads in umi_to_reads.items():
        centroid = assigned[umi]
        if centroid in centroid_reads:
            centroid_reads[centroid]['umi_group_size'] += len(reads)
        else:
            df = pd.DataFrame([reads[0]._asdict()])
            df['umi_group_size'] = len(reads)
            centroid_reads[centroid] = df
    umi_groups = pd.concat(centroid_reads.values(), ignore_index=True)
    return umi_groups

def umi_cluster_and_collapse(df):
    '''After umi clustering:
    collapse to unique primer:genome position
    and count how many umi groups support each unique position
    '''
    umi_grouped = cluster_umis(df)
    umi_grouped['read1_anchor'] = umi_grouped.apply(lambda row: row['read1_pos'] if row['alignStrand'] == 0 else row['read1_aend'], axis=1)
    umi_grouped_collapsed = umi_grouped[['primer', 'umiSeq', 'isIntendedSite', 'alignChrom', 'alignStrand', 'read1_anchor']].groupby(
        ['primer', 'isIntendedSite', 'alignChrom', 'alignStrand', 'read1_anchor']
    ).aggregate(len).reset_index()
    umi_grouped_collapsed.rename(columns={'umiSeq': 'umi_groups_cnt'}, inplace=True)
    return umi_grouped_collapsed

# ---- Genomic sequence extraction ----
def get_genomic_seq(row, ref_hg19):
    alignChrom = row.alignChrom
    alignStrand = row.alignStrand
    primer = row.primer
    read1_anchor = row.read1_anchor
    if alignStrand == 0:
        gen_seq = ref_hg19.fetch(alignChrom, read1_anchor-2-len(primer), read1_anchor+18)
    elif alignStrand == 1:
        gen_seq = ref_hg19.fetch(alignChrom, read1_anchor-18, read1_anchor+2+len(primer))
        gen_seq = Seq(gen_seq).reverse_complement()
    gen_seq = str(gen_seq).upper()
    return gen_seq

# Use ThreadPoolExecutor for parallel processing
def process_group(group):
    return umi_cluster_and_collapse(group)

def submit_groups(executor, aligns):
    futures = {}
    for _, group in aligns.groupby(['primer', 'alignChrom', 'alignStrand']):
        f = executor.submit(umi_cluster_and_collapse, group)
        futures[f] = None
    return futures

    
# ---- Main workflow ----
def main():
    args = parse_args()
    aligns_column_names = [
        'primerChrom', 'loc5', 'primerStrand', 'primer', 'umiSeq', 'isIntendedSite',
        'alignChrom', 'alignStrand', 'alignLocRand', 'alignLoc', 'readId',
        'read1_pos', 'read1_aend', 'cigar1', 'read2_pos', 'read2_aend', 'cigar2'
    ]
    
    print(f"Loading alignments from {args.input}...")
    aligns = pd.read_csv(args.input, sep="|", header=None, names=aligns_column_names)
    print(f"Loaded {len(aligns)} alignments.")
    
    print("Loading reference genome...")
    ref_hg19 = pysam.FastaFile(args.ref)
    
    print("Clustering UMIs and collapsing alignments...")
    #aligns_split = [group for _, group in aligns.groupby(['primer', 'alignChrom', 'alignStrand'])]

    num_threads = min(24, multiprocessing.cpu_count())
    aligns_clustered = []
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = submit_groups(executor, aligns)
        #[executor.submit(process_group, group) for group in aligns_split]
        for f in tqdm(as_completed(futures), total=len(futures)):
            aligns_clustered.append(f.result())
    aligns_umi_grouped = pd.concat(aligns_clustered, ignore_index=True, copy=False)
    
    print("Extracting genomic sequences...")
    aligns_umi_grouped['genomic_seq'] = aligns_umi_grouped.apply(lambda row: get_genomic_seq(row, ref_hg19), axis=1)
    aligns_umi_grouped.to_csv(args.output, index=False)
    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()
