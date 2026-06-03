import os
import sys
import argparse
import pandas as pd
from Bio.Seq import Seq

def rev_com(seq):
    return str(Seq(seq).reverse_complement())

def load_epcr07_hits(path, hits_columns):
    epcr07_hit_files = [
        os.path.join(root, f)
        for root, _, files in os.walk(path)
        for f in files
            if f.endswith('hits.txt')
    ]
    print(f"Combining {len(epcr07_hit_files)} hit files in {path}")
    epcr07_hits = pd.concat(
        [pd.read_csv(f, sep='|', header=None, names=hits_columns) for f in epcr07_hit_files],
        ignore_index=True
    )
    return epcr07_hits

def process_hits(epcr07_hits):
    epcr07_hits['isIntendedSite'] = epcr07_hits.apply(lambda row: 1 if row['ssw_align'] == ':'*len(str(row['primer'])) else 0, axis=1)
    epcr07_hits['alignChrom'] = epcr07_hits['chrom']
    epcr07_hits['alignStrand'] = epcr07_hits.apply(lambda row: row['strand'] if row['strand'] == 1 else 0, axis=1)
    epcr07_hits['read1_anchor'] = epcr07_hits.apply(lambda row: row['loc'] + 8 if row['alignStrand'] == 1 else row['loc'] - 7, axis=1)
    epcr07_hits['genomic_seq'] = epcr07_hits.apply(lambda row: rev_com(row['hit_seq']), axis=1)
    return epcr07_hits

def main():
    parser = argparse.ArgumentParser(description="Combine epcr07 hits from all hits.txt files in a directory and output a processed CSV.")
    parser.add_argument('--input_dir', required=True, help='Directory containing epcr07 hits.txt files')
    parser.add_argument('--output_csv', required=True, help='Output CSV file')
    args = parser.parse_args()

    hits_columns = [
        "primer",
        "hit_seq",
        "chrom",
        "loc",
        "strand",
        "ssw_score",
        "num_GC",
        "num_GT",
        "num_AT",
        "query_beg",
        "query_end",
        "subj_beg",
        "subj_end",
        "cigar",
        "dS",
        "dH",
        "dG",
        "Tm",
        "ssw_primer",
        "ssw_align",
        "ssw_genome",
        "ntthal_seq1",
        "ntthal_seq2",
        "ntthal_str1",
        "ntthal_str2",
    ]

    epcr07_hits = load_epcr07_hits(args.input_dir, hits_columns)
    epcr07_hits = process_hits(epcr07_hits)
    epcr07_hits[['primer', 'isIntendedSite', 'alignChrom', 'alignStrand', 'read1_anchor', 'genomic_seq']].to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()