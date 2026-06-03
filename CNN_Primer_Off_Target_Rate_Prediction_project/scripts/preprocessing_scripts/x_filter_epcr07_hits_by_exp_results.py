import pandas as pd
import argparse

def epcr_row_has_match_in_exp(row, exp_results):
    exp_match_by_pos = exp_results[(exp_results['alignChrom'] == row.alignChrom) &
                                   (exp_results['alignStrand'] == row.alignStrand) &
                                   (exp_results['read1_anchor'] > row.read1_anchor - 10) &
                                   (exp_results['read1_anchor'] < row.read1_anchor + 10)]
    return len(exp_match_by_pos) > 0

def main(exp_results_csv, epcr07_hits_csv, output_csv):
    exp_results = pd.read_csv(exp_results_csv, index_col=False)
    epcr07_hits = pd.read_csv(epcr07_hits_csv, index_col=False)

    # Filter out intended sites if column exists
    if 'isIntendedSite' in epcr07_hits.columns:
        epcr07_hits_filtered = epcr07_hits[epcr07_hits['isIntendedSite'] == 0]
    else:
        epcr07_hits_filtered = epcr07_hits.copy()

    # Filter out only primers that exist in exp_results
    exp_primers = set(exp_results['primer'])
    epcr07_hits_filtered = epcr07_hits_filtered[epcr07_hits_filtered['primer'].isin(exp_primers)]


    # Keep only rows without a matching experimental hit by position
    epcr07_hits_filtered = epcr07_hits_filtered[~epcr07_hits_filtered.apply(lambda row: epcr_row_has_match_in_exp(row, exp_results), axis=1)]

    # Add umi_groups_cnt = 0 column as these are alignments not found in exp_results
    epcr07_hits_filtered['umi_groups_cnt'] = 0

    # Output only columns present in exp_results
    epcr07_hits_filtered[exp_results.columns].to_csv(output_csv, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter epcr07 hits by experimental results and output matching hits in exp_results format.")
    parser.add_argument("--exp_results", required=True, help="Experimental results CSV file")
    parser.add_argument("--epcr07_hits", required=True, help="epcr07 combined hits CSV file")
    parser.add_argument("--output", required=True, help="Output CSV file")
    args = parser.parse_args()
    main(args.exp_results, args.epcr07_hits, args.output)
