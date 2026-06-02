# Data Quality Report — DHS-001Z_20220525_HA collapsed alignments

_Generated: 2026-05-31 17:57:28_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 10,745 |
| Groups with mixed labels      | 21 |
| Fraction of mixed groups      | 0.195440% |
| Rows in mixed groups          | 46 |
| Unique primers in mixed groups| 19 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...TTATCA---------` | 2 | 1 | 1 | chr21, chr7 | 11,038,979, 151,945,237 |
| 2 | `...TAAC----------A` | 2 | 1 | 1 | chr10, chr9 | 33,675,559, 89,725,107 |
| 3 | `...TGACA----------` | 2 | 1 | 1 | chr21, chr7 | 10,998,233, 151,904,411 |
| 4 | `...CACCA----------` | 2 | 1 | 1 | chr10, chr11 | 88,679,037, 121,233,180 |
| 5 | `...TCCAA----------` | 2 | 1 | 1 | chr1, chr7 | 148,873,514, 151,935,799 |
| 6 | `...AATGA----------` | 2 | 1 | 1 | chr7 | 6,018,302, 6,785,743 |
| 7 | `...TATCC----------` | 2 | 1 | 1 | chr17, chr7 | 16,097,871, 57,659,511 |
| 8 | `...AGACT----------` | 3 | 1 | 2 | chr1, chr2, chr7 | 91,886,498, 148,889,721, 151,919,626 |
| 9 | `...TT----------TTC` | 2 | 1 | 1 | chr19 | 9,014,741, 9,017,555 |
| 10 | `...CCTAC----------` | 2 | 1 | 1 | chr19 | 9,014,721, 9,017,535 |
| 11 | `...TCTG----------T` | 2 | 1 | 1 | chr17, chr20 | 16,068,357, 25,733,316 |
| 12 | `...CCTGG----------` | 2 | 1 | 1 | chr7 | 6,022,583, 6,781,354 |
| 13 | `...CAA----------TA` | 2 | 1 | 1 | chr15, chr22 | 20,489,366, 29,091,243 |
| 14 | `...CCTTC----------` | 3 | 1 | 2 | chr16, chr22 | 29,085,246, 32,369,662, 33,366,853 |
| 15 | `...GCTT----------G` | 2 | 1 | 1 | chr7 | 6,017,247, 6,786,807 |
| 16 | `...CCTTC----------` | 2 | 1 | 1 | chr2, chr7 | 91,888,498, 151,921,628 |
| 17 | `...GTGT----------G` | 3 | 1 | 2 | chr21, chr7, chrUn_gl000212 | 161, 11,020,915, 151,927,110 |
| 18 | `...CTTGC----------` | 2 | 1 | 1 | chr22, chr3 | 17,055,458, 178,938,877 |
| 19 | `...GTGTG----------` | 3 | 1 | 2 | chr1, chr2, chr7 | 91,888,526, 148,887,668, 151,921,656 |
| 20 | `...CCAA-----------` | 2 | 1 | 1 | chr1, chr7 | 148,873,513, 151,935,800 |
| 21 | `...CTGG-----------` | 2 | 1 | 1 | chr7 | 6,022,582, 6,781,355 |
