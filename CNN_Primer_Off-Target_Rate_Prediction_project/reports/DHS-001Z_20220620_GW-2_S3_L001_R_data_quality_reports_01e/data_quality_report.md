# Data Quality Report — DHS-001Z_20220620_GW-2_S3_L001_R collapsed alignments

_Generated: 2026-05-31 17:58:36_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 36,274 |
| Groups with mixed labels      | 28 |
| Fraction of mixed groups      | 0.077190% |
| Rows in mixed groups          | 63 |
| Unique primers in mixed groups| 20 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...GTTAT---------C` | 2 | 1 | 1 | chr21, chr7 | 11,038,980, 151,945,238 |
| 2 | `...TAATGA---------` | 2 | 1 | 1 | chr7 | 6,018,303, 6,785,742 |
| 3 | `...AAATCC---------` | 2 | 1 | 1 | chr21, chr7 | 11,058,210, 151,970,839 |
| 4 | `...TTATCA---------` | 2 | 1 | 1 | chr21, chr7 | 11,038,979, 151,945,237 |
| 5 | `...TAAC----------A` | 2 | 1 | 1 | chr10, chr9 | 33,675,559, 89,725,107 |
| 6 | `...TGACA----------` | 2 | 1 | 1 | chr21, chr7 | 10,998,233, 151,904,411 |
| 7 | `...CACCA----------` | 2 | 1 | 1 | chr10, chr11 | 88,679,037, 121,233,180 |
| 8 | `...TCCAA----------` | 2 | 1 | 1 | chr1, chr7 | 148,873,514, 151,935,799 |
| 9 | `...AATGA----------` | 2 | 1 | 1 | chr7 | 6,018,302, 6,785,743 |
| 10 | `...TATCC----------` | 2 | 1 | 1 | chr17, chr7 | 16,097,871, 57,659,511 |
| 11 | `...AGACT----------` | 3 | 1 | 2 | chr1, chr2, chr7 | 91,886,498, 148,889,721, 151,919,626 |
| 12 | `...TT----------TTC` | 2 | 1 | 1 | chr19 | 9,014,741, 9,017,555 |
| 13 | `...CCTAC----------` | 2 | 1 | 1 | chr19 | 9,014,721, 9,017,535 |
| 14 | `...TCTG----------T` | 2 | 1 | 1 | chr17, chr20 | 16,068,357, 25,733,316 |
| 15 | `...CCTGG----------` | 2 | 1 | 1 | chr7 | 6,022,583, 6,781,354 |
| 16 | `...CAA----------TA` | 2 | 1 | 1 | chr15, chr22 | 20,489,366, 29,091,243 |
| 17 | `...CCTTC----------` | 3 | 1 | 2 | chr16, chr22 | 29,085,246, 32,369,662, 33,366,853 |
| 18 | `...GCTT----------G` | 2 | 1 | 1 | chr7 | 6,017,247, 6,786,807 |
| 19 | `...CCTTC----------` | 4 | 1 | 3 | chr2, chr21, chr7, chrUn_gl000212 | 173,314, 11,015,471, 91,888,498, 151,921,628 |
| 20 | `...GTGT----------G` | 3 | 1 | 2 | chr21, chr7, chrUn_gl000212 | 161, 11,020,915, 151,927,110 |
| 21 | `...CTTGC----------` | 2 | 1 | 1 | chr22, chr3 | 17,055,458, 178,938,877 |
| 22 | `...GTGTG----------` | 3 | 1 | 2 | chr1, chr2, chr7 | 91,888,526, 148,887,668, 151,921,656 |
| 23 | `...CCAA-----------` | 2 | 1 | 1 | chr1, chr7 | 148,873,513, 151,935,800 |
| 24 | `...TTGC-----------` | 2 | 1 | 1 | chr22, chr3 | 17,055,457, 178,938,876 |
| 25 | `...TGA------------` | 2 | 1 | 1 | chr7 | 6,018,300, 6,785,745 |
| 26 | `...ATA------------` | 2 | 1 | 1 | chr15, chr22 | 20,489,368, 29,091,241 |
| 27 | `...------------TGC` | 2 | 1 | 1 | chr22, chr3 | 17,055,456, 178,938,875 |
| 28 | `...G--------------` | 3 | 1 | 2 | chr1, chr2, chr7 | 91,888,522, 148,887,672, 151,921,652 |
