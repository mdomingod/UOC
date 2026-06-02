# Data Quality Report — DHS-001Z_20220607_HA_S1_R collapsed alignments

_Generated: 2026-05-31 17:58:03_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 47,985 |
| Groups with mixed labels      | 36 |
| Fraction of mixed groups      | 0.075023% |
| Rows in mixed groups          | 81 |
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
| 2 | `...TTGACA---------` | 2 | 1 | 1 | chr21, chr7 | 10,998,232, 151,904,410 |
| 3 | `...ACACCA---------` | 2 | 1 | 1 | chr10, chr11 | 88,679,036, 121,233,181 |
| 4 | `...TAATGA---------` | 2 | 1 | 1 | chr7 | 6,018,303, 6,785,742 |
| 5 | `...CCAATA---------` | 2 | 1 | 1 | chr15, chr22 | 20,489,365, 29,091,244 |
| 6 | `...TCCTT---------C` | 3 | 1 | 2 | chr16, chr22 | 29,085,247, 32,369,663, 33,366,854 |
| 7 | `...AAATCC---------` | 2 | 1 | 1 | chr21, chr7 | 11,058,210, 151,970,839 |
| 8 | `...TTATCA---------` | 2 | 1 | 1 | chr21, chr7 | 11,038,979, 151,945,237 |
| 9 | `...TAAC----------A` | 2 | 1 | 1 | chr10, chr9 | 33,675,559, 89,725,107 |
| 10 | `...TGACA----------` | 2 | 1 | 1 | chr21, chr7 | 10,998,233, 151,904,411 |
| 11 | `...CACCA----------` | 2 | 1 | 1 | chr10, chr11 | 88,679,037, 121,233,180 |
| 12 | `...TCCAA----------` | 2 | 1 | 1 | chr1, chr7 | 148,873,514, 151,935,799 |
| 13 | `...AATGA----------` | 2 | 1 | 1 | chr7 | 6,018,302, 6,785,743 |
| 14 | `...TATCC----------` | 2 | 1 | 1 | chr17, chr7 | 16,097,871, 57,659,511 |
| 15 | `...AGACT----------` | 3 | 1 | 2 | chr1, chr2, chr7 | 91,886,498, 148,889,721, 151,919,626 |
| 16 | `...TT----------TTC` | 2 | 1 | 1 | chr19 | 9,014,741, 9,017,555 |
| 17 | `...CCTAC----------` | 2 | 1 | 1 | chr19 | 9,014,721, 9,017,535 |
| 18 | `...TCTG----------T` | 2 | 1 | 1 | chr17, chr20 | 16,068,357, 25,733,316 |
| 19 | `...CCTGG----------` | 2 | 1 | 1 | chr7 | 6,022,583, 6,781,354 |
| 20 | `...CAA----------TA` | 2 | 1 | 1 | chr15, chr22 | 20,489,366, 29,091,243 |
| 21 | `...CCTTC----------` | 3 | 1 | 2 | chr16, chr22 | 29,085,246, 32,369,662, 33,366,853 |
| 22 | `...GCTT----------G` | 2 | 1 | 1 | chr7 | 6,017,247, 6,786,807 |
| 23 | `...CCTTC----------` | 4 | 1 | 3 | chr2, chr21, chr7, chrUn_gl000212 | 173,314, 11,015,471, 91,888,498, 151,921,628 |
| 24 | `...GTGT----------G` | 3 | 1 | 2 | chr21, chr7, chrUn_gl000212 | 161, 11,020,915, 151,927,110 |
| 25 | `...CTTGC----------` | 2 | 1 | 1 | chr22, chr3 | 17,055,458, 178,938,877 |
| 26 | `...GTGTG----------` | 3 | 1 | 2 | chr1, chr2, chr7 | 91,888,526, 148,887,668, 151,921,656 |
| 27 | `...CCAA-----------` | 2 | 1 | 1 | chr1, chr7 | 148,873,513, 151,935,800 |
| 28 | `...CTAC-----------` | 2 | 1 | 1 | chr19 | 9,014,720, 9,017,534 |
| 29 | `...TTGC-----------` | 2 | 1 | 1 | chr22, chr3 | 17,055,457, 178,938,876 |
| 30 | `...ATCA-----------` | 2 | 1 | 1 | chr21, chr7 | 11,038,977, 151,945,235 |
| 31 | `...TGA------------` | 2 | 1 | 1 | chr7 | 6,018,300, 6,785,745 |
| 32 | `...ATA------------` | 2 | 1 | 1 | chr15, chr22 | 20,489,368, 29,091,241 |
| 33 | `...GTG------------` | 2 | 1 | 1 | chr21, chr7 | 11,020,913, 151,927,108 |
| 34 | `...GT------------G` | 2 | 1 | 1 | chr2, chr7 | 91,888,524, 151,921,654 |
| 35 | `...TC-------------` | 3 | 1 | 2 | chr16, chr22 | 29,085,243, 32,369,659, 33,366,850 |
| 36 | `...G--------------` | 3 | 1 | 2 | chr1, chr2, chr7 | 91,888,522, 148,887,672, 151,921,652 |
