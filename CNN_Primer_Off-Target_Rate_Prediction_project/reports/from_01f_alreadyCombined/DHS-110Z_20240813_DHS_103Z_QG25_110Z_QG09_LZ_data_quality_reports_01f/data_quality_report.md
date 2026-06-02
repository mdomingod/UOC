# Data Quality Report — DHS-110Z_20240813_DHS_103Z_QG25_110Z_QG09_LZ collapsed alignments

_Generated: 2026-05-31 18:15:51_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 118,069 |
| Groups with mixed labels      | 19 |
| Fraction of mixed groups      | 0.016092% |
| Rows in mixed groups          | 44 |
| Unique primers in mixed groups| 11 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...CTTCTG---------` | 2 | 10 | 0 | chr15, chr22 | 20,489,486, 29,091,123 |
| 2 | `...GAGTGC---------` | 2 | 10 | 0 | chr15, chr22 | 20,488,771, 29,091,838 |
| 3 | `...ATTTCT---------` | 3 | 100 | 0 | chr15, chr16, chr22 | 20,489,342, 29,091,267, 33,372,834 |
| 4 | `...CATA----------A` | 2 | 1 | 0 | chr10, chr22 | 29,083,784, 39,012,459 |
| 5 | `...AGTGC----------` | 2 | 10 | 0 | chr15, chr22 | 20,488,770, 29,091,839 |
| 6 | `...ATA----------CC` | 2 | 10 | 0 | chr15, chr22 | 20,488,868, 29,091,741 |
| 7 | `...----------TTTCT` | 3 | 100 | 0 | chr15, chr16, chr22 | 20,489,343, 29,091,266, 33,372,833 |
| 8 | `...GATCA----------` | 2 | 10 | 0 | chr16, chr22 | 29,091,181, 33,372,748 |
| 9 | `...CCTTC----------` | 3 | 100 | 0 | chr16, chr22 | 29,085,246, 32,369,662, 33,366,853 |
| 10 | `...AGTT----------G` | 3 | 100 | 0 | chr16, chr22 | 29,090,090, 32,375,001, 33,371,659 |
| 11 | `...CTCCA----------` | 3 | 100 | 0 | chr15, chr16, chr22 | 20,488,762, 29,091,847, 33,373,388 |
| 12 | `...AGTGG----------` | 2 | 10 | 0 | chr15, chr22 | 20,489,380, 29,091,229 |
| 13 | `...GTGC-----------` | 2 | 10 | 0 | chr15, chr22 | 20,488,769, 29,091,840 |
| 14 | `...TACC-----------` | 2 | 10 | 0 | chr15, chr22 | 20,488,867, 29,091,742 |
| 15 | `...TTCT-----------` | 3 | 100 | 0 | chr15, chr16, chr22 | 20,489,344, 29,091,265, 33,372,832 |
| 16 | `...TCC-----------A` | 2 | 10 | 0 | chr16, chr22 | 29,091,846, 33,373,387 |
| 17 | `...GTGG-----------` | 2 | 10 | 0 | chr15, chr22 | 20,489,381, 29,091,228 |
| 18 | `...C------------CA` | 2 | 10 | 0 | chr11, chr7 | 108,159,786, 53,854,941 |
| 19 | `...CCA------------` | 2 | 10 | 0 | chr16, chr22 | 29,091,845, 33,373,386 |
