# Data Quality Report — DHS-101_DHS_003Z_101Z_LNF06_20210616_JC collapsed alignments

_Generated: 2026-05-31 18:12:40_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 80,476 |
| Groups with mixed labels      | 23 |
| Fraction of mixed groups      | 0.028580% |
| Rows in mixed groups          | 46 |
| Unique primers in mixed groups| 14 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...GTATC---------C` | 2 | 1 | 0 | chr2, chr9 | 132,181,787, 80,409,458 |
| 2 | `...GGCCT---------C` | 2 | 1 | 0 | chr22, chr3 | 17,053,957, 178,937,377 |
| 3 | `...TGAAGC---------` | 2 | 1 | 0 | chr22, chr3 | 17,055,431, 178,938,850 |
| 4 | `...CTTTCT---------` | 2 | 1 | 0 | chr22, chr3 | 17,054,374, 178,937,797 |
| 5 | `...TAATTC---------` | 2 | 10 | 0 | chr22, chr3 | 17,052,950, 178,936,030 |
| 6 | `...AAGGG----------` | 2 | 10 | 0 | chr22, chr3 | 17,052,943, 178,936,023 |
| 7 | `...TATCC----------` | 2 | 1 | 0 | chr2, chr9 | 132,181,788, 80,409,459 |
| 8 | `...CCCAA----------` | 2 | 1 | 0 | chr22, chr3 | 17,054,397, 178,937,820 |
| 9 | `...GCCT----------C` | 2 | 1 | 0 | chr22, chr3 | 17,053,958, 178,937,378 |
| 10 | `...GAAGC----------` | 2 | 1 | 0 | chr22, chr3 | 17,055,432, 178,938,851 |
| 11 | `...TTGAC----------` | 2 | 1 | 0 | chr22, chr3 | 17,054,370, 178,937,793 |
| 12 | `...ATTC----------A` | 2 | 1 | 0 | chr22, chr3 | 17,055,451, 178,938,870 |
| 13 | `...CAAAC----------` | 2 | 1 | 0 | chr22, chr3 | 17,055,395, 178,938,814 |
| 14 | `...GAGCT----------` | 2 | 10 | 0 | chr22, chr3 | 17,052,974, 178,936,054 |
| 15 | `...TTTCT----------` | 2 | 1 | 0 | chr22, chr3 | 17,054,373, 178,937,796 |
| 16 | `...TGGGA----------` | 2 | 1 | 0 | chr22, chr3 | 17,054,441, 178,937,864 |
| 17 | `...TCTAG----------` | 2 | 10 | 0 | chr22, chr3 | 17,053,631, 178,937,053 |
| 18 | `...TTTAGC---------` | 2 | 1 | 0 | chr22, chr3 | 17,053,055, 178,936,136 |
| 19 | `...AATTC----------` | 2 | 10 | 0 | chr22, chr3 | 17,052,949, 178,936,029 |
| 20 | `...TTC-----------T` | 2 | 1 | 0 | chr22, chr3 | 17,054,372, 178,937,795 |
| 21 | `...TCT------------` | 2 | 1 | 0 | chr22, chr3 | 17,054,371, 178,937,794 |
| 22 | `...TAG------------` | 2 | 10 | 0 | chr22, chr3 | 17,053,629, 178,937,051 |
| 23 | `...TAGC-----------` | 2 | 1 | 0 | chr22, chr3 | 17,053,053, 178,936,134 |
