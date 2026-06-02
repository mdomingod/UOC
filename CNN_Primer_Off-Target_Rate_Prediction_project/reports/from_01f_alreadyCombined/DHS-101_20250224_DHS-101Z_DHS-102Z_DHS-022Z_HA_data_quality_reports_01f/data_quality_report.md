# Data Quality Report — DHS-101_20250224_DHS-101Z_DHS-102Z_DHS-022Z_HA collapsed alignments

_Generated: 2026-05-31 18:11:27_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 52,133 |
| Groups with mixed labels      | 22 |
| Fraction of mixed groups      | 0.042200% |
| Rows in mixed groups          | 44 |
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
| 3 | `...TGAAGC---------` | 2 | 10 | 0 | chr22, chr3 | 17,055,431, 178,938,850 |
| 4 | `...ATTGAC---------` | 2 | 1 | 0 | chr22, chr3 | 17,054,369, 178,937,792 |
| 5 | `...AAGGG----------` | 2 | 1 | 0 | chr22, chr3 | 17,052,943, 178,936,023 |
| 6 | `...TATCC----------` | 2 | 1 | 0 | chr2, chr9 | 132,181,788, 80,409,459 |
| 7 | `...CCCAA----------` | 2 | 10 | 0 | chr22, chr3 | 17,054,397, 178,937,820 |
| 8 | `...GCCT----------C` | 2 | 1 | 0 | chr22, chr3 | 17,053,958, 178,937,378 |
| 9 | `...GAAGC----------` | 2 | 10 | 0 | chr22, chr3 | 17,055,432, 178,938,851 |
| 10 | `...TTGAC----------` | 2 | 1 | 0 | chr22, chr3 | 17,054,370, 178,937,793 |
| 11 | `...ATTC----------A` | 2 | 10 | 0 | chr22, chr3 | 17,055,451, 178,938,870 |
| 12 | `...CAAAC----------` | 2 | 10 | 0 | chr22, chr3 | 17,055,395, 178,938,814 |
| 13 | `...GAGCT----------` | 2 | 10 | 0 | chr22, chr3 | 17,052,974, 178,936,054 |
| 14 | `...TTTCT----------` | 2 | 10 | 0 | chr22, chr3 | 17,054,373, 178,937,796 |
| 15 | `...TGGGA----------` | 2 | 1 | 0 | chr22, chr3 | 17,054,441, 178,937,864 |
| 16 | `...TCTAG----------` | 2 | 10 | 0 | chr22, chr3 | 17,053,631, 178,937,053 |
| 17 | `...TTTAGC---------` | 2 | 1 | 0 | chr22, chr3 | 17,053,055, 178,936,136 |
| 18 | `...AATTC----------` | 2 | 10 | 0 | chr22, chr3 | 17,052,949, 178,936,029 |
| 19 | `...TTC-----------T` | 2 | 10 | 0 | chr22, chr3 | 17,054,372, 178,937,795 |
| 20 | `...TCT------------` | 2 | 10 | 0 | chr22, chr3 | 17,054,371, 178,937,794 |
| 21 | `...TAGC-----------` | 2 | 1 | 0 | chr22, chr3 | 17,053,053, 178,936,134 |
| 22 | `...AG-------------` | 2 | 10 | 0 | chr22, chr3 | 17,053,628, 178,937,050 |
