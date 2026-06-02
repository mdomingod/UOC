# Data Quality Report — DHS-101_20231027_DHS_003Z_PX19_101Z_103Z_PX03_LZ collapsed alignments

_Generated: 2026-05-31 18:10:39_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 41,242 |
| Groups with mixed labels      | 17 |
| Fraction of mixed groups      | 0.041220% |
| Rows in mixed groups          | 34 |
| Unique primers in mixed groups| 14 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...GGCCT---------C` | 2 | 1 | 0 | chr22, chr3 | 17,053,957, 178,937,377 |
| 2 | `...AAGGG----------` | 2 | 10 | 0 | chr22, chr3 | 17,052,943, 178,936,023 |
| 3 | `...TATCC----------` | 2 | 10 | 0 | chr2, chr9 | 132,181,788, 80,409,459 |
| 4 | `...CCCAA----------` | 2 | 10 | 0 | chr22, chr3 | 17,054,397, 178,937,820 |
| 5 | `...GCCT----------C` | 2 | 1 | 0 | chr22, chr3 | 17,053,958, 178,937,378 |
| 6 | `...GAAGC----------` | 2 | 1 | 0 | chr22, chr3 | 17,055,432, 178,938,851 |
| 7 | `...TTGAC----------` | 2 | 1 | 0 | chr22, chr3 | 17,054,370, 178,937,793 |
| 8 | `...ATTC----------A` | 2 | 10 | 0 | chr22, chr3 | 17,055,451, 178,938,870 |
| 9 | `...CAAAC----------` | 2 | 10 | 0 | chr22, chr3 | 17,055,395, 178,938,814 |
| 10 | `...GAGCT----------` | 2 | 10 | 0 | chr22, chr3 | 17,052,974, 178,936,054 |
| 11 | `...TTTCT----------` | 2 | 10 | 0 | chr22, chr3 | 17,054,373, 178,937,796 |
| 12 | `...TGGGA----------` | 2 | 10 | 0 | chr22, chr3 | 17,054,441, 178,937,864 |
| 13 | `...TCTAG----------` | 2 | 10 | 0 | chr22, chr3 | 17,053,631, 178,937,053 |
| 14 | `...TTTAGC---------` | 2 | 10 | 0 | chr22, chr3 | 17,053,055, 178,936,136 |
| 15 | `...AATTC----------` | 2 | 10 | 0 | chr22, chr3 | 17,052,949, 178,936,029 |
| 16 | `...TTC-----------T` | 2 | 10 | 0 | chr22, chr3 | 17,054,372, 178,937,795 |
| 17 | `...AA-------------` | 2 | 10 | 0 | chr22, chr3 | 17,054,394, 178,937,817 |
