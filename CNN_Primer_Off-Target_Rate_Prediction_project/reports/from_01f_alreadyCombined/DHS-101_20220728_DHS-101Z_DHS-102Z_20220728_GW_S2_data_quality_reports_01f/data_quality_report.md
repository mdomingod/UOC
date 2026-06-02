# Data Quality Report — DHS-101_20220728_DHS-101Z_DHS-102Z_20220728_GW_S2 collapsed alignments

_Generated: 2026-05-31 18:09:56_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 60,318 |
| Groups with mixed labels      | 18 |
| Fraction of mixed groups      | 0.029842% |
| Rows in mixed groups          | 36 |
| Unique primers in mixed groups| 14 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...TAAGGG---------` | 2 | 1 | 0 | chr22, chr3 | 17,052,942, 178,936,022 |
| 2 | `...GGCCT---------C` | 2 | 1 | 0 | chr22, chr3 | 17,053,957, 178,937,377 |
| 3 | `...AAGGG----------` | 2 | 1 | 0 | chr22, chr3 | 17,052,943, 178,936,023 |
| 4 | `...TATCC----------` | 2 | 1 | 0 | chr2, chr9 | 132,181,788, 80,409,459 |
| 5 | `...CCCAA----------` | 2 | 1 | 0 | chr22, chr3 | 17,054,397, 178,937,820 |
| 6 | `...GCCT----------C` | 2 | 1 | 0 | chr22, chr3 | 17,053,958, 178,937,378 |
| 7 | `...GAAGC----------` | 2 | 1 | 0 | chr22, chr3 | 17,055,432, 178,938,851 |
| 8 | `...TTGAC----------` | 2 | 10 | 0 | chr22, chr3 | 17,054,370, 178,937,793 |
| 9 | `...ATTC----------A` | 2 | 10 | 0 | chr22, chr3 | 17,055,451, 178,938,870 |
| 10 | `...CAAAC----------` | 2 | 1 | 0 | chr22, chr3 | 17,055,395, 178,938,814 |
| 11 | `...GAGCT----------` | 2 | 10 | 0 | chr22, chr3 | 17,052,974, 178,936,054 |
| 12 | `...TTTCT----------` | 2 | 1 | 0 | chr22, chr3 | 17,054,373, 178,937,796 |
| 13 | `...TGGGA----------` | 2 | 10 | 0 | chr22, chr3 | 17,054,441, 178,937,864 |
| 14 | `...TCTAG----------` | 2 | 1 | 0 | chr22, chr3 | 17,053,631, 178,937,053 |
| 15 | `...TTTAGC---------` | 2 | 10 | 0 | chr22, chr3 | 17,053,055, 178,936,136 |
| 16 | `...AATTC----------` | 2 | 10 | 0 | chr22, chr3 | 17,052,949, 178,936,029 |
| 17 | `...TTC-----------T` | 2 | 1 | 0 | chr22, chr3 | 17,054,372, 178,937,795 |
| 18 | `...TAGC-----------` | 2 | 10 | 0 | chr22, chr3 | 17,053,053, 178,936,134 |
