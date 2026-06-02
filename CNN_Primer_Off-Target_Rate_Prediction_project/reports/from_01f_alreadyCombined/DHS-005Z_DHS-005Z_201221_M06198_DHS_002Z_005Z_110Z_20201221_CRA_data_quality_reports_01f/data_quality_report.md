# Data Quality Report — DHS-005Z_DHS-005Z_201221_M06198_DHS_002Z_005Z_110Z_20201221_CRA collapsed alignments

_Generated: 2026-05-31 18:02:20_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 98,401 |
| Groups with mixed labels      | 11 |
| Fraction of mixed groups      | 0.011179% |
| Rows in mixed groups          | 22 |
| Unique primers in mixed groups| 5 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...GTAACA---------` | 2 | 10 | 0 | chr10, chr9 | 33,675,560, 89,725,106 |
| 2 | `...CACA---------GC` | 2 | 10 | 0 | chr7 | 126,544,053, 65,082,687 |
| 3 | `...TAAC----------A` | 2 | 10 | 0 | chr10, chr9 | 33,675,559, 89,725,107 |
| 4 | `...TT----------TTC` | 2 | 1 | 0 | chr19 | 9,014,741, 9,017,555 |
| 5 | `...CCTAC----------` | 2 | 1 | 0 | chr19 | 9,014,721, 9,017,535 |
| 6 | `...ACAGC----------` | 2 | 10 | 0 | chr7 | 126,544,054, 65,082,686 |
| 7 | `...CTTGC----------` | 2 | 1 | 0 | chr22, chr3 | 17,055,458, 178,938,877 |
| 8 | `...CTAC-----------` | 2 | 1 | 0 | chr19 | 9,014,720, 9,017,534 |
| 9 | `...TTGC-----------` | 2 | 1 | 0 | chr22, chr3 | 17,055,457, 178,938,876 |
| 10 | `...A-------------C` | 2 | 1 | 0 | chr19 | 9,014,718, 9,017,532 |
| 11 | `...GC-------------` | 2 | 1 | 0 | chr22, chr3 | 17,055,455, 178,938,874 |
