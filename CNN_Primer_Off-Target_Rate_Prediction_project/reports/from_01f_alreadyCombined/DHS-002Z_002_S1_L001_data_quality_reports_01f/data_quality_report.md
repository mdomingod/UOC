# Data Quality Report — DHS-002Z_002_S1_L001 collapsed alignments

_Generated: 2026-05-31 17:58:50_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 61,407 |
| Groups with mixed labels      | 18 |
| Fraction of mixed groups      | 0.029313% |
| Rows in mixed groups          | 39 |
| Unique primers in mixed groups| 9 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...GTAACA---------` | 2 | 10 | 0 | chr10, chr9 | 33,675,560, 89,725,106 |
| 2 | `...TCCTT---------C` | 3 | 100 | 0 | chr16, chr22 | 29,085,247, 32,369,663, 33,366,854 |
| 3 | `...ATAGCC---------` | 2 | 10 | 0 | chr17, chr2 | 133,020,025, 45,219,702 |
| 4 | `...TAAC----------A` | 2 | 10 | 0 | chr10, chr9 | 33,675,559, 89,725,107 |
| 5 | `...CACCA----------` | 2 | 10 | 0 | chr10, chr11 | 121,233,180, 88,679,037 |
| 6 | `...AATGA----------` | 2 | 10 | 0 | chr7 | 6,018,302, 6,785,743 |
| 7 | `...CCTGG----------` | 2 | 1 | 0 | chr7 | 6,022,583, 6,781,354 |
| 8 | `...CAA----------TA` | 2 | 10 | 0 | chr15, chr22 | 20,489,366, 29,091,243 |
| 9 | `...CCTTC----------` | 3 | 100 | 0 | chr16, chr22 | 29,085,246, 32,369,662, 33,366,853 |
| 10 | `...GCTT----------G` | 2 | 10 | 0 | chr7 | 6,017,247, 6,786,807 |
| 11 | `...TAGC----------C` | 2 | 10 | 0 | chr17, chr2 | 133,020,026, 45,219,701 |
| 12 | `...CTTGC----------` | 2 | 1 | 0 | chr22, chr3 | 17,055,458, 178,938,877 |
| 13 | `...A-----------TGA` | 2 | 10 | 0 | chr7 | 6,018,301, 6,785,744 |
| 14 | `...CTGG-----------` | 2 | 1 | 0 | chr7 | 6,022,582, 6,781,355 |
| 15 | `...TTGC-----------` | 2 | 1 | 0 | chr22, chr3 | 17,055,457, 178,938,876 |
| 16 | `...TC-------------` | 2 | 10 | 0 | chr16, chr22 | 29,085,243, 33,366,850 |
| 17 | `...C--------------` | 3 | 100 | 0 | chr16, chr22 | 29,085,242, 32,369,658, 33,366,849 |
| 18 | `...C--------------` | 2 | 10 | 0 | chr17, chr2 | 133,020,030, 45,219,697 |
