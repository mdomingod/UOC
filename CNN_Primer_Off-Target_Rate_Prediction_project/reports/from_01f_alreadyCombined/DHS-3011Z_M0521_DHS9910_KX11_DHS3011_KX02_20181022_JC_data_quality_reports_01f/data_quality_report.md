# Data Quality Report — DHS-3011Z_M0521_DHS9910_KX11_DHS3011_KX02_20181022_JC collapsed alignments

_Generated: 2026-05-31 18:24:07_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 267,317 |
| Groups with mixed labels      | 44 |
| Fraction of mixed groups      | 0.016460% |
| Rows in mixed groups          | 94 |
| Unique primers in mixed groups| 37 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...ATAAT---------A` | 2 | 1 | 0 | chr3 | 10,030,360, 10,085,475 |
| 2 | `...GAGTAG---------` | 2 | 1 | 0 | chr1 | 155,186,216, 155,207,173 |
| 3 | `...TCCAG---------C` | 2 | 10 | 0 | chr20, chr7 | 117,188,958, 29,449,391 |
| 4 | `...TGGT---------AG` | 2 | 10 | 0 | chr3, chr5 | 195,712,617, 224,517 |
| 5 | `...TAGCC----------` | 2 | 1 | 0 | chr1 | 155,184,189, 155,204,816 |
| 6 | `...GTGTG----------` | 2 | 1 | 0 | chr2 | 152,437,979, 152,459,083 |
| 7 | `...TAATA----------` | 2 | 1 | 0 | chr3 | 10,030,361, 10,085,476 |
| 8 | `...AGTA----------G` | 2 | 1 | 0 | chr1 | 155,186,217, 155,207,174 |
| 9 | `...ATCAT----------` | 2 | 1 | 0 | chr6 | 31,975,083, 32,007,819 |
| 10 | `...TACAAG---------` | 2 | 1 | 0 | chr3 | 10,091,212, 11,912,925 |
| 11 | `...TAAC----------A` | 2 | 10 | 0 | chr10, chr9 | 33,675,559, 89,725,107 |
| 12 | `...TGCCA----------` | 2 | 1 | 0 | chr11, chr7 | 111,965,607, 135,129,593 |
| 13 | `...CAT----------TC` | 2 | 1 | 0 | chr1 | 155,184,472, 155,205,099 |
| 14 | `...TTCCT----------` | 2 | 1 | 0 | chr6 | 31,975,912, 32,008,647 |
| 15 | `...TCGTG----------` | 2 | 1 | 0 | chr6 | 31,976,003, 32,008,738 |
| 16 | `...GTTT----------C` | 2 | 1 | 0 | chr6 | 31,974,381, 32,007,117 |
| 17 | `...AAATC----------` | 2 | 1 | 0 | chr5 | 1,576,712, 251,164 |
| 18 | `...TCTTC----------` | 2 | 1 | 0 | chr6 | 31,974,787, 32,007,523 |
| 19 | `...GTTCC----------` | 2 | 1 | 0 | chr6 | 31,973,531, 32,006,265 |
| 20 | `...CAATC----------` | 2 | 1 | 0 | chr1 | 155,184,211, 155,204,838 |
| 21 | `...CCCAT----------` | 3 | 1 | 0 | chr2 | 152,440,192, 152,450,744, 152,461,297 |
| 22 | `...CACCA----------` | 2 | 10 | 0 | chr10, chr11 | 121,233,180, 88,679,037 |
| 23 | `...AAAGC----------` | 3 | 1 | 0 | chr1, chr11 | 111,965,751, 25,620,597, 25,725,748 |
| 24 | `...AATGA----------` | 2 | 10 | 0 | chr7 | 6,018,302, 6,785,743 |
| 25 | `...TCTCG----------` | 3 | 10 | 0 | chr2, chr22, chrX | 153,008,652, 16,871,073, 92,031,294 |
| 26 | `...CCTGG----------` | 2 | 1 | 0 | chr7 | 6,022,583, 6,781,354 |
| 27 | `...CAA----------TA` | 2 | 10 | 0 | chr15, chr22 | 20,489,366, 29,091,243 |
| 28 | `...GCCAA----------` | 2 | 1 | 0 | chr6 | 31,973,621, 32,006,355 |
| 29 | `...CTG----------TG` | 2 | 10 | 0 | chr2, chrX | 153,008,528, 92,031,170 |
| 30 | `...CCTTC----------` | 3 | 100 | 0 | chr16, chr22 | 29,085,246, 32,369,662, 33,366,853 |
| 31 | `...CCAGC----------` | 2 | 10 | 0 | chr20, chr7 | 117,188,957, 29,449,392 |
| 32 | `...GCTT----------G` | 2 | 10 | 0 | chr7 | 6,017,247, 6,786,807 |
| 33 | `...CATCT----------` | 2 | 1 | 0 | chr1 | 169,510,367, 169,510,637 |
| 34 | `...GCTGG----------` | 2 | 1 | 0 | chr6 | 31,976,200, 32,008,935 |
| 35 | `...TCCTC----------` | 2 | 1 | 0 | chr6 | 31,973,724, 32,006,458 |
| 36 | `...CAGG----------A` | 2 | 10 | 0 | chr20, chr7 | 117,188,805, 25,900,206 |
| 37 | `...GGATG----------` | 2 | 1 | 0 | chr1 | 155,187,107, 155,208,384 |
| 38 | `...GGTA----------G` | 4 | 100 | 0 | chr3, chr5 | 195,389,485, 195,712,616, 197,350,209, 224,518 |
| 39 | `...CTTGC----------` | 2 | 10 | 0 | chr22, chr3 | 17,055,458, 178,938,877 |
| 40 | `...CGATG----------` | 2 | 10 | 0 | chr5 | 1,593,242, 236,699 |
| 41 | `...GGTGT----------` | 2 | 1 | 0 | chr5 | 1,574,430, 254,487 |
| 42 | `...TTGC-----------` | 2 | 10 | 0 | chr22, chr3 | 17,055,457, 178,938,876 |
| 43 | `...GTGT-----------` | 2 | 1 | 0 | chr5 | 1,574,429, 254,488 |
| 44 | `...TA------------G` | 2 | 1 | 0 | chr1 | 155,186,219, 155,207,176 |
