# Data Quality Report — DHS-3011Z_M0512_TEPCR_3011Z_20181015_AZ-100227131 collapsed alignments

_Generated: 2026-05-31 18:19:48_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 261,524 |
| Groups with mixed labels      | 44 |
| Fraction of mixed groups      | 0.016824% |
| Rows in mixed groups          | 94 |
| Unique primers in mixed groups| 36 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...ACACCA---------` | 2 | 10 | 0 | chr10, chr11 | 121,233,181, 88,679,036 |
| 2 | `...CGGATG---------` | 2 | 1 | 0 | chr1 | 155,187,106, 155,208,383 |
| 3 | `...TGGT---------AG` | 2 | 1 | 0 | chr3, chr5 | 195,712,617, 224,517 |
| 4 | `...TAGCC----------` | 2 | 1 | 0 | chr1 | 155,184,189, 155,204,816 |
| 5 | `...GTGTG----------` | 2 | 1 | 0 | chr2 | 152,437,979, 152,459,083 |
| 6 | `...TAATA----------` | 2 | 1 | 0 | chr3 | 10,030,361, 10,085,476 |
| 7 | `...AGTA----------G` | 2 | 1 | 0 | chr1 | 155,186,217, 155,207,174 |
| 8 | `...ATCAT----------` | 2 | 1 | 0 | chr6 | 31,975,083, 32,007,819 |
| 9 | `...TACAAG---------` | 2 | 1 | 0 | chr3 | 10,091,212, 11,912,925 |
| 10 | `...TAAC----------A` | 2 | 10 | 0 | chr10, chr9 | 33,675,559, 89,725,107 |
| 11 | `...TGCCA----------` | 2 | 10 | 0 | chr11, chr7 | 111,965,607, 135,129,593 |
| 12 | `...CAT----------TC` | 2 | 1 | 0 | chr1 | 155,184,472, 155,205,099 |
| 13 | `...TCGTG----------` | 2 | 1 | 0 | chr6 | 31,976,003, 32,008,738 |
| 14 | `...GTTT----------C` | 2 | 1 | 0 | chr6 | 31,974,381, 32,007,117 |
| 15 | `...AAATC----------` | 2 | 1 | 0 | chr5 | 1,576,712, 251,164 |
| 16 | `...TCTTC----------` | 2 | 1 | 0 | chr6 | 31,974,787, 32,007,523 |
| 17 | `...GTTCC----------` | 2 | 1 | 0 | chr6 | 31,973,531, 32,006,265 |
| 18 | `...CAATC----------` | 2 | 1 | 0 | chr1 | 155,184,211, 155,204,838 |
| 19 | `...CCCAT----------` | 3 | 1 | 0 | chr2 | 152,440,192, 152,450,744, 152,461,297 |
| 20 | `...CACCA----------` | 2 | 10 | 0 | chr10, chr11 | 121,233,180, 88,679,037 |
| 21 | `...AAAGC----------` | 3 | 1 | 0 | chr1, chr11 | 111,965,751, 25,620,597, 25,725,748 |
| 22 | `...AATGA----------` | 2 | 10 | 0 | chr7 | 6,018,302, 6,785,743 |
| 23 | `...TCTCG----------` | 3 | 10 | 0 | chr2, chr22, chrX | 153,008,652, 16,871,073, 92,031,294 |
| 24 | `...CCTGG----------` | 2 | 10 | 0 | chr7 | 6,022,583, 6,781,354 |
| 25 | `...CAA----------TA` | 2 | 10 | 0 | chr15, chr22 | 20,489,366, 29,091,243 |
| 26 | `...GCCAA----------` | 2 | 1 | 0 | chr6 | 31,973,621, 32,006,355 |
| 27 | `...CTG----------TG` | 2 | 1 | 0 | chr2, chrX | 153,008,528, 92,031,170 |
| 28 | `...CCTTC----------` | 3 | 100 | 0 | chr16, chr22 | 29,085,246, 32,369,662, 33,366,853 |
| 29 | `...CCAGC----------` | 2 | 10 | 0 | chr20, chr7 | 117,188,957, 29,449,392 |
| 30 | `...GCTT----------G` | 2 | 1 | 0 | chr7 | 6,017,247, 6,786,807 |
| 31 | `...CATCT----------` | 2 | 1 | 0 | chr1 | 169,510,367, 169,510,637 |
| 32 | `...GCTGG----------` | 2 | 1 | 0 | chr6 | 31,976,200, 32,008,935 |
| 33 | `...TCCTC----------` | 2 | 1 | 0 | chr6 | 31,973,724, 32,006,458 |
| 34 | `...CAGG----------A` | 2 | 10 | 0 | chr20, chr7 | 117,188,805, 25,900,206 |
| 35 | `...GGATG----------` | 2 | 1 | 0 | chr1 | 155,187,107, 155,208,384 |
| 36 | `...GGTA----------G` | 4 | 1 | 0 | chr3, chr5 | 195,389,485, 195,712,616, 197,350,209, 224,518 |
| 37 | `...CTTGC----------` | 2 | 1 | 0 | chr22, chr3 | 17,055,458, 178,938,877 |
| 38 | `...CGATG----------` | 2 | 1 | 0 | chr5 | 1,593,242, 236,699 |
| 39 | `...GGTGT----------` | 2 | 1 | 0 | chr5 | 1,574,430, 254,487 |
| 40 | `...GCCA-----------` | 2 | 10 | 0 | chr11, chr7 | 111,965,606, 135,129,592 |
| 41 | `...ATT-----------C` | 2 | 1 | 0 | chr1 | 155,184,471, 155,205,098 |
| 42 | `...A-----------TGA` | 2 | 10 | 0 | chr7 | 6,018,301, 6,785,744 |
| 43 | `...TTGC-----------` | 2 | 1 | 0 | chr22, chr3 | 17,055,457, 178,938,876 |
| 44 | `...GTGT-----------` | 2 | 1 | 0 | chr5 | 1,574,429, 254,488 |
