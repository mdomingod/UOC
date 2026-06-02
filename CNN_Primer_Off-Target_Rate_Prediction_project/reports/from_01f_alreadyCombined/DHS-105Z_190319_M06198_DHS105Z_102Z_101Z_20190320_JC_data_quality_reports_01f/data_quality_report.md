# Data Quality Report — DHS-105Z_190319_M06198_DHS105Z_102Z_101Z_20190320_JC collapsed alignments

_Generated: 2026-05-31 18:13:55_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 8,142 |
| Groups with mixed labels      | 36 |
| Fraction of mixed groups      | 0.442152% |
| Rows in mixed groups          | 75 |
| Unique primers in mixed groups| 34 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...GCCAA---------C` | 2 | 1 | 0 | chr1, chrM | 568,033, 7,483 |
| 2 | `...TTCTC----------` | 2 | 1 | 0 | chr1, chrM | 4,347, 564,896 |
| 3 | `...ATCGG----------` | 3 | 1 | 0 | chr11, chr5, chrM | 10,531,173, 1,234, 79,947,304 |
| 4 | `...TTCGA----------` | 2 | 1 | 0 | chr1, chrM | 567,893, 7,343 |
| 5 | `...AAATG----------` | 2 | 1 | 0 | chr1, chrM | 567,261, 6,711 |
| 6 | `...AATGG----------` | 2 | 1 | 0 | chr1, chrM | 5,661, 566,210 |
| 7 | `...TCA----------CG` | 2 | 1 | 0 | chr1, chrM | 4,648, 565,197 |
| 8 | `...TTAG----------C` | 2 | 1 | 0 | chrM, chrX | 125,606,618, 786 |
| 9 | `...CAAC----------A` | 2 | 1 | 0 | chrM, chrX | 2,965, 55,208,097 |
| 10 | `...CTATT----------` | 2 | 1 | 0 | chr1, chrM | 567,645, 7,095 |
| 11 | `...TCAC----------A` | 2 | 1 | 0 | chr1, chrM | 5,275, 565,824 |
| 12 | `...CCAAC----------` | 2 | 1 | 0 | chr1, chrM | 568,034, 7,484 |
| 13 | `...ATA----------TG` | 2 | 1 | 0 | chr1, chrM | 567,442, 6,892 |
| 14 | `...TCAGG----------` | 2 | 1 | 0 | chr1, chrM | 570,196, 9,648 |
| 15 | `...GA----------GGG` | 2 | 1 | 0 | chr1, chrM | 568,839, 8,291 |
| 16 | `...ATAG----------G` | 2 | 1 | 0 | chr1, chrM | 4,390, 564,939 |
| 17 | `...CTTAG----------` | 2 | 1 | 0 | chr1, chrM | 566,256, 5,707 |
| 18 | `...TAAGC----------` | 2 | 1 | 0 | chr1, chrM | 566,246, 5,697 |
| 19 | `...AGAT----------G` | 2 | 1 | 0 | chr1, chrM | 4,501, 565,050 |
| 20 | `...ATATG----------` | 2 | 1 | 0 | chr1, chrM | 567,369, 6,819 |
| 21 | `...CAATG----------` | 2 | 1 | 0 | chr5, chrM | 14,698, 93,905,195 |
| 22 | `...AAAA----------G` | 2 | 1 | 0 | chr1, chrM | 4,275, 564,824 |
| 23 | `...TAAGG----------` | 2 | 1 | 0 | chr1, chrM | 569,402, 8,854 |
| 24 | `...TAGGC----------` | 2 | 1 | 0 | chr1, chrM | 566,535, 5,986 |
| 25 | `...TAAAA----------` | 2 | 1 | 0 | chr1, chrM | 4,558, 565,107 |
| 26 | `...TATGG----------` | 2 | 1 | 0 | chr1, chrM | 567,927, 7,377 |
| 27 | `...TACA----------C` | 3 | 1 | 0 | chr11, chr5, chrM | 10,530,932, 1,475, 79,947,063 |
| 28 | `...CATG----------T` | 2 | 1 | 0 | chr1, chrM | 5,189, 565,738 |
| 29 | `...AATCG----------` | 2 | 1 | 0 | chr1, chrM | 4,149, 564,698 |
| 30 | `...TAATC----------` | 2 | 1 | 0 | chr1, chrM | 569,442, 8,894 |
| 31 | `...GCGAA----------` | 2 | 1 | 0 | chr1, chrM | 566,464, 5,915 |
| 32 | `...CTTA----------A` | 2 | 1 | 0 | chr1, chrM | 569,692, 9,144 |
| 33 | `...TGAGC----------` | 2 | 1 | 0 | chr1, chrM | 569,396, 8,848 |
| 34 | `...CAGGG----------` | 2 | 1 | 0 | chr1, chrM | 4,067, 564,616 |
| 35 | `...GTGT----------G` | 3 | 1 | 0 | chr1, chr17, chrM | 51,183,693, 567,968, 7,418 |
| 36 | `...G------------GC` | 2 | 1 | 0 | chr1, chrM | 566,537, 5,988 |
