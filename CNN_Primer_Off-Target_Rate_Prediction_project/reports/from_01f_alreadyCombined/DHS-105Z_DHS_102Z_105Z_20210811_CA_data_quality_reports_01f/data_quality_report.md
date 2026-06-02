# Data Quality Report — DHS-105Z_DHS_102Z_105Z_20210811_CA collapsed alignments

_Generated: 2026-05-31 18:14:04_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 9,443 |
| Groups with mixed labels      | 49 |
| Fraction of mixed groups      | 0.518903% |
| Rows in mixed groups          | 101 |
| Unique primers in mixed groups| 34 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...GATCG---------G` | 2 | 1 | 0 | chr5, chrM | 1,235, 79,947,303 |
| 2 | `...CAAATG---------` | 2 | 1 | 0 | chr1, chrM | 567,262, 6,712 |
| 3 | `...GCCAA---------C` | 2 | 1 | 0 | chr1, chrM | 568,033, 7,483 |
| 4 | `...AATATG---------` | 2 | 1 | 0 | chr1, chrM | 567,441, 6,891 |
| 5 | `...GATAGG---------` | 2 | 1 | 0 | chr1, chrM | 4,391, 564,940 |
| 6 | `...GCTTAG---------` | 2 | 1 | 0 | chr1, chrM | 566,257, 5,708 |
| 7 | `...CTAAG---------C` | 2 | 1 | 0 | chr1, chrM | 566,245, 5,696 |
| 8 | `...CCTTAA---------` | 2 | 1 | 0 | chr1, chrM | 569,691, 9,143 |
| 9 | `...ATGAGC---------` | 2 | 1 | 0 | chr1, chrM | 569,395, 8,847 |
| 10 | `...TTCTC----------` | 2 | 1 | 0 | chr1, chrM | 4,347, 564,896 |
| 11 | `...ATCGG----------` | 3 | 1 | 0 | chr11, chr5, chrM | 10,531,173, 1,234, 79,947,304 |
| 12 | `...TTCGA----------` | 2 | 1 | 0 | chr1, chrM | 567,893, 7,343 |
| 13 | `...AAATG----------` | 2 | 1 | 0 | chr1, chrM | 567,261, 6,711 |
| 14 | `...AATGG----------` | 2 | 1 | 0 | chr1, chrM | 5,661, 566,210 |
| 15 | `...TCA----------CG` | 2 | 1 | 0 | chr1, chrM | 4,648, 565,197 |
| 16 | `...TTAG----------C` | 2 | 1 | 0 | chrM, chrX | 125,606,618, 786 |
| 17 | `...CAAC----------A` | 2 | 1 | 0 | chrM, chrX | 2,965, 55,208,097 |
| 18 | `...CTATT----------` | 2 | 1 | 0 | chr1, chrM | 567,645, 7,095 |
| 19 | `...TCAC----------A` | 2 | 1 | 0 | chr1, chrM | 5,275, 565,824 |
| 20 | `...CCAAC----------` | 2 | 1 | 0 | chr1, chrM | 568,034, 7,484 |
| 21 | `...ATA----------TG` | 2 | 1 | 0 | chr1, chrM | 567,442, 6,892 |
| 22 | `...TCAGG----------` | 2 | 1 | 0 | chr1, chrM | 570,196, 9,648 |
| 23 | `...GA----------GGG` | 2 | 1 | 0 | chr1, chrM | 568,839, 8,291 |
| 24 | `...ATAG----------G` | 2 | 1 | 0 | chr1, chrM | 4,390, 564,939 |
| 25 | `...CTTAG----------` | 2 | 1 | 0 | chr1, chrM | 566,256, 5,707 |
| 26 | `...TAAGC----------` | 2 | 1 | 0 | chr1, chrM | 566,246, 5,697 |
| 27 | `...AGAT----------G` | 2 | 1 | 0 | chr1, chrM | 4,501, 565,050 |
| 28 | `...ATATG----------` | 2 | 1 | 0 | chr1, chrM | 567,369, 6,819 |
| 29 | `...CAATG----------` | 2 | 1 | 0 | chr5, chrM | 14,698, 93,905,195 |
| 30 | `...AAAA----------G` | 2 | 1 | 0 | chr1, chrM | 4,275, 564,824 |
| 31 | `...TAAGG----------` | 2 | 1 | 0 | chr1, chrM | 569,402, 8,854 |
| 32 | `...TAGGC----------` | 2 | 1 | 0 | chr1, chrM | 566,535, 5,986 |
| 33 | `...TAAAA----------` | 2 | 1 | 0 | chr1, chrM | 4,558, 565,107 |
| 34 | `...TATGG----------` | 2 | 1 | 0 | chr1, chrM | 567,927, 7,377 |
| 35 | `...TACA----------C` | 3 | 1 | 0 | chr11, chr5, chrM | 10,530,932, 1,475, 79,947,063 |
| 36 | `...CATG----------T` | 2 | 1 | 0 | chr1, chrM | 5,189, 565,738 |
| 37 | `...AATCG----------` | 2 | 1 | 0 | chr1, chrM | 4,149, 564,698 |
| 38 | `...TAATC----------` | 2 | 1 | 0 | chr1, chrM | 569,442, 8,894 |
| 39 | `...GCGAA----------` | 2 | 1 | 0 | chr1, chrM | 566,464, 5,915 |
| 40 | `...CTTA----------A` | 2 | 1 | 0 | chr1, chrM | 569,692, 9,144 |
| 41 | `...TGAGC----------` | 2 | 1 | 0 | chr1, chrM | 569,396, 8,848 |
| 42 | `...CAGGG----------` | 2 | 1 | 0 | chr1, chrM | 4,067, 564,616 |
| 43 | `...GTGT----------G` | 3 | 1 | 0 | chr1, chr17, chrM | 51,183,693, 567,968, 7,418 |
| 44 | `...ATGG-----------` | 2 | 1 | 0 | chr1, chrM | 5,662, 566,211 |
| 45 | `...A-----------GGG` | 2 | 1 | 0 | chr1, chrM | 568,838, 8,290 |
| 46 | `...AAA-----------A` | 2 | 1 | 0 | chr1, chrM | 4,557, 565,106 |
| 47 | `...AAA------------` | 2 | 1 | 0 | chr1, chrM | 4,556, 565,105 |
| 48 | `...AG-------------` | 2 | 1 | 0 | chr1, chrM | 4,278, 564,827 |
| 49 | `...---------------` | 2 | 1 | 0 | chr1, chrM | 568,834, 8,286 |
