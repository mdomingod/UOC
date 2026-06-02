# Data Quality Report — DHS-105Z_M05215_DHS003Z_TEPCR_Buffer_4_Lot_Comparison_Test_105Z_20181008_AZ-99172080 collapsed alignments

_Generated: 2026-05-31 18:14:10_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 6,329 |
| Groups with mixed labels      | 21 |
| Fraction of mixed groups      | 0.331806% |
| Rows in mixed groups          | 43 |
| Unique primers in mixed groups| 21 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...TTCTC----------` | 2 | 1 | 0 | chr1, chrM | 4,347, 564,896 |
| 2 | `...ATCGG----------` | 2 | 1 | 0 | chr11, chrM | 10,531,173, 1,234 |
| 3 | `...AATGG----------` | 2 | 1 | 0 | chr1, chrM | 5,661, 566,210 |
| 4 | `...TCA----------CG` | 2 | 1 | 0 | chr1, chrM | 4,648, 565,197 |
| 5 | `...CTATT----------` | 2 | 1 | 0 | chr1, chrM | 567,645, 7,095 |
| 6 | `...CCAAC----------` | 2 | 1 | 0 | chr1, chrM | 568,034, 7,484 |
| 7 | `...ATA----------TG` | 2 | 1 | 0 | chr1, chrM | 567,442, 6,892 |
| 8 | `...GA----------GGG` | 2 | 1 | 0 | chr1, chrM | 568,839, 8,291 |
| 9 | `...CTTAG----------` | 2 | 1 | 0 | chr1, chrM | 566,256, 5,707 |
| 10 | `...AGAT----------G` | 2 | 1 | 0 | chr1, chrM | 4,501, 565,050 |
| 11 | `...ATATG----------` | 2 | 1 | 0 | chr1, chrM | 567,369, 6,819 |
| 12 | `...CAATG----------` | 2 | 1 | 0 | chr5, chrM | 14,698, 93,905,195 |
| 13 | `...TAAGG----------` | 2 | 1 | 0 | chr1, chrM | 569,402, 8,854 |
| 14 | `...TAGGC----------` | 2 | 1 | 0 | chr1, chrM | 566,535, 5,986 |
| 15 | `...TACA----------C` | 3 | 1 | 0 | chr11, chr5, chrM | 10,530,932, 1,475, 79,947,063 |
| 16 | `...CATG----------T` | 2 | 1 | 0 | chr1, chrM | 5,189, 565,738 |
| 17 | `...AATCG----------` | 2 | 1 | 0 | chr1, chrM | 4,149, 564,698 |
| 18 | `...GCGAA----------` | 2 | 1 | 0 | chr1, chrM | 566,464, 5,915 |
| 19 | `...TGAGC----------` | 2 | 1 | 0 | chr1, chrM | 569,396, 8,848 |
| 20 | `...CAGGG----------` | 2 | 1 | 0 | chr1, chrM | 4,067, 564,616 |
| 21 | `...GTGT----------G` | 2 | 1 | 0 | chr1, chrM | 567,968, 7,418 |
