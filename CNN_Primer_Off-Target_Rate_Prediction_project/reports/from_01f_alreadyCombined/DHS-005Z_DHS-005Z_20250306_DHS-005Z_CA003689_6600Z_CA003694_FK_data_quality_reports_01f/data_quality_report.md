# Data Quality Report — DHS-005Z_DHS-005Z_20250306_DHS-005Z_CA003689_6600Z_CA003694_FK collapsed alignments

_Generated: 2026-05-31 18:04:03_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 96,656 |
| Groups with mixed labels      | 5 |
| Fraction of mixed groups      | 0.005173% |
| Rows in mixed groups          | 10 |
| Unique primers in mixed groups| 5 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...TAAC----------A` | 2 | 10 | 0 | chr10, chr9 | 33,675,559, 89,725,107 |
| 2 | `...TT----------TTC` | 2 | 1 | 0 | chr19 | 9,014,741, 9,017,555 |
| 3 | `...CCTAC----------` | 2 | 1 | 0 | chr19 | 9,014,721, 9,017,535 |
| 4 | `...ACAGC----------` | 2 | 10 | 0 | chr7 | 126,544,054, 65,082,686 |
| 5 | `...CTTGC----------` | 2 | 10 | 0 | chr22, chr3 | 17,055,458, 178,938,877 |
