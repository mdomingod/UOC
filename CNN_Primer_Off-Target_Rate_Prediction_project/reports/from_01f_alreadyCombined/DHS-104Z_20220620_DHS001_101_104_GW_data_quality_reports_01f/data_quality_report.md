# Data Quality Report — DHS-104Z_20220620_DHS001_101_104_GW collapsed alignments

_Generated: 2026-05-31 18:13:39_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 6,841 |
| Groups with mixed labels      | 2 |
| Fraction of mixed groups      | 0.029235% |
| Rows in mixed groups          | 4 |
| Unique primers in mixed groups| 1 |

### Interpretation

A small number of alignment patterns occur at multiple genomic locations, 
where one location is the designed target and others are off-target hits. 
These groups are resolved using `max(isIntendedSite)` during collapse, 
labeling the alignment pattern as on-target if any contributing site was designed as such.

### Per-group breakdown

| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |
|---|---|---:|---:|---:|---|---|
| 1 | `...TCATG----------` | 2 | 1 | 0 | chr10 | 96,535,206, 96,702,008 |
| 2 | `...ATG------------` | 2 | 1 | 0 | chr10 | 96,535,204, 96,702,006 |
