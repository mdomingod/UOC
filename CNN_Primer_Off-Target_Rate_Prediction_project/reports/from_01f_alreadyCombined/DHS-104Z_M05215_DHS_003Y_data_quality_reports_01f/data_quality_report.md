# Data Quality Report — DHS-104Z_M05215_DHS_003Y collapsed alignments

_Generated: 2026-05-31 18:13:47_

## Mixed on-target / off-target label analysis

After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report 
identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.

### Summary

| Metric | Value |
|---|---|
| Total alignment groups        | 8,017 |
| Groups with mixed labels      | 3 |
| Fraction of mixed groups      | 0.037420% |
| Rows in mixed groups          | 6 |
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
| 3 | `...--------------G` | 2 | 1 | 0 | chr10 | 96,535,202, 96,702,004 |
