# scripts_marta/data_quality_report.py
from pathlib import Path
from datetime import datetime
import pandas as pd


def compute_mixed_labels_report(combined_data: pd.DataFrame) -> dict:
    """
    Run the mixed-label diagnostic on the combined_data dataframe.
    Returns a dict of computed statistics — no I/O, no printing.
    """
    # Find groups with mixed isIntendedSite labels
    mixed = (
        combined_data
        .groupby(['prm_aligned', 'gsq_aligned'])['isIntendedSite']
        .agg(['min', 'max'])
        .reset_index()
    )
    mixed_groups = mixed[mixed['min'] != mixed['max']]
    n_mixed = len(mixed_groups)

    n_total_groups = mixed['min'].size

    # Pull the original rows for the mixed groups
    mixed_keys = list(zip(mixed_groups['prm_aligned'], mixed_groups['gsq_aligned']))
    mixed_rows = combined_data[
        combined_data.set_index(['prm_aligned', 'gsq_aligned']).index.isin(mixed_keys)
    ].copy()

    # Build the per-group breakdown as a list of dicts
    per_group = []
    for (prm, gsq), grp in mixed_rows.groupby(['prm_aligned', 'gsq_aligned']):
        per_group.append({
            'primer_tail':         prm[-15:],
            'n_rows':              len(grp),
            'n_on_target':         int(grp['isIntendedSite'].sum()),
            'n_off_target':        int((~grp['isIntendedSite'].astype(bool)).sum()),
            'positions':           sorted(grp['read1_anchor'].unique().tolist()),
            'chromosomes':         sorted(grp['alignChrom'].unique().tolist()),
            'n_unique_chromosomes': grp['alignChrom'].nunique(),
        })

    return {
        'n_total_groups':         int(n_total_groups),
        'n_mixed_groups':         int(n_mixed),
        'fraction_mixed':         float(n_mixed / n_total_groups) if n_total_groups else 0.0,
        'n_rows_in_mixed_groups': int(len(mixed_rows)),
        'n_unique_primers':       int(mixed_rows['primer'].nunique()) if n_mixed else 0,
        'per_group':              per_group,
    }


def render_markdown_report(report: dict, dataset_name: str = "Combined dataset") -> str:
    """
    Convert the structured report dict into a Markdown string.
    """
    md = []
    md.append(f"# Data Quality Report — {dataset_name}")
    md.append(f"")
    md.append(f"_Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}_")
    md.append(f"")
    md.append(f"## Mixed on-target / off-target label analysis")
    md.append(f"")
    md.append(f"After collapsing alignments by `(prm_aligned, gsq_aligned)`, this report ")
    md.append(f"identifies groups whose contributing rows have inconsistent `isIntendedSite` labels.")
    md.append(f"")

    # Summary table
    md.append(f"### Summary")
    md.append(f"")
    md.append(f"| Metric | Value |")
    md.append(f"|---|---|")
    md.append(f"| Total alignment groups        | {report['n_total_groups']:,} |")
    md.append(f"| Groups with mixed labels      | {report['n_mixed_groups']:,} |")
    md.append(f"| Fraction of mixed groups      | {report['fraction_mixed']:.6%} |")
    md.append(f"| Rows in mixed groups          | {report['n_rows_in_mixed_groups']:,} |")
    md.append(f"| Unique primers in mixed groups| {report['n_unique_primers']:,} |")
    md.append(f"")

    # Interpretation
    md.append(f"### Interpretation")
    md.append(f"")
    if report['n_mixed_groups'] == 0:
        md.append(f"No mixed groups detected. The on-target labeling is consistent across all collapsed alignments.")
    else:
        md.append(f"A small number of alignment patterns occur at multiple genomic locations, ")
        md.append(f"where one location is the designed target and others are off-target hits. ")
        md.append(f"These groups are resolved using `max(isIntendedSite)` during collapse, ")
        md.append(f"labeling the alignment pattern as on-target if any contributing site was designed as such.")
    md.append(f"")

    # Per-group details
    if report['per_group']:
        md.append(f"### Per-group breakdown")
        md.append(f"")
        md.append(f"| # | Primer tail | Rows | On-target | Off-target | Chromosomes | Positions |")
        md.append(f"|---|---|---:|---:|---:|---|---|")
        for i, g in enumerate(report['per_group'], 1):
            chroms = ", ".join(g['chromosomes'])
            positions = ", ".join(f"{int(p):,}" for p in g['positions'])
            md.append(
                f"| {i} | `...{g['primer_tail']}` "
                f"| {g['n_rows']} "
                f"| {g['n_on_target']} "
                f"| {g['n_off_target']} "
                f"| {chroms} "
                f"| {positions} |"
            )
        md.append(f"")

    return "\n".join(md)


def save_report(report: dict, output_dir: Path, dataset_name: str = "Combined dataset"):
    """
    Save the report in both Markdown (for humans) and JSON (for machines).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Markdown
    md_path = output_dir / "data_quality_report.md"
    md_path.write_text(render_markdown_report(report, dataset_name))

    # JSON — for programmatic access later
    import json
    json_path = output_dir / "data_quality_report.json"
    json_path.write_text(json.dumps(report, indent=2, default=str))

    print(f"Saved: {md_path}")
    print(f"Saved: {json_path}")