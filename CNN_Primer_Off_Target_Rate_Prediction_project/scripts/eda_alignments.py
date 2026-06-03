#!/usr/bin/env python3
"""
EDA starter pack for the raw SPE-primer off-target data: {run_name}_alignments.csv

Expected schema (one row per primer-target pair):
    primer | isIntendedSite | alignChrom | alignStrand | read1_anchor | umi_groups_cnt | genomic_seq

The script answers the questions that matter for this thesis:
  - How is the regression target (umi_groups_cnt) distributed, and how do
    on-target and off-target sites differ?  -> the training signal
  - How severe is the zero/low-count imbalance?               -> modelling difficulty
  - Do simple features (GC, length) explain amplification?    -> motivates the CNN
  - What is the panel-level off-target rate?                   -> the quantity H2 predicts

Usage:
    python eda_alignments.py RUN_alignments.csv               # one run
    python eda_alignments.py path/to/folder/                  # all *_alignments.csv in a folder

Outputs:
    ./eda_outputs/*.png   (figures)
    ./eda_outputs/eda_summary.csv   (metrics table)
Dependencies: pandas, numpy, matplotlib
"""

import sys
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = Path("eda_outputs_epcr07+exp01")   # !!! edit this to avoid overwriting previous EDA outputs
ON_COLOR, OFF_COLOR = "#2563eb", "#ea580c"   # on-target blue, off-target orange
COLS = ["primer", "isIntendedSite", "alignChrom", "alignStrand",
        "read1_anchor", "umi_groups_cnt", "genomic_seq"]


# --------------------------------------------------------------------------- #
# Load & derive
# --------------------------------------------------------------------------- #
def load(path: str) -> pd.DataFrame:
    p = Path(path)
    files = sorted(glob.glob(str(p / "*_alignments.csv"))) if p.is_dir() else [str(p)] #"*_alignments_from_data.csv" for only alignments_from_data
    if not files:
        sys.exit(f"No *_alignments.csv found in {path}")
    frames = []
    for f in files:
        df = pd.read_csv(f)
        df["run"] = Path(f).name.replace("_alignments.csv", "")
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)

    missing = [c for c in COLS if c not in df.columns]
    if missing:
        sys.exit(f"Missing expected columns: {missing}")

    # derived columns
    df["umi_groups_cnt"] = pd.to_numeric(df["umi_groups_cnt"], errors="coerce").fillna(0).astype(int)
    df["site_class"] = np.where(df["isIntendedSite"] == 1, "on-target", "off-target")
    df["is_zero"]    = df["umi_groups_cnt"] == 0
    df["primer_len"] = df["primer"].str.len()
    df["gc"]         = df["primer"].str.upper().str.count("[GC]") / df["primer_len"]
    df["genomic_len"] = df["genomic_seq"].str.len()
    return df


# --------------------------------------------------------------------------- #
# Summary metrics
# --------------------------------------------------------------------------- #
def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for run, g in [("ALL", df)] + list(df.groupby("run")):
        on  = g[g.site_class == "on-target"]
        off = g[g.site_class == "off-target"]
        umi_on, umi_off = on.umi_groups_cnt.sum(), off.umi_groups_cnt.sum()
        rate = umi_off / (umi_on + umi_off) if (umi_on + umi_off) else np.nan
        rows.append({
            "run": run,
            "n_sites": len(g),
            "n_primers": g.primer.nunique(),
            "n_on_target": len(on),
            "n_off_target": len(off),
            "pct_off_target_sites": 100 * len(off) / len(g),
            "pct_zero_umi": 100 * g.is_zero.mean(),
            "umi_mean_on": on.umi_groups_cnt.mean(),
            "umi_mean_off": off.umi_groups_cnt.mean(),
            "umi_max": g.umi_groups_cnt.max(),
            "panel_offtarget_rate_pct": 100 * rate,   # UMI-weighted: the H2 quantity
        })
    return pd.DataFrame(rows).set_index("run").round(3)


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
def _save(fig, name, subdir=None ):
    target = OUTDIR / subdir if subdir else OUTDIR
    OUTDIR.mkdir(exist_ok=True)
    target.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(target / name, dpi=150)
    plt.close(fig)
    print("saved", target / name)


import textwrap

def plot_class_balance(df, per_fig=4, ncols=2):
    """On-target vs off-target counts, one subplot per panel (run),
    paginated across several figures (default 4 panels per figure, 2x2).
    24 panels -> 6 images named 01_class_balance_p1.png ... _p6.png."""
    runs = sorted(df["run"].unique())
    order, colors = ["off-target", "on-target"], [OFF_COLOR, ON_COLOR]
    nrows = int(np.ceil(per_fig / ncols))

    # split the run list into pages of `per_fig`
    pages = [runs[i:i + per_fig] for i in range(0, len(runs), per_fig)]

    for p, page in enumerate(pages, start=1):
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(4.0 * ncols, 3.6 * nrows),
                                 squeeze=False, sharey=True)
        for ax, run in zip(axes.ravel(), page):
            g = df[df["run"] == run]
            counts = g["site_class"].value_counts().reindex(order, fill_value=0)
            ax.bar(order, counts.values, color=colors)
            for i, v in enumerate(counts.values):
                pct = 100 * v / len(g) if len(g) else 0
                ax.text(i, v, f"{v:,}\n({pct:.1f}%)", ha="center", va="bottom", fontsize=8)
            ax.set_title("\n".join(textwrap.wrap((run).replace("_alignments_from_data.csv", ""), 22)), fontsize=8)
            ax.margins(y=0.18)
        for ax in axes.ravel():            # blank any unused cells on the last page
            if not ax.has_data():
                ax.axis("off")
        axes[0, 0].set_ylabel("number of sites")
        fig.suptitle(f"On-target vs off-target sites per panel "
                     f"(page {p} of {len(pages)})", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        balance_out = Path("01_class_balance_fig")
        _save(fig, subdir=balance_out, name=f"class_balance_p{p}.png")


def plot_umi_distribution(df, per_fig=4, ncols=2):
    """log1p UMI-count distribution by site class, one subplot per panel (run),
    paginated (default 4 panels per figure, 2x2).
    24 panels -> 6 images named 02_umi_distribution_p1.png ... _p6.png."""
    runs = sorted(df["run"].unique())
    nrows = int(np.ceil(per_fig / ncols))

    # shared bins across ALL panels so pages are comparable
    bins = np.linspace(0, np.log1p(df.umi_groups_cnt.max() + 1), 40)
    pages = [runs[i:i + per_fig] for i in range(0, len(runs), per_fig)]

    for p, page in enumerate(pages, start=1):
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(4.5 * ncols, 3.6 * nrows),
                                 squeeze=False, sharex=True)
        for ax, run in zip(axes.ravel(), page):
            g = df[df["run"] == run]
            for cls, color in [("off-target", OFF_COLOR), ("on-target", ON_COLOR)]: # !!!!
                vals = np.log1p(g.loc[g.site_class == cls, "umi_groups_cnt"]) # !!!!
                ax.hist(vals, bins=bins, alpha=0.6, label=cls, color=color) # !!!!
            ax.set_yscale("log") # !!!!
            ax.set_title("\n".join(textwrap.wrap((run).replace("_alignments_from_data.csv", ""), 22)), fontsize=8) # !!!!
            ax.set_xlabel("log(1 + umi_groups_cnt)") # !!!!
            ax.legend(fontsize=7)
        for ax in axes.ravel():            # blank unused cells on the last page
            if not ax.has_data():
                ax.axis("off")
        axes[0, 0].set_ylabel("sites (log)") # <<< 
        fig.suptitle(f"UMI-count distribution by site class " # <<<
                     f"(page {p} of {len(pages)})", fontsize=12) # <<<
        fig.tight_layout(rect=[0, 0, 1, 0.93]) # <<<
        umiDis_out = Path("02_umi_distribution_fig")
        _save(fig, subdir=umiDis_out, name=f"umi_distribution_p{p}.png")


def trim_title(text: str, max_len: int = 30, min_display: int = 10) -> str:
    """
    Trim text to max_len with ellipsis.
    Always shows at least min_display characters before the ellipsis.
    
    Examples with max_len=15, min_display=10:
      "Short text"           → "Short text"       (under max_len, no trim)
      "A longer title here"  → "A longer ti..."   (trimmed at 15)
      "AB"                   → "AB"               (under min_display, no trim)
    """
    if len(text) <= max_len:
        return text
    # Always show at least min_display characters before "..."
    cut = max(min_display, max_len - 3)   # -3 to account for "..."
    return text[:cut] + "..."

def plot_features_vs_target(df, per_fig: int = 4):
    """
    Per-run primer length and GC content distributions.
    A weak relationship between these features and UMI count
    is the argument for learning from raw sequence (CNN).
    """
    runs  = sorted(df["run"].unique())
    pages = [runs[i : i + per_fig] for i in range(0, len(runs), per_fig)]
    ncols = 2   # col 0 = primer length, col 1 = GC content

    for p, page in enumerate(pages, start=1):
        nrows = len(page)          # one row per run on this page
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(4.5 * ncols, 3.6 * nrows),
            squeeze=False,         # always 2-D axes array
        )

        for r, run in enumerate(page):
            g = df[df["run"] == run]     # per-run slice — used below

            # Column 0 — primer length
            axes[r, 0].hist(g["primer_len"].dropna(), bins=30, color="#0d9488")
            axes[r, 0].set_title(f"{trim_title((run).replace("_alignments_from_data.csv", ""))} — \nPrimer length")
            axes[r, 0].set_xlabel("bp")
            axes[r, 0].set_ylabel("sites")

            # Column 1 — GC content
            axes[r, 1].hist(g["gc"].dropna(), bins=30, color="#0d9488")
            axes[r, 1].set_title(f"{trim_title((run).replace("_alignments_from_data.csv", ""))} — \nGC content")
            axes[r, 1].set_xlabel("GC fraction")
            axes[r, 1].set_ylabel("sites")

        # Blank any unused cells on the last page
        for ax in axes.ravel():
            if not ax.has_data():
                ax.axis("off")

        fig.suptitle(
            f"Feature distributions by run "
            f"(page {p} of {len(pages)})",
            fontsize=12,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.93])

        feat_out = Path("05_features_vs_target_fig")
        _save(fig, subdir=feat_out, name=f"features_vs_target_p{p}.png")


def plot_genomic_distribution(df):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    off = df[df.site_class == "off-target"]
    chrom = off.alignChrom.astype(str).value_counts().head(25)
    axes[0].bar(chrom.index, chrom.values, color=OFF_COLOR)
    axes[0].set_title("Off-target sites in all panels per chromosome")
    axes[0].set_ylabel("off-target sites"); axes[0].tick_params(axis="x", rotation=90, labelsize=7)
    strand = df.groupby(["site_class", "alignStrand"]).size().unstack(fill_value=0)
    strand.plot(kind="bar", ax=axes[1], color=["#1d4ed8", "#60a5fa"])
    axes[1].set_title("Strand balance"); axes[1].set_xlabel(""); axes[1].tick_params(axis="x", rotation=0)
    axes[1].legend(title="alignStrand")
    _save(fig, "06_genomic_distribution.png")


def plot_panel_offtarget_rate(summary):
    summary.index = summary.index.str.replace("_alignments_from_data.csv", "")
    runs = summary.drop(index="ALL", errors="ignore")
    if len(runs) < 1:
        return
    fig, ax = plt.subplots(figsize=(max(5, 0.5* len(runs)), 10))
    ax.bar(runs.index, runs.panel_offtarget_rate_pct, color=OFF_COLOR)
    #ax.barh(runs.index, runs.panel_offtarget_rate_pct, color=OFF_COLOR)
    ax.set_ylabel("off-target rate (%)") 
    ax.set_title("Panel-level off-target rate (UMI-weighted)")
    ax.set_xticks(range(len(runs)))                  # fix tick positions first
    ax.set_xticklabels(runs.index, rotation=45, ha="right", fontsize=8)
    #ax.tick_params(axis="x", rotation=90, labelsize=8)
    fig.subplots_adjust(bottom=0.30)     # reserve space for rotated labels
    _save(fig, name="07_panel_offtarget_rate.png")

def plot_n_primers_by_run(summary):
    summary.index = summary.index.str.replace("_alignments_from_data.csv", "")
    runs = summary.drop(index="ALL", errors="ignore")
    if len(runs) < 1:
        return
    fig, ax = plt.subplots(figsize=(max(5, 0.5* len(runs)), 10))
    ax.bar(runs.index, runs.n_primers, color=OFF_COLOR)
    #ax.barh(runs.index, runs.panel_offtarget_rate_pct, color=OFF_COLOR)
    ax.set_ylabel("unique primers") 
    ax.set_title("Number of unique primers per run")
    ax.set_xticks(range(len(runs)))                  # fix tick positions first
    ax.set_xticklabels(runs.index, rotation=45, ha="right", fontsize=8)
    #ax.tick_params(axis="x", rotation=90, labelsize=8)
    fig.subplots_adjust(bottom=0.30)     # reserve space for rotated labels
    _save(fig, name="08_n_primers_by_run.png")


# --------------------------------------------------------------------------- #
def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    df = load(sys.argv[1])
    print(f"Loaded {len(df):,} sites from {df.run.nunique()} run(s).")

    summary = summarize(df)
    OUTDIR.mkdir(exist_ok=True)
    summary.to_csv( OUTDIR / "eda_summary.csv")
    print("\n=== Summary metrics ===")
    print(summary.to_string())

    plot_class_balance(df)
    plot_umi_distribution(df)
    plot_features_vs_target(df)
    plot_genomic_distribution(df)
    plot_panel_offtarget_rate(summary)
    plot_n_primers_by_run(summary)
    print("\nDone. Figures and eda_summary.csv are in ./eda_outputs/")


if __name__ == "__main__":
    main()
