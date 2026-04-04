#!/usr/bin/env python3
"""
07p2_plot_figure_3A_chip_annotation_diverging.py

H3K27me3 ChIP-seq annotation diverging bar chart for Figure 3A.
Mirrors the manuscript diverging-bar style used for the current ATAC Figure 4A panel.
  left  = Lost in patient  (less H3K27me3)
  right = Gained in patient (more H3K27me3)
Split by genomic annotation class.
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from _repo_paths import PROJECT_ROOT
IN_TSV = os.path.join(PROJECT_ROOT, "outputs", "figures", "tables",
                      "Figure_3A_chip_annotation_groups.tsv")
OUT_PDF      = os.path.join(PROJECT_ROOT, "outputs", "figures", "main",
                            "python_plots", "Figure_3A_chip_annotation_diverging.pdf")
OUT_PDF_MAIN = os.path.join(PROJECT_ROOT, "outputs", "figures", "main",
                            "Figure_3A_chip_annotation_diverging.pdf")
OUT_PNG      = os.path.join(PROJECT_ROOT, "outputs", "figures", "main",
                            "python_plots", "Figure_3A_chip_annotation_diverging.png")

DIR_LEFT  = "Lost in patient"      # less H3K27me3
DIR_RIGHT = "Gained in patient"    # more H3K27me3
C_LEFT    = "#2f6bff"              # blue  — consistent with ATAC figure
C_RIGHT   = "#d1495b"              # red   — consistent with ATAC figure

# Promoter first; Distal replaces Intergenic vs the ATAC version; Other excluded (n=1)
GROUP_ORDER = ["Promoter", "Distal", "Intron", "Exon", "UTR"]


def load_counts(path):
    counts = {
        DIR_LEFT:  {g: 0 for g in GROUP_ORDER},
        DIR_RIGHT: {g: 0 for g in GROUP_ORDER},
    }
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            direction = row.get("direction", "").strip()
            group     = row.get("annotation_group", "").strip()
            if direction not in counts or group not in counts[direction]:
                continue
            try:
                counts[direction][group] = int(row["n_regions"])
            except (TypeError, ValueError, KeyError):
                continue
    return counts


def main():
    os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)

    counts     = load_counts(IN_TSV)
    left_vals  = [-counts[DIR_LEFT][g]  for g in GROUP_ORDER]
    right_vals = [ counts[DIR_RIGHT][g] for g in GROUP_ORDER]
    y_pos      = list(range(len(GROUP_ORDER)))

    max_abs  = max(max(abs(v) for v in left_vals), max(abs(v) for v in right_vals))
    text_pad = max_abs * 0.03

    fig, ax = plt.subplots(figsize=(6.0, 4.15))

    ax.barh(y_pos, left_vals,  color=C_LEFT,  alpha=0.90, height=0.62)
    ax.barh(y_pos, right_vals, color=C_RIGHT, alpha=0.90, height=0.62)

    ax.axvline(0, color="#9e9e9e", lw=0.9)
    ax.grid(axis="x", color="#ebebeb", linewidth=0.7)
    ax.set_axisbelow(True)

    # Count labels
    for y, val in zip(y_pos, left_vals):
        if abs(val) > 0:
            ax.text(val - text_pad, y, f"{abs(val):,}",
                    ha="right", va="center", fontsize=8.5, color="#333")
    for y, val in zip(y_pos, right_vals):
        if val > 0:
            ax.text(val + text_pad, y, f"{val:,}",
                    ha="left", va="center", fontsize=8.5, color="#333")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(GROUP_ORDER, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("Differential H3K27me3 peaks", fontsize=10)
    ax.tick_params(axis="x", labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_xlim(-max_abs * 1.38, max_abs * 1.38)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{abs(int(x)):,}"))

    ax.text(0.02, 1.04, "Less H3K27me3 in patient",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9.5, color=C_LEFT, fontweight="bold")
    ax.text(0.98, 1.04, "More H3K27me3 in patient",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9.5, color=C_RIGHT, fontweight="bold")
    ax.text(0.5, -0.17, "H3K27me3 differential peaks by genomic annotation (FDR < 0.05)",
            transform=ax.transAxes, ha="center", va="top",
            fontsize=7.5, color="#777")

    fig.tight_layout()
    fig.savefig(OUT_PDF,      bbox_inches="tight")
    fig.savefig(OUT_PDF_MAIN, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_PDF}")
    print(f"Saved: {OUT_PDF_MAIN}")
    print(f"Saved: {OUT_PNG}")

    # Summary to terminal
    print("\nH3K27me3 diverging bar — counts:")
    print(f"{'Category':<12} {'Lost':>8} {'Gained':>8}")
    for g, lv, rv in zip(GROUP_ORDER, left_vals, right_vals):
        print(f"{g:<12} {abs(lv):>8,} {rv:>8,}")


if __name__ == "__main__":
    main()
