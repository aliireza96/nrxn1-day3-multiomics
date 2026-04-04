#!/usr/bin/env python3
"""
07o_plot_figure_3B_rna_chip_scatter.py

Standalone manuscript scatter panel for the Figure 3 concordance slot:
- x-axis: RNA log2 fold change (patient / control)
- y-axis: H3K27me3 log2 fold change (patient / control)
- points: all shared RNA-H3K27me3 overlap genes, with concordant groups highlighted
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _repo_paths import PROJECT_ROOT


IN_TSV = os.path.join(PROJECT_ROOT, "outputs", "integration", "tables", "rna_chip_main_overlap.tsv")
OUT_PDF = os.path.join(PROJECT_ROOT, "outputs", "figures", "main", "Figure_3B_rna_chip_scatter.pdf")
OUT_PNG = os.path.join(PROJECT_ROOT, "outputs", "figures", "main", "Figure_3B_rna_chip_scatter.png")

C_GAIN = "#d1495b"   # H3K27me3 gained / RNA down
C_LOST = "#2f6bff"   # H3K27me3 lost / RNA up
C_NON = "#b8b8b8"    # non-concordant overlap genes


def load_rows(path):
    gained = []
    lost = []
    non_concordant = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                rna = float(row["rna_log2fc"])
                chip = float(row["chip_log2fc"])
            except (TypeError, ValueError, KeyError):
                continue

            if chip > 0 and rna < 0:
                gained.append((rna, chip))
            elif chip < 0 and rna > 0:
                lost.append((rna, chip))
            else:
                non_concordant.append((rna, chip))
    return gained, lost, non_concordant


def main():
    gained, lost, non_concordant = load_rows(IN_TSV)
    all_pts = gained + lost + non_concordant
    if not all_pts:
        raise SystemExit("No concordant RNA-H3K27me3 overlap genes found.")

    fig, ax = plt.subplots(figsize=(4.8, 4.8))

    if non_concordant:
        ax.scatter(
            [p[0] for p in non_concordant],
            [p[1] for p in non_concordant],
            s=14,
            c=C_NON,
            alpha=0.55,
            linewidths=0,
            label=f"Non-concordant overlap (n={len(non_concordant)})",
            zorder=1,
        )
    if gained:
        gained_xs = [p[0] for p in gained]
        gained_ys = [p[1] for p in gained]
        ax.scatter(
            gained_xs,
            gained_ys,
            s=34,
            c="white",
            alpha=0.9,
            linewidths=0,
            zorder=2,
        )
        ax.scatter(
            gained_xs,
            gained_ys,
            s=22,
            c=C_GAIN,
            alpha=0.92,
            edgecolors="#8f2230",
            linewidths=0.35,
            label=f"H3K27me3 gained / RNA down (n={len(gained)})",
            zorder=3,
        )
    if lost:
        ax.scatter(
            [p[0] for p in lost],
            [p[1] for p in lost],
            s=18,
            c=C_LOST,
            alpha=0.75,
            linewidths=0,
            label=f"H3K27me3 lost / RNA up (n={len(lost)})",
            zorder=2,
        )

    ax.axhline(0, color="#999", lw=0.7, ls="--", zorder=0)
    ax.axvline(0, color="#999", lw=0.7, ls="--", zorder=0)

    ax.set_xlim(-10, 10)
    ax.set_ylim(-4, 4)
    ax.set_xticks([-10, -5, 0, 5, 10])
    ax.set_yticks([-4, -2, 0, 2, 4])
    ax.set_xlabel("RNA log$_2$ fold change (patient / control)", fontsize=11)
    ax.set_ylabel("H3K27me3 log$_2$ fold change (patient / control)", fontsize=11)
    ax.tick_params(labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.text(
        0.02, 0.98,
        f"RNA-H3K27me3 overlap genes: n={len(all_pts)}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        color="#444",
    )
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        frameon=True,
        edgecolor="#ccc",
        fontsize=9,
        markerscale=1.0,
    )

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.24)
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_PDF}")
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
