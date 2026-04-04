#!/usr/bin/env python3
"""
07g2_plot_figure_4A_promoter_dars.py

Promoter-focused Figure 4A panel for Nature journal submission.

Input:
  outputs/atac_seq/atac_main_local_controls/tables/promoter_dars_standardized.tsv

Outputs:
  outputs/figures/main/Figure_4A_promoter_dars.pdf
  outputs/figures/main/Figure_4A_promoter_dars.png
  outputs/figures/main/python_plots/Figure_4A_promoter_dars.pdf
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from _repo_paths import PROJECT_ROOT, OUTPUTS_DIR

IN_TSV = OUTPUTS_DIR / "atac_seq" / "atac_main_local_controls" / "tables" / "promoter_dars_standardized.tsv"
OUT_MAIN_PDF = OUTPUTS_DIR / "figures" / "main" / "Figure_4A_promoter_dars.pdf"
OUT_MAIN_PNG = OUTPUTS_DIR / "figures" / "main" / "Figure_4A_promoter_dars.png"
OUT_PYTHON_PDF = OUTPUTS_DIR / "figures" / "main" / "python_plots" / "Figure_4A_promoter_dars.pdf"

PADJ_THRESH = 0.05
COLOR_MORE = "#d1495b"
COLOR_LESS = "#2f6bff"
CAPTION = "Promoter DARs (HOMER promoter-TSS annotation); 3,716 total."


def load_counts(path):
    total = 0
    more_open = 0
    less_open = 0

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                lfc = float(row["log2FoldChange"])
                padj_raw = row.get("padj", "")
                padj = float(padj_raw) if padj_raw not in ("", "NA", "nan") else 1.0
            except (KeyError, ValueError):
                continue
            if padj >= PADJ_THRESH:
                continue
            total += 1
            if lfc > 0:
                more_open += 1
            elif lfc < 0:
                less_open += 1

    return total, more_open, less_open


def style_matplotlib():
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
    plt.rcParams["font.size"] = 8
    plt.rcParams["axes.labelsize"] = 8
    plt.rcParams["xtick.labelsize"] = 8
    plt.rcParams["ytick.labelsize"] = 8


def plot_counts(total, more_open, less_open):
    labels = ["More open\nin patient", "Less open\nin patient"]
    values = [more_open, less_open]
    colors = [COLOR_MORE, COLOR_LESS]
    y_pos = [0.0, 0.52]

    fig, ax = plt.subplots(figsize=(3.0, 2.8))

    max_val = max(values) if values else 1
    outside_pad = max_val * 0.02

    ax.barh(y_pos, values, color=colors, height=0.22, edgecolor="none")

    for y, value, color in zip(y_pos, values, colors):
        ax.text(
            value + outside_pad,
            y,
            f"{value:,}",
            ha="left",
            va="center",
            fontsize=8,
            color=color,
            clip_on=False,
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Number of promoter DARs (padj < 0.05)", fontsize=8)
    ax.tick_params(axis="y", length=0, labelsize=8)
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)

    ax.set_xlim(0, max_val * 1.18)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):,}" if x >= 0 else ""))
    ax.set_ylim(0.72, -0.20)

    fig.text(0.5, 0.02, CAPTION, ha="center", va="bottom", fontsize=7, color="#666666")
    fig.tight_layout(rect=(0, 0.09, 1, 1))

    return fig


def main():
    style_matplotlib()
    total, more_open, less_open = load_counts(IN_TSV)
    fig = plot_counts(total, more_open, less_open)

    OUT_MAIN_PDF.parent.mkdir(parents=True, exist_ok=True)
    OUT_PYTHON_PDF.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(OUT_MAIN_PDF, bbox_inches="tight")
    fig.savefig(OUT_MAIN_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PYTHON_PDF, bbox_inches="tight")
    plt.close(fig)

    print(f"Total promoter DARs: {total:,}")
    print(f"More open in patient: {more_open:,}")
    print(f"Less open in patient: {less_open:,}")
    print(f"Saved: {OUT_MAIN_PDF}")
    print(f"Saved: {OUT_MAIN_PNG}")
    print(f"Saved: {OUT_PYTHON_PDF}")


if __name__ == "__main__":
    main()
