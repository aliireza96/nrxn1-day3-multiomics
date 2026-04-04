#!/usr/bin/env python3
"""
Figure 2F: Chromatin regulatory machinery is disrupted at expression and
splicing levels in NRXN1alpha-null cells. Entry point is unbiased topGO
enrichment of spliced genes in chromatin remodeling (GO:0006338, 179 genes,
p=6.5e-8, Fig 2E). This panel zooms into mechanistically relevant factors:
JARID2 (PRC2 targeting cofactor, down) predicts H3K27me3 redistribution;
EZH2 dominant transcript (up) suggests maintained catalytic capacity; CBX7
(Polycomb domain reader, down) suggests impaired domain maintenance.
Filled circles = gene also has significant alternative splicing events (rMATS).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from _repo_paths import OUTPUTS_DIR


FIG_W = 6.5
FIG_H = 5.0
C_UP = "#d1495b"
C_DOWN = "#2f6bff"
C_REF = "#999999"
C_GROUP = "#777777"


GENE_ROWS = [
    {
        "gene": "JARID2",
        "label": "JARID2",
        "group": "Polycomb / PRC2",
        "log2fc": -0.984,
        "padj": 4.67e-7,
        "layer": "gene",
    },
    {
        "gene": "EZH1",
        "label": "EZH1",
        "group": "Polycomb / PRC2",
        "log2fc": -0.552,
        "padj": 1.54e-4,
        "layer": "gene",
    },
    {
        "gene": "EZH2",
        "label": "EZH2*",
        "group": "Polycomb / PRC2",
        "log2fc": 0.909,
        "padj": 3.92e-2,
        "layer": "dominant transcript",
    },
    {
        "gene": "CBX7",
        "label": "CBX7",
        "group": "Polycomb / PRC2",
        "log2fc": -1.334,
        "padj": 6.40e-4,
        "layer": "gene",
    },
    {
        "gene": "RBBP7",
        "label": "RBBP7",
        "group": "Polycomb / PRC2",
        "log2fc": 0.499,
        "padj": 5.44e-6,
        "layer": "gene",
    },
]


def find_project_root() -> Path:
    here = Path(__file__).resolve()
    for parent in [here.parent, *here.parents]:
        if (parent / ".nrxn1_project_root").exists():
            return parent
    raise RuntimeError("Could not find project root sentinel .nrxn1_project_root")

PROJECT_ROOT = find_project_root()
OUT_MAIN_PDF = OUTPUTS_DIR / "figures" / "main" / "Figure_2F_prc2_splicing.pdf"
OUT_MAIN_PNG = OUTPUTS_DIR / "figures" / "main" / "Figure_2F_prc2_splicing.png"
OUT_PY_PDF = OUTPUTS_DIR / "figures" / "main" / "python_plots" / "Figure_2F_chromatin_remodeler_panel.pdf"


def main() -> None:
    plt.rcParams["font.family"] = "DejaVu Sans"

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    y_positions = list(range(len(GENE_ROWS)))

    for y, row in zip(y_positions, GENE_ROWS):
        color = C_UP if row["log2fc"] > 0 else C_DOWN
        ax.hlines(
            y,
            0,
            row["log2fc"],
            color=color,
            linewidth=2,
            alpha=0.7,
            zorder=1,
        )
        ax.scatter(
            row["log2fc"],
            y,
            s=90,
            marker="o",
            facecolor=color,
            edgecolor="none",
            zorder=3,
        )

    ax.axvline(0, color=C_REF, linewidth=0.8, linestyle="--", zorder=0)

    ax.set_xlim(-2.0, 2.0)
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.set_ylim(-0.7, len(GENE_ROWS) - 0.3)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([row["label"] for row in GENE_ROWS], fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("log2 fold change (patient vs. control)", fontsize=10)
    ax.tick_params(axis="y", length=0, pad=4)
    ax.tick_params(axis="x", labelsize=9)

    for label in ax.get_yticklabels():
        label.set_fontstyle("italic")

    ax.text(
        0.02,
        1.03,
        "Decreased in patient",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        color=C_DOWN,
        fontweight="bold",
    )
    ax.text(
        0.98,
        1.03,
        "Increased in patient",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color=C_UP,
        fontweight="bold",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.22, right=0.80, top=0.84, bottom=0.27)
    fig.text(
        0.22,
        0.105,
        "All genes shown: FDR < 0.05",
        ha="left",
        va="bottom",
        fontsize=8,
        color="#555555",
    )
    fig.text(
        0.22,
        0.030,
        "* EZH2: dominant transcript (ENST00000460911, log2FC +0.91);\n"
        "  gene-level log2FC +0.28, padj 0.026",
        ha="left",
        va="bottom",
        fontsize=7.5,
        color="#777777",
    )

    OUT_MAIN_PDF.parent.mkdir(parents=True, exist_ok=True)
    OUT_PY_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_MAIN_PDF, bbox_inches="tight")
    fig.savefig(OUT_MAIN_PNG, dpi=200, bbox_inches="tight")
    fig.savefig(OUT_PY_PDF, bbox_inches="tight")
    plt.close(fig)

    print("gene\tgroup\tlog2FC\tpadj\tlayer")
    for row in GENE_ROWS:
        print(
            f"{row['gene']}\t{row['group']}\t{row['log2fc']:+.3f}\t"
            f"{row['padj']:.2e}\t{row['layer']}"
        )
    print(f"Saved: {OUT_MAIN_PDF}")
    print(f"Saved: {OUT_MAIN_PNG}")
    print(f"Saved: {OUT_PY_PDF}")


if __name__ == "__main__":
    main()
