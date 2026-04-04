#!/usr/bin/env python3
"""
07s_plot_figure_4C_rna_atac_go.py
Figure 4C — GO enrichment of RNA-ATAC concordant genes (topGO, BP).

Curated term selection: biologically relevant terms grouped thematically,
matching the style of Figures 2B and 3D.

Inputs (written by 08_go_enrichment.R):
  outputs/go_enrichment/tables/rna_atac_main_up_topgo_bp_full.tsv
  outputs/go_enrichment/tables/rna_atac_main_down_topgo_bp_full.tsv

Output:
  outputs/figures/main/Figure_4C_rna_atac_go.pdf
"""

import os
import csv
import math
import textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from _repo_paths import PROJECT_ROOT

# ── Paths ──────────────────────────────────────────────────────────────────────
TOPGO_DIR = os.path.join(PROJECT_ROOT, "outputs", "go_enrichment", "tables")
FIG_MAIN  = os.path.join(PROJECT_ROOT, "outputs", "figures", "main")

IN_UP   = os.path.join(TOPGO_DIR, "rna_atac_main_up_topgo_bp_full.tsv")
IN_DOWN = os.path.join(TOPGO_DIR, "rna_atac_main_down_topgo_bp_full.tsv")
OUT_PDF = os.path.join(FIG_MAIN,  "Figure_4C_rna_atac_go.pdf")
OUT_PNG = os.path.join(FIG_MAIN,  "Figure_4C_rna_atac_go.png")

# ── Palette ────────────────────────────────────────────────────────────────────
PALETTE = {
    "Increased in patient": "#D55E00",
    "Decreased in patient": "#0072B2",
}

# ── Curated term selection ─────────────────────────────────────────────────────
# Hand-picked biologically relevant terms, ordered within each group from
# least to most significant (plotted bottom-to-top within group).
# Groups: UP terms (increased in patient) and DOWN terms (decreased in patient).

UP_TERMS = [
    # Neural development / guidance (least significant first)
    "nervous system development",
    "positive regulation of axonogenesis",
    "semaphorin-plexin signaling pathway",
    # Mesenchymal / skeletal
    "BMP signaling pathway",
    "collagen fibril organization",
    "chondrocyte development",
]

DOWN_TERMS = [
    # Synaptic / neural
    "excitatory chemical synaptic transmission",
    "action potential",
    # ECM / cell adhesion / signaling
    "stem cell proliferation",
    "positive regulation of cell junction assembly",
    "canonical Wnt signaling pathway",
    "extracellular matrix organization",
]

# ── Constants ──────────────────────────────────────────────────────────────────
BAR_ALPHA  = 0.28
BAR_HEIGHT = 0.55
FONT_TICK  = 10.0
FONT_ANNOT = 9.0
WRAP_WIDTH = 32     # chars before wrapping label to second line
P05        = -math.log10(0.05)
ROW_SPACING = 1.12


def dot_size(count):
    return 36 + count * 8


def representative_counts(counts):
    unique = sorted(set(int(c) for c in counts))
    if len(unique) <= 4:
        return unique
    idx = [0, len(unique) // 3, (2 * len(unique)) // 3, len(unique) - 1]
    values = [unique[i] for i in idx]
    deduped = []
    for value in values:
        if value not in deduped:
            deduped.append(value)
    return deduped


def load_topgo(path):
    """Return dict of Term -> row for all rows in a topGO full table."""
    table = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                pval_raw = (row.get("weight01Fisher_num") or
                            row.get("weight01Fisher", "1")).lstrip("<").strip()
                table[row["Term"]] = {
                    "term":  row["Term"],
                    "pval":  float(pval_raw),
                    "score": -math.log10(max(float(pval_raw), 1e-30)),
                    "count": int(row.get("Significant", 0)),
                }
            except (ValueError, KeyError):
                continue
    return table


def wrap_label(term, width=WRAP_WIDTH):
    """Wrap a long term label to at most 2 lines."""
    lines = textwrap.wrap(term, width=width, break_long_words=False)
    return "\n".join(lines[:2])


def plot_panel(ax, rows, title, color, x_max, show_threshold_label):
    y_pos = np.arange(len(rows)) * ROW_SPACING
    labels = [wrap_label(r["term"]) for r in rows]
    scores = [r["score"] for r in rows]
    counts = [r["count"] for r in rows]
    dot_sizes = [dot_size(c) for c in counts]

    ax.barh(y_pos, scores, color=color, alpha=BAR_ALPHA, height=BAR_HEIGHT, zorder=1)
    ax.scatter(
        scores,
        y_pos,
        color=color,
        s=dot_sizes,
        zorder=3,
        linewidths=0.8,
        edgecolors="#333333",
        alpha=1.0,
    )

    ax.axvline(P05, color="#999999", linewidth=0.9, linestyle="--", zorder=2)
    if show_threshold_label:
        ax.text(
            P05 + 0.04,
            -0.08,
            "p = 0.05",
            transform=ax.get_xaxis_transform(),
            va="top",
            ha="left",
            fontsize=FONT_ANNOT - 1,
            color="#666666",
            clip_on=False,
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=FONT_TICK, linespacing=0.9)
    ax.set_xlim(0, x_max)
    ax.set_ylim(-0.5, y_pos[-1] + 0.45)
    ax.tick_params(axis="x", labelsize=10)
    ax.tick_params(axis="y", pad=4)
    ax.set_title(title, fontsize=12, fontweight="bold", color=color, loc="left", pad=10)

    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)


def main():
    for path in (IN_UP, IN_DOWN):
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"Missing topGO result: {path}\n"
                "Run 08_topgo_rna_chip_go_rebuild.R first."
            )

    up_table   = load_topgo(IN_UP)
    down_table = load_topgo(IN_DOWN)

    # Build ordered rows: DOWN (bottom) then UP (top), each in list order
    rows = []
    for term in DOWN_TERMS:
        if term in down_table:
            rows.append({**down_table[term], "group": "Decreased in patient"})
        else:
            print(f"WARNING: DOWN term not found: '{term}'")

    n_down = len(rows)

    for term in UP_TERMS:
        if term in up_table:
            rows.append({**up_table[term], "group": "Increased in patient"})
        else:
            print(f"WARNING: UP term not found: '{term}'")

    down_rows = [r for r in rows if r["group"] == "Decreased in patient"]
    up_rows = [r for r in rows if r["group"] == "Increased in patient"]
    all_counts = [r["count"] for r in rows]
    all_scores = [r["score"] for r in rows]
    x_max = min(4.5, max(all_scores) * 1.30)

    fig_h = max(5.3, len(rows) * 0.58 + 1.4)
    fig = plt.figure(figsize=(7.3, fig_h))
    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[max(len(up_rows), 1), max(len(down_rows), 1)],
        hspace=0.26,
    )
    ax_up = fig.add_subplot(gs[0])
    ax_dn = fig.add_subplot(gs[1], sharex=ax_up)

    plot_panel(
        ax_up,
        up_rows,
        "Increased in patient",
        PALETTE["Increased in patient"],
        x_max,
        show_threshold_label=True,
    )
    plot_panel(
        ax_dn,
        down_rows,
        "Decreased in patient",
        PALETTE["Decreased in patient"],
        x_max,
        show_threshold_label=False,
    )
    ax_up.tick_params(axis="x", labelbottom=False)
    ax_dn.set_xlabel(r"$-\log_{10}$(adjusted p-value)", fontsize=11)

    size_handles = [
        Line2D(
            [0], [0], marker="o", color="none",
            markerfacecolor="#666666",
            markeredgecolor="#333333",
            markersize=math.sqrt(dot_size(c)),
            label=f"n = {c}",
        )
        for c in representative_counts(all_counts)
    ]
    fig.legend(
        handles=size_handles,
        title="Gene count",
        title_fontsize=9.5,
        fontsize=9.5,
        frameon=True,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.03),
        ncol=len(size_handles),
    )

    fig.subplots_adjust(left=0.50, right=0.92, bottom=0.16, top=0.95)

    os.makedirs(FIG_MAIN, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, dpi=300)
    fig.savefig(OUT_PNG, bbox_inches="tight", dpi=150)
    plt.close(fig)

    import warnings; warnings.filterwarnings("ignore")
    try:
        from pypdf import PdfReader
        p = PdfReader(OUT_PDF).pages[0]
        w = float(p.mediabox.width) / 72
        h = float(p.mediabox.height) / 72
    except ImportError:
        from PyPDF2 import PdfFileReader
        with open(OUT_PDF, "rb") as fh:
            p = PdfFileReader(fh).getPage(0)
            w = (float(p.mediaBox.getUpperRight_x()) - float(p.mediaBox.getLowerLeft_x())) / 72
            h = (float(p.mediaBox.getUpperRight_y()) - float(p.mediaBox.getLowerLeft_y())) / 72
    print(f"Saved: {OUT_PDF}")
    print(f"Size:  {w:.3f} x {h:.3f} in")


if __name__ == "__main__":
    main()
