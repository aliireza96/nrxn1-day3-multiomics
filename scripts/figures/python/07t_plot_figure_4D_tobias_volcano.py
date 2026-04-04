#!/usr/bin/env python3
"""
07t_plot_figure_4D_tobias_volcano.py
Figure 4D — TOBIAS TF footprinting volcano (q95 threshold, curated labels).

Labels chosen based on:
  - Lam et al. 2019 (same patient, neural fate shift toward radial glia/astroglia)
  - Project context: BMP/TGF-β axis, EMT module (Fig 3E), action potential deficits

More accessible in patient (YAP/glial/aberrant programs):
  TEAD2, TEAD4, MEF2C, SOX9, SOX13, ONECUT1, MEIS1

Less accessible in patient (lost neural identity/function):
  KLF5, SP4, ZEB1, SP3, KLF4, Pou5f1::Sox2, EN1

Input:
  outputs/figures/tables/Figure_4B_tobias_threshold_q95.tsv

Output:
  outputs/figures/main/Figure_4D_tobias_volcano.pdf
"""

import os
import csv
import warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from _repo_paths import PROJECT_ROOT

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
FIG_TABLES = os.path.join(PROJECT_ROOT, "outputs", "figures", "tables")
FIG_MAIN   = os.path.join(PROJECT_ROOT, "outputs", "figures", "main")

IN_TSV  = os.path.join(FIG_TABLES, "Figure_4B_tobias_threshold_q95.tsv")
OUT_PDF = os.path.join(FIG_MAIN,   "Figure_4D_tobias_volcano.pdf")
OUT_PNG = os.path.join(FIG_MAIN,   "Figure_4D_tobias_volcano.png")

# ── Curated labels ─────────────────────────────────────────────────────────────
LABEL_MORE = {
    "TEAD2",      # Hippo/YAP target; radial glia & mesenchymal fate
    "TEAD4",      # Hippo/YAP; highest binding change among TEADs
    "MEF2C",      # Excitatory/inhibitory neuron maturation; altered subtype
    "SOX9",       # Master gliogenesis/radial-glia TF; Lam: more glia in patient
    "SOX13",      # Aberrant non-neural SOX activity
    "ONECUT1",    # Fate-switching TF inappropriate in neural progenitors
    "MEIS1",      # Brain regionalization; altered developmental program
}

LABEL_LESS = {
    "KLF5",           # Top significance; progenitor maintenance
    "SP4",            # Neuronal gene expression; links to AP deficits (Lam et al.)
    "ZEB1",           # EMT master regulator; connects to Fig 3E EMT module
    "SP3",            # Same family as SP4; neuronal program
    "KLF4",           # Progenitor/stemness maintenance
    "Pou5f1::Sox2",   # Neural stem identity composite; Lam: altered NES identity
    "EN1",            # Midbrain/hindbrain homeobox; lost neural patterning
}

# ── Colours ────────────────────────────────────────────────────────────────────
COL_MORE = "#d1495b"   # red (more accessible)
COL_LESS = "#2f6bff"   # blue (less accessible)
COL_NS   = "#c7c7c7"   # grey (NS)

# ── Display name overrides (cleaner axis labels) ───────────────────────────────
DISPLAY_NAME = {
    "Pou5f1::Sox2": "POU5F1::SOX2",
}


def main():
    rows = []
    with open(IN_TSV, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append({
                "name":    r["name"],
                "x":       float(r["binding_change"]),
                "y":       float(r["neg_log10_p"]),
                "status":  r["sig_status"],
            })

    # Assign colours and label flags
    for r in rows:
        if r["status"] == "More accessible in patient":
            r["color"] = COL_MORE
            r["alpha"] = 0.72
            r["size"]  = 8
        elif r["status"] == "Less accessible in patient":
            r["color"] = COL_LESS
            r["alpha"] = 0.72
            r["size"]  = 8
        else:
            r["color"] = COL_NS
            r["alpha"] = 0.45
            r["size"]  = 6

        r["label"] = (
            (r["name"] in LABEL_MORE and r["status"] == "More accessible in patient") or
            (r["name"] in LABEL_LESS and r["status"] == "Less accessible in patient")
        )

    xs = np.array([r["x"] for r in rows])
    ys = np.array([r["y"] for r in rows])

    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    # NS points first (background)
    ns  = [r for r in rows if r["status"] == "NS"]
    ax.scatter([r["x"] for r in ns], [r["y"] for r in ns],
               c=COL_NS, s=6, alpha=0.40, linewidths=0, zorder=1)

    # Significant points
    for status, col in [("Less accessible in patient", COL_LESS),
                        ("More accessible in patient", COL_MORE)]:
        pts = [r for r in rows if r["status"] == status]
        ax.scatter([r["x"] for r in pts], [r["y"] for r in pts],
                   c=col, s=8, alpha=0.72, linewidths=0, zorder=2)

    # Vertical zero line
    ax.axvline(0, color="#aaaaaa", linewidth=0.5, zorder=0)

    # ── Labels with manual repel ──────────────────────────────────────────────
    # Deduplicate labeled list: if a TF name appears more than once, keep highest y
    seen = {}
    for r in rows:
        if r["label"]:
            if r["name"] not in seen or r["y"] > seen[r["name"]]["y"]:
                seen[r["name"]] = r
    labeled = list(seen.values())

    # Manual per-label text positions (xytext) to avoid overlapping.
    # Coordinates in data units: (xytext_x, xytext_y).
    MANUAL_XY = {
        # More accessible (right side)
        "TEAD2":   (0.282, 183),
        "TEAD4":   (0.282, 174),
        "MEF2C":   (0.103, 181),
        "ONECUT1": (0.098, 162),
        "SOX13":   (0.210, 142),
        "SOX9":    (0.213, 128),
        "MEIS1":   (0.096, 127),
        # Less accessible (left side)
        "KLF5":        (-0.255, 195),
        "SP4":         (-0.148, 193),
        "ZEB1":        (-0.278, 183),
        "SP3":         (-0.268, 170),
        "Pou5f1::Sox2": (-0.278, 157),
        "KLF4":        (-0.138, 165),
        "EN1":         (-0.288, 138),
    }

    # Try to use adjustText if available, otherwise use manual positions
    try:
        from adjustText import adjust_text
        texts = []
        for r in labeled:
            col = COL_MORE if r["status"] == "More accessible in patient" else COL_LESS
            display = DISPLAY_NAME.get(r["name"], r["name"])
            t = ax.text(r["x"], r["y"], display,
                        fontsize=7.5, color=col, fontweight="bold",
                        ha="center", va="bottom",
                        bbox=dict(boxstyle="round,pad=0.12", fc="white",
                                  ec=col, lw=0.5, alpha=0.88))
            texts.append(t)
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
                    expand_points=(1.3, 1.5), expand_text=(1.2, 1.3),
                    force_text=(0.4, 0.6), force_points=(0.2, 0.4),
                    only_move={"points": "y", "texts": "xy"})
    except ImportError:
        # Manual curated-position labeling with leader lines
        for r in labeled:
            col = COL_MORE if r["status"] == "More accessible in patient" else COL_LESS
            display = DISPLAY_NAME.get(r["name"], r["name"])
            if r["name"] in MANUAL_XY:
                tx, ty = MANUAL_XY[r["name"]]
            else:
                dx = 0.010 if r["x"] > 0 else -0.010
                tx, ty = r["x"] + dx, r["y"] + 4
            ha = "left" if tx >= r["x"] else "right"
            ax.annotate(
                display,
                xy=(r["x"], r["y"]),
                xytext=(tx, ty),
                fontsize=7.5, color=col, fontweight="bold",
                ha=ha, va="center",
                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5,
                                shrinkA=0, shrinkB=3),
                bbox=dict(boxstyle="round,pad=0.12", fc="white",
                          ec=col, lw=0.5, alpha=0.88),
            )

    # Axes labels and formatting
    ax.set_xlabel("Bound-score change (patient \u2212 control)", fontsize=11)
    ax.set_ylabel("\u2212log\u2081\u2080(p-value)", fontsize=11)
    ax.tick_params(labelsize=9.5)

    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)

    # Legend
    from matplotlib.lines import Line2D
    leg_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_LESS,
               markersize=6, label="Less accessible in patient"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_MORE,
               markersize=6, label="More accessible in patient"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_NS,
               markersize=6, label="NS"),
    ]
    ax.legend(handles=leg_handles, fontsize=8.5, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)

    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.20, top=0.97)

    os.makedirs(FIG_MAIN, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, dpi=300)
    fig.savefig(OUT_PNG, bbox_inches="tight", dpi=150)
    plt.close(fig)

    from pypdf import PdfReader
    p = PdfReader(OUT_PDF).pages[0]
    w = float(p.mediabox.width) / 72
    h = float(p.mediabox.height) / 72
    print(f"Saved: {OUT_PDF}")
    print(f"Size:  {w:.3f} x {h:.3f} in")


if __name__ == "__main__":
    main()
