#!/usr/bin/env python3
"""
Figure 1A — Study design schematic  (v2 — modern styling)
NRXN1alpha iPSC multi-omic profiling at Day 3 of neural induction

Outputs:
  outputs/figures/main/Figure_1A_study_design.pdf
  outputs/figures/main/Figure_1A_study_design.png
  outputs/figures/main/python_plots/Figure_1A_study_design.pdf
"""

import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse, FancyBboxPatch, FancyArrowPatch
from matplotlib.path import Path
import matplotlib.patheffects as pe
from _repo_paths import PROJECT_ROOT, OUTPUTS_DIR


FIG_MAIN = OUTPUTS_DIR / "figures" / "main"
FIG_PY = FIG_MAIN / "python_plots"
FIG_MAIN.mkdir(parents=True, exist_ok=True)
FIG_PY.mkdir(parents=True, exist_ok=True)

# ── palette ───────────────────────────────────────────────────────────────────
TEAL      = "#2E9EB0"
TEAL_LT   = "#A8D8E0"
PURPLE    = "#6B4A96"
PURPLE_LT = "#C8B8E0"
GRAY      = "#3A3A3A"
MIDGRAY   = "#888888"
LGRAY     = "#CCCCCC"
PANEL_BG  = "#EAF5F8"     # light blue for expression panel
BG        = "white"

# ── figure ────────────────────────────────────────────────────────────────────
fig_w, fig_h = 9.0, 9.8
fig = plt.figure(figsize=(fig_w, fig_h), facecolor=BG)
ax  = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
ax.set_facecolor(BG)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Expression panel with shaded background   (y 0.66 – 0.97)
# ═══════════════════════════════════════════════════════════════════════════════

PNL_X0, PNL_X1 = 0.12, 0.93
PNL_Y0, PNL_Y1 = 0.66, 0.97

panel_bg = FancyBboxPatch((PNL_X0, PNL_Y0), PNL_X1 - PNL_X0, PNL_Y1 - PNL_Y0,
                          boxstyle="round,pad=0.005",
                          facecolor=PANEL_BG, edgecolor=TEAL_LT, lw=1.0, zorder=1)
ax.add_patch(panel_bg)

# Curve domain (x fraction, within panel)
CX0, CX1 = PNL_X0 + 0.02, PNL_X1 - 0.02
cx = np.linspace(CX0, CX1, 600)
t  = (cx - CX0) / (CX1 - CX0)   # 0 → 1

# ── control curve: two clear peaks ──────────────────────────────────────────
# Peak 1 near t=0.20, trough at t=0.42, Peak 2 near t=0.72, gentle plateau
ctrl_baseline = 0.795
amp = 0.055
peak1 = amp * 1.40 * np.exp(-((t - 0.20)**2) / 0.008)
trough= -amp * 0.60 * np.exp(-((t - 0.42)**2) / 0.012)
peak2 = amp * 1.25 * np.exp(-((t - 0.72)**2) / 0.020)
ctrl_vals = ctrl_baseline + peak1 + trough + peak2

# ── patient curve: flat baseline ─────────────────────────────────────────────
pat_baseline = 0.710
pat_vals = np.full_like(cx, pat_baseline)

# Day 3 x-position: at the first peak
day3_x = CX0 + 0.20 * (CX1 - CX0)

# Fill under control curve (very faint)
ax.fill_between(cx, ctrl_baseline, ctrl_vals, color=TEAL, alpha=0.10, zorder=2)

# Plot curves
ax.plot(cx, ctrl_vals, color=TEAL,   lw=2.4, solid_capstyle="round", zorder=3)
ax.plot(cx, pat_vals,  color=PURPLE, lw=2.4, solid_capstyle="round", zorder=3)

# ── Bracket spanning the two peaks ──────────────────────────────────────────
# Peak positions
pk1_t, pk2_t = 0.20, 0.72
pk1_x = CX0 + pk1_t * (CX1 - CX0)
pk2_x = CX0 + pk2_t * (CX1 - CX0)
brk_top = 0.952

for px in [pk1_x, pk2_x]:
    ax.plot([px, px], [brk_top - 0.012, brk_top], color=TEAL, lw=1.3)
ax.plot([pk1_x, pk2_x], [brk_top, brk_top], color=TEAL, lw=1.3)
ax.text((pk1_x + pk2_x) / 2, brk_top + 0.006,
        "NRXN1α expression", ha="center", va="bottom",
        fontsize=8.0, color=TEAL, fontstyle="italic", fontweight="bold")

# ── Labels on the right of each curve ────────────────────────────────────────
ax.text(PNL_X1 - 0.005, ctrl_baseline + 0.004,
        "Control (NRXN1α +/+)", ha="right", va="bottom",
        fontsize=8.2, color=TEAL, fontweight="bold")
ax.text(PNL_X1 - 0.005, pat_baseline - 0.004,
        "Patient (NRXN1α −/−)", ha="right", va="top",
        fontsize=8.2, color=PURPLE, fontweight="bold")

# ── "NRXN1α expression" y-axis label ─────────────────────────────────────────
ax.text(PNL_X0 + 0.012, (ctrl_baseline + pat_baseline) / 2,
        "NRXN1α\nexpression", ha="center", va="center",
        fontsize=7.5, color=GRAY, rotation=90, linespacing=1.3,
        fontstyle="italic")

# ── Day 3 dashed line through expression panel ───────────────────────────────
ax.plot([day3_x, day3_x], [PNL_Y0, brk_top - 0.018],
        color=GRAY, lw=1.1, linestyle="--", zorder=2, alpha=0.55)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Timeline   (y ≈ 0.560 – 0.610)
# ═══════════════════════════════════════════════════════════════════════════════

TL_Y    = 0.590
TL_X0   = 0.08
TL_X1   = 0.93

# Main timeline arrow
ax.annotate("", xy=(TL_X1, TL_Y), xytext=(TL_X0, TL_Y),
            arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=1.8,
                            mutation_scale=12), zorder=4)

# Stage positions (fraction along timeline)
STAGE_XS = {
    "iPSC":                  TL_X0 + 0.02,
    "Day 3":                 day3_x,
    "Neural Stem\nCells":    TL_X0 + 0.55 * (TL_X1 - TL_X0),
    "Neural\nDifferentiation": TL_X0 + 0.88 * (TL_X1 - TL_X0),
}

for label, sx in STAGE_XS.items():
    if label == "Day 3":
        ax.plot([sx, sx], [TL_Y - 0.015, TL_Y + 0.015],
                color=GRAY, lw=2.0, zorder=5)
    else:
        ax.plot(sx, TL_Y, "o", color="white", ms=9, zorder=5,
                mec=GRAY, mew=1.8)
    y_off = -0.038 if "\n" in label else -0.028
    ax.text(sx, TL_Y + y_off, label, ha="center", va="top",
            fontsize=8.0, color=GRAY, linespacing=1.3)

# Day 3 label above timeline
ax.text(day3_x, TL_Y + 0.022, "Day 3", ha="center", va="bottom",
        fontsize=8.5, color=GRAY, fontweight="bold")

# Dashed continuation between NES and ND
nes_x = STAGE_XS["Neural Stem\nCells"]
nd_x  = STAGE_XS["Neural\nDifferentiation"]
ax.plot([nes_x + 0.025, nd_x - 0.025], [TL_Y, TL_Y],
        color=GRAY, lw=1.8, linestyle=":", zorder=3)

# "Neural Induction" span bracket below
bkt_y  = TL_Y - 0.065
bkt_x0 = STAGE_XS["iPSC"]
bkt_x1 = nes_x
for bx in [bkt_x0, bkt_x1]:
    ax.plot([bx, bx], [TL_Y - 0.005, bkt_y], color=LGRAY, lw=1.0)
ax.plot([bkt_x0, bkt_x1], [bkt_y, bkt_y], color=LGRAY, lw=1.0)
ax.text((bkt_x0 + bkt_x1) / 2, bkt_y - 0.007,
        "Neural Induction", ha="center", va="top",
        fontsize=7.0, color=MIDGRAY, fontstyle="italic")

# Downward arrow from Day 3 into sample section
ax.annotate("", xy=(day3_x, 0.490), xytext=(day3_x, TL_Y - 0.014),
            arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=1.4,
                            mutation_scale=10))


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Sample groups   (y ≈ 0.33 – 0.52)
# ═══════════════════════════════════════════════════════════════════════════════

BRK_Y   = 0.485
CTRL_CX = 0.210
PAT_CX  = 0.565

# Horizontal split bracket
for bx in [CTRL_CX, PAT_CX]:
    ax.plot([bx, bx], [BRK_Y, BRK_Y - 0.022], color=GRAY, lw=1.3)
ax.plot([CTRL_CX, PAT_CX], [BRK_Y, BRK_Y], color=GRAY, lw=1.3)

DISH_RX  = 0.038
DISH_RY  = 0.025
DISH_SEP = 0.105   # centre-to-centre, well clear of 2*DISH_RX=0.076
DISH_Y   = 0.408
LABEL_Y  = 0.474
SUBLAB_Y = 0.456

def draw_dish(cx, cy, color, n_dots=7):
    dish = Ellipse((cx, cy), width=DISH_RX*2, height=DISH_RY*2,
                   facecolor=color + "22", edgecolor=color,
                   lw=1.6, zorder=4)
    ax.add_patch(dish)
    # inner ring (lid outline)
    inner = Ellipse((cx, cy), width=DISH_RX*1.55, height=DISH_RY*1.55,
                    facecolor="none", edgecolor=color + "60", lw=0.7, zorder=4)
    ax.add_patch(inner)
    rng = np.random.default_rng(abs(hash((round(cx,4), round(cy,4)))) % (2**31))
    placed = 0
    for _ in range(200):
        if placed >= n_dots:
            break
        ddx = rng.uniform(-DISH_RX*0.65, DISH_RX*0.65)
        ddy = rng.uniform(-DISH_RY*0.65, DISH_RY*0.65)
        if (ddx/(DISH_RX*0.65))**2 + (ddy/(DISH_RY*0.65))**2 <= 1:
            ax.plot(cx + ddx, cy + ddy, "o",
                    color=color, ms=2.5, alpha=0.85, zorder=5)
            placed += 1

# Arrows from bracket to group labels
for bx in [CTRL_CX, PAT_CX]:
    ax.annotate("", xy=(bx, LABEL_Y + 0.015), xytext=(bx, BRK_Y - 0.023),
                arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=1.2,
                                mutation_scale=9))

# ── Control ───────────────────────────────────────────────────────────────────
ax.text(CTRL_CX, LABEL_Y,   "Control lines",           ha="center", va="center",
        fontsize=8.8, color=TEAL,  fontweight="bold")
ax.text(CTRL_CX, SUBLAB_Y,  "Ctrl-7, Ctrl-10, Ctrl-14", ha="center", va="center",
        fontsize=7.6, color=GRAY)

for i, (dlabel, ndots) in enumerate([("Low-density", 5), ("High-density", 15)]):
    dx = CTRL_CX - DISH_SEP/2 + i*DISH_SEP
    draw_dish(dx, DISH_Y, TEAL, n_dots=ndots)
    ax.text(dx, DISH_Y - DISH_RY - 0.015, dlabel, ha="center", va="top",
            fontsize=7.0, color=GRAY)

# ── Patient ───────────────────────────────────────────────────────────────────
ax.text(PAT_CX, LABEL_Y,   "Patient clones",  ha="center", va="center",
        fontsize=8.8, color=PURPLE, fontweight="bold")
ax.text(PAT_CX, SUBLAB_Y,  "PI, PII, PIII",   ha="center", va="center",
        fontsize=7.6, color=GRAY)

for i, (dlabel, ndots) in enumerate([("Low-density", 5), ("High-density", 15)]):
    dx = PAT_CX - DISH_SEP/2 + i*DISH_SEP
    draw_dish(dx, DISH_Y, PURPLE, n_dots=ndots)
    ax.text(dx, DISH_Y - DISH_RY - 0.015, dlabel, ha="center", va="top",
            fontsize=7.0, color=GRAY)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — Three assay boxes   (y ≈ 0.195 – 0.315)
# ═══════════════════════════════════════════════════════════════════════════════

ASSAY_CX  = (CTRL_CX + PAT_CX) / 2
ASSAY_Y   = 0.248
BOT_LABEL = DISH_Y - DISH_RY - 0.015 - 0.018   # bottom of density labels

# Connecting bracket from both groups into centre
brk2_y = BOT_LABEL - 0.018
for bx in [CTRL_CX, PAT_CX]:
    ax.plot([bx, bx], [BOT_LABEL, brk2_y], color=GRAY, lw=1.0)
ax.plot([CTRL_CX, PAT_CX], [brk2_y, brk2_y], color=GRAY, lw=1.0)

ax.annotate("", xy=(ASSAY_CX, ASSAY_Y + 0.086),
            xytext=(ASSAY_CX, brk2_y),
            arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=1.3,
                            mutation_scale=10))

assay_defs = [
    ("RNA-seq",           "Transcriptome",      "#2980B9", "#D6EAF8"),
    ("ATAC-seq",          "Chromatin\nopenness", "#1A9E5C", "#D5F5E3"),
    ("H3K27me3\nChIP-seq","Repressive\nepigenome","#7D3C98","#E8DAEF"),
]

BOX_W, BOX_H = 0.148, 0.080
BOX_SEP = 0.160
box_xs = [ASSAY_CX - BOX_SEP, ASSAY_CX, ASSAY_CX + BOX_SEP]

for (name, subtitle, col, bgcol), bx in zip(assay_defs, box_xs):
    r = FancyBboxPatch((bx - BOX_W/2, ASSAY_Y - BOX_H/2), BOX_W, BOX_H,
                       boxstyle="round,pad=0.010",
                       facecolor=bgcol, edgecolor=col, lw=1.5, zorder=3)
    ax.add_patch(r)
    n_lines = name.count("\n") + 1
    name_y = ASSAY_Y + 0.012 if n_lines == 1 else ASSAY_Y + 0.020
    ax.text(bx, name_y, name, ha="center", va="center",
            fontsize=8.2, color=col, fontweight="bold",
            linespacing=1.25, zorder=4)
    sub_y = ASSAY_Y - 0.018 if n_lines == 1 else ASSAY_Y - 0.020
    ax.text(bx, sub_y, subtitle, ha="center", va="center",
            fontsize=6.8, color=GRAY, linespacing=1.2, zorder=4)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Bioinformatics funnel   (y ≈ 0.02 – 0.175)
# ═══════════════════════════════════════════════════════════════════════════════

BIO_CX  = ASSAY_CX
BIO_TOP = ASSAY_Y - BOX_H/2 - 0.010

ax.annotate("", xy=(BIO_CX, 0.155), xytext=(BIO_CX, BIO_TOP),
            arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=1.3,
                            mutation_scale=10))

BW, BH = 0.310, 0.058
BIO_Y = 0.126
rb = FancyBboxPatch((BIO_CX - BW/2, BIO_Y - BH/2), BW, BH,
                    boxstyle="round,pad=0.010",
                    facecolor="#F4F4F4", edgecolor=MIDGRAY, lw=1.4, zorder=3)
ax.add_patch(rb)
ax.text(BIO_CX, BIO_Y, "Bioinformatics integration",
        ha="center", va="center",
        fontsize=8.8, color=GRAY, fontweight="bold", zorder=4)

ax.annotate("", xy=(BIO_CX, 0.065), xytext=(BIO_CX, BIO_Y - BH/2 - 0.004),
            arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=1.3,
                            mutation_scale=10))

OW, OH = 0.480, 0.052
OUT_Y = 0.040
ro = FancyBboxPatch((BIO_CX - OW/2, OUT_Y - OH/2), OW, OH,
                    boxstyle="round,pad=0.010",
                    facecolor="#FEF9E7", edgecolor="#D4AC0D", lw=1.5, zorder=3)
ax.add_patch(ro)
ax.text(BIO_CX, OUT_Y, "Dysregulated genes and signaling pathways",
        ha="center", va="center",
        fontsize=8.4, color="#7D6608", fontweight="bold", zorder=4)


# ═══════════════════════════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════════════════════════
out_pdf = FIG_MAIN / "Figure_1A_study_design.pdf"
out_png = FIG_MAIN / "Figure_1A_study_design.png"
out_py_pdf = FIG_PY / "Figure_1A_study_design.pdf"
fig.savefig(out_pdf, bbox_inches="tight", dpi=300, facecolor=BG)
fig.savefig(out_png, bbox_inches="tight", dpi=300, facecolor=BG)
fig.savefig(out_py_pdf, bbox_inches="tight", dpi=300, facecolor=BG)
plt.close(fig)
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
print(f"Saved: {out_py_pdf}")
