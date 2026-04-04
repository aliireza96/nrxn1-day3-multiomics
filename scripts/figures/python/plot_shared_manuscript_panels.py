#!/usr/bin/env python3
"""
plot_shared_manuscript_panels.py

Consolidated driver for the remaining Python manuscript panels that still share
one plotting entrypoint.
Saves to: outputs/figures/main/

Figures generated
-----------------
Figure_2C_rna_go_directional_categories.pdf
    Main RNA GO directional panel.

Figure_2D_splicing_summary_bar.pdf
    Main splicing summary panel.

Figure_3E_mechanistic_heatmap.pdf
    Curated mechanistic heatmap panel.

Standalone manuscript panels now live in their own dedicated scripts
(for example Figure 3B, Figure 4A-4E, and Figure 5B) and are intentionally
not generated from this shared file anymore.
"""

import argparse
import os
import csv
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.colorbar import ColorbarBase
from _repo_paths import PROJECT_ROOT

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
FIG_TABLES  = os.path.join(PROJECT_ROOT, "outputs", "figures", "tables")
FIG_OUT     = os.path.join(PROJECT_ROOT, "outputs", "figures", "main")
os.makedirs(FIG_OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# Shared palette (consistent: red = up/gained/more open, blue = down/lost/less open)
# ---------------------------------------------------------------------------
C_RED    = "#C0392B"
C_BLUE   = "#2980B9"
C_GRAY   = "#C8C8C8"

THEME_PALETTE = {
    "Spliceosome / snRNP":        "#D55E00",
    "Cell adhesion / ECM":        "#0072B2",
    "Signaling":                  "#CC79A7",
    "Neural / regionalization":   "#009E73",
    "Calcium / signaling":        "#E69F00",
    "Development":                "#999999",
}

def neg_log10(p_iter):
    return [-math.log10(max(pv, 1e-320)) for pv in p_iter]


# ===========================================================================
# Figure 2B  — RNA GO directional lollipop  (p=0.05 label position fix)
# ===========================================================================
UP_DATA = [
    ("embryonic digestive tract morphogenesis",         5,  1.1e-4, "Development"),
    ("Notch signaling pathway",                        16,  1.5e-5, "Signaling"),
    ("synaptic membrane adhesion",                      6,  2.4e-4, "Cell adhesion / ECM"),
    ("collagen fibril organization",                   11,  4.7e-6, "Cell adhesion / ECM"),
    ("homophilic cell adhesion via\nmembrane adhesion molecules",
                                                       17,  6.0e-6, "Cell adhesion / ECM"),
    ("spliceosomal tri-snRNP complex assembly",         6,  7.1e-5, "Spliceosome / snRNP"),
    ("formation of quadruple SL/U4/U5/U6 snRNP",       5,  2.0e-6, "Spliceosome / snRNP"),
]
DN_DATA = [
    ("skeletal system development",                    50,  2.5e-5, "Development"),
    ("positive regulation of MAPK cascade",            38,  5.3e-4, "Calcium / signaling"),
    ("cellular response to calcium ion",               13,  9.7e-5, "Calcium / signaling"),
    ("neural crest cell migration",                    10,  5.0e-4, "Neural / regionalization"),
    ("facial nerve structural organization",            4,  4.8e-4, "Neural / regionalization"),
    ("neuron remodeling",                               6,  3.4e-5, "Neural / regionalization"),
    ("rhombomere development",                          5,  4.1e-6, "Neural / regionalization"),
]

BAR_ALPHA  = 0.28
BAR_HEIGHT = 0.62
FONT_ANNOT = 12
TEXT_PAD   = 0.05
P05        = -math.log10(0.05)
COLOR_UP   = "#B85C00"
COLOR_DN   = "#1A5F8A"
SIZE_LEGEND_COUNTS = [5, 15, 30, 50]


def go_dot_size(n_genes):
    return 40 + n_genes * 6


def _go_panel(ax, data, direction, x_max):
    terms   = [d[0] for d in data]
    n_genes = [d[1] for d in data]
    cats    = [d[3] for d in data]
    colors  = [THEME_PALETTE[c] for c in cats]
    scores  = neg_log10([d[2] for d in data])
    dot_sizes = [go_dot_size(ng) for ng in n_genes]
    n       = len(data)
    y_pos   = np.arange(n)

    ax.barh(y_pos, scores, color=colors, alpha=BAR_ALPHA, height=BAR_HEIGHT, zorder=1)
    ax.scatter(scores, y_pos, s=dot_sizes, c=colors,
               zorder=3, linewidths=0.8, edgecolors="#333333", alpha=1.0)
    for i in range(1, n):
        if cats[i] != cats[i - 1]:
            ax.axhline(i - 0.5, color="#bbb", lw=0.7, ls="--", zorder=0)

    # Anchor the threshold label to the top plot edge rather than inside the data area.
    ax.axvline(P05, color="#888", lw=0.8, ls=":", zorder=0)
    ax.text(
        P05 + 0.04,
        1.01,
        "p\u202f=\u202f0.05",
        transform=ax.get_xaxis_transform(),
        ha="left",
        va="bottom",
        fontsize=11,
        color="#666666",
        clip_on=False,
    )

    dir_label = "\u25b2 Increased in patient" if direction == "up" else "\u25bc Decreased in patient"
    dir_color = COLOR_UP if direction == "up" else COLOR_DN
    ax.set_title(dir_label, fontsize=14, fontweight="bold", color=dir_color,
                 loc="left", pad=18)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(terms, fontsize=13)
    ax.set_xlim(0, x_max)
    ax.set_ylim(-0.6, n - 0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", labelsize=12)


def plot_fig2b():
    all_scores = neg_log10([d[2] for d in UP_DATA + DN_DATA])
    x_max = max(all_scores) * 1.28
    fig = plt.figure(figsize=(9.0, 8.4))   # taller to give room above each panel title
    gs  = gridspec.GridSpec(2, 1, height_ratios=[len(UP_DATA), len(DN_DATA)], hspace=0.60)
    ax_up = fig.add_subplot(gs[0])
    ax_dn = fig.add_subplot(gs[1])
    _go_panel(ax_up, UP_DATA, "up", x_max)
    _go_panel(ax_dn, DN_DATA, "dn", x_max)
    ax_dn.set_xlabel(r"$-\log_{10}$ (p-value, topGO weight01)", fontsize=13)

    all_cats = sorted(set(d[3] for d in UP_DATA + DN_DATA))
    category_handles = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=THEME_PALETTE[c], markersize=10,
               markeredgecolor="#333333", label=c)
        for c in all_cats
    ]
    size_handles = [
        Line2D(
            [0], [0],
            marker="o",
            color="w",
            markerfacecolor="#666666",
            markeredgecolor="white",
            markersize=math.sqrt(go_dot_size(n)),
            label=f"n = {n}",
        )
        for n in SIZE_LEGEND_COUNTS
    ]
    fig.subplots_adjust(bottom=0.31)
    fig.legend(
        handles=category_handles,
        title="Category",
        title_fontsize=12,
        fontsize=12,
        ncol=3,
        loc="lower center",
        bbox_to_anchor=(0.43, 0.015),
        frameon=True,
        edgecolor="#ccc",
    )
    fig.legend(
        handles=size_handles,
        title="Gene count",
        title_fontsize=12,
        fontsize=11,
        ncol=4,
        loc="lower center",
        bbox_to_anchor=(0.60, 0.135),
        frameon=True,
        edgecolor="#ccc",
    )
    out_pdf = os.path.join(FIG_OUT, "Figure_2C_rna_go_directional_categories.pdf")
    out_png = os.path.join(FIG_OUT, "Figure_2C_rna_go_directional_categories.png")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_pdf}")
    print(f"  Saved: {out_png}")


# ===========================================================================
# Figure 2D  — Splicing summary: diverging bar, red=increased, blue=decreased
# ===========================================================================
SIG_TSV = os.path.join(ROOT_DIR, "results", "splicing", "tables",
                        "significant_events_jcec.tsv")
ET_ORDER = ["SE", "RI", "A3SS", "A5SS", "MXE"]
ET_LABELS_2D = {
    "SE":   "Skipped exon (SE)",
    "RI":   "Retained intron (RI)",
    "A3SS": "Alt. 3\u2032 splice site (A3SS)",
    "A5SS": "Alt. 5\u2032 splice site (A5SS)",
    "MXE":  "Mutually exclusive exons (MXE)",
}


def _load_splicing(path):
    counts = {et: {"inc": 0, "dec": 0} for et in ET_ORDER}
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            et = row.get("event_type", "").strip()
            if et not in counts:
                continue
            try:
                ild = float(row["IncLevelDifference"])
                if ild > 0:
                    counts[et]["inc"] += 1
                else:
                    counts[et]["dec"] += 1
            except (ValueError, KeyError):
                pass
    return counts


def plot_fig2d():
    counts = _load_splicing(SIG_TSV)
    et_sorted = sorted(ET_ORDER, key=lambda e: counts[e]["inc"] + counts[e]["dec"],
                       reverse=True)
    n        = len(et_sorted)
    y_pos    = np.arange(n)
    inc_vals = [counts[et]["inc"] for et in et_sorted]
    dec_vals = [counts[et]["dec"] for et in et_sorted]
    totals   = [counts[et]["inc"] + counts[et]["dec"] for et in et_sorted]
    labels   = [ET_LABELS_2D[et] for et in et_sorted]
    max_side = max(max(inc_vals), max(dec_vals))
    x_right  =  max_side * 1.32
    x_left   = -max_side * 0.72

    # Build figure with extra bottom space for legend
    fig, ax = plt.subplots(figsize=(8.8, 5.2))   # taller to clear title
    fig.subplots_adjust(bottom=0.25, top=0.82)

    bar_h = 0.44
    ax.barh(y_pos, inc_vals, height=bar_h, color=C_RED,  alpha=0.85)
    ax.barh(y_pos, [-v for v in dec_vals], height=bar_h, color=C_BLUE, alpha=0.85)

    # Count labels at bar ends
    for i, iv in enumerate(inc_vals):
        ax.text(iv + max_side * 0.012, i, f"{iv:,}",
                va="center", ha="left", fontsize=12, color="#333")
    for i, dv in enumerate(dec_vals):
        ax.text(-dv - max_side * 0.012, i, f"{dv:,}",
                va="center", ha="right", fontsize=12, color="#333")
    # Total count — placed ABOVE each bar pair (at y = i + bar_h*0.72) to avoid
    # horizontal collision with the per-direction count labels
    for i, tot in enumerate(totals):
        ax.text(0, i + bar_h * 0.72, f"n = {tot:,}",
                va="bottom", ha="center", fontsize=11, color="#666", style="italic")

    ax.set_xlim(x_left, x_right)
    ax.set_ylim(-0.6, n - 0.1)   # slight top margin for "n = N" labels above top bar
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=13)
    ax.axvline(0, color="#999", lw=0.8)

    max_tick = int(round(max_side / 500)) * 500
    x_ticks  = list(range(-max_tick, max_tick + 1, 500))
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([str(abs(v)) for v in x_ticks], fontsize=12)
    ax.set_xlabel("Number of significant splicing events (FDR < 0.05)", fontsize=13)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend placed below axes, outside the plot area
    legend_handles = [
        Patch(color=C_RED,  label="Increased inclusion in patient"),
        Patch(color=C_BLUE, label="Decreased inclusion in patient"),
    ]
    fig.legend(handles=legend_handles, fontsize=12, ncol=2,
               loc="lower center", bbox_to_anchor=(0.5, 0.02),
               frameon=True, edgecolor="#ccc")

    out = os.path.join(FIG_OUT, "Figure_2D_splicing_summary_bar.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


# ===========================================================================
# Figure 3C mechanistic heatmap
# — identical RdBu_r scale; EMT split into Epithelial + Mesenchymal
# ===========================================================================
MECH_TSV = os.path.join(FIG_TABLES, "Figure_3C_mechanistic_heatmap_source.tsv")

# Module order as requested: neural fate → TGF-β/BMP suppressors → epithelial → mesenchymal
MODULE_ORDER = [
    "Neural fate markers",
    "TGF-\u03b2 / BMP suppressors",
    "Epithelial markers",
    "Mesenchymal markers",
]
MODULE_LABELS = {
    "Neural fate markers": "Neural fate\nmarkers",
    "TGF-\u03b2 / BMP suppressors": "TGF-\u03b2 / BMP\nsuppressors",
    "Epithelial markers": "Epithelial\nmarkers",
    "Mesenchymal markers": "Mesenchymal\nmarkers",
}

# EMT genes split by direction (epithelial = down, mesenchymal = up)
EMT_EPITHELIAL  = {"CLDN7", "CDH1"}
EMT_MESENCHYMAL = {"TWIST2", "SNAI2", "CDH2"}
MECH_EXCLUDE = {"TBXT", "TDGF1"}

# Shared colormap: identical for both assay columns
MECH_CMAP    = "RdBu_r"
MECH_VMIN    = -5.5
MECH_VCENTER = 0.0
MECH_VMAX    = 5.5

COL_LABELS_MECH = {
    # Keep headers short so the compact panel stays legible in a 3-panel row.
    "rna_log2fc":  "RNA",
    "chip_log2fc": "H3K27\nme3",
}
COL_ORDER_MECH = ["rna_log2fc", "chip_log2fc"]


def _load_mech(path):
    rows = []
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            gene = row["gene_name"]
            if gene in MECH_EXCLUDE:
                continue

            chip_raw = row.get("chip_log2fc", "NA")
            chip_val = float(chip_raw) if chip_raw not in ("", "NA", "nan") else None

            # Reassign module based on EMT split
            orig_mod = row["module"]
            if orig_mod == "EMT / Mesenchymal":
                if gene in EMT_EPITHELIAL:
                    module = "Epithelial markers"
                else:
                    module = "Mesenchymal markers"
            elif orig_mod == "TGF-beta / SMAD":
                module = "TGF-\u03b2 / BMP suppressors"
            else:
                module = "Neural fate markers"

            rows.append({
                "gene":   gene,
                "module": module,
                "tier":   row["evidence_tier"],
                "rna":    float(row["rna_log2fc"]),
                "chip":   chip_val,
            })
    return rows


def plot_fig3c_mech():
    rows = _load_mech(MECH_TSV)

    # Order genes by module (as in MODULE_ORDER), then by RNA within each module.
    ordered = []
    module_spans = {}
    for mod in MODULE_ORDER:
        mod_rows = [r for r in rows if r["module"] == mod]
        mod_rows.sort(key=lambda r: r["rna"])
        if not mod_rows:
            continue
        start_y = len(ordered)
        ordered.extend(mod_rows)
        end_y = len(ordered) - 1
        module_spans[mod] = [start_y, end_y]

    n_genes = len(ordered)
    n_cols  = len(COL_ORDER_MECH)

    # Layout — compact, but keep enough room for readable horizontal module labels.
    cell_w    = 0.52   # widen heatmap body so the assembled panel reads more cleanly
    cell_h    = 0.26   # shorter rows
    label_w   = 0.50   # gene name column (genes are ≤6 chars)
    mod_lbl_w = 0.92   # extra room for manuscript-style module labels
    top_h     = 0.46   # column-header labels + top padding
    cbar_h    = 0.22   # colorbar height in inches
    bot_h     = cbar_h + 0.24   # colorbar + bottom margin

    total_w = label_w + n_cols * cell_w + mod_lbl_w + 0.16
    total_h = n_genes * cell_h + top_h + bot_h

    fig = plt.figure(figsize=(total_w, total_h))

    left_m   = label_w / total_w
    hm_w_f   = (n_cols * cell_w) / total_w
    hm_h_f   = (n_genes * cell_h) / total_h
    bot_m    = bot_h / total_h

    ax_hm = fig.add_axes([left_m, bot_m, hm_w_f, hm_h_f])

    cmap_obj = plt.get_cmap(MECH_CMAP)
    norm     = mcolors.TwoSlopeNorm(vmin=MECH_VMIN, vcenter=MECH_VCENTER, vmax=MECH_VMAX)
    na_color = "#D5D8DC"

    for ci, col in enumerate(COL_ORDER_MECH):
        vals = [r["rna"] if col == "rna_log2fc" else r["chip"] for r in ordered]
        for ri, v in enumerate(vals):
            fc = na_color if v is None else cmap_obj(norm(v))
            ax_hm.add_patch(Rectangle((ci - 0.5, ri - 0.5), 1, 1,
                                      facecolor=fc, edgecolor="white", lw=0.5))

    ax_hm.set_xlim(-0.5, n_cols - 0.5)
    ax_hm.set_ylim(-0.5, n_genes - 0.5)
    ax_hm.invert_yaxis()
    ax_hm.set_xticks(range(n_cols))
    ax_hm.set_xticklabels([COL_LABELS_MECH[c] for c in COL_ORDER_MECH],
                           fontsize=8.5, fontweight="bold")
    ax_hm.xaxis.set_ticks_position("top")
    ax_hm.xaxis.set_label_position("top")
    ax_hm.tick_params(axis="x", pad=2)
    ax_hm.set_yticks(range(n_genes))
    ax_hm.set_yticklabels([r["gene"] for r in ordered], fontsize=9, style="italic")
    ax_hm.tick_params(axis="both", length=0)

    # Module label x: left edge inside the label column to the right of the heatmap.
    label_x = left_m + hm_w_f + 0.05 / total_w

    for mod in MODULE_ORDER:
        if mod not in module_spans:
            continue
        start, end = module_spans[mod]
        if start > 0:
            ax_hm.axhline(start - 0.5, color="#222222", lw=1.6, zorder=10)

        s_top = bot_m + hm_h_f - (start * cell_h) / total_h
        s_ht  = ((end - start + 1) * cell_h) / total_h
        s_center = s_top - s_ht / 2.0

        # Wrapped horizontal label aligned to the grouped rows.
        fig.text(label_x, s_center, MODULE_LABELS[mod],
                 va="center", ha="left", fontsize=7.5, color="#222222",
                 fontweight="bold", linespacing=0.95,
                 transform=fig.transFigure)

    # Shared colorbar (same scale for both columns)
    cbar_bot = 0.12 / total_h   # small bottom margin below colorbar
    cbar_ht  = cbar_h / total_h
    ax_cb = fig.add_axes([left_m, cbar_bot, hm_w_f, cbar_ht])
    cb = ColorbarBase(ax_cb, cmap=cmap_obj, norm=norm, orientation="horizontal")
    cb.ax.tick_params(labelsize=7.5)
    cb.set_label("log\u2082FC", fontsize=7.5)
    cb.set_ticks([-5, 0, 5])

    # NA indicator next to colorbar
    na_ax = fig.add_axes([left_m + hm_w_f + 0.02,
                          cbar_bot + cbar_ht * 0.15,
                          0.025, cbar_ht * 0.7])
    na_ax.add_patch(Rectangle((0, 0), 1, 1, facecolor=na_color, edgecolor="#aaa", lw=0.5))
    na_ax.axis("off")
    fig.text(left_m + hm_w_f + 0.055, cbar_bot + cbar_ht * 0.5,
             "NA", va="center", ha="left", fontsize=7.5, color="#555",
             transform=fig.transFigure)

    out = os.path.join(FIG_OUT, "Figure_3E_mechanistic_heatmap.pdf")
    out_copy = os.path.join(FIG_OUT, "python_plots", "Figure_3E_mechanistic_heatmap.pdf")
    os.makedirs(os.path.dirname(out_copy), exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out_copy, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    print(f"  Saved: {out_copy}")


# ===========================================================================
# Figure selection / main
# ===========================================================================
def _outputs_fig2b():
    return [
        os.path.join(FIG_OUT, "Figure_2C_rna_go_directional_categories.pdf"),
        os.path.join(FIG_OUT, "Figure_2C_rna_go_directional_categories.png"),
    ]


def _outputs_fig2d():
    return [os.path.join(FIG_OUT, "Figure_2D_splicing_summary_bar.pdf")]


def _outputs_fig3e():
    return [
        os.path.join(FIG_OUT, "Figure_3E_mechanistic_heatmap.pdf"),
        os.path.join(FIG_OUT, "python_plots", "Figure_3E_mechanistic_heatmap.pdf"),
    ]


FIGURE_SPECS = [
    {
        "key": "2B",
        "title": "Figure 2B — RNA GO directional (p=0.05 label fix)",
        "func": plot_fig2b,
        "outputs": _outputs_fig2b,
    },
    {
        "key": "2D",
        "title": "Figure 2D — Splicing diverging bar (title + legend fix)",
        "func": plot_fig2d,
        "outputs": _outputs_fig2d,
    },
    {
        "key": "3E",
        "title": "Figure 3E — Mechanistic heatmap",
        "func": plot_fig3c_mech,
        "outputs": _outputs_fig3e,
    },
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate Python-owned manuscript figures. By default, "
            "existing outputs are treated as locked and skipped."
        )
    )
    parser.add_argument(
        "--figures",
        nargs="+",
        choices=[spec["key"] for spec in FIGURE_SPECS],
        help="Only generate the specified figure panels (for example: --figures 2B 3E).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerate selected figures even if their outputs already exist.",
    )
    return parser.parse_args()


def run_selected_figures(selected_keys=None, force=False):
    print("Generating Python-owned manuscript figures ...")
    print()

    selected = set(selected_keys) if selected_keys else None
    specs = [spec for spec in FIGURE_SPECS if selected is None or spec["key"] in selected]

    for idx, spec in enumerate(specs, start=1):
        print(f"[{idx}/{len(specs)}] {spec['title']}")
        missing = [path for path in spec["outputs"]() if not os.path.exists(path)]
        if not force and not missing:
            print("  Skipped: outputs already exist (use --force to regenerate)")
            continue
        spec["func"]()

    print()
    print("Done. Outputs in:", FIG_OUT)


if __name__ == "__main__":
    args = parse_args()
    run_selected_figures(selected_keys=args.figures, force=args.force)
