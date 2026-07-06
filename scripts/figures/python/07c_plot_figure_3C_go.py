#!/usr/bin/env python3
"""
07c_plot_figure_3C_go.py
Redesigned Figure 3C (left panel) — H3K27me3 GO enrichment (clusterProfiler)
Genes with gained H3K27me3 AND decreased expression in NRXN1alpha-null cells.

Input table (written by 07_plot_figures.R):
  outputs/figures/tables/Figure_3C_chip_go_bp.tsv

Output:
  outputs/figures/main/Figure_3C_chip_go.pdf
  outputs/figures/main/Figure_3C_chip_go.png


Design notes:
  - Top 12 terms by p.adjust from clusterProfiler output
  - Two thematic groups: Synaptic/neuronal (orange) and Mesenchymal/neural crest (green)
  - Colors avoid blue/red (reserved for volcano up/down)
  - Gene counts encoded as dot size with a compact size legend
  - Group separator line; p=0.05 threshold dotted line
  - Category legend below axes
"""

import os
import csv
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from _repo_paths import PROJECT_ROOT, OUTPUTS_DIR

# ── Paths ──────────────────────────────────────────────────────────────────────
FIG_TABLES = os.path.join(PROJECT_ROOT, "outputs", "figures", "tables")
FIG_MAIN   = os.path.join(PROJECT_ROOT, "outputs", "figures", "main")

IN_TSV   = os.path.join(FIG_TABLES, "Figure_3C_chip_go_bp.tsv")
OUT_PDF  = os.path.join(FIG_MAIN,   "Figure_3C_chip_go.pdf")
OUT_PNG  = os.path.join(FIG_MAIN,   "Figure_3C_chip_go.png")

# ── Thematic palette ───────────────────────────────────────────────────────────
# Both colors match the thematic palette used in Figure 2B
PALETTE = {
    "Synaptic / neuronal":            "#D55E00",
    "Neural crest / developmental":   "#009E73",
}

# ── Group assignment ───────────────────────────────────────────────────────────
# Terms are assigned to thematic groups based on biological content.
# Terms not listed here are classified as "Synaptic / neuronal" by default
# (the majority of terms reflect synaptic/neuronal biology).
# NOTE: this panel avoids the label "mesenchymal": the developmental terms here
# are driven by neural-crest / early-lineage genes (SOX10, RET, semaphorins,
# FGF19, LEF1) that GO also tags as "mesenchyme" because the neural crest is an
# embryonic source of ectomesenchyme. The up-regulated mesenchymal EMT effectors
# (TWIST2, SNAI2, CDH2) are a distinct, oppositely-directed program shown in
# Figure 3D, so "mesenchymal" is reserved for that panel to avoid confusion.
DEVELOPMENTAL_TERMS = {
    "neural crest cell development",
    "neural crest cell differentiation",
    "neural crest cell migration",
    "stem cell differentiation",
    "stem cell development",
    "mesenchymal cell differentiation",
    "mesenchyme development",
    "fibroblast growth factor receptor signaling pathway",
    "extracellular matrix organization",
    "extracellular structure organization",
}

# Relabel a redundant representative to the most biologically accurate term
# among the terms that collapse into it (its genes are canonical neural crest).
RELABEL = {
    "mesenchymal cell differentiation": "neural crest cell development",
    "mesenchyme development":           "neural crest cell development",
}

# Gene-overlap (Jaccard) threshold for collapsing semantically redundant GO
# terms driven by overlapping gene sets (REVIGO / clusterProfiler::simplify
# style). Terms sharing >= this fraction of genes with an already-kept, more
# significant term are treated as redundant and dropped.
JACCARD_CUTOFF = 0.5
PADJ_CUTOFF    = 0.05

# ── Constants ──────────────────────────────────────────────────────────────────
TOP_N      = 12
BAR_ALPHA  = 0.28
BAR_HEIGHT = 0.55
P05        = -math.log10(0.05)
FIGURE_TITLE = "Genes with increased H3K27me3 and decreased expression"
ROW_SPACING = 1.18


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


def _jaccard(a, b):
    a, b = set(a), set(b)
    union = a | b
    return len(a & b) / len(union) if union else 0.0


def read_chip_go_tsv(path, top_n=12):
    """Read clusterProfiler GO table and return a non-redundant set of terms.

    Steps: (1) keep terms with p.adjust < PADJ_CUTOFF; (2) sort by significance;
    (3) greedily drop any term whose driver-gene set overlaps an already-kept,
    more significant term by >= JACCARD_CUTOFF (semantic redundancy reduction);
    (4) relabel collapsed representatives to their most accurate term; (5) return
    the top_n most significant survivors.
    """
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                padj = float(row["p.adjust"])
                if padj >= PADJ_CUTOFF:
                    continue
                rows.append({
                    "term":  row["Description"],
                    "padj":  padj,
                    "count": int(row["Count"]),
                    "score": -math.log10(padj),
                    "genes": frozenset(str(row["geneID"]).split("/")),
                })
            except (ValueError, KeyError):
                continue
    rows.sort(key=lambda r: r["padj"])

    kept = []
    for r in rows:
        if any(_jaccard(r["genes"], k["genes"]) >= JACCARD_CUTOFF for k in kept):
            continue
        r["term"] = RELABEL.get(r["term"].lower(), r["term"])
        kept.append(r)
    return kept[:top_n]


def assign_group(term):
    """Assign thematic group based on term text."""
    term_lower = term.lower()
    for dev_term in DEVELOPMENTAL_TERMS:
        if dev_term in term_lower:
            return "Neural crest / developmental"
    return "Synaptic / neuronal"


def main():
    rows = read_chip_go_tsv(IN_TSV, top_n=TOP_N)

    # Assign groups and sort: Synaptic first (top of chart), Mesenchymal below
    for r in rows:
        r["group"] = assign_group(r["term"])

    # Sort within each group by score descending, then stack groups:
    # Mesenchymal at bottom (low y), Synaptic at top (high y) for visual grouping
    synaptic   = sorted([r for r in rows if r["group"] == "Synaptic / neuronal"],
                        key=lambda r: r["score"])
    mesenchymal = sorted([r for r in rows if r["group"] == "Neural crest / developmental"],
                         key=lambda r: r["score"])

    # Plot order: developmental at bottom, separator, synaptic at top
    ordered = mesenchymal + synaptic
    n = len(ordered)

    terms   = [r["term"] for r in ordered]
    scores  = [r["score"] for r in ordered]
    counts  = [r["count"] for r in ordered]
    colors  = [PALETTE[r["group"]] for r in ordered]
    dot_sizes = [dot_size(c) for c in counts]
    y_pos   = np.arange(n) * ROW_SPACING

    x_max = max(scores) * 1.25

    fig, ax = plt.subplots(figsize=(7.8, 6.1))

    # bars
    ax.barh(y_pos, scores, color=colors, alpha=BAR_ALPHA,
            height=BAR_HEIGHT, zorder=1)
    # dots
    ax.scatter(scores, y_pos, s=dot_sizes, c=colors,
               zorder=3, linewidths=0.8, edgecolors="#333333", alpha=1.0)

    # group separator line between mesenchymal and synaptic blocks
    sep_idx = len(mesenchymal)
    if 0 < sep_idx < n:
        ax.axhline((sep_idx - 0.5) * ROW_SPACING, color="#bbb", linewidth=0.8, linestyle="-", zorder=0)

    # p=0.05 threshold
    ax.axvline(P05, color="#888", linewidth=0.8, linestyle=":", zorder=0)
    ax.text(
        P05 + 0.03,
        1.01,
        "p\u202f=\u202f0.05",
        transform=ax.get_xaxis_transform(),
        ha="left",
        va="bottom",
        fontsize=10,
        color="#666666",
        clip_on=False,
    )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(terms, fontsize=11.5)
    ax.set_xlim(0, x_max)
    ax.set_ylim(-0.6, y_pos[-1] + 0.4)
    ax.set_xlabel(r"$-\log_{10}$ (adjusted p-value)", fontsize=11)
    fig.suptitle(FIGURE_TITLE, fontsize=12.5, fontweight="bold", y=0.975)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", labelsize=10)
    ax.tick_params(axis="y", pad=5)

    # category legend below axes
    legend_handles = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=PALETTE[g], markersize=10,
               markeredgecolor="#333333", label=g)
        for g in PALETTE
    ]
    size_handles = [
        Line2D(
            [0], [0], marker="o", color="w",
            markerfacecolor="#666666",
            markeredgecolor="#333333",
            markersize=math.sqrt(dot_size(c)),
            label=f"n = {c}",
        )
        for c in representative_counts(counts)
    ]
    fig.subplots_adjust(left=0.49, right=0.98, bottom=0.33, top=0.87)
    fig.legend(
        handles=legend_handles,
        title="Category",
        title_fontsize=10,
        fontsize=9.5,
        ncol=1,
        loc="lower left",
        bbox_to_anchor=(0.49, 0.045),
        frameon=True,
        edgecolor="#ccc",
    )
    fig.legend(
        handles=size_handles,
        title="Gene count",
        title_fontsize=10,
        fontsize=9.5,
        ncol=1,
        loc="lower right",
        bbox_to_anchor=(0.98, 0.045),
        frameon=True,
        edgecolor="#ccc",
    )

    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_PDF}")
    print(f"Saved: {OUT_PNG}")
    print()
    print("NOTE: Regenerate the assembled panel after this script:")
    print("  cd manuscript_reanalysis/assembly_work")
    print("  pdflatex -interaction=nonstopmode Figure_3C_assembly.tex")
    print("  cp Figure_3C_assembly.pdf ../outputs/figures/main/Figure_3C_panel.pdf")


if __name__ == "__main__":
    main()
