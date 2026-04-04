#!/usr/bin/env python3
"""
07i_plot_figure_4B_rna_atac_scatter.py
Final manuscript Figure 4B - RNA vs ATAC promoter log2FC concordance scatter.
Uses the canonical Python version with concordance colors, without gene labels.

Input:
  outputs/integration/tables/rna_atac_main_overlap.tsv

Output:
  outputs/figures/main/Figure_4B_rna_atac_scatter.pdf
  outputs/figures/main/Figure_4B_rna_atac_scatter.png
"""

import os, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _repo_paths import PROJECT_ROOT, OUTPUTS_DIR
try:
    from scipy import stats as _scipy_stats
    def linregress(x, y):
        return _scipy_stats.linregress(x, y)
except ImportError:
    # Pure-numpy fallback
    def linregress(x, y):
        n = len(x)
        mx, my = x.mean(), y.mean()
        slope = ((x - mx) * (y - my)).sum() / ((x - mx)**2).sum()
        intercept = my - slope * mx
        y_hat = slope * x + intercept
        ss_res = ((y - y_hat)**2).sum()
        ss_tot = ((y - my)**2).sum()
        r_val = (((x - mx) * (y - my)).sum() /
                 (np.sqrt(((x - mx)**2).sum()) * np.sqrt(((y - my)**2).sum())))
        # p-value from t-distribution approximation
        t_stat = r_val * np.sqrt((n - 2) / (1 - r_val**2 + 1e-30))
        from math import erfc, sqrt
        p_val = 2 * (1 - 0.5 * erfc(-abs(t_stat) / sqrt(2)))
        return slope, intercept, r_val, p_val, None

ATAC_TSV = os.path.join(PROJECT_ROOT, "outputs", "integration", "tables", "rna_atac_main_overlap.tsv")
OUT_PDF  = os.path.join(PROJECT_ROOT, "outputs", "figures", "main", "Figure_4B_rna_atac_scatter.pdf")
OUT_PNG  = os.path.join(PROJECT_ROOT, "outputs", "figures", "main", "Figure_4B_rna_atac_scatter.png")
X_LIM = (-15, 15)

# Concordant: both down or both up
C_BOTH_DOWN = "#2980B9"   # blue  — both decreased
C_BOTH_UP   = "#C0392B"   # red   — both increased
C_DISCORD   = "#BBBBB0"   # gray  — discordant

def load_data(path):
    rows = []
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            try:
                rows.append({
                    "gene": row["gene_name"],
                    "rna":  float(row["rna_log2fc"]),
                    "atac": float(row["atac_log2fc"]),
                })
            except (ValueError, KeyError):
                pass
    return rows


def main():
    rows = load_data(ATAC_TSV)

    rna_vals  = np.array([r["rna"]  for r in rows])
    atac_vals = np.array([r["atac"] for r in rows])

    # Correlation
    slope, intercept, r_val, p_val, _ = linregress(rna_vals, atac_vals)
    # Classify quadrant
    def classify(r):
        if r["rna"] < 0 and r["atac"] < 0:
            return "both_down"
        elif r["rna"] > 0 and r["atac"] > 0:
            return "both_up"
        return "discord"

    for r in rows:
        r["cls"] = classify(r)

    both_down = [r for r in rows if r["cls"] == "both_down"]
    both_up   = [r for r in rows if r["cls"] == "both_up"]
    discord   = [r for r in rows if r["cls"] == "discord"]

    fig, ax = plt.subplots(figsize=(5.8, 5.4))

    ax.scatter([r["rna"] for r in discord],   [r["atac"] for r in discord],
               s=12, c=C_DISCORD, alpha=0.40, linewidths=0, zorder=1)
    ax.scatter([r["rna"] for r in both_down], [r["atac"] for r in both_down],
               s=14, c=C_BOTH_DOWN, alpha=0.60, linewidths=0, zorder=2)
    ax.scatter([r["rna"] for r in both_up],   [r["atac"] for r in both_up],
               s=14, c=C_BOTH_UP, alpha=0.60, linewidths=0, zorder=2)

    # Regression line
    x_line = np.array([X_LIM[0], X_LIM[1]])
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, color="#444", lw=1.2, zorder=3)

    # Correlation stats
    ax.text(0.97, 0.05,
            f"r = {r_val:.2f}\np = {p_val:.1e}\nn = {len(rows):,}",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.5, color="#333",
            bbox=dict(fc="white", ec="#ccc", pad=3))

    # Quadrant counts
    ax.text(0.03, 0.97, f"Both \u2193: {len(both_down):,}", transform=ax.transAxes,
            ha="left", va="top", fontsize=8, color=C_BOTH_DOWN, fontweight="bold")
    ax.text(0.97, 0.97, f"Both \u2191: {len(both_up):,}", transform=ax.transAxes,
            ha="right", va="top", fontsize=8, color=C_BOTH_UP, fontweight="bold")

    # Reference lines
    ax.axhline(0, color="#aaa", lw=0.6, ls="--")
    ax.axvline(0, color="#aaa", lw=0.6, ls="--")

    ax.set_xlabel("RNA log$_2$ fold change (patient / control)", fontsize=10)
    ax.set_ylabel("ATAC log$_2$ fold change (patient / control)", fontsize=10)
    ax.set_xlim(*X_LIM)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)

    fig.tight_layout(pad=0.9)
    fig.savefig(OUT_PDF, bbox_inches="tight", dpi=300)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"r={r_val:.3f}, p={p_val:.2e}, n={len(rows)}")
    print(f"Saved: {OUT_PDF}")
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
