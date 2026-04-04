#!/usr/bin/env python3
"""
07n_plot_figure_2E_splicing_go.py

Generate the main Figure 2E GO Biological Process panel for genes with
significant splicing events using the splicing topGO output table.
"""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from _repo_paths import OUTPUTS_DIR


TOPGO_PATH = OUTPUTS_DIR / "splicing_go" / "tables" / "splicing_topgo_bp_full.tsv"
OUT_PDF = OUTPUTS_DIR / "figures" / "main" / "Figure_2E_splicing_go.pdf"
OUT_PNG = OUTPUTS_DIR / "figures" / "main" / "Figure_2E_splicing_go.png"

CURATED_TERMS = [
    ("Wnt signaling pathway", "Transcription / signaling"),
    ("neural precursor cell proliferation", "Transcription / signaling"),
    ("positive regulation of DNA-templated transcription", "Transcription / signaling"),
    ("RNA splicing", "Spliceosome / RNA splicing"),
    ("regulation of alternative mRNA splicing, via spliceosome", "Spliceosome / RNA splicing"),
    ("chromatin organization", "Chromatin organization"),
    ("chromatin remodeling", "Chromatin organization"),
]

PALETTE = {
    "Transcription / signaling": "#009E73",
    "Spliceosome / RNA splicing": "#E69F00",
    "Chromatin organization": "#CC79A7",
}

BAR_ALPHA = 0.28
BAR_HEIGHT = 0.55
ROW_SPACING = 1.12
P05 = -math.log10(0.05)
FIGURE_TITLE = "Genes with significant splicing events"
WRAP_WIDTH = 30


def dot_size(count: float) -> float:
    return 28 + count * 0.55


def load_rows(path: Path) -> list[dict[str, object]]:
    import csv

    rows_by_term: dict[str, dict[str, object]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                rows_by_term[row["Term"]] = {
                    "term": row["Term"],
                    "count": int(row["Significant"]),
                    "pvalue": float(row["weight01Fisher_num"]),
                    "score": float(row["enrichment_score"]),
                }
            except (KeyError, ValueError):
                continue

    rows: list[dict[str, object]] = []
    missing: list[str] = []
    for term, category in CURATED_TERMS:
        if term not in rows_by_term:
            missing.append(term)
            continue
        rows.append({**rows_by_term[term], "category": category})

    if missing:
        raise RuntimeError(
            "Missing expected splicing topGO terms: " + ", ".join(missing)
        )

    return rows


def wrap_label(text: str, width: int = WRAP_WIDTH) -> str:
    lines = textwrap.wrap(text, width=width, break_long_words=False)
    return "\n".join(lines[:2])


def representative_counts(counts: np.ndarray) -> list[int]:
    if counts.size == 0:
        return []

    max_count = float(np.max(counts))
    if max_count <= 25:
        step = 5
    elif max_count <= 100:
        step = 10
    elif max_count <= 250:
        step = 25
    else:
        step = 50

    raw_values = np.linspace(float(np.min(counts)), max_count, 4)
    rounded = [max(step, int(step * round(value / step))) for value in raw_values]

    deduped: list[int] = []
    for value in rounded:
        if value not in deduped:
            deduped.append(value)
    return deduped


def main() -> None:
    rows = load_rows(TOPGO_PATH)
    if not rows:
        raise RuntimeError("No curated terms found in splicing topGO table.")

    descriptions = [wrap_label(str(r["term"])) for r in rows]
    counts = np.array([float(r["count"]) for r in rows], dtype=float)
    scores = np.array([float(r["score"]) for r in rows], dtype=float)
    categories = [str(r["category"]) for r in rows]
    colors = [PALETTE[category] for category in categories]
    dot_sizes = [dot_size(count) for count in counts]
    y = np.arange(len(rows)) * ROW_SPACING

    plt.rcParams["font.family"] = "DejaVu Sans"
    fig, ax = plt.subplots(figsize=(7.8, 5.4))
    ax.barh(y, scores, color=colors, alpha=BAR_ALPHA, height=BAR_HEIGHT, zorder=1)
    ax.scatter(
        scores,
        y,
        s=dot_sizes,
        c=colors,
        edgecolors="#333333",
        linewidths=0.8,
        zorder=3,
        alpha=1.0,
    )
    ax.axvline(P05, color="#888", linewidth=0.8, linestyle=":", zorder=0)
    ax.text(
        P05 + 0.03,
        1.01,
        "p = 0.05",
        transform=ax.get_xaxis_transform(),
        ha="left",
        va="bottom",
        fontsize=10,
        color="#666666",
        clip_on=False,
    )
    for idx in range(1, len(rows)):
        if categories[idx] != categories[idx - 1]:
            ax.axhline((idx - 0.5) * ROW_SPACING, color="#bbbbbb", linewidth=0.8, linestyle="-", zorder=0)

    ax.set_yticks(y)
    ax.set_yticklabels(descriptions, fontsize=11)
    ax.set_xlabel(r"$-\log_{10}$ (weight01 Fisher)", fontsize=11)
    fig.suptitle(FIGURE_TITLE, fontsize=12.5, fontweight="bold", y=0.975)
    ax.tick_params(axis="x", labelsize=10)
    ax.tick_params(axis="y", length=0, pad=5)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, scores.max() * 1.18)
    ax.set_ylim(-0.6, y[-1] + 0.45)

    category_handles = [
        Line2D(
            [0], [0],
            marker="o",
            color="w",
            markerfacecolor=PALETTE[category],
            markeredgecolor="#333333",
            markeredgewidth=0.8,
            markersize=9,
            label=category,
        )
        for category in ["Transcription / signaling", "Spliceosome / RNA splicing", "Chromatin organization"]
    ]
    size_handles = [
        Line2D(
            [0], [0],
            marker="o",
            color="w",
            markerfacecolor="#666666",
            markeredgecolor="#333333",
            markeredgewidth=0.8,
            markersize=math.sqrt(dot_size(count)),
            label=f"n = {count}",
        )
        for count in representative_counts(counts)
    ]
    fig.legend(
        category_handles,
        [handle.get_label() for handle in category_handles],
        title="Category",
        title_fontsize=10,
        fontsize=9.5,
        frameon=True,
        loc="lower center",
        bbox_to_anchor=(0.31, 0.035),
        ncol=2,
    )
    fig.legend(
        size_handles,
        [handle.get_label() for handle in size_handles],
        title="Gene count",
        title_fontsize=10,
        fontsize=9.5,
        frameon=True,
        loc="lower center",
        bbox_to_anchor=(0.82, 0.005),
        ncol=1,
    )

    fig.subplots_adjust(left=0.44, right=0.95, top=0.87, bottom=0.30)
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_PDF}")
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
