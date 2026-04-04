#!/usr/bin/env python3
"""
plot_figure_5B_integration_heatmap.py

Generate the locked manuscript Figure 5B integration heatmap from the curated
display-order table and the triple-overlap matrix.

Inputs:
  data/processed_inputs/integration/updated_heatmap_curated_table.csv
  outputs/integration/tables/triple_overlap_main_heatmap_matrix.tsv

Outputs:
  outputs/figures/main/Figure_5B_integration_heatmap.pdf
  outputs/figures/main/Figure_5B_integration_heatmap.png
  outputs/figures/tables/Figure_5B_integration_heatmap_source_data.tsv
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from _repo_paths import PROJECT_ROOT, OUTPUTS_DIR


CURATED_PATH = (
    PROJECT_ROOT
    / "data"
    / "processed_inputs"
    / "integration"
    / "updated_heatmap_curated_table.csv"
)
MATRIX_PATH = OUTPUTS_DIR / "integration" / "tables" / "triple_overlap_main_heatmap_matrix.tsv"
OUT_PDF = OUTPUTS_DIR / "figures" / "main" / "Figure_5B_integration_heatmap.pdf"
OUT_PNG = OUTPUTS_DIR / "figures" / "main" / "Figure_5B_integration_heatmap.png"
OUT_TABLE = OUTPUTS_DIR / "figures" / "tables" / "Figure_5B_integration_heatmap_source_data.tsv"

ROW_ORDER = [
    ("rna_log2fc", "RNA-seq"),
    ("atac_log2fc", "ATAC-seq"),
    ("chip_log2fc", "H3K27me3"),
]

DISPLAY_LABEL = {
    "Synaptic": "Synaptic",
    "Neural identity": "Neural identity",
    "Epithelial": "Epithelial",
    "Developmental signaling": "Developmental\nsignaling",
    "Mesenchymal": "Mesenchymal",
}

EXPECTED_CATEGORY_ORDER = [
    "Synaptic",
    "Neural identity",
    "Epithelial",
    "Developmental signaling",
    "Mesenchymal",
]


def load_curated_rows(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("status") != "keep":
                continue
            rows.append(row)
    rows.sort(key=lambda row: (int(row["category_order"]), row["gene_name"]))
    return rows


def load_matrix(path: Path) -> dict[str, dict[str, float]]:
    rows: dict[str, dict[str, float]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["gene_name"]] = {
                "rna_log2fc": float(row["rna_log2fc"]),
                "atac_log2fc": float(row["atac_log2fc"]),
                "chip_log2fc": float(row["chip_log2fc"]),
            }
    return rows


def signed_scale_row(values: np.ndarray) -> np.ndarray:
    scaled = np.zeros_like(values, dtype=float)
    neg_mask = values < 0
    pos_mask = values > 0

    if np.any(neg_mask):
        neg_min = np.min(values[neg_mask])
        if neg_min < 0:
            scaled[neg_mask] = values[neg_mask] / abs(neg_min)
    if np.any(pos_mask):
        pos_max = np.max(values[pos_mask])
        if pos_max > 0:
            scaled[pos_mask] = values[pos_mask] / pos_max
    return scaled


def build_order(curated_rows: list[dict[str, str]], matrix: dict[str, dict[str, float]]):
    ordered_rows = [row for row in curated_rows if row["gene_name"] in matrix]
    categories = [row["visual_category"] for row in ordered_rows]
    ordered_categories = []
    for category in categories:
        if category not in ordered_categories:
            ordered_categories.append(category)

    if ordered_categories != EXPECTED_CATEGORY_ORDER:
        raise RuntimeError(
            f"Unexpected category order: {ordered_categories}. "
            f"Expected {EXPECTED_CATEGORY_ORDER}."
        )

    spans: list[tuple[str, int, int]] = []
    cursor = 0
    for category in EXPECTED_CATEGORY_ORDER:
        genes = [row["gene_name"] for row in ordered_rows if row["visual_category"] == category]
        if not genes:
            continue
        start = cursor
        cursor += len(genes)
        spans.append((category, start, cursor - 1))

    return ordered_rows, spans


def write_source_table(rows: list[dict[str, str]], matrix: dict[str, dict[str, float]], scaled: np.ndarray) -> None:
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TABLE.open("w", newline="") as handle:
        fieldnames = [
            "gene_order",
            "category_order",
            "visual_category",
            "gene_name",
            "status",
            "reason",
            "rna_log2fc",
            "rna_scaled",
            "atac_log2fc",
            "atac_scaled",
            "chip_log2fc",
            "chip_scaled",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for idx, row in enumerate(rows):
            gene = row["gene_name"]
            writer.writerow(
                {
                    "gene_order": idx + 1,
                    "category_order": row["category_order"],
                    "visual_category": row["visual_category"],
                    "gene_name": gene,
                    "status": row["status"],
                    "reason": row["reason"],
                    "rna_log2fc": matrix[gene]["rna_log2fc"],
                    "rna_scaled": scaled[0, idx],
                    "atac_log2fc": matrix[gene]["atac_log2fc"],
                    "atac_scaled": scaled[1, idx],
                    "chip_log2fc": matrix[gene]["chip_log2fc"],
                    "chip_scaled": scaled[2, idx],
                }
            )


def main() -> None:
    curated_rows = load_curated_rows(CURATED_PATH)
    matrix_rows = load_matrix(MATRIX_PATH)
    ordered_rows, spans = build_order(curated_rows, matrix_rows)
    gene_order = [row["gene_name"] for row in ordered_rows]

    if len(gene_order) != 20:
        raise RuntimeError(f"Expected 20 curated genes, found {len(gene_order)}.")

    raw_matrix = np.array(
        [[matrix_rows[gene][key] for gene in gene_order] for key, _ in ROW_ORDER],
        dtype=float,
    )
    scaled_matrix = np.vstack([signed_scale_row(raw_matrix[row_idx, :]) for row_idx in range(raw_matrix.shape[0])])
    write_source_table(ordered_rows, matrix_rows, scaled_matrix)

    cmap = LinearSegmentedColormap.from_list("signed_fc", ["#2f6bff", "#ffffff", "#d1495b"])

    fig = plt.figure(figsize=(7.4, 3.2))
    label_ax = fig.add_axes([0.08, 0.70, 0.78, 0.09])
    heat_ax = fig.add_axes([0.08, 0.32, 0.78, 0.36])
    cbar_ax = fig.add_axes([0.88, 0.32, 0.025, 0.36])

    im = heat_ax.imshow(
        scaled_matrix,
        aspect="auto",
        interpolation="none",
        cmap=cmap,
        vmin=-1.0,
        vmax=1.0,
    )

    heat_ax.set_yticks(np.arange(len(ROW_ORDER)))
    heat_ax.set_yticklabels([label for _, label in ROW_ORDER], fontsize=9)
    heat_ax.set_xticks(np.arange(len(gene_order)))
    heat_ax.set_xticklabels(gene_order, rotation=58, ha="right", rotation_mode="anchor", fontsize=7)
    heat_ax.tick_params(axis="x", length=0, pad=2)
    heat_ax.tick_params(axis="y", length=0)
    for spine in heat_ax.spines.values():
        spine.set_visible(False)

    label_ax.set_xlim(-0.5, len(gene_order) - 0.5)
    label_ax.set_ylim(0, 1)
    label_ax.axis("off")

    for category, start, end in spans:
        if end < len(gene_order) - 1:
            heat_ax.axvline(end + 0.5, color="black", linewidth=2.2)
            label_ax.axvline(end + 0.5, color="black", linewidth=2.2, ymin=0.10, ymax=0.82)
        label_ax.text(
            start + (end - start) / 2.0,
            0.22,
            DISPLAY_LABEL.get(category, category),
            ha="center",
            va="bottom",
            fontsize=7.2,
            color="black",
            fontweight="normal",
            linespacing=0.95,
        )

    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label("Normalized log2FC\n(per-assay signed scaling)", fontsize=8)

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {OUT_PDF}")
    print(f"Saved: {OUT_PNG}")
    print(f"Saved: {OUT_TABLE}")


if __name__ == "__main__":
    main()
