#!/usr/bin/env python3
"""
plot_figure_S3_qpcr_validation.py

Generate Supplementary Figure S3 RT-qPCR validation panels.

Inputs:
  data/processed_inputs/qpcr/all_genes_expression.txt
  data/processed_inputs/qpcr/qPCR_stats_summary.txt

Outputs:
  outputs/figures/supp/Figure_S3A_SMAD7_qPCR.pdf
  outputs/figures/supp/Figure_S3B_Ecadherin_qPCR.pdf
"""

import csv
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _repo_paths import PROJECT_ROOT, OUTPUTS_DIR
EXPR_PATH = PROJECT_ROOT / "data" / "processed_inputs" / "qpcr" / "all_genes_expression.txt"
STATS_PATH = PROJECT_ROOT / "data" / "processed_inputs" / "qpcr" / "qPCR_stats_summary.txt"
OUT_DIR = OUTPUTS_DIR / "figures" / "supp"

GENE_CONFIG = {
    "SMAD7": {
        "display_name": "SMAD7",
        "output_name": "Figure_S3A_SMAD7_qPCR.pdf",
        "fold_change_label": "4.6\u00d7 \u2193",
        "p_label": "p = 0.09\u2020",
    },
    "E-cadherin": {
        "display_name": "E-cadherin",
        "output_name": "Figure_S3B_Ecadherin_qPCR.pdf",
        "fold_change_label": "1.6\u00d7 \u2193",
        "p_label": "p = 0.19\u2020",
    },
}

CTRL_KEYS = ["Ctrl_indiv1", "Ctrl_indiv2", "Ctrl_indiv3"]
DEL_KEYS = ["Del_clone1", "Del_clone2", "Del_clone3"]
CTRL_X = 0.0
DEL_X = 1.0
JITTER_OFFSETS = (-0.06, 0.0, 0.06)
CTRL_FACE = "#222222"
CTRL_EDGE = "#222222"
DEL_FACE = "white"
DEL_EDGE = "#222222"
SUMMARY_COLOR = "#444444"
FOOTNOTE = "\u2020 Exploratory; patient group = 3 clones from one individual"
GROUP_LABELS = ["Control", "NRXN1\u03b1-del"]


def style_matplotlib():
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
    plt.rcParams["font.size"] = 8
    plt.rcParams["axes.labelsize"] = 8
    plt.rcParams["axes.titlesize"] = 9
    plt.rcParams["xtick.labelsize"] = 8
    plt.rcParams["ytick.labelsize"] = 8


def load_expression(path):
    out = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row["Gene"]
            if gene not in GENE_CONFIG:
                continue
            out[gene] = {
                "control": [float(row[key]) for key in CTRL_KEYS],
                "deletion": [float(row[key]) for key in DEL_KEYS],
            }
    return out


def load_stats(path):
    header = None
    rows = {}
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("Gene\t"):
                header = line.split("\t")
                continue
            if header is None or line.startswith("LEGEND"):
                continue
            parts = line.split("\t")
            if len(parts) != len(header):
                continue
            row = dict(zip(header, parts))
            gene = row["Gene"]
            if gene in GENE_CONFIG:
                rows[gene] = row
    return rows


def compute_summary(values):
    mean_val = statistics.mean(values)
    sd_val = statistics.stdev(values)
    return mean_val, sd_val


def validate_against_stats(expression_data, stats_data):
    for gene, values in expression_data.items():
        stats_row = stats_data[gene]
        ctrl_mean, ctrl_sd = compute_summary(values["control"])
        del_mean, del_sd = compute_summary(values["deletion"])
        fc = del_mean / ctrl_mean

        if abs(ctrl_mean - float(stats_row["Ctrl_mean"])) > 0.11:
            raise ValueError(f"{gene}: control mean does not match summary table")
        if abs(ctrl_sd - float(stats_row["Ctrl_SD"])) > 0.11:
            raise ValueError(f"{gene}: control SD does not match summary table")
        if abs(del_mean - float(stats_row["Del_mean"])) > 0.11:
            raise ValueError(f"{gene}: deletion mean does not match summary table")
        if abs(del_sd - float(stats_row["Del_SD"])) > 0.11:
            raise ValueError(f"{gene}: deletion SD does not match summary table")
        if abs(fc - float(stats_row["FC_del_over_ctrl"])) > 0.001:
            raise ValueError(f"{gene}: fold change does not match summary table")


def draw_mean_sd(ax, x_pos, values):
    mean_val, sd_val = compute_summary(values)
    cap_half = 0.09
    mean_half = 0.13
    ax.vlines(x_pos, mean_val - sd_val, mean_val + sd_val, color=SUMMARY_COLOR, lw=1.3, zorder=2)
    ax.hlines([mean_val - sd_val, mean_val + sd_val], x_pos - cap_half, x_pos + cap_half,
              color=SUMMARY_COLOR, lw=1.1, zorder=2)
    ax.hlines(mean_val, x_pos - mean_half, x_pos + mean_half, color=SUMMARY_COLOR, lw=2.0, zorder=3)


def plot_gene(gene, values, config, output_path):
    ctrl_vals = values["control"]
    del_vals = values["deletion"]
    ctrl_mean, ctrl_sd = compute_summary(ctrl_vals)
    del_mean, del_sd = compute_summary(del_vals)

    fig, ax = plt.subplots(figsize=(3.0, 3.2))
    fig.subplots_adjust(left=0.20, right=0.97, top=0.85, bottom=0.28)

    for offset, value in zip(JITTER_OFFSETS, ctrl_vals):
        ax.scatter(
            CTRL_X + offset,
            value,
            s=34,
            facecolor=CTRL_FACE,
            edgecolor=CTRL_EDGE,
            linewidth=0.8,
            zorder=4,
        )

    for offset, value in zip(JITTER_OFFSETS, del_vals):
        ax.scatter(
            DEL_X + offset,
            value,
            s=34,
            facecolor=DEL_FACE,
            edgecolor=DEL_EDGE,
            linewidth=1.0,
            zorder=4,
        )

    draw_mean_sd(ax, CTRL_X, ctrl_vals)
    draw_mean_sd(ax, DEL_X, del_vals)

    ymax = max(max(ctrl_vals), max(del_vals), ctrl_mean + ctrl_sd, del_mean + del_sd) * 1.30
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(0, ymax)
    ax.set_xticks([CTRL_X, DEL_X])
    ax.set_xticklabels(GROUP_LABELS)
    ax.set_ylabel("Relative expression")
    ax.set_title(config["display_name"], fontweight="bold", pad=8)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(axis="x", length=0)

    annotation = f'{config["fold_change_label"]}\n{config["p_label"]}'
    ax.text(
        0.98,
        0.96,
        annotation,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        color="#222222",
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "none", "alpha": 0.9},
    )

    fig.text(0.5, 0.045, FOOTNOTE, ha="center", va="center", fontsize=6.7, color="#555555")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main():
    style_matplotlib()
    expression_data = load_expression(EXPR_PATH)
    stats_data = load_stats(STATS_PATH)
    validate_against_stats(expression_data, stats_data)

    for gene, config in GENE_CONFIG.items():
        out_path = OUT_DIR / config["output_name"]
        plot_gene(gene, expression_data[gene], config, out_path)
        print(f"Saved: {out_path}")

        ctrl_mean, ctrl_sd = compute_summary(expression_data[gene]["control"])
        del_mean, del_sd = compute_summary(expression_data[gene]["deletion"])
        print(
            f"{gene}: control mean={ctrl_mean:.1f}, control SD={ctrl_sd:.1f}, "
            f"deletion mean={del_mean:.1f}, deletion SD={del_sd:.1f}"
        )


if __name__ == "__main__":
    main()
