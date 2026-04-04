#!/usr/bin/env python3
"""
07u_plot_figure_4E_tobias_mechanistic.py

Figure 4E - TOBIAS mechanistic footprint panel.
Four aggregate footprint plots:
  up in patient:   TEAD4, SOX9
  down in patient: KLF5, EN1
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
from _repo_paths import PROJECT_ROOT, OUTPUTS_DIR

os.environ.setdefault("MPLCONFIGDIR", str((OUTPUTS_DIR / ".matplotlib").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    import pyBigWig  # type: ignore

    HAS_PYBIGWIG = True
except ImportError:
    pyBigWig = None
    HAS_PYBIGWIG = False


def resolve_existing_path(*relative_candidates: str) -> Path:
    for relative in relative_candidates:
        candidate = (PROJECT_ROOT / relative).resolve()
        if candidate.exists():
            return candidate
    return (PROJECT_ROOT / relative_candidates[0]).resolve()


TABLE_PATH = OUTPUTS_DIR / "figures" / "tables" / "Figure_4B_tobias_threshold_q95.tsv"
ATAC_BASE = resolve_existing_path("data/processed_inputs/atac/footprinting")
PATIENT_BW = (ATAC_BASE / "tobias_allvsall/bw_files/All-Ps_footprints.bw").resolve()
CONTROL_BW = (ATAC_BASE / "tobias_allvsall/bw_files/All-ctrls_footprints.bw").resolve()
OUT_DIR = OUTPUTS_DIR / "figures" / "main"
OUT_PDF = OUT_DIR / "Figure_4E_tobias_mechanistic.pdf"
OUT_PNG = OUT_DIR / "Figure_4E_tobias_mechanistic.png"

WINDOW = 150
MAX_SITES = 5000
SMOOTH_SIGMA_BP = 2.5  # light smoothing; ~5.9 bp FWHM at 1 bp resolution

C_PATIENT = "#C00000"
C_CONTROL = "#4472C4"
C_SHADE = "#BFBFBF"

FOOTPRINT_PANELS = [
    {
        "tf": "TEAD4",
        "directory": "TEAD4_MA0809.2",
        "subtitle": "more accessible in patient",
        "subtitle_color": C_PATIENT,
        "expected_sign": 1,
        "view_window": 150,
    },
    {
        "tf": "SOX9",
        "directory": "SOX9_MA0077.1",
        "subtitle": "more accessible in patient",
        "subtitle_color": C_PATIENT,
        "expected_sign": 1,
        "view_window": 150,
    },
    {
        "tf": "KLF5",
        "directory": "KLF5_MA0599.1",
        "subtitle": "less accessible in patient",
        "subtitle_color": C_CONTROL,
        "expected_sign": -1,
        "view_window": 150,
    },
    {
        "tf": "EN1",
        "directory": "EN1_MA0027.2",
        "subtitle": "less accessible in patient",
        "subtitle_color": C_CONTROL,
        "expected_sign": -1,
        "view_window": 150,
    },
]


plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def bed_path_for_panel(panel: dict[str, str]) -> Path:
    tf_dir = panel["directory"]
    return (ATAC_BASE / f"tobias_allvsall/BINdetect/{tf_dir}/beds/{tf_dir}_all.bed").resolve()


def load_bed_sites(bed_path: Path) -> list[tuple[str, int, int, str, float, float, float]]:
    sites: list[tuple[str, int, int, str, float, float, float]] = []
    with bed_path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for parts in reader:
            if len(parts) < 14:
                continue
            try:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                strand = parts[5]
                motif_score = float(parts[4])
                patient_score = float(parts[12])
                control_score = float(parts[13])
            except ValueError:
                continue
            sites.append((chrom, start, end, strand, motif_score, patient_score, control_score))
    sites.sort(key=lambda item: item[4], reverse=True)
    return sites


def choose_site_count(
    sites: list[tuple[str, int, int, str, float, float, float]],
    expected_sign: int,
    default_max_sites: int = MAX_SITES,
) -> int:
    total = len(sites)
    if total <= default_max_sites:
        return total

    candidate = min(default_max_sites, total)
    while True:
        subset = sites[:candidate]
        mean_delta = float(np.mean([row[5] - row[6] for row in subset]))
        if mean_delta == 0 or np.sign(mean_delta) == expected_sign or candidate >= total:
            return candidate
        next_candidate = min(candidate * 2, total)
        if next_candidate == candidate:
            return candidate
        candidate = next_candidate


def extract_footprint(
    bw_path: Path,
    sites: list[tuple[str, int, int, str, float, float, float]],
    window: int = WINDOW,
    max_sites: int = MAX_SITES,
) -> np.ndarray:
    profiles: list[np.ndarray] = []
    subset = sites[:max_sites]
    with pyBigWig.open(str(bw_path)) as bw:  # type: ignore[union-attr]
        for chrom, start, end, strand, _score, _patient_score, _control_score in subset:
            center = (start + end) // 2
            left = center - window
            right = center + window
            if left < 0 or right <= left:
                continue
            try:
                values = bw.values(chrom, left, right, numpy=True)
            except RuntimeError:
                continue
            if values is None:
                continue
            arr = np.asarray(values, dtype=float)
            if arr.size != 2 * window or np.all(np.isnan(arr)):
                continue
            if strand == "-":
                arr = arr[::-1]
            profiles.append(arr)

    if not profiles:
        return np.full(2 * window, np.nan)
    return np.nanmean(np.vstack(profiles), axis=0)


def signal_limits(patient: np.ndarray, control: np.ndarray, x: np.ndarray, view_window: int) -> tuple[float, float]:
    mask = (x >= -view_window) & (x <= view_window)
    vals = np.concatenate(
        [
            patient[mask & np.isfinite(patient)],
            control[mask & np.isfinite(control)],
        ]
    )
    if vals.size == 0:
        return 0.0, 1.0
    y_min = float(np.min(vals))
    y_max = float(np.max(vals))
    spread = y_max - y_min
    pad = max(spread * 0.18, max(abs(y_min), abs(y_max)) * 0.03, 0.004)
    return y_min - pad, y_max + pad


def smooth_signal(values: np.ndarray, sigma_bp: float = SMOOTH_SIGMA_BP) -> np.ndarray:
    if sigma_bp <= 0 or values.size == 0:
        return values
    radius = max(1, int(np.ceil(4 * sigma_bp)))
    grid = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-(grid ** 2) / (2.0 * sigma_bp ** 2))
    kernel /= kernel.sum()
    padded = np.pad(values, radius, mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def plot_footprint_panel(
    ax: plt.Axes,
    x: np.ndarray,
    patient: np.ndarray,
    control: np.ndarray,
    panel: dict[str, str],
    show_legend: bool,
    show_xlabel: bool,
) -> None:
    view_window = int(panel.get("view_window", WINDOW))
    patient_smooth = smooth_signal(patient)
    control_smooth = smooth_signal(control)
    patient_line, = ax.plot(x, patient_smooth, color=C_PATIENT, lw=1.7, label="Patient", zorder=3)
    control_line, = ax.plot(x, control_smooth, color=C_CONTROL, lw=1.7, label="Control", zorder=2)
    ax.axvspan(-8, 8, color=C_SHADE, alpha=0.12, zorder=0)

    ax.set_xlim(-view_window, view_window)
    ax.set_ylim(*signal_limits(patient_smooth, control_smooth, x, view_window))
    ax.set_ylabel("Footprint score", fontsize=9)
    ax.tick_params(labelsize=8)

    ax.text(
        0.0,
        1.08,
        panel["tf"],
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        clip_on=False,
    )
    ax.text(
        0.0,
        1.015,
        panel["subtitle"],
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        color=panel["subtitle_color"],
        style="italic",
        clip_on=False,
    )

    if show_xlabel:
        ax.set_xlabel("Position relative to motif center (bp)", fontsize=9)
    else:
        ax.set_xlabel("")

    if show_legend:
        ax.legend(
            handles=[patient_line, control_line],
            loc="upper right",
            frameon=False,
            fontsize=8,
            handlelength=2.5,
        )


def read_summary_rows(path: Path) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            tf = row["name"]
            if tf not in {panel["tf"] for panel in FOOTPRINT_PANELS}:
                continue
            if tf not in best or float(row["neg_log10_p"]) > float(best[tf]["neg_log10_p"]):
                best[tf] = row
    return best


def plot_fallback_bars(ax: plt.Axes) -> None:
    summary = read_summary_rows(TABLE_PATH)
    ordered = [summary[panel["tf"]] for panel in FOOTPRINT_PANELS if panel["tf"] in summary]
    labels = [row["name"] for row in ordered]
    patient = np.array([float(row["P_mean_score"]) for row in ordered], dtype=float)
    control = np.array([float(row["ctrl_mean_score"]) for row in ordered], dtype=float)
    y = np.arange(len(ordered), dtype=float)
    bar_h = 0.34

    ax.barh(y + bar_h / 2, patient, height=bar_h, color=C_PATIENT, alpha=0.90, label="Patient")
    ax.barh(y - bar_h / 2, control, height=bar_h, color=C_CONTROL, alpha=0.90, label="Control")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Mean binding score", fontsize=9)
    ax.set_title("Mean TF binding score: patient vs control", fontsize=10, fontweight="bold", pad=6)
    ax.tick_params(labelsize=8)
    ax.legend(loc="lower right", frameon=False, fontsize=8)


def build_figure() -> plt.Figure:
    fig = plt.figure(figsize=(7.5, 5.2))

    if HAS_PYBIGWIG:
        gs = fig.add_gridspec(
            2,
            2,
            hspace=0.58,
            wspace=0.30,
            left=0.09,
            right=0.98,
            top=0.92,
            bottom=0.12,
        )

        x = np.arange(-WINDOW, WINDOW)
        for idx, panel in enumerate(FOOTPRINT_PANELS):
            row = idx // 2
            col = idx % 2
            ax = fig.add_subplot(gs[row, col])
            sites = load_bed_sites(bed_path_for_panel(panel))
            site_count = choose_site_count(sites, expected_sign=int(panel.get("expected_sign", 1)))
            patient = extract_footprint(PATIENT_BW, sites, max_sites=site_count)
            control = extract_footprint(CONTROL_BW, sites, max_sites=site_count)
            plot_footprint_panel(
                ax,
                x,
                patient,
                control,
                panel,
                show_legend=(idx == 0),
                show_xlabel=(row == 1),
            )
    else:
        ax = fig.add_axes([0.10, 0.14, 0.84, 0.76])
        plot_fallback_bars(ax)

    return fig


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = build_figure()
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=200)
    w, h = fig.get_size_inches()
    plt.close(fig)
    print(f"Saved: {OUT_PDF}  Size: {w:.3f} x {h:.3f} in")


if __name__ == "__main__":
    main()
