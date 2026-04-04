#!/usr/bin/env python3
"""
Assemble main and supplementary figures into complete PDFs, preserving source
panels as vector PDFs where possible.

Dependencies: pypdf, reportlab

Inputs:
  outputs/figures/main/*.pdf
  outputs/figures/supp/*.pdf

Outputs:
  outputs/figures/main/Figure_<n>_complete.pdf
  outputs/figures/supp/Figure_S<n>_complete.pdf
"""

import io
from copy import copy
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

from pypdf import PdfReader, PdfWriter, Transformation
from pypdf._page import PageObject
from pypdf.generic import RectangleObject
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent.parent
FIG_MAIN = ROOT_DIR / "outputs" / "figures" / "main"
FIG_SUPP = ROOT_DIR / "outputs" / "figures" / "supp"
FIG_SUPP_RNA = FIG_SUPP / "female_only_rna"
FIG_SUPP_MULTI = FIG_SUPP / "female_only_multiomic"
FIG_SUPP_QC = FIG_SUPP / "qc_panels"


def _page_dims_in(path, crop=None):
    with path.open("rb") as fh:
        reader = PdfReader(fh)
        page = reader.pages[0]
        llx = float(page.mediabox.left)
        lly = float(page.mediabox.bottom)
        urx = float(page.mediabox.right)
        ury = float(page.mediabox.top)

    crop = crop or {}
    llx += 72 * crop.get("left", 0.0)
    lly += 72 * crop.get("bottom", 0.0)
    urx -= 72 * crop.get("right", 0.0)
    ury -= 72 * crop.get("top", 0.0)
    return (urx - llx) / 72.0, (ury - lly) / 72.0


def _load_cropped_page(path, crop=None):
    data = path.read_bytes()
    reader = PdfReader(io.BytesIO(data))
    page = copy(reader.pages[0])

    llx = float(page.mediabox.left)
    lly = float(page.mediabox.bottom)
    urx = float(page.mediabox.right)
    ury = float(page.mediabox.top)

    crop = crop or {}
    llx += 72 * crop.get("left", 0.0)
    lly += 72 * crop.get("bottom", 0.0)
    urx -= 72 * crop.get("right", 0.0)
    ury -= 72 * crop.get("top", 0.0)

    page.mediabox = RectangleObject((llx, lly, urx, ury))
    page.cropbox = RectangleObject((llx, lly, urx, ury))
    return page, (urx - llx) / 72.0, (ury - lly) / 72.0


def _fit_within(src_w_in, src_h_in, max_w_in, max_h_in):
    scale = min(max_w_in / src_w_in, max_h_in / src_h_in)
    return src_w_in * scale, src_h_in * scale


def _make_overlay(page_w_pt, page_h_pt, masks, labels):
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=(page_w_pt, page_h_pt))
    c.setFillColorRGB(1, 1, 1)
    for mask in masks:
        c.rect(mask["x"] * 72.0, mask["y"] * 72.0,
               mask["w"] * 72.0, mask["h"] * 72.0,
               stroke=0, fill=1)
    c.setFont("Helvetica-Bold", 14)
    c.setFillColorRGB(0, 0, 0)
    for label in labels:
        c.drawString(label["x"] * 72.0, label["y"] * 72.0, label["text"])
    c.save()
    buf.seek(0)
    return PdfReader(buf).pages[0]


def _assemble(page_w_in, page_h_in, panels, out_path):
    page = PageObject.create_blank_page(width=page_w_in * 72.0, height=page_h_in * 72.0)
    labels = []
    masks = []

    for panel in panels:
        src_page, src_w_in, src_h_in = _load_cropped_page(panel["path"], panel.get("crop"))

        if "width" in panel:
            scale = panel["width"] / src_w_in
            placed_w = panel["width"]
            placed_h = src_h_in * scale
        else:
            scale = panel["height"] / src_h_in
            placed_h = panel["height"]
            placed_w = src_w_in * scale

        x_in = panel["x"]
        y_in = panel["top"] - placed_h if "top" in panel else panel["y"]

        page.merge_transformed_page(
            src_page,
            Transformation((scale, 0, 0, scale, x_in * 72.0, y_in * 72.0)),
        )

        if panel.get("mask_top"):
            mask_above = panel.get("mask_above", 0.0)
            masks.append({
                "x": x_in,
                "y": y_in + placed_h - panel["mask_top"],
                "w": placed_w,
                "h": panel["mask_top"] + mask_above,
            })
        labels.append(
            {
                "text": panel["label"],
                "x": panel.get("label_x", x_in),
                "y": panel["label_y"],
            }
        )

    overlay = _make_overlay(page_w_in * 72.0, page_h_in * 72.0, masks, labels)
    page.merge_page(overlay)

    writer = PdfWriter()
    writer.add_page(page)
    try:
        writer.compress_identical_objects(remove_identicals=True, remove_orphans=True)
    except Exception:
        pass
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as fh:
        writer.write(fh)
    print(f"Saved: {out_path}")


def _assemble_pages(page_specs, out_path):
    writer = PdfWriter()
    for spec in page_specs:
        page_w_in = spec["page_w"]
        page_h_in = spec["page_h"]
        page = PageObject.create_blank_page(width=page_w_in * 72.0, height=page_h_in * 72.0)
        labels = []
        masks = []

        for panel in spec["panels"]:
            src_page, src_w_in, src_h_in = _load_cropped_page(panel["path"], panel.get("crop"))

            if "width" in panel:
                scale = panel["width"] / src_w_in
                placed_w = panel["width"]
                placed_h = src_h_in * scale
            else:
                scale = panel["height"] / src_h_in
                placed_h = panel["height"]
                placed_w = src_w_in * scale

            x_in = panel["x"]
            y_in = panel["top"] - placed_h if "top" in panel else panel["y"]

            page.merge_transformed_page(
                src_page,
                Transformation((scale, 0, 0, scale, x_in * 72.0, y_in * 72.0)),
            )

            if panel.get("mask_top"):
                mask_above = panel.get("mask_above", 0.0)
                masks.append({
                    "x": x_in,
                    "y": y_in + placed_h - panel["mask_top"],
                    "w": placed_w,
                    "h": panel["mask_top"] + mask_above,
                })

            labels.append(
                {
                    "text": panel["label"],
                    "x": panel.get("label_x", x_in),
                    "y": panel["label_y"],
                }
            )

        overlay = _make_overlay(page_w_in * 72.0, page_h_in * 72.0, masks, labels)
        page.merge_page(overlay)
        writer.add_page(page)

    try:
        writer.compress_identical_objects(remove_identicals=True, remove_orphans=True)
    except Exception:
        pass

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as fh:
        writer.write(fh)
    print(f"Saved: {out_path}")


def _render_pdf_thumbnail(pdf_path, png_path, width_px):
    png_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["sips", "-s", "format", "png", str(pdf_path), "--out", str(png_path)],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if not png_path.exists():
        raise FileNotFoundError(f"Could not render preview for {pdf_path}")


def _assemble_image_pages(page_specs, out_path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    c = None
    for idx, page_spec in enumerate(page_specs):
        page_w_pt = page_spec["page_w"] * 72.0
        page_h_pt = page_spec["page_h"] * 72.0
        if c is None:
            c = canvas.Canvas(str(out_path), pagesize=(page_w_pt, page_h_pt))
        else:
            c.setPageSize((page_w_pt, page_h_pt))

        c.setFillColorRGB(1, 1, 1)
        c.rect(0, 0, page_w_pt, page_h_pt, stroke=0, fill=1)
        c.setFillColorRGB(0, 0, 0)
        c.setFont("Helvetica-Bold", 14)

        for panel in page_spec["panels"]:
            c.drawImage(
                ImageReader(str(panel["path"])),
                panel["x"] * 72.0,
                panel["y"] * 72.0,
                width=panel["width"] * 72.0,
                height=panel["height"] * 72.0,
                preserveAspectRatio=False,
                mask="auto",
            )
            c.drawString(panel["x"] * 72.0, panel["label_y"] * 72.0, panel["label"])

        if idx < len(page_specs) - 1:
            c.showPage()

    c.save()
    print(f"Saved: {out_path}")


def build_figure1_complete():
    page_w = 7.50
    page_h = 7.10
    margin = 0.18
    gutter = 0.18
    label_gap = 0.08
    label_headroom = 0.20

    top_cell_w = (page_w - 2 * margin - gutter) / 2.0
    top_cell_h = 3.35
    bottom_cell_w = page_w - 2 * margin
    bottom_cell_h = page_h - 2 * margin - gutter - top_cell_h

    top_row_y = page_h - margin - top_cell_h
    bottom_row_y = margin

    sources = {
        "A": FIG_MAIN / "Figure_1A_study_design.pdf",
        "B": FIG_MAIN / "Figure_1B_nrxn1_isoforms.pdf",
        "C": FIG_MAIN / "Figure_1C_density_pca_threepanel.pdf",
    }

    panels = []
    for label, cell_x, cell_y, cell_w, cell_h in (
        ("A", margin, top_row_y, top_cell_w, top_cell_h),
        ("B", margin + top_cell_w + gutter, top_row_y, top_cell_w, top_cell_h),
        ("C", margin, bottom_row_y, bottom_cell_w, bottom_cell_h),
    ):
        src_w, src_h = _page_dims_in(sources[label])
        placed_w, placed_h = _fit_within(src_w, src_h, cell_w, cell_h - label_headroom)
        panel_x = cell_x + (cell_w - placed_w) / 2.0
        panel_y = cell_y + (cell_h - label_headroom - placed_h) / 2.0
        panels.append(
            {
                "label": label,
                "path": sources[label],
                "x": panel_x,
                "y": panel_y,
                "width": placed_w,
                "height": placed_h,
                "label_y": panel_y + placed_h + label_gap,
            }
        )

    _assemble(page_w, page_h, panels, FIG_MAIN / "Figure_1_complete.pdf")


def build_figure2_complete():
    page_w = 7.50
    page_h = 10.05
    left = 0.18
    gap = 0.14

    row1_top = 9.15
    row2_top = 5.90
    row3_top = 2.90
    label_gap = 0.10

    panels = [
        {
            "label": "A",
            "path": FIG_MAIN / "Figure_2A_rna_heatmap_all_genes.pdf",
            "x": left,
            "top": row1_top,
            "width": 2.58,
            "label_y": row1_top + label_gap,
        },
        {
            "label": "B",
            "path": FIG_MAIN / "Figure_2B_rna_volcano.pdf",
            "x": left + 2.58 + gap,
            "top": row1_top,
            "width": 4.40,
            "label_y": row1_top + label_gap,
        },
        {
            "label": "C",
            "path": FIG_MAIN / "Figure_2C_rna_go_directional_categories.pdf",
            "x": left,
            "top": row2_top,
            "width": 4.18,
            "label_y": row2_top + label_gap,
        },
        {
            "label": "D",
            "path": FIG_MAIN / "Figure_2D_splicing_summary_bar.pdf",
            "x": left + 4.18 + gap,
            "top": 4.20,
            "width": 2.82,
            "label_y": 4.30,
        },
        {
            "label": "E",
            "path": FIG_MAIN / "Figure_2E_splicing_go.pdf",
            "x": left,
            "top": row3_top,
            "width": 4.10,
            "label_y": row3_top + label_gap,
        },
        {
            "label": "F",
            "path": FIG_MAIN / "Figure_2F_prc2_splicing.pdf",
            "x": left + 4.10 + gap,
            "top": row3_top,
            "width": 3.08,
            "crop": {"top": 0.12},
            "label_y": row3_top + label_gap,
        },
    ]

    _assemble(page_w, page_h, panels, FIG_MAIN / "Figure_2_complete.pdf")


def build_figure3_complete():
    page_w = 7.50
    page_h = 6.90
    gutter = 0.18
    row1_top = 6.35
    row2_top = 3.20
    label_gap = 0.10

    paths = {
        "A": FIG_MAIN / "Figure_3A_chip_annotation_diverging.pdf",
        "B": FIG_MAIN / "Figure_3B_rna_chip_scatter.pdf",
        "C": FIG_MAIN / "Figure_3C_chip_go.pdf",
        "D": FIG_MAIN / "Figure_3D_mechanistic_heatmap.pdf",
    }
    crops = {key: {} for key in paths}

    heights = {
        "A": 2.05,
        "B": 2.05,
        "C": 2.50,
        "D": 2.95,
    }
    widths = {}
    for key in ("A", "B", "C", "D"):
        w_in, h_in = _page_dims_in(paths[key], crops[key])
        widths[key] = heights[key] * (w_in / h_in)

    row1_total = widths["A"] + widths["B"] + gutter
    row2_total = widths["C"] + widths["D"] + gutter

    row1_left = (page_w - row1_total) / 2.0
    row2_left = (page_w - row2_total) / 2.0

    x_a = row1_left
    x_b = x_a + widths["A"] + gutter
    x_c = row2_left
    x_d = x_c + widths["C"] + gutter

    panels = [
        {
            "label": "A",
            "path": paths["A"],
            "crop": crops["A"],
            "x": x_a,
            "top": row1_top,
            "height": heights["A"],
            "label_y": row1_top + label_gap,
        },
        {
            "label": "B",
            "path": paths["B"],
            "crop": crops["B"],
            "x": x_b,
            "top": row1_top,
            "height": heights["B"],
            "label_y": row1_top + label_gap,
        },
        {
            "label": "C",
            "path": paths["C"],
            "crop": crops["C"],
            "x": x_c,
            "top": row2_top,
            "height": heights["C"],
            "label_y": row2_top + label_gap,
        },
        {
            "label": "D",
            "path": paths["D"],
            "crop": crops["D"],
            "x": x_d,
            "top": row2_top,
            "height": heights["D"],
            "label_y": row2_top + label_gap,
        },
    ]

    _assemble(page_w, page_h, panels, FIG_MAIN / "Figure_3_complete.pdf")


def build_figure4_complete():
    page_w = 7.50
    page_h = 5.60
    gutter = 0.18
    row1_top = 5.05
    row2_top = 2.55
    label_gap = 0.10

    paths = {
        "A": FIG_MAIN / "Figure_4A_promoter_dars.pdf",
        "B": FIG_MAIN / "Figure_4B_rna_atac_scatter.pdf",
        "C": FIG_MAIN / "Figure_4C_rna_atac_go.pdf",
        "D": FIG_MAIN / "Figure_4D_tobias_volcano.pdf",
        "E": FIG_MAIN / "Figure_4E_tobias_mechanistic.pdf",
    }
    crops = {key: {} for key in paths}

    heights = {
        "A": 1.55,
        "B": 1.92,
        "C": 2.12,
        "D": 2.18,
        "E": 2.18,
    }
    widths = {}
    for key in paths:
        w_in, h_in = _page_dims_in(paths[key], crops[key])
        widths[key] = heights[key] * (w_in / h_in)

    row1_total = widths["A"] + widths["B"] + widths["C"] + 2 * gutter
    row2_total = widths["D"] + widths["E"] + gutter
    row1_left = (page_w - row1_total) / 2.0
    row2_left = (page_w - row2_total) / 2.0

    x_a = row1_left
    x_b = x_a + widths["A"] + gutter
    x_c = x_b + widths["B"] + gutter
    x_d = row2_left
    x_e = x_d + widths["D"] + gutter

    panels = [
        {
            "label": "A",
            "path": paths["A"],
            "crop": crops["A"],
            "x": x_a,
            "top": row1_top,
            "height": heights["A"],
            "label_y": row1_top + label_gap,
        },
        {
            "label": "B",
            "path": paths["B"],
            "crop": crops["B"],
            "x": x_b,
            "top": row1_top,
            "height": heights["B"],
            "label_y": row1_top + label_gap,
        },
        {
            "label": "C",
            "path": paths["C"],
            "crop": crops["C"],
            "x": x_c,
            "top": row1_top,
            "height": heights["C"],
            "label_y": row1_top + label_gap,
        },
        {
            "label": "D",
            "path": paths["D"],
            "crop": crops["D"],
            "x": x_d,
            "top": row2_top,
            "height": heights["D"],
            "label_y": row2_top + label_gap,
        },
        {
            "label": "E",
            "path": paths["E"],
            "crop": crops["E"],
            "x": x_e,
            "top": row2_top,
            "height": heights["E"],
            "label_y": row2_top + label_gap,
        },
    ]

    _assemble(page_w, page_h, panels, FIG_MAIN / "Figure_4_complete.pdf")


def build_figure5_complete():
    page_w = 8.30
    gutter = 0.18
    label_gap = 0.08
    page_h = 7.25

    sources = {
        "A": FIG_MAIN / "Figure_5A_triple_overlap_inclusive.pdf",
        "B": FIG_MAIN / "Figure_5B_integration_heatmap.pdf",
        "C": FIG_MAIN / "Figure_5C_schematic_integrative_model.pdf",
    }

    ratio_a = _page_dims_in(sources["A"])[0] / _page_dims_in(sources["A"])[1]
    ratio_b = _page_dims_in(sources["B"])[0] / _page_dims_in(sources["B"])[1]

    height_a = 1.55
    height_b = 1.35
    width_a = height_a * ratio_a
    width_b = height_b * ratio_b
    row1_total = width_a + width_b + gutter
    row1_left = (page_w - row1_total) / 2.0
    row1_top = 6.10

    width_c = 7.70
    row2_top = 4.35

    panels = [
        {
            "label": "A",
            "path": sources["A"],
            "x": row1_left,
            "top": row1_top,
            "height": height_a,
            "label_y": row1_top + label_gap,
        },
        {
            "label": "B",
            "path": sources["B"],
            "x": row1_left + width_a + gutter,
            "top": row1_top,
            "height": height_b,
            "label_y": row1_top + label_gap,
        },
        {
            "label": "C",
            "path": sources["C"],
            "x": (page_w - width_c) / 2.0,
            "top": row2_top,
            "width": width_c,
            "label_y": row2_top + label_gap,
        },
    ]

    _assemble(page_w, page_h, panels, FIG_MAIN / "Figure_5_complete.pdf")


def build_figure_s1_qc():
    out_pdf = FIG_SUPP / "Figure_S1_qc.pdf"
    page_w = 11.20
    page_h = 8.95

    panels = [
        {
            "label": "A",
            "path": FIG_SUPP_QC / "Figure_S1A_rna_sample_correlation.pdf",
            "x": 0.28,
            "top": 8.14,
            "width": 5.25,
            "label_x": 0.16,
            "label_y": 7.98,
        },
        {
            "label": "B",
            "path": FIG_SUPP_QC / "Figure_S1B_chip_frip_scores.pdf",
            "x": 5.72,
            "top": 8.14,
            "width": 5.18,
            "label_x": 5.60,
            "label_y": 7.98,
        },
        {
            "label": "C",
            "path": FIG_SUPP_QC / "Figure_S1C_atac_fragment_composition.pdf",
            "x": 0.28,
            "top": 3.45,
            "width": 3.46,
            "label_x": 0.16,
            "label_y": 3.29,
        },
        {
            "label": "D",
            "path": FIG_SUPP_QC / "Figure_S1D_atac_tss_enrichment_scores.pdf",
            "x": 3.87,
            "top": 3.45,
            "width": 3.34,
            "label_x": 3.75,
            "label_y": 3.29,
        },
        {
            "label": "E",
            "path": FIG_SUPP_QC / "Figure_S1E_atac_frip_scores.pdf",
            "x": 7.47,
            "top": 3.45,
            "width": 3.34,
            "label_x": 7.35,
            "label_y": 3.29,
        },
    ]

    _assemble(page_w, page_h, panels, out_pdf)


def build_figure_s2_female_only():
    page_w = 12.00
    page_h = 9.45
    gutter = 0.28
    row1_top = 8.90
    row2_top = 4.15
    label_gap = 0.10

    sources = {
        "A": FIG_SUPP_RNA / "Figure_2C_rna_heatmap_all_genes_female_only.pdf",
        "B": FIG_SUPP_RNA / "Figure_2A_rna_volcano_female_only.pdf",
        "C": FIG_SUPP_RNA / "Figure_2B_rna_go_directional_female_only.pdf",
        "D": FIG_SUPP_MULTI / "Figure_4B_tobias_volcano_female_only.pdf",
    }

    crops = {
        "A": {"left": 0.18, "right": 0.48, "top": 0.04, "bottom": 0.42},
        "B": {"left": 0.10, "right": 0.05, "top": 0.03, "bottom": 0.12},
        "C": {"left": 0.42, "right": 0.22, "top": 0.05, "bottom": 0.16},
        "D": {"left": 0.10, "right": 0.05, "top": 0.03, "bottom": 0.12},
    }
    panels = [
        {
            "label": "A",
            "path": sources["A"],
            "crop": crops["A"],
            "x": 0.38,
            "top": row1_top,
            "width": 3.95,
            "label_y": row1_top + label_gap,
        },
        {
            "label": "B",
            "path": sources["B"],
            "crop": crops["B"],
            "x": 0.38 + 3.95 + gutter,
            "top": row1_top,
            "width": 5.95,
            "label_y": row1_top + label_gap,
        },
        {
            "label": "C",
            "path": sources["C"],
            "crop": crops["C"],
            "x": 0.58,
            "top": row2_top,
            "width": 4.00,
            "label_x": 0.46,
            "label_y": row2_top + 0.06,
        },
        {
            "label": "D",
            "path": sources["D"],
            "crop": crops["D"],
            "x": 0.58 + 4.00 + gutter,
            "top": row2_top,
            "width": 5.45,
            "label_x": 4.76,
            "label_y": row2_top + 0.06,
        },
    ]

    out_pdf = FIG_SUPP / "Figure_S2_female_only_replication.pdf"
    _assemble(page_w, page_h, panels, out_pdf)


def build_figure1c_threepanel_pca():
    page_w = 12.6
    page_h = 4.6
    margin = 0.28
    gutter = 0.20
    label_gap = 0.08
    label_headroom = 0.20

    cell_w = (page_w - 2 * margin - 2 * gutter) / 3.0
    cell_h = page_h - 2 * margin
    y0 = margin

    sources = {
        "A": FIG_MAIN / "archive" / "Figure_1D_rna_pca.pdf",
        "B": FIG_MAIN / "archive" / "Figure_1D_atac_pca.pdf",
        "C": FIG_MAIN / "archive" / "Figure_1D_chip_pca.pdf",
    }

    panels = []
    for idx, label in enumerate(("A", "B", "C")):
        path = sources[label]
        src_w, src_h = _page_dims_in(path)
        placed_w, placed_h = _fit_within(src_w, src_h, cell_w, cell_h - label_headroom)
        cell_x = margin + idx * (cell_w + gutter)
        panel_x = cell_x + (cell_w - placed_w) / 2.0
        panel_y = y0 + (cell_h - label_headroom - placed_h) / 2.0
        panels.append(
            {
                "label": label,
                "path": path,
                "x": panel_x,
                "y": panel_y,
                "width": placed_w,
                "height": placed_h,
                "label_y": panel_y + placed_h + label_gap,
            }
        )

    out_pdf = FIG_MAIN / "archive" / "Figure_1C_density_pca_threepanel_from_archived_panels.pdf"
    _assemble(page_w, page_h, panels, out_pdf)


def build_figure_s3_qpcr():
    page_w = 7.50
    page_h = 3.85
    margin = 0.18
    gutter = 0.18
    label_gap = 0.08
    label_headroom = 0.22

    cell_w = (page_w - 2 * margin - gutter) / 2.0
    cell_h = page_h - 2 * margin
    row_y = margin

    sources = {
        "A": FIG_SUPP / "Figure_S3A_SMAD7_qPCR.pdf",
        "B": FIG_SUPP / "Figure_S3B_Ecadherin_qPCR.pdf",
    }

    panels = []
    for label, cell_x in (("A", margin), ("B", margin + cell_w + gutter)):
        src_w, src_h = _page_dims_in(sources[label])
        placed_w, placed_h = _fit_within(src_w, src_h, cell_w, cell_h - label_headroom)
        panel_x = cell_x + (cell_w - placed_w) / 2.0
        panel_y = row_y + (cell_h - label_headroom - placed_h) / 2.0
        panels.append(
            {
                "label": label,
                "path": sources[label],
                "x": panel_x,
                "y": panel_y,
                "width": placed_w,
                "height": placed_h,
                "label_y": panel_y + placed_h + label_gap,
            }
        )

    out_pdf = FIG_SUPP / "Figure_S3_complete.pdf"
    _assemble(page_w, page_h, panels, out_pdf)

if __name__ == "__main__":
    build_figure1_complete()
    build_figure2_complete()
    build_figure3_complete()
    build_figure4_complete()
    build_figure5_complete()
