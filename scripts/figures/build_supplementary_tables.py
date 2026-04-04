#!/usr/bin/env python3
"""
Build manuscript-facing supplementary table workbooks in the final modality-
based ST1-ST6 layout without relying on external spreadsheet libraries.

Inputs:
  outputs/ from the assay-level rebuild scripts
  data/processed_inputs/qc/
  data/processed_inputs/qpcr/

Outputs:
  outputs/manuscript_supplementary_tables/tables/
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile


def find_project_root() -> Path:
    here = Path(__file__).resolve()
    for parent in [here.parent, *here.parents]:
        if (parent / ".nrxn1_project_root").exists():
            return parent
    raise RuntimeError("Could not find project root sentinel .nrxn1_project_root")


PROJECT_ROOT = find_project_root()
RESULTS_DIR = PROJECT_ROOT / "outputs"
OUT_DIR = RESULTS_DIR / "manuscript_supplementary_tables" / "tables"


@dataclass(frozen=True)
class SheetSpec:
    name: str
    source: Path
    delimiter: str = "\t"


@dataclass(frozen=True)
class WorkbookSpec:
    output_name: str
    sheets: tuple[SheetSpec, ...]


WORKBOOKS: tuple[WorkbookSpec, ...] = (
    WorkbookSpec(
        output_name="ST1_Sample_Metadata.xlsx",
        sheets=(
            SheetSpec(
                "Samples",
                PROJECT_ROOT / "data" / "processed_inputs" / "qc" / "Table_S1_metadata.csv",
                delimiter=",",
            ),
        ),
    ),
    WorkbookSpec(
        output_name="ST2_RNAseq_Main.xlsx",
        sheets=(
            SheetSpec(
                "Gene_DE",
                RESULTS_DIR
                / "rna_seq"
                / "rna_main_all_controls"
                / "tables"
                / "rna_main_all_controls_gene_results_sig_lfc.tsv",
            ),
            SheetSpec(
                "Transcript_DE",
                RESULTS_DIR
                / "rna_seq"
                / "rna_main_all_controls"
                / "tables"
                / "rna_main_all_controls_aggregated_transcript_results_sig.tsv",
            ),
            SheetSpec(
                "Splicing_Events",
                RESULTS_DIR / "splicing" / "tables" / "significant_events_jcec.tsv",
            ),
            SheetSpec(
                "Splicing_GO",
                RESULTS_DIR / "splicing_go" / "tables" / "splicing_topgo_bp_full.tsv",
            ),
            SheetSpec(
                "Density_Sensitivity",
                RESULTS_DIR / "figures" / "tables" / "ST6_density_sensitivity_summary.tsv",
            ),
            SheetSpec(
                "RNA_GO_Up",
                RESULTS_DIR / "go_enrichment" / "tables" / "rna_main_up_topgo_bp_full.tsv",
            ),
            SheetSpec(
                "RNA_GO_Down",
                RESULTS_DIR / "go_enrichment" / "tables" / "rna_main_down_topgo_bp_full.tsv",
            ),
        ),
    ),
    WorkbookSpec(
        output_name="ST3_ATAC_Main.xlsx",
        sheets=(
            SheetSpec(
                "All_DARs",
                RESULTS_DIR
                / "atac_seq"
                / "atac_main_local_controls"
                / "tables"
                / "differential_accessibility_standardized.tsv",
            ),
            SheetSpec(
                "Promoter_DARs",
                RESULTS_DIR
                / "atac_seq"
                / "atac_main_local_controls"
                / "tables"
                / "promoter_dars_standardized.tsv",
            ),
            SheetSpec(
                "TOBIAS",
                RESULTS_DIR / "figures" / "tables" / "Figure_4B_tobias_threshold_q95.tsv",
            ),
            SheetSpec(
                "RNA_ATAC_GO_Up",
                RESULTS_DIR / "go_enrichment" / "tables" / "rna_atac_main_up_topgo_bp_full.tsv",
            ),
            SheetSpec(
                "RNA_ATAC_GO_Down",
                RESULTS_DIR / "go_enrichment" / "tables" / "rna_atac_main_down_topgo_bp_full.tsv",
            ),
        ),
    ),
    WorkbookSpec(
        output_name="ST4_H3K27me3_Main.xlsx",
        sheets=(
            SheetSpec(
                "All_Regions",
                RESULTS_DIR
                / "chip_seq"
                / "chip_main_all_controls"
                / "tables"
                / "differential_regions_standardized.tsv",
            ),
            SheetSpec(
                "Promoter_Regions",
                RESULTS_DIR
                / "chip_seq"
                / "chip_main_all_controls"
                / "tables"
                / "promoter_regions_standardized.tsv",
            ),
            SheetSpec(
                "RNA_H3K27me3_GO",
                RESULTS_DIR / "go_enrichment" / "tables" / "rna_chip_main_topgo_bp_full.tsv",
            ),
        ),
    ),
    WorkbookSpec(
        output_name="ST5_Integration_Main.xlsx",
        sheets=(
            SheetSpec(
                "RNA_ATAC",
                RESULTS_DIR / "integration" / "tables" / "rna_atac_main_overlap.tsv",
            ),
            SheetSpec(
                "RNA_H3K27me3",
                RESULTS_DIR / "integration" / "tables" / "rna_chip_main_overlap.tsv",
            ),
            SheetSpec(
                "Triple_Overlap",
                RESULTS_DIR / "integration" / "tables" / "triple_overlap_main.tsv",
            ),
            SheetSpec(
                "Heatmap_Matrix",
                RESULTS_DIR / "integration" / "tables" / "triple_overlap_main_heatmap_matrix.tsv",
            ),
            SheetSpec(
                "Triple_Overlap_GO",
                RESULTS_DIR / "go_enrichment" / "tables" / "triple_main_topgo_bp_full.tsv",
            ),
        ),
    ),
    WorkbookSpec(
        output_name="ST6_qPCR_Validation.xlsx",
        sheets=(
            SheetSpec(
                "qPCR_Stats",
                PROJECT_ROOT / "data" / "processed_inputs" / "qpcr" / "qPCR_stats_summary.txt",
            ),
        ),
    ),
)


def col_letter(index: int) -> str:
    letters: list[str] = []
    while index:
        index, rem = divmod(index - 1, 26)
        letters.append(chr(65 + rem))
    return "".join(reversed(letters))


def read_table(path: Path, delimiter: str) -> tuple[list[str], list[list[str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        rows = list(reader)
    if not rows:
        return [], []
    return rows[0], rows[1:]


def maybe_numeric(value: str) -> str | None:
    stripped = value.strip()
    if stripped == "":
        return None
    lowered = stripped.lower()
    if lowered in {"nan", "inf", "-inf"}:
        return None
    try:
        float(stripped)
    except ValueError:
        return None
    return stripped


def cell_xml(ref: str, value: str, header: bool = False) -> str:
    style = ' s="1"' if header else ""
    numeric = maybe_numeric(value)
    if numeric is not None:
        return f'<c r="{ref}"{style}><v>{numeric}</v></c>'
    return f'<c r="{ref}" t="inlineStr"{style}><is><t>{escape(value)}</t></is></c>'


def build_sheet_xml(columns: list[str], rows: list[list[str]]) -> str:
    last_col = col_letter(max(1, len(columns)))
    last_row = max(1, len(rows) + 1)

    preview_rows = rows[:250]
    col_widths: list[int] = []
    for idx, col_name in enumerate(columns):
        max_len = len(str(col_name))
        for row in preview_rows:
            if idx < len(row):
                max_len = max(max_len, len(str(row[idx])))
        col_widths.append(max(12, min(32, max_len + 2)))

    xml_rows: list[str] = []
    header_cells = [
        cell_xml(f"{col_letter(idx)}1", col_name, header=True)
        for idx, col_name in enumerate(columns, start=1)
    ]
    xml_rows.append(f'<row r="1">{"".join(header_cells)}</row>')

    for row_num, row in enumerate(rows, start=2):
        row_cells = []
        for idx, value in enumerate(row, start=1):
            row_cells.append(cell_xml(f"{col_letter(idx)}{row_num}", value))
        xml_rows.append(f'<row r="{row_num}">{"".join(row_cells)}</row>')

    col_defs = [
        f'<col min="{idx}" max="{idx}" width="{width}" customWidth="1"/>'
        for idx, width in enumerate(col_widths, start=1)
    ]

    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<dimension ref="A1:{last_col}{last_row}"/>'
        '<sheetViews><sheetView workbookViewId="0">'
        '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
        '</sheetView></sheetViews>'
        '<sheetFormatPr defaultRowHeight="15"/>'
        f'<cols>{"".join(col_defs)}</cols>'
        f'<sheetData>{"".join(xml_rows)}</sheetData>'
        '</worksheet>'
    )


def build_workbook(spec: WorkbookSpec) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUT_DIR / spec.output_name

    sheet_payloads: list[tuple[str, str]] = []
    for sheet in spec.sheets:
        columns, rows = read_table(sheet.source, sheet.delimiter)
        sheet_payloads.append((sheet.name, build_sheet_xml(columns, rows)))

    content_overrides = [
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
        '<Override PartName="/xl/styles.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
        '<Override PartName="/docProps/core.xml" '
        'ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>',
        '<Override PartName="/docProps/app.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>',
    ]
    for idx in range(1, len(sheet_payloads) + 1):
        content_overrides.append(
            f'<Override PartName="/xl/worksheets/sheet{idx}.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        )

    workbook_sheets = []
    workbook_rels = []
    for idx, (name, _) in enumerate(sheet_payloads, start=1):
        workbook_sheets.append(
            f'<sheet name="{escape(name)}" sheetId="{idx}" r:id="rId{idx}"/>'
        )
        workbook_rels.append(
            f'<Relationship Id="rId{idx}" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
            f'Target="worksheets/sheet{idx}.xml"/>'
        )
    workbook_rels.append(
        f'<Relationship Id="rId{len(sheet_payloads) + 1}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/>'
    )

    with ZipFile(output_path, "w", compression=ZIP_DEFLATED) as zf:
        zf.writestr(
            "[Content_Types].xml",
            """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  """
            + "\n  ".join(content_overrides)
            + "\n</Types>",
        )
        zf.writestr(
            "_rels/.rels",
            """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>""",
        )
        zf.writestr(
            "docProps/core.xml",
            f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"
 xmlns:dc="http://purl.org/dc/elements/1.1/"
 xmlns:dcterms="http://purl.org/dc/terms/"
 xmlns:dcmitype="http://purl.org/dc/dcmitype/"
 xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:title>{escape(spec.output_name.replace('.xlsx', ''))}</dc:title>
  <dc:creator>Codex</dc:creator>
</cp:coreProperties>""",
        )
        zf.writestr(
            "docProps/app.xml",
            """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"
 xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>Codex</Application>
</Properties>""",
        )
        zf.writestr(
            "xl/workbook.xml",
            """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    """
            + "\n    ".join(workbook_sheets)
            + "\n  </sheets>\n</workbook>",
        )
        zf.writestr(
            "xl/_rels/workbook.xml.rels",
            """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  """
            + "\n  ".join(workbook_rels)
            + "\n</Relationships>",
        )
        zf.writestr(
            "xl/styles.xml",
            """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="2">
    <font><sz val="11"/><name val="Calibri"/></font>
    <font><b/><sz val="11"/><name val="Calibri"/></font>
  </fonts>
  <fills count="2">
    <fill><patternFill patternType="none"/></fill>
    <fill><patternFill patternType="gray125"/></fill>
  </fills>
  <borders count="1">
    <border><left/><right/><top/><bottom/><diagonal/></border>
  </borders>
  <cellStyleXfs count="1">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0"/>
  </cellStyleXfs>
  <cellXfs count="2">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
    <xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>
  </cellXfs>
  <cellStyles count="1">
    <cellStyle name="Normal" xfId="0" builtinId="0"/>
  </cellStyles>
</styleSheet>""",
        )
        for idx, (_, xml) in enumerate(sheet_payloads, start=1):
            zf.writestr(f"xl/worksheets/sheet{idx}.xml", xml)

    print(f"Saved: {output_path}")


def main() -> None:
    for workbook in WORKBOOKS:
        build_workbook(workbook)


if __name__ == "__main__":
    main()
